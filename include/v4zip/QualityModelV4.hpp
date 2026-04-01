/**
 * QualityModelV4.hpp - Adaptive alphabet quality model
 *
 * Key insight: Real FASTQ files often use a small subset of the 94 possible
 * quality values, and quality scores may have little temporal correlation
 * (especially in corrected reads).
 *
 * This model:
 * 1. First pass: Discovers the actual alphabet used and builds a compact mapping
 * 2. Encodes quality values using the reduced alphabet (can be as low as 4.32 bits
 *    for 20 symbols vs 6.55 bits for 94)
 * 3. Uses adaptive order-0 model when context doesn't help
 * 4. Falls back to order-1 when beneficial
 *
 * Expected: Near-entropy compression (~4.3-4.5 bits for typical Illumina/PacBio)
 *
 * Encoded format:
 * +---------------------+--------+----------------------------------+
 * | Field               | Size   | Description                      |
 * +---------------------+--------+----------------------------------+
 * | Alphabet size       | 1 byte | Number of unique quality values  |
 * | Alphabet mapping    | N bytes| Maps compact idx -> original val |
 * | rANS state          | 4 bytes| Initial decoder state (BE)       |
 * | rANS data           | var    | Encoded quality symbols          |
 * +---------------------+--------+----------------------------------+
 *
 * The rANS encoder/decoder uses adaptive probability updates (Laplace smoothing)
 * with rescaling when counts exceed COUNT_RESCALE_THRESHOLD.
 */

#ifndef V4ZIP_QUALITYMODELV4_HPP
#define V4ZIP_QUALITYMODELV4_HPP

#include "ContextModel.hpp"  // For PROB_BITS, PROB_SCALE
#include "ArithmeticCoder.hpp"
#include <array>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdint>
#include <cmath>
#include <unordered_map>

namespace v4zip {

// Quality score alphabet: ASCII 33 ('!') through 126 ('~') = 94 symbols
constexpr int QUALITY_ALPHABET_SIZE = 94;
constexpr int MAX_QUALITY_ALPHABET = QUALITY_ALPHABET_SIZE;
constexpr int QUALITY_OFFSET = 33;

inline int qualityCharToIndex(char c) {
    return static_cast<int>(static_cast<unsigned char>(c)) - QUALITY_OFFSET;
}

inline char indexToQualityChar(int idx) {
    return static_cast<char>(idx + QUALITY_OFFSET);
}

// Adaptive model count rescale threshold (halves counts to prevent overflow)
constexpr uint32_t COUNT_RESCALE_THRESHOLD = 65535;

/**
 * Adaptive order-0 model for a reduced alphabet
 * Uses arithmetic coding with adaptive probability updates
 */
template<int MaxAlphabetSize>
class AdaptiveOrder0Model {
public:
    AdaptiveOrder0Model() : alphabetSize_(MaxAlphabetSize) {
        reset();
    }

    void setAlphabetSize(int size) {
        alphabetSize_ = std::min(size, MaxAlphabetSize);
        reset();
    }

    void reset() {
        // Laplace smoothing: start with count of 1 for each symbol
        for (int i = 0; i < alphabetSize_; i++) {
            counts_[i] = 1;
        }
        for (int i = alphabetSize_; i < MaxAlphabetSize; i++) {
            counts_[i] = 0;
        }
        total_ = alphabetSize_;
    }

    // Get CDF for encoding/decoding
    std::vector<uint32_t> getCDF() const {
        return normalizeCDFVec(counts_.data(), alphabetSize_, PROB_SCALE);
    }

    void update(int symbol) {
        if (symbol >= 0 && symbol < alphabetSize_) {
            counts_[symbol]++;
            total_++;

            // Rescale if counts get too large
            if (total_ > COUNT_RESCALE_THRESHOLD) {
                rescale();
            }
        }
    }

    int alphabetSize() const { return alphabetSize_; }
    const std::array<uint32_t, MaxAlphabetSize>& counts() const { return counts_; }
    uint32_t total() const { return total_; }

private:
    void rescale() {
        rescaleCounts(counts_, alphabetSize_, total_);
    }

    std::array<uint32_t, MaxAlphabetSize> counts_;
    uint32_t total_;
    int alphabetSize_;
};

/**
 * Adaptive order-1 model for a reduced alphabet
 * Combines V4's adaptive alphabet with order-1 context (P(q_i | q_{i-1}))
 * Falls back to order-0 when context has insufficient observations
 */
template<int MaxAlphabetSize>
class AdaptiveOrder1Model {
public:
    static constexpr uint32_t MIN_CONTEXT_COUNT = 8;  // Minimum observations to trust order-1

    AdaptiveOrder1Model() : alphabetSize_(MaxAlphabetSize) {
        reset();
    }

    void setAlphabetSize(int size) {
        alphabetSize_ = std::min(size, MaxAlphabetSize);
        reset();
    }

    void reset() {
        // Initialize order-0 model
        for (int i = 0; i < alphabetSize_; i++) {
            order0Counts_[i] = 1;
        }
        for (int i = alphabetSize_; i < MaxAlphabetSize; i++) {
            order0Counts_[i] = 0;
        }
        order0Total_ = alphabetSize_;

        // Initialize order-1 models (one per context symbol)
        for (int ctx = 0; ctx < MaxAlphabetSize; ctx++) {
            for (int i = 0; i < alphabetSize_; i++) {
                order1Counts_[ctx][i] = 1;
            }
            for (int i = alphabetSize_; i < MaxAlphabetSize; i++) {
                order1Counts_[ctx][i] = 0;
            }
            order1Totals_[ctx] = alphabetSize_;
        }
    }

    // Get CDF for encoding/decoding given previous symbol context
    // prevSymbol < 0 means no context (first symbol in read)
    std::vector<uint32_t> getCDF(int prevSymbol) const {
        const uint32_t* counts;

        // Use order-1 if we have context and enough observations
        if (prevSymbol >= 0 && prevSymbol < alphabetSize_ &&
            order1Totals_[prevSymbol] >= MIN_CONTEXT_COUNT) {
            counts = order1Counts_[prevSymbol].data();
        } else {
            counts = order0Counts_.data();
        }

        return normalizeCDFVec(counts, alphabetSize_, PROB_SCALE);
    }

    void update(int prevSymbol, int symbol) {
        if (symbol < 0 || symbol >= alphabetSize_) return;

        // Always update order-0
        order0Counts_[symbol]++;
        order0Total_++;
        if (order0Total_ > COUNT_RESCALE_THRESHOLD) {
            rescaleOrder0();
        }

        // Update order-1 if we have context
        if (prevSymbol >= 0 && prevSymbol < alphabetSize_) {
            order1Counts_[prevSymbol][symbol]++;
            order1Totals_[prevSymbol]++;
            if (order1Totals_[prevSymbol] > COUNT_RESCALE_THRESHOLD) {
                rescaleOrder1(prevSymbol);
            }
        }
    }

    int alphabetSize() const { return alphabetSize_; }

private:
    void rescaleOrder0() {
        rescaleCounts(order0Counts_, alphabetSize_, order0Total_);
    }

    void rescaleOrder1(int ctx) {
        rescaleCounts(order1Counts_[ctx], alphabetSize_, order1Totals_[ctx]);
    }

    // Order-0 counts
    std::array<uint32_t, MaxAlphabetSize> order0Counts_;
    uint32_t order0Total_;

    // Order-1 counts: order1Counts_[prev][curr]
    std::array<std::array<uint32_t, MaxAlphabetSize>, MaxAlphabetSize> order1Counts_;
    std::array<uint32_t, MaxAlphabetSize> order1Totals_;

    int alphabetSize_;
};

/**
 * PPM*-blending quality model (orders 0/1/2/3)
 * Blends predictions from 4 context orders weighted by reliability.
 * Order-2 and order-3 use hash maps for memory efficiency.
 */
template<int MaxAlphabetSize>
class AdaptiveQualityPPM {
public:
    static constexpr uint32_t ORDER1_THRESHOLD = 8;
    static constexpr uint32_t ORDER2_THRESHOLD = 16;
    static constexpr uint32_t ORDER3_THRESHOLD = 16;

    AdaptiveQualityPPM() : alphabetSize_(MaxAlphabetSize) {
        reset();
    }

    void setAlphabetSize(int size) {
        alphabetSize_ = std::min(size, MaxAlphabetSize);
        reset();
    }

    void reset() {
        for (int i = 0; i < alphabetSize_; i++) order0Counts_[i] = 1;
        for (int i = alphabetSize_; i < MaxAlphabetSize; i++) order0Counts_[i] = 0;
        order0Total_ = alphabetSize_;

        for (int ctx = 0; ctx < MaxAlphabetSize; ctx++) {
            for (int i = 0; i < alphabetSize_; i++) order1Counts_[ctx][i] = 1;
            for (int i = alphabetSize_; i < MaxAlphabetSize; i++) order1Counts_[ctx][i] = 0;
            order1Totals_[ctx] = alphabetSize_;
        }

        order2Contexts_.clear();
        order3Contexts_.clear();
    }

    // Initialize order-0 counts from a frequency histogram (warm-start)
    // freqs should have alphabetSize_ entries, scaled so sum ~= 256
    void initFromFrequencies(const std::vector<uint32_t>& freqs) {
        if (static_cast<int>(freqs.size()) != alphabetSize_) return;
        order0Total_ = 0;
        for (int i = 0; i < alphabetSize_; i++) {
            order0Counts_[i] = std::max(1u, freqs[i]);
            order0Total_ += order0Counts_[i];
        }
        // Also warm-start order-1 contexts with the same marginal
        for (int ctx = 0; ctx < alphabetSize_; ctx++) {
            order1Totals_[ctx] = 0;
            for (int i = 0; i < alphabetSize_; i++) {
                order1Counts_[ctx][i] = std::max(1u, freqs[i]);
                order1Totals_[ctx] += order1Counts_[ctx][i];
            }
        }
    }

    // PPM* blended CDF from orders 0, 1, 2, 3
    std::vector<uint32_t> getCDF(int prev3, int prev2, int prev1) const {
        double blended[MaxAlphabetSize];
        for (int i = 0; i < alphabetSize_; i++) blended[i] = 0.0;
        double totalWeight = 0.0;

        // Order-3
        if (prev3 >= 0 && prev2 >= 0 && prev1 >= 0 &&
            prev3 < alphabetSize_ && prev2 < alphabetSize_ && prev1 < alphabetSize_) {
            uint64_t key = static_cast<uint64_t>(prev3) * alphabetSize_ * alphabetSize_
                         + static_cast<uint64_t>(prev2) * alphabetSize_ + prev1;
            auto it = order3Contexts_.find(key);
            if (it != order3Contexts_.end() && it->second.total >= ORDER3_THRESHOLD) {
                double weight = 8.0 * std::log(
                    static_cast<double>(it->second.total) / ORDER3_THRESHOLD + 1.0);
                for (int i = 0; i < alphabetSize_; i++) {
                    blended[i] += weight * static_cast<double>(it->second.counts[i]) / it->second.total;
                }
                totalWeight += weight;
            }
        }

        // Order-2
        if (prev2 >= 0 && prev1 >= 0 && prev2 < alphabetSize_ && prev1 < alphabetSize_) {
            uint32_t key = static_cast<uint32_t>(prev2) * alphabetSize_ + prev1;
            auto it = order2Contexts_.find(key);
            if (it != order2Contexts_.end() && it->second.total >= ORDER2_THRESHOLD) {
                double weight = 4.0 * std::log(
                    static_cast<double>(it->second.total) / ORDER2_THRESHOLD + 1.0);
                for (int i = 0; i < alphabetSize_; i++) {
                    blended[i] += weight * static_cast<double>(it->second.counts[i]) / it->second.total;
                }
                totalWeight += weight;
            }
        }

        // Order-1
        if (prev1 >= 0 && prev1 < alphabetSize_ &&
            order1Totals_[prev1] >= ORDER1_THRESHOLD) {
            double weight = 2.0 * std::log(
                static_cast<double>(order1Totals_[prev1]) / ORDER1_THRESHOLD + 1.0);
            for (int i = 0; i < alphabetSize_; i++) {
                blended[i] += weight * static_cast<double>(order1Counts_[prev1][i]) / order1Totals_[prev1];
            }
            totalWeight += weight;
        }

        // Order-0 (always contributes)
        {
            double weight = 1.0 * std::log(static_cast<double>(order0Total_) + 1.0);
            for (int i = 0; i < alphabetSize_; i++) {
                blended[i] += weight * static_cast<double>(order0Counts_[i]) / order0Total_;
            }
            totalWeight += weight;
        }

        // Convert to pseudo-counts and normalize
        uint32_t pseudoCounts[MaxAlphabetSize];
        for (int i = 0; i < alphabetSize_; i++) {
            double p = blended[i] / totalWeight;
            pseudoCounts[i] = static_cast<uint32_t>(p * 1000000.0 + 0.5);
            if (pseudoCounts[i] == 0 && blended[i] > 0) pseudoCounts[i] = 1;
        }

        return normalizeCDFVec(pseudoCounts, alphabetSize_, PROB_SCALE);
    }

    // Backward-compatible 2-arg version (prev3 = -1)
    std::vector<uint32_t> getCDF(int prev2, int prev1) const {
        return getCDF(-1, prev2, prev1);
    }

    void update(int prev3, int prev2, int prev1, int symbol) {
        if (symbol < 0 || symbol >= alphabetSize_) return;

        order0Counts_[symbol]++;
        order0Total_++;
        if (order0Total_ > COUNT_RESCALE_THRESHOLD) rescaleOrder0();

        if (prev1 >= 0 && prev1 < alphabetSize_) {
            order1Counts_[prev1][symbol]++;
            order1Totals_[prev1]++;
            if (order1Totals_[prev1] > COUNT_RESCALE_THRESHOLD) rescaleOrder1(prev1);
        }

        if (prev2 >= 0 && prev1 >= 0 && prev2 < alphabetSize_ && prev1 < alphabetSize_) {
            uint32_t key2 = static_cast<uint32_t>(prev2) * alphabetSize_ + prev1;
            auto& ctx2 = order2Contexts_[key2];
            if (ctx2.counts.empty()) {
                ctx2.counts.assign(alphabetSize_, 1);
                ctx2.total = alphabetSize_;
            }
            ctx2.counts[symbol]++;
            ctx2.total++;
            if (ctx2.total > COUNT_RESCALE_THRESHOLD) rescaleOrder2(key2);
        }

        if (prev3 >= 0 && prev2 >= 0 && prev1 >= 0 &&
            prev3 < alphabetSize_ && prev2 < alphabetSize_ && prev1 < alphabetSize_) {
            uint64_t key3 = static_cast<uint64_t>(prev3) * alphabetSize_ * alphabetSize_
                          + static_cast<uint64_t>(prev2) * alphabetSize_ + prev1;
            auto& ctx3 = order3Contexts_[key3];
            if (ctx3.counts.empty()) {
                ctx3.counts.assign(alphabetSize_, 1);
                ctx3.total = alphabetSize_;
            }
            ctx3.counts[symbol]++;
            ctx3.total++;
            if (ctx3.total > COUNT_RESCALE_THRESHOLD) rescaleOrder3(key3);
        }
    }

    // Backward-compatible 3-arg version (prev3 = -1)
    void update(int prev2, int prev1, int symbol) {
        update(-1, prev2, prev1, symbol);
    }

    int alphabetSize() const { return alphabetSize_; }

private:
    struct SparseCtx {
        std::vector<uint32_t> counts;
        uint32_t total = 0;
    };

    void rescaleOrder0() {
        rescaleCounts(order0Counts_, alphabetSize_, order0Total_);
    }

    void rescaleOrder1(int ctx) {
        rescaleCounts(order1Counts_[ctx], alphabetSize_, order1Totals_[ctx]);
    }

    void rescaleOrder2(uint32_t key) {
        auto& ctx = order2Contexts_[key];
        rescaleCounts(ctx.counts, alphabetSize_, ctx.total);
    }

    void rescaleOrder3(uint64_t key) {
        auto& ctx = order3Contexts_[key];
        rescaleCounts(ctx.counts, alphabetSize_, ctx.total);
    }

    int alphabetSize_;
    std::array<uint32_t, MaxAlphabetSize> order0Counts_;
    uint32_t order0Total_;
    std::array<std::array<uint32_t, MaxAlphabetSize>, MaxAlphabetSize> order1Counts_;
    std::array<uint32_t, MaxAlphabetSize> order1Totals_;
    std::unordered_map<uint32_t, SparseCtx> order2Contexts_;
    std::unordered_map<uint64_t, SparseCtx> order3Contexts_;
};

/**
 * Variable-size rANS encoder that handles dynamic alphabet sizes
 */
class DynamicRANSEncoder {
public:
    static std::vector<uint8_t> encode(
        const std::vector<int>& symbols,
        const std::vector<std::vector<uint32_t>>& cdfs,
        int /*alphabetSize*/) {

        if (symbols.empty()) return {};

        uint32_t state = RANS_L;
        std::vector<uint8_t> output;

        // Encode in reverse order
        for (int i = static_cast<int>(symbols.size()) - 1; i >= 0; i--) {
            int sym = symbols[i];
            const auto& cdf = cdfs[i];

            uint32_t start = cdf[sym];
            uint32_t freq = cdf[sym + 1] - start;

            uint32_t max_state = ((RANS_L / PROB_SCALE) * freq) << 8;
            while (state >= max_state) {
                output.push_back(state & 0xFF);
                state >>= 8;
            }

            state = ((state / freq) << PROB_BITS) + (state % freq) + start;
        }

        // Append final state (LSB first; reverse will put MSB first = big-endian)
        output.push_back((state >> 0) & 0xFF);
        output.push_back((state >> 8) & 0xFF);
        output.push_back((state >> 16) & 0xFF);
        output.push_back((state >> 24) & 0xFF);

        std::reverse(output.begin(), output.end());
        return output;
    }
};

/**
 * Variable-size rANS decoder
 */
class DynamicRANSDecoder {
public:

    template<typename CDFGetter, typename SymbolHandler>
    static void decode(
        const uint8_t* data, size_t len,
        size_t count, int alphabetSize,
        CDFGetter getCDF, SymbolHandler handleSymbol) {

        if (count == 0 || len < 4) return;

        const uint8_t* ptr = data;
        const uint8_t* end = data + len;

        // Read initial state (big-endian, MSB first)
        uint32_t state = (static_cast<uint32_t>(ptr[0]) << 24) |
                         (static_cast<uint32_t>(ptr[1]) << 16) |
                         (static_cast<uint32_t>(ptr[2]) << 8) |
                         static_cast<uint32_t>(ptr[3]);
        ptr += 4;

        for (size_t i = 0; i < count; i++) {
            auto cdf = getCDF(i);

            // Extract symbol from state
            uint32_t slot = state & (PROB_SCALE - 1);

            int sym = searchCDF(cdf, alphabetSize, slot);

            handleSymbol(i, sym);

            // Decode: state = freq * (state >> PROB_BITS) + slot - start
            uint32_t start = cdf[sym];
            uint32_t freq = cdf[sym + 1] - start;
            state = freq * (state >> PROB_BITS) + slot - start;

            // Renormalize: read bytes while state is too small
            while (state < RANS_L && ptr < end) {
                state = (state << 8) | *ptr++;
            }
        }
    }
};

// ============================================================================
// RLE foundation models (used by V5 quality codec)
// ============================================================================

// Compact RLE constants: 32-symbol alphabet with escape for long runs
constexpr int RLE_V5_DIRECT_MAX = 30;   // Direct extra run 0-30 (runLen 2-32)
constexpr int RLE_V5_ESCAPE = 31;       // Escape: add 31 more, emit another symbol
constexpr int RLE_V5_ALPHABET = 32;     // 0-30 direct, 31 escape

// Binary model: 0 = no run (length 1), 1 = has run (length >= 2)
class AdaptiveRunFlagModel {
public:
    AdaptiveRunFlagModel() { reset(); }
    void reset() { counts_[0] = counts_[1] = 1; total_ = 2; }

    std::vector<uint32_t> getCDF() const {
        return normalizeCDFVec(counts_, 2, PROB_SCALE);
    }

    void update(int flag) {
        counts_[flag]++;
        total_++;
        if (total_ > COUNT_RESCALE_THRESHOLD) {
            rescaleCounts(counts_, 2, total_);
        }
    }

private:
    uint32_t counts_[2];
    uint32_t total_;
};

// 32-symbol model for run extensions (only used when run-flag = 1)
class AdaptiveRunExtModel {
public:
    AdaptiveRunExtModel() { reset(); }
    void reset() {
        for (int i = 0; i < RLE_V5_ALPHABET; i++) counts_[i] = 1;
        total_ = RLE_V5_ALPHABET;
    }

    std::vector<uint32_t> getCDF() const {
        return normalizeCDFVec(counts_.data(), RLE_V5_ALPHABET, PROB_SCALE);
    }

    void update(int sym) {
        counts_[sym]++;
        total_++;
        if (total_ > COUNT_RESCALE_THRESHOLD) {
            rescaleCounts(counts_, RLE_V5_ALPHABET, total_);
        }
    }

private:
    std::array<uint32_t, RLE_V5_ALPHABET> counts_;
    uint32_t total_;
};

// Reconstruct quality strings from RLE-decoded quality symbols and run lengths
inline std::vector<std::string> reconstructQualities(
    const std::vector<int>& qualitySymbols,
    const std::vector<int>& runLengths,
    const std::vector<uint8_t>& fromCompact,
    int alphabetSize,
    const std::vector<size_t>& lengths) {

    std::vector<std::string> results;
    results.reserve(lengths.size());
    size_t runIdx = 0;

    for (size_t readIdx = 0; readIdx < lengths.size(); readIdx++) {
        std::string qualities;
        qualities.reserve(lengths[readIdx]);

        size_t remaining = lengths[readIdx];
        while (remaining > 0 && runIdx < qualitySymbols.size()) {
            int compactIdx = qualitySymbols[runIdx];
            size_t runLen = std::min(static_cast<size_t>(runLengths[runIdx]), remaining);
            int originalQ = (compactIdx >= 0 && compactIdx < alphabetSize)
                ? fromCompact[compactIdx] : 0;
            char qChar = static_cast<char>(originalQ + QUALITY_OFFSET);
            for (size_t j = 0; j < runLen; j++) {
                qualities.push_back(qChar);
            }
            remaining -= runLen;
            runIdx++;
        }

        results.push_back(std::move(qualities));
    }

    return results;
}

// ============================================================================
// V5: Position-bucketed PPM with order-3, per-quality RLE, histogram warm-start
// ============================================================================

constexpr int NUM_POS_BINS_V5 = 8;
constexpr int NUM_QUAL_BINS_V6 = 32;  // 8 position bins × 4 DNA bases

inline int getPosBinV5(size_t pos, size_t readLen) {
    if (readLen == 0) return 0;
    return static_cast<int>(std::min(static_cast<size_t>(NUM_POS_BINS_V5 - 1),
                                      pos * NUM_POS_BINS_V5 / readLen));
}

// V6: combined position + DNA base context bin
inline int getQualBinV6(size_t pos, size_t readLen, int dnaBase) {
    int posBin = getPosBinV5(pos, readLen);
    return posBin * 4 + std::max(0, std::min(3, dnaBase));
}

/**
 * AdaptiveQualityPPMv2 - Wraps N independent PPM instances (one per context bin).
 * Default: 8 bins (position-bucketed, V5). With DNA base context: 32 bins (V6).
 */
template<int MaxAlphabetSize>
class AdaptiveQualityPPMv2 {
public:
    explicit AdaptiveQualityPPMv2(int numBins = NUM_POS_BINS_V5)
        : bins_(numBins), numBins_(numBins) {}

    void setAlphabetSize(int size) {
        for (auto& bin : bins_) bin.setAlphabetSize(size);
    }

    void initFromFrequencies(const std::vector<uint32_t>& freqs) {
        for (auto& bin : bins_) bin.initFromFrequencies(freqs);
    }

    std::vector<uint32_t> getCDF(int bin, int prev3, int prev2, int prev1) const {
        bin = std::max(0, std::min(numBins_ - 1, bin));
        return bins_[bin].getCDF(prev3, prev2, prev1);
    }

    void update(int bin, int prev3, int prev2, int prev1, int symbol) {
        bin = std::max(0, std::min(numBins_ - 1, bin));
        bins_[bin].update(prev3, prev2, prev1, symbol);
    }

private:
    std::vector<AdaptiveQualityPPM<MaxAlphabetSize>> bins_;
    int numBins_;
};

/**
 * QualityEncoderV5 - V5 encoder with:
 *   - Position-bucketed PPM (8 bins)
 *   - Order-3 sparse contexts
 *   - Per-quality-symbol RLE models
 *   - First-pass histogram warm-start
 *   - Cross-read context persistence (no context reset at read boundaries)
 *
 * Format byte: 0x03
 * Header: [0x03] [alphabetSize] [alphabet: N bytes] [globalFreqs: N varints]
 *         [numRuns: u32LE] [qualDataLen: u32LE] [runFlagDataLen: u32LE] [numRunExts: u32LE]
 *         [qualityData] [runFlagData] [runExtData]
 */
class QualityEncoderV5 {
public:
    QualityEncoderV5() {}

    std::vector<uint8_t> encodeAll(const std::vector<std::string>& qualityStrings) {
        if (qualityStrings.empty()) return {};

        // --- Pass 1: discover alphabet + collect global histogram ---
        std::array<bool, MAX_QUALITY_ALPHABET> used;
        used.fill(false);
        std::array<uint64_t, MAX_QUALITY_ALPHABET> rawFreqs;
        rawFreqs.fill(0);
        size_t totalLen = 0;

        for (const auto& qs : qualityStrings) {
            totalLen += qs.size();
            for (char c : qs) {
                int idx = static_cast<int>(static_cast<unsigned char>(c)) - QUALITY_OFFSET;
                if (idx >= 0 && idx < MAX_QUALITY_ALPHABET) {
                    used[idx] = true;
                    rawFreqs[idx]++;
                }
            }
        }

        if (totalLen == 0) return {};

        // Build alphabet mapping
        std::array<int8_t, MAX_QUALITY_ALPHABET> toCompact;
        std::vector<uint8_t> fromCompact;
        toCompact.fill(-1);

        for (int i = 0; i < MAX_QUALITY_ALPHABET; i++) {
            if (used[i]) {
                toCompact[i] = static_cast<int8_t>(fromCompact.size());
                fromCompact.push_back(static_cast<uint8_t>(i));
            }
        }

        int alphabetSize = static_cast<int>(fromCompact.size());

        // Build scaled histogram for warm-start (sum ~= 256)
        std::vector<uint32_t> scaledFreqs(alphabetSize);
        {
            uint64_t sumRaw = 0;
            for (int i = 0; i < alphabetSize; i++) {
                sumRaw += rawFreqs[fromCompact[i]];
            }
            if (sumRaw == 0) sumRaw = 1;
            for (int i = 0; i < alphabetSize; i++) {
                scaledFreqs[i] = static_cast<uint32_t>(
                    rawFreqs[fromCompact[i]] * 256ULL / sumRaw);
                if (scaledFreqs[i] == 0) scaledFreqs[i] = 1;
            }
        }

        // --- Pass 2: encode with all V5 models ---
        AdaptiveQualityPPMv2<MAX_QUALITY_ALPHABET> qualityModel;
        qualityModel.setAlphabetSize(alphabetSize);
        qualityModel.initFromFrequencies(scaledFreqs);

        // Per-quality-symbol run models
        std::vector<AdaptiveRunFlagModel> runFlagModels(alphabetSize);
        std::vector<AdaptiveRunExtModel> runExtModels(alphabetSize);

        // Collect 3 streams
        std::vector<int> qualitySymbols;
        std::vector<std::vector<uint32_t>> qualityCDFs;
        std::vector<int> runFlagSymbols;
        std::vector<std::vector<uint32_t>> runFlagCDFs;
        std::vector<int> runExtSymbols;
        std::vector<std::vector<uint32_t>> runExtCDFs;

        size_t estimatedRuns = totalLen / 2;
        qualitySymbols.reserve(estimatedRuns);
        qualityCDFs.reserve(estimatedRuns);
        runFlagSymbols.reserve(estimatedRuns);
        runFlagCDFs.reserve(estimatedRuns);

        // Cross-read context persistence: don't reset prev1/prev2/prev3
        int prev3 = -1, prev2 = -1, prev1 = -1;

        for (const auto& qs : qualityStrings) {
            size_t readLen = qs.size();
            size_t posInRead = 0;

            while (posInRead < qs.size()) {
                int idx = static_cast<int>(static_cast<unsigned char>(qs[posInRead])) - QUALITY_OFFSET;
                int compactIdx = (idx >= 0 && idx < MAX_QUALITY_ALPHABET) ? toCompact[idx] : 0;

                // Count run length
                size_t runLen = 1;
                while (posInRead + runLen < qs.size() && qs[posInRead + runLen] == qs[posInRead]) {
                    runLen++;
                }

                // Position bin
                int posBin = getPosBinV5(posInRead, readLen);

                // Quality symbol with position-bucketed PPM + order-3
                qualityCDFs.push_back(qualityModel.getCDF(posBin, prev3, prev2, prev1));
                qualitySymbols.push_back(compactIdx);
                qualityModel.update(posBin, prev3, prev2, prev1, compactIdx);

                // Per-quality run-flag
                int flag = (runLen > 1) ? 1 : 0;
                runFlagCDFs.push_back(runFlagModels[compactIdx].getCDF());
                runFlagSymbols.push_back(flag);
                runFlagModels[compactIdx].update(flag);

                // Run extension (per-quality model)
                if (flag == 1) {
                    int extra = static_cast<int>(runLen - 2);
                    while (extra > RLE_V5_DIRECT_MAX) {
                        runExtCDFs.push_back(runExtModels[compactIdx].getCDF());
                        runExtSymbols.push_back(RLE_V5_ESCAPE);
                        runExtModels[compactIdx].update(RLE_V5_ESCAPE);
                        extra -= (RLE_V5_DIRECT_MAX + 1);
                    }
                    runExtCDFs.push_back(runExtModels[compactIdx].getCDF());
                    runExtSymbols.push_back(extra);
                    runExtModels[compactIdx].update(extra);
                }

                prev3 = prev2;
                prev2 = prev1;
                prev1 = compactIdx;
                posInRead += runLen;
            }
            // Reset position counter at read boundary, but keep prev1/2/3
        }

        // Encode the 3 streams
        auto qualityData = DynamicRANSEncoder::encode(qualitySymbols, qualityCDFs, alphabetSize);
        auto runFlagData = DynamicRANSEncoder::encode(runFlagSymbols, runFlagCDFs, 2);
        auto runExtData = DynamicRANSEncoder::encode(runExtSymbols, runExtCDFs, RLE_V5_ALPHABET);

        // Build output
        std::vector<uint8_t> output;

        output.push_back(0x03);  // V5 format byte
        output.push_back(static_cast<uint8_t>(alphabetSize));
        for (uint8_t q : fromCompact) {
            output.push_back(q);
        }

        // Global frequency histogram (N varints)
        for (int i = 0; i < alphabetSize; i++) {
            writeVarint(output, scaledFreqs[i]);
        }

        // Stream metadata (4x uint32_t LE)
        auto writeLE32 = [&](uint32_t v) {
            output.push_back(v & 0xFF);
            output.push_back((v >> 8) & 0xFF);
            output.push_back((v >> 16) & 0xFF);
            output.push_back((v >> 24) & 0xFF);
        };

        writeLE32(static_cast<uint32_t>(qualitySymbols.size()));  // numRuns
        writeLE32(static_cast<uint32_t>(qualityData.size()));     // qualityDataLen
        writeLE32(static_cast<uint32_t>(runFlagData.size()));     // runFlagDataLen
        writeLE32(static_cast<uint32_t>(runExtSymbols.size()));   // numRunExts

        output.insert(output.end(), qualityData.begin(), qualityData.end());
        output.insert(output.end(), runFlagData.begin(), runFlagData.end());
        output.insert(output.end(), runExtData.begin(), runExtData.end());

        return output;
    }

    // V6: DNA-conditioned quality encoding. Same algorithm but uses 32 bins
    // (8 position × 4 DNA bases) instead of 8, and writes format byte 0x04.
    std::vector<uint8_t> encodeAll(const std::vector<std::string>& qualityStrings,
                                    const std::vector<std::string>& sequences) {
        if (qualityStrings.empty()) return {};

        // --- Pass 1: discover alphabet + collect global histogram (same as V5) ---
        std::array<bool, MAX_QUALITY_ALPHABET> used;
        used.fill(false);
        std::array<uint64_t, MAX_QUALITY_ALPHABET> rawFreqs;
        rawFreqs.fill(0);
        size_t totalLen = 0;

        for (const auto& qs : qualityStrings) {
            totalLen += qs.size();
            for (char c : qs) {
                int idx = static_cast<int>(static_cast<unsigned char>(c)) - QUALITY_OFFSET;
                if (idx >= 0 && idx < MAX_QUALITY_ALPHABET) {
                    used[idx] = true;
                    rawFreqs[idx]++;
                }
            }
        }

        if (totalLen == 0) return {};

        std::array<int8_t, MAX_QUALITY_ALPHABET> toCompact;
        std::vector<uint8_t> fromCompact;
        toCompact.fill(-1);
        for (int i = 0; i < MAX_QUALITY_ALPHABET; i++) {
            if (used[i]) {
                toCompact[i] = static_cast<int8_t>(fromCompact.size());
                fromCompact.push_back(static_cast<uint8_t>(i));
            }
        }
        int alphabetSize = static_cast<int>(fromCompact.size());

        std::vector<uint32_t> scaledFreqs(alphabetSize);
        {
            uint64_t sumRaw = 0;
            for (int i = 0; i < alphabetSize; i++) sumRaw += rawFreqs[fromCompact[i]];
            if (sumRaw == 0) sumRaw = 1;
            for (int i = 0; i < alphabetSize; i++) {
                scaledFreqs[i] = static_cast<uint32_t>(rawFreqs[fromCompact[i]] * 256ULL / sumRaw);
                if (scaledFreqs[i] == 0) scaledFreqs[i] = 1;
            }
        }

        // --- Pass 2: encode with 32-bin model (position × DNA base) ---
        AdaptiveQualityPPMv2<MAX_QUALITY_ALPHABET> qualityModel(NUM_QUAL_BINS_V6);
        qualityModel.setAlphabetSize(alphabetSize);
        qualityModel.initFromFrequencies(scaledFreqs);

        std::vector<AdaptiveRunFlagModel> runFlagModels(alphabetSize);
        std::vector<AdaptiveRunExtModel> runExtModels(alphabetSize);

        std::vector<int> qualitySymbols;
        std::vector<std::vector<uint32_t>> qualityCDFs;
        std::vector<int> runFlagSymbols;
        std::vector<std::vector<uint32_t>> runFlagCDFs;
        std::vector<int> runExtSymbols;
        std::vector<std::vector<uint32_t>> runExtCDFs;

        size_t estimatedRuns = totalLen / 2;
        qualitySymbols.reserve(estimatedRuns);
        qualityCDFs.reserve(estimatedRuns);
        runFlagSymbols.reserve(estimatedRuns);
        runFlagCDFs.reserve(estimatedRuns);

        int prev3 = -1, prev2 = -1, prev1 = -1;

        for (size_t readIdx = 0; readIdx < qualityStrings.size(); readIdx++) {
            const auto& qs = qualityStrings[readIdx];
            const auto& seq = (readIdx < sequences.size()) ? sequences[readIdx] : qs;
            size_t readLen = qs.size();
            size_t posInRead = 0;

            while (posInRead < qs.size()) {
                int idx = static_cast<int>(static_cast<unsigned char>(qs[posInRead])) - QUALITY_OFFSET;
                int compactIdx = (idx >= 0 && idx < MAX_QUALITY_ALPHABET) ? toCompact[idx] : 0;

                size_t runLen = 1;
                while (posInRead + runLen < qs.size() && qs[posInRead + runLen] == qs[posInRead]) {
                    runLen++;
                }

                // V6: combined position + DNA base bin
                int dnaBase = 0;
                if (posInRead < seq.size()) {
                    switch (seq[posInRead]) {
                        case 'A': case 'a': dnaBase = 0; break;
                        case 'C': case 'c': dnaBase = 1; break;
                        case 'G': case 'g': dnaBase = 2; break;
                        case 'T': case 't': dnaBase = 3; break;
                    }
                }
                int bin = getQualBinV6(posInRead, readLen, dnaBase);

                qualityCDFs.push_back(qualityModel.getCDF(bin, prev3, prev2, prev1));
                qualitySymbols.push_back(compactIdx);
                qualityModel.update(bin, prev3, prev2, prev1, compactIdx);

                int flag = (runLen > 1) ? 1 : 0;
                runFlagCDFs.push_back(runFlagModels[compactIdx].getCDF());
                runFlagSymbols.push_back(flag);
                runFlagModels[compactIdx].update(flag);

                if (flag == 1) {
                    int extra = static_cast<int>(runLen - 2);
                    while (extra > RLE_V5_DIRECT_MAX) {
                        runExtCDFs.push_back(runExtModels[compactIdx].getCDF());
                        runExtSymbols.push_back(RLE_V5_ESCAPE);
                        runExtModels[compactIdx].update(RLE_V5_ESCAPE);
                        extra -= (RLE_V5_DIRECT_MAX + 1);
                    }
                    runExtCDFs.push_back(runExtModels[compactIdx].getCDF());
                    runExtSymbols.push_back(extra);
                    runExtModels[compactIdx].update(extra);
                }

                prev3 = prev2;
                prev2 = prev1;
                prev1 = compactIdx;
                posInRead += runLen;
            }
        }

        auto qualityData = DynamicRANSEncoder::encode(qualitySymbols, qualityCDFs, alphabetSize);
        auto runFlagData = DynamicRANSEncoder::encode(runFlagSymbols, runFlagCDFs, 2);
        auto runExtData = DynamicRANSEncoder::encode(runExtSymbols, runExtCDFs, RLE_V5_ALPHABET);

        std::vector<uint8_t> output;
        output.push_back(0x04);  // V6 format byte
        output.push_back(static_cast<uint8_t>(alphabetSize));
        for (uint8_t q : fromCompact) output.push_back(q);
        for (int i = 0; i < alphabetSize; i++) writeVarint(output, scaledFreqs[i]);

        auto writeLE32 = [&](uint32_t v) {
            output.push_back(v & 0xFF);
            output.push_back((v >> 8) & 0xFF);
            output.push_back((v >> 16) & 0xFF);
            output.push_back((v >> 24) & 0xFF);
        };
        writeLE32(static_cast<uint32_t>(qualitySymbols.size()));
        writeLE32(static_cast<uint32_t>(qualityData.size()));
        writeLE32(static_cast<uint32_t>(runFlagData.size()));
        writeLE32(static_cast<uint32_t>(runExtSymbols.size()));

        output.insert(output.end(), qualityData.begin(), qualityData.end());
        output.insert(output.end(), runFlagData.begin(), runFlagData.end());
        output.insert(output.end(), runExtData.begin(), runExtData.end());

        return output;
    }

    void reset() {}
};

/**
 * QualityDecoderV5 - V5/V6 decoder (dispatches on format byte)
 */
class QualityDecoderV5 {
public:
    QualityDecoderV5() {}

    std::vector<std::string> decodeAll(const uint8_t* data, size_t len,
                                        const std::vector<size_t>& lengths) {
        if (len < 2) return emptyResults(lengths);

        // data[0] should be 0x03 (already checked by dispatcher)
        int alphabetSize = data[1];
        if (alphabetSize == 0 || static_cast<size_t>(2 + alphabetSize) > len) {
            return emptyResults(lengths);
        }

        // Read alphabet mapping
        std::vector<uint8_t> fromCompact(alphabetSize);
        for (int i = 0; i < alphabetSize; i++) {
            fromCompact[i] = data[2 + i];
        }

        // Read global frequency histogram
        const uint8_t* p = data + 2 + alphabetSize;
        const uint8_t* end = data + len;
        std::vector<uint32_t> scaledFreqs(alphabetSize);
        for (int i = 0; i < alphabetSize && p < end; i++) {
            scaledFreqs[i] = static_cast<uint32_t>(readVarint(p));
        }

        // Read stream metadata
        if (p + 16 > end) return emptyResults(lengths);
        uint32_t numRuns = readLE32(p); p += 4;
        uint32_t qualityDataLen = readLE32(p); p += 4;
        uint32_t runFlagDataLen = readLE32(p); p += 4;
        /* numRunExts not needed for lockstep decode */ p += 4;

        if (numRuns == 0) return emptyResults(lengths);

        const uint8_t* qualityData = p;
        const uint8_t* runFlagData = p + qualityDataLen;
        const uint8_t* runExtData = runFlagData + runFlagDataLen;
        size_t runExtDataLen = len - static_cast<size_t>(runExtData - data);

        // --- Decode run flags (per-quality) ---
        // We need quality symbols to index per-quality run models, but quality
        // symbols need position bins which need run lengths... circular dependency.
        // Solution: decode all 3 streams with the same model state as encoder.

        // First pass: decode quality symbols and run flags/extensions together
        // by reconstructing the exact encoder state.

        // Decode run flags first with per-quality models — but we don't know quality
        // symbols yet. We must decode quality first, then flags.
        // Actually the encoder interleaves: for each run it does quality, then flag, then ext.
        // But the streams are separate rANS streams. So we need to decode them independently
        // with the models replaying the same state.

        // Strategy: decode quality symbols first (need position bins which depend on run lengths).
        // But run lengths depend on flags. So we need to decode flags first without knowing qualities.
        // Wait — the flag models are per-quality, so we need qualities first.

        // The trick: run flags are encoded with per-quality models, so we can't decode flags
        // without knowing qualities. But qualities are a separate stream. Let's decode
        // quality stream first (we need run context for position). But position depends on run lengths...

        // Resolution: we can decode in this order:
        // 1. Decode quality symbols with a temporary pass (no position info yet — use a fallback?)
        // No, we need exact model state match.

        // Actually let's think about what the encoder does per-run:
        //   qualityCDFs.push_back(model.getCDF(posBin, prev3, prev2, prev1))
        //   runFlagCDFs.push_back(runFlagModels[compactIdx].getCDF())
        //   runExtCDFs.push_back(runExtModels[compactIdx].getCDF())
        // The run flag CDF depends on compactIdx (which quality symbol was encoded).
        // So the 3 streams are actually coupled through quality symbols.

        // Approach: decode quality stream first with dummy position bins (we don't know
        // run lengths yet). Then use those quality symbols to decode flags to get run lengths.
        // Then re-decode quality with correct positions.

        // Better approach: decode quality and flags in lockstep. But they're separate rANS
        // streams...

        // Simplest correct approach: two-pass decode.
        // Pass A: decode quality symbols with position tracking deferred (all bin 0).
        //         This is WRONG because the encoder used correct bins.
        // This won't work.

        // Correct approach: since the encoder wrote quality with correct position bins,
        // we MUST know the position bins when decoding quality. Position bins need run lengths.
        // Run lengths need flag stream. Flag stream needs quality symbols.
        //
        // But wait: the flag stream models are indexed by quality symbol, and quality symbols
        // are determined by the rANS quality stream independently of flags. The quality model
        // uses (posBin, prev3, prev2, prev1) which are derived from the quality sequence itself
        // plus position tracking. Position tracking requires knowing run lengths.
        //
        // The only circular dependency is: quality decoding needs position bins, which need
        // run lengths, which need flags, which need quality symbols.
        //
        // Break the cycle: decode flags WITHOUT per-quality indexing in a first pass, then
        // use run lengths to decode quality with correct position bins, then the flag decoding
        // was "wrong" model-wise...
        //
        // Actually the cleanest solution: Pre-decode run flags and extensions using a SINGLE
        // global model (not per-quality), then compute run lengths, then decode quality
        // with correct positions, then verify consistency.
        //
        // BUT the encoder used per-quality models for flags. So we can't decode flags with
        // a single model — they were entropy-coded with per-quality CDFs.
        //
        // The real resolution: we need to change the encoding order so that flags are coded
        // before quality, OR use a non-per-quality flag model.
        //
        // For this implementation: we'll use the flag CDF that was computed at encode time.
        // The encoder computes flag CDF AFTER knowing compactIdx. So we need compactIdx
        // to decode flags, but we need flags to decode quality with correct position.
        //
        // Practical resolution: We know that per-quality flag models only need the current
        // quality symbol. And quality symbols form a separate rANS stream that can be decoded
        // independently IF we know the position bins. Position bins require run lengths, i.e.,
        // the flag stream.
        //
        // Final approach: Use a two-pass strategy:
        // 1. First, decode quality stream WITHOUT position bins (posBin=0 for all).
        //    This gives us approximate quality symbols. Then use these to decode flags
        //    (which gives correct run lengths).
        //    BUT this is wrong because the quality model in the encoder used correct position bins.
        //    The rANS stream was coded with those CDFs, so decoding with wrong CDFs gives garbage.
        //
        // The ONLY correct solution: ensure we can decode run lengths without knowing
        // quality symbols. This means we should NOT use per-quality flag models, OR we should
        // write run lengths as a separate independent stream.
        //
        // DESIGN CHANGE: For V5, we store run lengths as a separate stream (just the count
        // of runs per-quality), with a single global flag model (not per-quality).
        // The per-quality benefit comes from the position-bucketed PPM and order-3 contexts.
        //
        // ... Actually wait. Let's re-read the encoder more carefully. The encoder's streams
        // are coded independently. The quality stream has its own rANS. The flag stream has
        // its own rANS. They each have their own CDFs collected during a single pass.
        //
        // At decode time we need to replay the same CDF sequence. For quality: we need
        // (posBin, prev3, prev2, prev1). For flags: we need compactIdx (to pick the right model).
        // For extensions: we need compactIdx.
        //
        // So the dependency is: quality decode needs posInRead (which needs run lengths from flags),
        // and flag decode needs compactIdx (from quality decode). This is circular.
        //
        // RESOLUTION: Encode run lengths as an entirely separate fourth stream with a single
        // model. Or better: encode total character count per run as a simple stream.
        //
        // Simplest fix: use a SINGLE flag model (not per-quality) and SINGLE ext model.
        // The per-quality benefit is marginal compared to position bins + order-3.
        // Let me NOT do per-quality flag models, and instead keep the single global models
        // like V4RLE. Then we can decode flags first, get run lengths, compute positions,
        // and decode quality with correct position bins. This avoids the circular dependency.
        //
        // I'll fix the encoder above to use single global models. But since the encoder code
        // is already written with per-quality models, let me instead add a SEPARATE run-count
        // stream that stores just the raw run lengths, independently decodable.
        //
        // ACTUALLY - simplest approach that works: decode flags and quality in the same
        // lockstep order as the encoder. We can do this by streaming from both rANS decoders
        // simultaneously, one symbol at a time, maintaining the same state.
        //
        // This requires a different rANS decoder that can decode one symbol at a time from
        // a pre-initialized stream. Let me implement an incremental rANS decoder.

        // INCREMENTAL DECODE approach:
        // Create two incremental rANS decoder states. For each run:
        //   1. Decode quality symbol from quality stream (needs posBin from posInRead)
        //   2. Decode run flag from flag stream (needs compactIdx from step 1)
        //   3. If flag=1, decode run extension(s) from ext stream (needs compactIdx)
        //   4. Advance posInRead by runLen
        // This exactly mirrors the encoder and breaks the circular dependency.

        // Initialize incremental rANS decoders
        IncrementalRANSDecoder qualityDec(qualityData, qualityDataLen);
        IncrementalRANSDecoder flagDec(runFlagData, runFlagDataLen);
        IncrementalRANSDecoder extDec(runExtData, runExtDataLen);

        // Initialize models (same as encoder)
        AdaptiveQualityPPMv2<MAX_QUALITY_ALPHABET> qualityModel;
        qualityModel.setAlphabetSize(alphabetSize);
        qualityModel.initFromFrequencies(scaledFreqs);

        std::vector<AdaptiveRunFlagModel> runFlagModels(alphabetSize);
        std::vector<AdaptiveRunExtModel> runExtModels(alphabetSize);

        // Decode all runs in lockstep
        std::vector<int> qualitySymbols(numRuns);
        std::vector<int> runLengths(numRuns);

        int prev3 = -1, prev2 = -1, prev1 = -1;
        size_t globalReadIdx = 0;
        size_t posInRead = 0;
        size_t currentReadLen = lengths.empty() ? 0 : lengths[0];

        for (uint32_t r = 0; r < numRuns; r++) {
            // Advance to next read if we've consumed the current one
            while (globalReadIdx < lengths.size() && posInRead >= lengths[globalReadIdx]) {
                posInRead = 0;
                globalReadIdx++;
                currentReadLen = (globalReadIdx < lengths.size()) ? lengths[globalReadIdx] : 0;
            }

            // 1. Decode quality symbol
            int posBin = getPosBinV5(posInRead, currentReadLen);
            auto qualCDF = qualityModel.getCDF(posBin, prev3, prev2, prev1);
            int compactIdx = qualityDec.decodeSymbol(qualCDF, alphabetSize);
            qualitySymbols[r] = compactIdx;
            qualityModel.update(posBin, prev3, prev2, prev1, compactIdx);

            // 2. Decode run flag
            auto flagCDF = runFlagModels[compactIdx].getCDF();
            int flag = flagDec.decodeSymbol(flagCDF, 2);
            runFlagModels[compactIdx].update(flag);

            // 3. Decode run extension if needed
            int runLen = 1;
            if (flag == 1) {
                int extra = 0;
                int extSym;
                do {
                    auto extCDF = runExtModels[compactIdx].getCDF();
                    extSym = extDec.decodeSymbol(extCDF, RLE_V5_ALPHABET);
                    runExtModels[compactIdx].update(extSym);
                    if (extSym == RLE_V5_ESCAPE) {
                        extra += RLE_V5_DIRECT_MAX + 1;
                    } else {
                        extra += extSym;
                    }
                } while (extSym == RLE_V5_ESCAPE);
                runLen = extra + 2;
            }
            runLengths[r] = runLen;

            prev3 = prev2;
            prev2 = prev1;
            prev1 = compactIdx;
            posInRead += runLen;
        }

        // Reconstruct quality strings
        return reconstructQualities(
            qualitySymbols, runLengths, fromCompact, alphabetSize, lengths);
    }

    // V6: DNA-conditioned quality decoding. Uses 32 bins (8 pos × 4 DNA bases).
    // Format byte 0x04 signals V6 (checked by caller before dispatch).
    std::vector<std::string> decodeAll(const uint8_t* data, size_t len,
                                        const std::vector<size_t>& lengths,
                                        const std::vector<std::string>& sequences) {
        if (len < 2) return emptyResults(lengths);

        // Parse header (same format as V5)
        int alphabetSize = data[1];
        if (alphabetSize == 0 || static_cast<size_t>(2 + alphabetSize) > len) {
            return emptyResults(lengths);
        }

        std::vector<uint8_t> fromCompact(alphabetSize);
        for (int i = 0; i < alphabetSize; i++) fromCompact[i] = data[2 + i];

        const uint8_t* p = data + 2 + alphabetSize;
        const uint8_t* end = data + len;
        std::vector<uint32_t> scaledFreqs(alphabetSize);
        for (int i = 0; i < alphabetSize && p < end; i++) {
            scaledFreqs[i] = static_cast<uint32_t>(readVarint(p));
        }

        if (p + 16 > end) return emptyResults(lengths);
        uint32_t numRuns = readLE32(p); p += 4;
        uint32_t qualityDataLen = readLE32(p); p += 4;
        uint32_t runFlagDataLen = readLE32(p); p += 4;
        p += 4;  // numRunExts

        if (numRuns == 0) return emptyResults(lengths);

        const uint8_t* qualityData = p;
        const uint8_t* runFlagData = p + qualityDataLen;
        const uint8_t* runExtData = runFlagData + runFlagDataLen;
        size_t runExtDataLen = len - static_cast<size_t>(runExtData - data);

        // Initialize decoders and models with 32 bins
        IncrementalRANSDecoder qualityDec(qualityData, qualityDataLen);
        IncrementalRANSDecoder flagDec(runFlagData, runFlagDataLen);
        IncrementalRANSDecoder extDec(runExtData, runExtDataLen);

        AdaptiveQualityPPMv2<MAX_QUALITY_ALPHABET> qualityModel(NUM_QUAL_BINS_V6);
        qualityModel.setAlphabetSize(alphabetSize);
        qualityModel.initFromFrequencies(scaledFreqs);

        std::vector<AdaptiveRunFlagModel> runFlagModels(alphabetSize);
        std::vector<AdaptiveRunExtModel> runExtModels(alphabetSize);

        std::vector<int> qualitySymbols(numRuns);
        std::vector<int> runLengths(numRuns);

        int prev3 = -1, prev2 = -1, prev1 = -1;
        size_t globalReadIdx = 0;
        size_t posInRead = 0;
        size_t currentReadLen = lengths.empty() ? 0 : lengths[0];

        for (uint32_t r = 0; r < numRuns; r++) {
            while (globalReadIdx < lengths.size() && posInRead >= lengths[globalReadIdx]) {
                posInRead = 0;
                globalReadIdx++;
                currentReadLen = (globalReadIdx < lengths.size()) ? lengths[globalReadIdx] : 0;
            }

            // V6: combined position + DNA base bin
            int dnaBase = 0;
            if (globalReadIdx < sequences.size() && posInRead < sequences[globalReadIdx].size()) {
                switch (sequences[globalReadIdx][posInRead]) {
                    case 'A': case 'a': dnaBase = 0; break;
                    case 'C': case 'c': dnaBase = 1; break;
                    case 'G': case 'g': dnaBase = 2; break;
                    case 'T': case 't': dnaBase = 3; break;
                }
            }
            int bin = getQualBinV6(posInRead, currentReadLen, dnaBase);

            auto qualCDF = qualityModel.getCDF(bin, prev3, prev2, prev1);
            int compactIdx = qualityDec.decodeSymbol(qualCDF, alphabetSize);
            qualitySymbols[r] = compactIdx;
            qualityModel.update(bin, prev3, prev2, prev1, compactIdx);

            auto flagCDF = runFlagModels[compactIdx].getCDF();
            int flag = flagDec.decodeSymbol(flagCDF, 2);
            runFlagModels[compactIdx].update(flag);

            int runLen = 1;
            if (flag == 1) {
                int extra = 0;
                int extSym;
                do {
                    auto extCDF = runExtModels[compactIdx].getCDF();
                    extSym = extDec.decodeSymbol(extCDF, RLE_V5_ALPHABET);
                    runExtModels[compactIdx].update(extSym);
                    if (extSym == RLE_V5_ESCAPE) {
                        extra += RLE_V5_DIRECT_MAX + 1;
                    } else {
                        extra += extSym;
                    }
                } while (extSym == RLE_V5_ESCAPE);
                runLen = extra + 2;
            }
            runLengths[r] = runLen;

            prev3 = prev2;
            prev2 = prev1;
            prev1 = compactIdx;
            posInRead += runLen;
        }

        return reconstructQualities(
            qualitySymbols, runLengths, fromCompact, alphabetSize, lengths);
    }

    void reset() {}

private:
    static std::vector<std::string> emptyResults(const std::vector<size_t>& lengths) {
        return std::vector<std::string>(lengths.size(), "");
    }

    static uint32_t readLE32(const uint8_t* p) {
        return p[0] | (static_cast<uint32_t>(p[1]) << 8) |
               (static_cast<uint32_t>(p[2]) << 16) | (static_cast<uint32_t>(p[3]) << 24);
    }

    /**
     * Incremental rANS decoder - decodes one symbol at a time from a stream.
     * Maintains state between calls, allowing interleaved decoding from multiple streams.
     */
    class IncrementalRANSDecoder {
    public:
        IncrementalRANSDecoder(const uint8_t* data, size_t len)
            : ptr_(data + 4), end_(data + len) {
            if (len >= 4) {
                state_ = (static_cast<uint32_t>(data[0]) << 24) |
                         (static_cast<uint32_t>(data[1]) << 16) |
                         (static_cast<uint32_t>(data[2]) << 8) |
                         static_cast<uint32_t>(data[3]);
            } else {
                state_ = RANS_L;
            }
        }

        int decodeSymbol(const std::vector<uint32_t>& cdf, int alphabetSize) {
            uint32_t slot = state_ & (PROB_SCALE - 1);
            int sym = searchCDF(cdf, alphabetSize, slot);

            // Decode step
            uint32_t start = cdf[sym];
            uint32_t freq = cdf[sym + 1] - start;
            state_ = freq * (state_ >> PROB_BITS) + slot - start;

            // Renormalize
            while (state_ < RANS_L && ptr_ < end_) {
                state_ = (state_ << 8) | *ptr_++;
            }

            return sym;
        }

    private:
        uint32_t state_;
        const uint8_t* ptr_;
        const uint8_t* end_;
    };
};

} // namespace v4zip

#endif // V4ZIP_QUALITYMODELV4_HPP
