/**
 * CompressorV4.hpp - DNA sequence compressor with adaptive context modeling
 *
 * Fully lossless: preserves FASTA headers and non-DNA characters.
 * Uses online learning where encoder and decoder update identically —
 * no model storage needed in the compressed stream.
 *
 * Key components:
 * - PPM* context model (up to order 16) with V₄ orbit pooling
 * - TwoLayerMixer: neural network blending of per-order predictions
 * - AdaptiveProbabilityMap: context-dependent probability refinement
 * - SubstitutionalContext: edit-tolerant contexts for repeat regions
 * - MatchModel: PAQ-style longest-match predictor
 * - V₄-symmetric enhancements: canonical match, coset features,
 *   orbit-stratified forgetting, symmetric edit classes, soft equivariant blending
 * - Multi-record FASTA with RLE-encoded ambiguous bases
 * - Interleaved rANS entropy coding
 */

#ifndef V4ZIP_COMPRESSOR_V4_HPP
#define V4ZIP_COMPRESSOR_V4_HPP

#include "Alphabet.hpp"
#include "V4Group.hpp"
#include "ArithmeticCoder.hpp"
#include "BinaryFormat.hpp"
#include "FastaIO.hpp"
#include "IdentifierCodec.hpp"
#include "Types.hpp"
#include "Varint.hpp"
#include <string>
#include <vector>
#include <array>
#include <unordered_map>
#include <cstdint>
#include <cmath>
#include <cstring>
#include <algorithm>
#include <numeric>
#include <cassert>

namespace v4zip {

constexpr uint64_t FIBONACCI_HASH_CONSTANT = 11400714819323198485ULL;  // 2^64 / phi
constexpr uint32_t MAGIC_V4 = 0x34524756;  // "VGR4"
constexpr uint8_t VERSION_V4_MAJOR = 4;
constexpr uint8_t VERSION_V4_MINOR = 9;  // V9: memory-bounded chunked FASTA encoding

// Non-DNA run record (position, length, character)
struct NonDNARun {
    uint64_t position;
    uint64_t length;
    char character;
};

// RLE-encoded ambiguous run (V5+): position, length, character
struct AmbiguousRun {
    uint64_t position;
    uint64_t length;
    char character;
};

// Convert per-position ambiguous list to RLE runs
static inline std::vector<AmbiguousRun> ambiguousToRuns(
        const std::vector<std::pair<size_t, char>>& ambiguous) {
    std::vector<AmbiguousRun> runs;
    if (ambiguous.empty()) return runs;
    AmbiguousRun cur{ambiguous[0].first, 1, ambiguous[0].second};
    for (size_t i = 1; i < ambiguous.size(); i++) {
        if (ambiguous[i].second == cur.character &&
            ambiguous[i].first == cur.position + cur.length) {
            cur.length++;
        } else {
            runs.push_back(cur);
            cur = {ambiguous[i].first, 1, ambiguous[i].second};
        }
    }
    runs.push_back(cur);
    return runs;
}

// ============================================================================
// Stretch/Squash lookup tables (replaces exp/log in hot path)
// ============================================================================

struct StretchSquashLUT {
    // stretch: [0..4095] → int16_t logit (fixed-point, scale 64)
    // Maps probability p/4096 to log(p/(1-p)) * 64
    static constexpr int TABLE_SIZE = 4096;
    static constexpr int LOGIT_SCALE = 64;  // Fixed-point scale for logits
    static constexpr int SQUASH_TABLE_SIZE = 8192;  // 2 * TABLE_SIZE for signed range

    int16_t stretchTable[TABLE_SIZE];
    uint16_t squashTable[SQUASH_TABLE_SIZE];  // Maps logit*64 ∈ [-4096..4095] → prob ∈ [0..4095]

    StretchSquashLUT() {
        // Build stretch table: prob index [0..4095] → logit
        for (int i = 0; i < TABLE_SIZE; i++) {
            double p = (i + 0.5) / TABLE_SIZE;
            p = std::max(1e-6, std::min(1.0 - 1e-6, p));
            double logit = std::log(p / (1.0 - p));
            logit = std::max(-64.0, std::min(64.0, logit));
            stretchTable[i] = static_cast<int16_t>(logit * LOGIT_SCALE);
        }
        // Build squash table: logit index [0..SQUASH_TABLE_SIZE-1] maps [-4096..4095] → prob
        for (int i = 0; i < SQUASH_TABLE_SIZE; i++) {
            double logit = (i - TABLE_SIZE) / static_cast<double>(LOGIT_SCALE);
            double p = 1.0 / (1.0 + std::exp(-logit));
            squashTable[i] = static_cast<uint16_t>(p * (TABLE_SIZE - 1) + 0.5);
        }
    }

    int16_t stretch(uint16_t prob4096) const {
        return stretchTable[std::min(static_cast<int>(prob4096), TABLE_SIZE - 1)];
    }

    uint16_t squash(int32_t logit_scaled) const {
        int idx = static_cast<int>(logit_scaled) + 4096;
        idx = std::max(0, std::min(8191, idx));
        return squashTable[idx];
    }

    static const StretchSquashLUT& instance() {
        static StretchSquashLUT lut;
        return lut;
    }
};

// ============================================================================
// Compact count type: uint16_t[5] with forgetting factor
// ============================================================================

using CompactCounts = std::array<uint16_t, 5>;

static constexpr uint16_t FORGET_THRESHOLD = 1024;

// Orbit-stratified forgetting: smaller orbits accumulate pooled counts faster,
// so they need earlier rescaling to prevent count saturation.
// orbitSz 4 (default): standard threshold. orbitSz 2: half threshold (2× pooling rate).
inline void applyForgetting(CompactCounts& c, uint8_t orbitSz = 4) {
    uint16_t threshold = (orbitSz <= 2) ? (FORGET_THRESHOLD / 2) : FORGET_THRESHOLD;
    if (c[4] > threshold) {
        for (int i = 0; i < 5; i++) {
            c[i] = std::max(static_cast<uint16_t>(1), static_cast<uint16_t>(c[i] >> 1));
        }
        // Recompute total
        c[4] = c[0] + c[1] + c[2] + c[3];
    }
}

// ============================================================================
// FlatHashMap: open-addressing with Fibonacci hashing for context counts
// ============================================================================

template<typename Value>
class FlatHashMap {
public:
    struct Entry {
        uint64_t key;
        Value value;
        uint8_t occupied;  // 0 = empty, 1 = occupied
    };

    FlatHashMap() : size_(0), shift_(64), mask_(0), maxCapacity_(0) {}

    explicit FlatHashMap(size_t initialCapacity) : size_(0), maxCapacity_(0) {
        size_t cap = 16;
        while (cap < initialCapacity) cap <<= 1;
        table_.resize(cap);
        mask_ = cap - 1;
        shift_ = 64 - __builtin_ctzll(cap);
        for (auto& e : table_) e.occupied = 0;
    }

    // Set maximum entry count. When exceeded, low-count entries are evicted.
    // Value type must have operator[](4) returning the total count for eviction to work.
    void setMaxCapacity(size_t maxCap) { maxCapacity_ = maxCap; }

    Value* find(uint64_t key) {
        if (table_.empty()) return nullptr;
        size_t idx = fibHash(key);
        for (size_t probe = 0; probe <= mask_; probe++) {
            size_t i = (idx + probe) & mask_;
            if (!table_[i].occupied) return nullptr;
            if (table_[i].key == key) return &table_[i].value;
        }
        return nullptr;
    }

    const Value* find(uint64_t key) const {
        if (table_.empty()) return nullptr;
        size_t idx = fibHash(key);
        for (size_t probe = 0; probe <= mask_; probe++) {
            size_t i = (idx + probe) & mask_;
            if (!table_[i].occupied) return nullptr;
            if (table_[i].key == key) return &table_[i].value;
        }
        return nullptr;
    }

    Value& operator[](uint64_t key) {
        if (table_.empty()) grow();
        // Check load factor — grow or evict
        if (size_ * 10 > mask_ * 7) {
            if (maxCapacity_ > 0 && size_ >= maxCapacity_) {
                evictLowCounts();
            } else {
                grow();
            }
        }

        size_t idx = fibHash(key);
        for (size_t probe = 0; ; probe++) {
            size_t i = (idx + probe) & mask_;
            if (!table_[i].occupied) {
                table_[i].key = key;
                table_[i].occupied = 1;
                table_[i].value = Value{};
                size_++;
                return table_[i].value;
            }
            if (table_[i].key == key) {
                return table_[i].value;
            }
        }
    }

    // Prefetch the cache line for a given key (non-blocking)
    void prefetch(uint64_t key) const {
        if (!table_.empty()) {
            size_t idx = fibHash(key);
            __builtin_prefetch(&table_[idx & mask_], 0, 1);
        }
    }

    void clear() {
        for (auto& e : table_) e.occupied = 0;
        size_ = 0;
    }

    size_t size() const { return size_; }

    void reserve(size_t n) {
        size_t needed = n * 10 / 7 + 16;  // Account for 70% load factor
        size_t cap = 16;
        while (cap < needed) cap <<= 1;
        if (cap > table_.size()) {
            auto old = std::move(table_);
            table_.resize(cap);
            mask_ = cap - 1;
            shift_ = 64 - __builtin_ctzll(cap);
            for (auto& e : table_) e.occupied = 0;
            size_ = 0;
            for (auto& e : old) {
                if (e.occupied) {
                    (*this)[e.key] = e.value;
                }
            }
        }
    }

private:
    size_t fibHash(uint64_t key) const {
        // Fibonacci hashing: multiply by golden ratio constant, take high bits
        return static_cast<size_t>((key * FIBONACCI_HASH_CONSTANT) >> shift_);
    }

    void grow() {
        size_t newCap = table_.empty() ? 256 : table_.size() * 2;
        // Respect max capacity: don't grow table beyond what maxCapacity needs
        if (maxCapacity_ > 0) {
            size_t maxTableCap = maxCapacity_ * 10 / 7 + 16;
            size_t maxPow2 = 16;
            while (maxPow2 < maxTableCap) maxPow2 <<= 1;
            newCap = std::min(newCap, maxPow2);
        }
        rebuild(newCap, 0);  // threshold=0 means keep all
    }

    // Evict entries with total count (value[4]) at or below threshold.
    // Rebuilds the table to eliminate tombstone issues with open addressing.
    void evictLowCounts() {
        // Start with threshold that should remove ~25-50% of entries
        uint16_t threshold = 8;  // Above Laplace prior (4)
        rebuild(table_.size(), threshold);
        // If still too full, increase threshold
        while (size_ * 10 > mask_ * 7 && threshold < 256) {
            threshold *= 2;
            rebuild(table_.size(), threshold);
        }
    }

    void rebuild(size_t newCap, uint16_t minTotalCount) {
        auto old = std::move(table_);
        table_.resize(newCap);
        mask_ = newCap - 1;
        shift_ = 64 - __builtin_ctzll(newCap);
        for (auto& e : table_) e.occupied = 0;
        size_ = 0;
        for (auto& e : old) {
            if (e.occupied && e.value[4] > minTotalCount) {
                (*this)[e.key] = e.value;
            }
        }
    }

    std::vector<Entry> table_;
    size_t size_;
    int shift_;
    size_t mask_;
    size_t maxCapacity_;  // 0 = unbounded
};

// ============================================================================
// Two-layer neural mixer with LUT, per-weight adaptive LR
// ============================================================================

class TwoLayerMixer {
public:
    static constexpr int HIDDEN_SIZE = 12;
    static constexpr double BASE_LR = 0.03;

    explicit TwoLayerMixer(int numInputs = 0)
        : numInputs_(numInputs) {
        if (numInputs > 0) {
            // Layer 1: numInputs*4 → HIDDEN_SIZE (shared across bases)
            w1_.assign(HIDDEN_SIZE, std::vector<double>(numInputs, 0.0));
            b1_.assign(HIDDEN_SIZE, 0.0);
            // Initialize with small random-ish weights
            for (int h = 0; h < HIDDEN_SIZE; h++) {
                for (int i = 0; i < numInputs; i++) {
                    // Deterministic pseudo-random init based on position
                    w1_[h][i] = 0.1 * std::sin(h * 13.0 + i * 7.0 + 1.0);
                }
                b1_[h] = 0.01 * (h - HIDDEN_SIZE / 2);
            }
            // Layer 2: HIDDEN_SIZE → 4 (per base output)
            w2_.assign(4, std::vector<double>(HIDDEN_SIZE, 0.0));
            b2_.assign(4, 0.0);
            for (int s = 0; s < 4; s++) {
                for (int h = 0; h < HIDDEN_SIZE; h++) {
                    w2_[s][h] = 0.1 * std::sin(s * 17.0 + h * 11.0 + 3.0);
                }
            }
            // Per-weight squared gradient accumulators (AdaGrad)
            g1_.assign(HIDDEN_SIZE, std::vector<double>(numInputs, 1.0));
            gb1_.assign(HIDDEN_SIZE, 1.0);
            g2_.assign(4, std::vector<double>(HIDDEN_SIZE, 1.0));
            gb2_.assign(4, 1.0);
            // Auxiliary head: direct hidden → 4-class prediction (GLN-style local loss)
            wAux_.assign(4, std::vector<double>(HIDDEN_SIZE, 0.0));
            bAux_.assign(4, 0.0);
            gAux_.assign(4, std::vector<double>(HIDDEN_SIZE, 1.0));
            gbAux_.assign(4, 1.0);
            for (int s = 0; s < 4; s++) {
                for (int h = 0; h < HIDDEN_SIZE; h++) {
                    wAux_[s][h] = 0.1 * std::sin(s * 19.0 + h * 13.0 + 5.0);
                }
            }
        }
    }

    // Accept both vector and array via pointer+size
    std::array<double, 4> predict(
            const std::array<double, 4>* modelProbs, int numProbs) const {
        const auto& lut = StretchSquashLUT::instance();
        const int n = std::min(numInputs_, numProbs);

        // Cache all per-model per-base stretch values (eliminates redundant computation in skip)
        stretchCache_.resize(n * 4);
        for (int k = 0; k < n; k++) {
            for (int s = 0; s < 4; s++) {
                uint16_t p = static_cast<uint16_t>(modelProbs[k][s] * 4095.0 + 0.5);
                stretchCache_[k * 4 + s] = lut.stretch(p);
            }
        }

        // Compute mean stretched input per model from cache
        stretched_.resize(numInputs_);
        for (int k = 0; k < n; k++) {
            int32_t sum = stretchCache_[k*4] + stretchCache_[k*4+1]
                        + stretchCache_[k*4+2] + stretchCache_[k*4+3];
            stretched_[k] = sum / (4.0 * StretchSquashLUT::LOGIT_SCALE);
        }

        // Layer 1: hidden = sigmoid(W1 * stretched + b1)
        // Uses squash LUT instead of exp() for sigmoid
        hidden_.resize(HIDDEN_SIZE);
        for (int h = 0; h < HIDDEN_SIZE; h++) {
            double z = b1_[h];
            for (int k = 0; k < numInputs_; k++) {
                z += w1_[h][k] * stretched_[k];
            }
            // Squash LUT: maps logit*LOGIT_SCALE → probability in [0,4095]
            int32_t z_scaled = static_cast<int32_t>(z * StretchSquashLUT::LOGIT_SCALE);
            hidden_[h] = lut.squash(z_scaled) * (1.0 / 4095.0);
        }

        // Layer 2: logits = W2 * hidden + b2 + skip connection
        std::array<double, 4> logits = {0.0, 0.0, 0.0, 0.0};
        for (int s = 0; s < 4; s++) {
            double z = b2_[s];
            for (int h = 0; h < HIDDEN_SIZE; h++) {
                z += w2_[s][h] * hidden_[h];
            }
            // Direct skip connection from cached stretch values
            for (int k = 0; k < n; k++) {
                z += stretchCache_[k * 4 + s] * (1.0 / StretchSquashLUT::LOGIT_SCALE);
            }
            logits[s] = z;
        }

        // Softmax via squash LUT: for each base, compute P(base) from logit differences
        // Exact softmax for 4 classes using pairwise squash
        double maxLogit = *std::max_element(logits.begin(), logits.end());
        std::array<double, 4> probs;
        double sum = 0.0;
        for (int s = 0; s < 4; s++) {
            double diff = logits[s] - maxLogit;
            // Use squash LUT for exp(diff)/(1+exp(diff)), then recover exp(diff)
            // For softmax we need exp(logit_s - max). squash gives sigmoid.
            // exp(x) = sigmoid(x) / (1 - sigmoid(x)) for x <= 0
            int32_t d_scaled = static_cast<int32_t>(diff * StretchSquashLUT::LOGIT_SCALE);
            uint16_t sig = lut.squash(d_scaled);
            // Map to [0,1] and compute exp(diff) = sig/(4095-sig), clamped
            double sigD = std::max(1.0, static_cast<double>(sig));
            double denomD = std::max(1.0, 4095.0 - sigD);
            probs[s] = sigD / denomD;
            sum += probs[s];
        }
        for (int s = 0; s < 4; s++) {
            probs[s] /= sum;
            probs[s] = std::max(1e-6, std::min(1.0 - 1e-6, probs[s]));
        }
        return probs;
    }

    void update(const std::array<double, 4>* /*modelProbs*/,
                const std::array<double, 4>& currentProbs, int target) {
        if (numInputs_ == 0) return;

        // Compute output error (softmax cross-entropy gradient)
        std::array<double, 4> dOut;
        for (int s = 0; s < 4; s++) {
            dOut[s] = currentProbs[s] - (s == target ? 1.0 : 0.0);
        }

        // Backprop through layer 2: dW2, db2, dHidden
        std::array<double, HIDDEN_SIZE> dHidden = {};
        for (int s = 0; s < 4; s++) {
            for (int h = 0; h < HIDDEN_SIZE; h++) {
                double grad = dOut[s] * hidden_[h];
                g2_[s][h] += grad * grad;
                double lr = BASE_LR * fastInvSqrt(g2_[s][h] + 1e-8);
                w2_[s][h] -= lr * grad;
                dHidden[h] += dOut[s] * w2_[s][h];
            }
            double grad_b = dOut[s];
            gb2_[s] += grad_b * grad_b;
            double lr_b = BASE_LR * fastInvSqrt(gb2_[s] + 1e-8);
            b2_[s] -= lr_b * grad_b;
        }

        // GLN-style auxiliary loss: direct hidden → output prediction
        // Gives Layer 1 stronger gradient signal without backprop attenuation
        std::array<double, 4> auxLogits = {0.0, 0.0, 0.0, 0.0};
        for (int s = 0; s < 4; s++) {
            double z = bAux_[s];
            for (int h = 0; h < HIDDEN_SIZE; h++) {
                z += wAux_[s][h] * hidden_[h];
            }
            auxLogits[s] = z;
        }
        // Softmax for auxiliary prediction
        double maxAux = *std::max_element(auxLogits.begin(), auxLogits.end());
        std::array<double, 4> auxProbs;
        double auxSum = 0.0;
        for (int s = 0; s < 4; s++) {
            auxProbs[s] = std::exp(std::min(20.0, auxLogits[s] - maxAux));
            auxSum += auxProbs[s];
        }
        for (int s = 0; s < 4; s++) auxProbs[s] /= auxSum;

        // Auxiliary gradient and weight update
        std::array<double, HIDDEN_SIZE> dHiddenAux = {};
        for (int s = 0; s < 4; s++) {
            double dAux = auxProbs[s] - (s == target ? 1.0 : 0.0);
            for (int h = 0; h < HIDDEN_SIZE; h++) {
                double grad = dAux * hidden_[h];
                gAux_[s][h] += grad * grad;
                double lr = BASE_LR * fastInvSqrt(gAux_[s][h] + 1e-8);
                wAux_[s][h] -= lr * grad;
                dHiddenAux[h] += dAux * wAux_[s][h];
            }
            double grad_b = dAux;
            gbAux_[s] += grad_b * grad_b;
            double lr_b = BASE_LR * fastInvSqrt(gbAux_[s] + 1e-8);
            bAux_[s] -= lr_b * grad_b;
        }

        // Combine backprop + auxiliary gradients for Layer 1
        // alpha=0.5: equal weight to backprop and direct local signal
        constexpr double ALPHA = 0.5;
        for (int h = 0; h < HIDDEN_SIZE; h++) {
            double bp = dHidden[h] * hidden_[h] * (1.0 - hidden_[h]);  // backprop through sigmoid
            double aux = dHiddenAux[h] * hidden_[h] * (1.0 - hidden_[h]);  // aux through sigmoid
            dHidden[h] = ALPHA * bp + (1.0 - ALPHA) * aux;
        }

        // Update layer 1 weights with combined gradient
        for (int h = 0; h < HIDDEN_SIZE; h++) {
            for (int k = 0; k < numInputs_; k++) {
                double grad = dHidden[h] * stretched_[k];
                g1_[h][k] += grad * grad;
                double lr = BASE_LR * fastInvSqrt(g1_[h][k] + 1e-8);
                w1_[h][k] -= lr * grad;
            }
            double grad_b = dHidden[h];
            gb1_[h] += grad_b * grad_b;
            double lr_b = BASE_LR * fastInvSqrt(gb1_[h] + 1e-8);
            b1_[h] -= lr_b * grad_b;
        }
    }

private:
    // Fast reciprocal square root (1 Newton-Raphson iteration, ~3-5x faster than 1/sqrt)
    static double fastInvSqrt(double x) {
        double half = 0.5 * x;
        uint64_t bits;
        std::memcpy(&bits, &x, sizeof(bits));
        static constexpr uint64_t FAST_INV_SQRT_MAGIC = 0x5FE6EB50C7B537A9ULL;
        bits = FAST_INV_SQRT_MAGIC - (bits >> 1);
        double y;
        std::memcpy(&y, &bits, sizeof(y));
        y *= 1.5 - half * y * y;  // Newton-Raphson
        return y;
    }

    int numInputs_;
    // Layer 1
    std::vector<std::vector<double>> w1_;
    std::vector<double> b1_;
    std::vector<std::vector<double>> g1_;  // AdaGrad accumulators
    std::vector<double> gb1_;
    // Layer 2
    std::vector<std::vector<double>> w2_;
    std::vector<double> b2_;
    std::vector<std::vector<double>> g2_;
    std::vector<double> gb2_;
    // Auxiliary head (GLN-style local loss for Layer 1)
    std::vector<std::vector<double>> wAux_;
    std::vector<double> bAux_;
    std::vector<std::vector<double>> gAux_;
    std::vector<double> gbAux_;
    // Scratch (mutable for const predict)
    mutable std::vector<double> stretched_;
    mutable std::vector<double> hidden_;
    mutable std::vector<int16_t> stretchCache_;  // Cached per-model per-base stretch values
};

// ============================================================================
// Adaptive Probability Map (APM) — PAQ-style post-mixer refinement
// ============================================================================

class AdaptiveProbabilityMap {
public:
    static constexpr int NUM_BUCKETS = 32;

    explicit AdaptiveProbabilityMap(int numCtx = 16) : numCtx_(numCtx) {
        table_.resize(numCtx_ * 4 * NUM_BUCKETS);
        counts_.resize(numCtx_ * 4 * NUM_BUCKETS, 0);
        for (int c = 0; c < numCtx_; c++) {
            for (int b = 0; b < 4; b++) {
                for (int bucket = 0; bucket < NUM_BUCKETS; bucket++) {
                    double p = (bucket + 0.5) / NUM_BUCKETS;
                    table_[idx(c, b, bucket)] = static_cast<uint16_t>(p * UINT16_MAX);
                }
            }
        }
    }

    std::array<double, 4> refine(const std::array<double, 4>& probs, int ctx) const {
        int c = ctx % numCtx_;
        std::array<double, 4> refined;
        double sum = 0.0;
        for (int b = 0; b < 4; b++) {
            int bucket = static_cast<int>(probs[b] * (NUM_BUCKETS - 1) + 0.5);
            bucket = std::max(0, std::min(NUM_BUCKETS - 1, bucket));
            refined[b] = table_[idx(c, b, bucket)] / static_cast<double>(UINT16_MAX);
            refined[b] = std::max(1e-6, refined[b]);
            sum += refined[b];
        }
        for (int b = 0; b < 4; b++) {
            refined[b] /= sum;
        }
        return refined;
    }

    void update(const std::array<double, 4>& probs, int ctx, int actual) {
        int c = ctx % numCtx_;
        for (int b = 0; b < 4; b++) {
            int bucket = static_cast<int>(probs[b] * (NUM_BUCKETS - 1) + 0.5);
            bucket = std::max(0, std::min(NUM_BUCKETS - 1, bucket));
            uint16_t target = (b == actual) ? UINT16_MAX : 0;
            int count = std::min(static_cast<int>(counts_[idx(c, b, bucket)]), 255);
            int rate = count + 4;
            int diff = static_cast<int>(target) - static_cast<int>(table_[idx(c, b, bucket)]);
            table_[idx(c, b, bucket)] += static_cast<int16_t>(
                (diff + (diff > 0 ? rate / 2 : -rate / 2)) / rate);
            if (counts_[idx(c, b, bucket)] < 255) counts_[idx(c, b, bucket)]++;
        }
    }

private:
    int idx(int c, int b, int bucket) const {
        return (c * 4 + b) * NUM_BUCKETS + bucket;
    }
    int numCtx_;
    mutable std::vector<uint16_t> table_;
    std::vector<uint8_t> counts_;
};

// ============================================================================
// Substitutional context: tracks edited context for mismatch tolerance
// ============================================================================

struct SubstitutionalContext {
    uint64_t editedCtx = 0;
    uint64_t editMask = 0;
    uint64_t editClassMask = 0;  // V₄ edit class: 1=transversion, 0=complement/no-edit
    int order;
    int threshold;
    uint64_t ctxMask;
    bool v4EditClass = false;  // Enable V₄ edit classification (modelVersion >= 5)

    SubstitutionalContext() : order(0), threshold(0), ctxMask(0) {}

    SubstitutionalContext(int ord, int thresh)
        : editedCtx(0), editMask(0), editClassMask(0), order(ord), threshold(thresh),
          ctxMask((1ULL << (2 * ord)) - 1) {}

    void reset() {
        editedCtx = 0;
        editMask = 0;
        editClassMask = 0;
    }

    void update(int predictedBase, int actualBase) {
        bool isEdit = (predictedBase != actualBase);
        int baseToUse = isEdit ? predictedBase : actualBase;
        editedCtx = ((editedCtx << 2) | baseToUse) & ctxMask;
        editMask = (editMask << 1) | (isEdit ? 1ULL : 0ULL);
        // V₄ edit classification: complement edits (A↔T, C↔G) are within
        // the same V₄ orbit pair; transversions (A↔C, A↔G, etc.) cross orbits.
        // XOR==3 identifies complement pairs: A(0)^T(3)=3, C(1)^G(2)=3.
        bool isTransversion = isEdit && ((predictedBase ^ actualBase) != 3);
        editClassMask = (editClassMask << 1) | (isTransversion ? 1ULL : 0ULL);
    }

    bool isActive() const {
        uint64_t window = (1ULL << order) - 1;
        if (v4EditClass) {
            // V₄-aware: only count transversion edits (complement edits are V₄-trivial)
            int transversions = __builtin_popcountll(editClassMask & editMask & window);
            return transversions <= threshold;
        }
        int edits = __builtin_popcountll(editMask & window);
        return edits <= threshold;
    }
};

// ============================================================================
// Match model: PAQ-style longest-match prediction for repetitive DNA
// ============================================================================

class MatchModel {
public:
    static constexpr int TABLE_BITS = 20;  // 1M entries
    static constexpr int TABLE_SIZE = 1 << TABLE_BITS;
    static constexpr uint64_t TABLE_MASK = TABLE_SIZE - 1;
    static constexpr int K_MATCHES = 3;    // Track top-K matches per bucket
    // Max history: 128M bases = 32 MB packed. Prevents unbounded growth for large files.
    // Well within int32_t range (2^27 < 2^31) so stored positions never overflow.
    static constexpr size_t MAX_HISTORY_BASES = 128 * 1024 * 1024;

    // matchOrder > 0 enables V₄ canonicalization of context before hashing.
    // multiMatch enables K-entry buckets for multiple simultaneous matches.
    explicit MatchModel(int matchOrder = 0, bool multiMatch = false)
        : matchOrder_(matchOrder), multiMatch_(multiMatch) {
        if (multiMatch_) {
            multitable_.assign(TABLE_SIZE, MatchBucket{});
        } else {
            table_.assign(TABLE_SIZE, MatchEntry{-1, false});
        }
        history_.reserve(1 << 16);  // Packed: 4 bases/byte, so 64K bytes = 256K bases
        for (int m = 0; m < K_MATCHES; m++) {
            matchPos_[m] = -1;
            matchLen_[m] = 0;
            matchComplemented_[m] = false;
        }
    }

    // Get prediction for match slot m (0 = primary, 1/2 = secondary matches)
    std::array<double, 4> predict(uint64_t /*context*/, int slot = 0) const {
        static const std::array<double, 4> UNIFORM = {0.25, 0.25, 0.25, 0.25};

        int m = std::min(slot, K_MATCHES - 1);
        if (matchPos_[m] >= 0 && matchPos_[m] < static_cast<int>(historyLen_) && matchLen_[m] >= 2) {
            int predicted = getBase(matchPos_[m]);
            if (matchComplemented_[m]) predicted = 3 - predicted;
            double conf = std::min(0.90, 0.25 + matchLen_[m] * 0.05);
            std::array<double, 4> probs;
            double rest = (1.0 - conf) / 3.0;
            for (int b = 0; b < 4; b++) {
                probs[b] = (b == predicted) ? conf : rest;
            }
            return probs;
        }

        return UNIFORM;
    }

    // Number of active match slots (matchLen >= 2)
    int numActiveMatches() const {
        int count = 0;
        int limit = multiMatch_ ? K_MATCHES : 1;
        for (int m = 0; m < limit; m++) {
            if (matchPos_[m] >= 0 && matchPos_[m] < static_cast<int>(historyLen_) && matchLen_[m] >= 2) {
                count++;
            }
        }
        return count;
    }

    int matchLen(int slot = 0) const {
        return matchLen_[std::min(slot, K_MATCHES - 1)];
    }

    bool isMultiMatch() const { return multiMatch_; }

    // Update match state with actual symbol
    void update(uint64_t context, int base) {
        uint64_t canonCtx = context;
        bool currentCompl = false;
        if (matchOrder_ > 0) {
            uint64_t mask = (1ULL << (2 * matchOrder_)) - 1;
            auto ctxT = canonicalWithTransform(context & mask, matchOrder_);
            canonCtx = ctxT.first;
            currentCompl = (ctxT.second == Transform::C || ctxT.second == Transform::CR);
        }

        if (multiMatch_) {
            updateMulti(canonCtx, currentCompl, base);
        } else {
            updateSingle(canonCtx, currentCompl, base);
        }

        // Store current position
        uint32_t hash = hashContext(canonCtx);
        if (multiMatch_) {
            multitable_[hash].push(MatchEntry{static_cast<int32_t>(historyLen_), currentCompl});
        } else {
            table_[hash] = MatchEntry{static_cast<int32_t>(historyLen_), currentCompl};
        }
        pushBase(base);
    }

    void reset() {
        if (multiMatch_) {
            for (auto& b : multitable_) b = MatchBucket{};
        } else {
            for (auto& e : table_) { e.position = -1; e.complemented = false; }
        }
        history_.clear();
        historyLen_ = 0;
        for (int m = 0; m < K_MATCHES; m++) {
            matchPos_[m] = -1;
            matchLen_[m] = 0;
            matchComplemented_[m] = false;
        }
    }

private:
    struct MatchEntry {
        int32_t position = -1;
        bool complemented = false;
    };

    struct MatchBucket {
        MatchEntry entries[K_MATCHES];
        uint8_t count = 0;

        void push(const MatchEntry& e) {
            // Shift entries down, insert at front (most recent first)
            for (int i = K_MATCHES - 1; i > 0; i--) {
                entries[i] = entries[i - 1];
            }
            entries[0] = e;
            if (count < K_MATCHES) count++;
        }
    };

    static uint32_t hashContext(uint64_t ctx) {
        return static_cast<uint32_t>((ctx * FIBONACCI_HASH_CONSTANT) >> (64 - TABLE_BITS));
    }

    void updateSingle(uint64_t canonCtx, bool currentCompl, int base) {
        if (matchPos_[0] >= 0 && matchPos_[0] < static_cast<int>(historyLen_)) {
            int expected = getBase(matchPos_[0]);
            if (matchComplemented_[0]) expected = 3 - expected;
            if (expected == base) {
                matchLen_[0]++;
                matchPos_[0]++;
            } else {
                matchPos_[0] = -1;
                matchLen_[0] = 0;
                uint32_t hash = hashContext(canonCtx);
                const auto& entry = table_[hash];
                if (entry.position >= 0 && entry.position < static_cast<int>(historyLen_)) {
                    matchPos_[0] = entry.position;
                    matchComplemented_[0] = (currentCompl != entry.complemented);
                }
            }
        } else {
            uint32_t hash = hashContext(canonCtx);
            const auto& entry = table_[hash];
            if (entry.position >= 0 && entry.position < static_cast<int>(historyLen_)) {
                bool compl_ = (currentCompl != entry.complemented);
                int storedBase = getBase(entry.position);
                int expected = compl_ ? (3 - storedBase) : storedBase;
                if (expected == base) {
                    matchPos_[0] = entry.position + 1;
                    matchLen_[0] = 1;
                    matchComplemented_[0] = compl_;
                }
            }
        }
    }

    void updateMulti(uint64_t canonCtx, bool currentCompl, int base) {
        // Update each active match slot
        for (int m = 0; m < K_MATCHES; m++) {
            if (matchPos_[m] >= 0 && matchPos_[m] < static_cast<int>(historyLen_)) {
                int expected = getBase(matchPos_[m]);
                if (matchComplemented_[m]) expected = 3 - expected;
                if (expected == base) {
                    matchLen_[m]++;
                    matchPos_[m]++;
                } else {
                    matchPos_[m] = -1;
                    matchLen_[m] = 0;
                }
            }
        }

        // Count active matches
        int active = 0;
        for (int m = 0; m < K_MATCHES; m++) {
            if (matchPos_[m] >= 0) active++;
        }

        // If we have empty slots, try to fill them from the bucket
        if (active < K_MATCHES) {
            uint32_t hash = hashContext(canonCtx);
            const auto& bucket = multitable_[hash];
            for (int e = 0; e < bucket.count && active < K_MATCHES; e++) {
                const auto& entry = bucket.entries[e];
                if (entry.position < 0 || entry.position >= static_cast<int>(historyLen_)) continue;

                // Check this entry isn't already tracked
                bool alreadyTracked = false;
                for (int m = 0; m < K_MATCHES; m++) {
                    if (matchPos_[m] == entry.position) { alreadyTracked = true; break; }
                }
                if (alreadyTracked) continue;

                bool compl_ = (currentCompl != entry.complemented);
                int storedBase = getBase(entry.position);
                int expected = compl_ ? (3 - storedBase) : storedBase;
                if (expected == base) {
                    // Find first empty slot
                    for (int m = 0; m < K_MATCHES; m++) {
                        if (matchPos_[m] < 0) {
                            matchPos_[m] = entry.position + 1;
                            matchLen_[m] = 1;
                            matchComplemented_[m] = compl_;
                            active++;
                            break;
                        }
                    }
                }
            }
        }
    }

    // Packed 2-bit history access (4 bases per byte, 16x memory reduction)
    int getBase(size_t pos) const {
        return (history_[pos >> 2] >> (2 * (pos & 3))) & 3;
    }

    void pushBase(int base) {
        // Bound history to prevent unbounded memory growth on large files.
        // When limit is reached, clear history and invalidate active matches.
        // Hash table entries become stale (rejected by position validity checks).
        // Encoder and decoder hit this at the same position, staying synchronized.
        if (historyLen_ >= MAX_HISTORY_BASES) {
            history_.clear();
            historyLen_ = 0;
            for (int m = 0; m < K_MATCHES; m++) {
                matchPos_[m] = -1;
                matchLen_[m] = 0;
                matchComplemented_[m] = false;
            }
        }
        if ((historyLen_ & 3) == 0) history_.push_back(0);
        history_.back() |= static_cast<uint8_t>(base << (2 * (historyLen_ & 3)));
        historyLen_++;
    }

    int matchOrder_;
    bool multiMatch_;
    mutable int matchPos_[K_MATCHES] = {-1, -1, -1};
    mutable int matchLen_[K_MATCHES] = {0, 0, 0};
    mutable bool matchComplemented_[K_MATCHES] = {false, false, false};
    std::vector<MatchEntry> table_;
    std::vector<MatchBucket> multitable_;
    std::vector<uint8_t> history_;  // 2-bit packed: 4 bases per byte
    size_t historyLen_ = 0;         // Number of bases stored
};

// ============================================================================
// CopyModel: long-range V₄-canonical copy model for repetitive DNA
// ============================================================================

class CopyModel {
public:
    static constexpr int HASH_BITS = 20;           // 1M entries (was 4M; reduced for memory)
    static constexpr int HASH_SIZE = 1 << HASH_BITS;
    static constexpr uint64_t HASH_MASK = HASH_SIZE - 1;
    static constexpr int NUM_CHANNELS = 4;
    static constexpr int ANCHOR_K = 16;           // k-mer size for re-anchoring
    static constexpr int POSITIONS_PER_BUCKET = 2; // Stored positions per hash entry (was 4)
    static constexpr int WINDOW_SIZE = 32;         // Mismatch tolerance window
    static constexpr int MAX_MISMATCHES = 3;       // Per window
    // Max history: 128M bases = 32 MB packed. Same limit as MatchModel.
    static constexpr size_t MAX_HISTORY_BASES = 128 * 1024 * 1024;

    struct CopyChannel {
        int64_t sourcePos = -1;
        int matchLen = 0;
        int mismatches = 0;
        uint64_t mismatchWindow = 0;  // Sliding window bitmask of recent mismatches
        Transform v4Transform = Transform::I;
        bool active = false;
    };

    struct HashEntry {
        int32_t positions[POSITIONS_PER_BUCKET];
        bool complemented[POSITIONS_PER_BUCKET];
        uint8_t count = 0;

        HashEntry() {
            for (int i = 0; i < POSITIONS_PER_BUCKET; i++) {
                positions[i] = -1;
                complemented[i] = false;
            }
        }

        void push(int32_t pos, bool compl_) {
            // Ring buffer: shift down, insert at front
            for (int i = POSITIONS_PER_BUCKET - 1; i > 0; i--) {
                positions[i] = positions[i - 1];
                complemented[i] = complemented[i - 1];
            }
            positions[0] = pos;
            complemented[0] = compl_;
            if (count < POSITIONS_PER_BUCKET) count++;
        }
    };

    static constexpr int SHORT_ANCHOR_K = 8;       // Short anchor for diverged repeats
    static constexpr int SHORT_HASH_BITS = 18;     // 256K entries for short anchor
    static constexpr int SHORT_HASH_SIZE = 1 << SHORT_HASH_BITS;

    explicit CopyModel(int anchorOrder = ANCHOR_K)
        : anchorOrder_(anchorOrder) {
        hashTable_.assign(HASH_SIZE, HashEntry{});
        shortHashTable_.assign(SHORT_HASH_SIZE, HashEntry{});
        history_.reserve(1 << 16);
    }

    // Get prediction from a specific channel
    std::array<double, 4> predict(int channel) const {
        static const std::array<double, 4> UNIFORM = {0.25, 0.25, 0.25, 0.25};
        int ch = std::min(channel, NUM_CHANNELS - 1);
        const auto& c = channels_[ch];
        if (!c.active || c.sourcePos < 0 || c.sourcePos >= static_cast<int64_t>(historyLen_))
            return UNIFORM;

        int predicted = getBase(static_cast<size_t>(c.sourcePos));
        // Apply V₄ transform: complement if C or CR
        if (c.v4Transform == Transform::C || c.v4Transform == Transform::CR)
            predicted = 3 - predicted;

        // Confidence scales with match length, decays with mismatches
        double conf = std::min(0.90, 0.25 + 0.03 * c.matchLen);
        // Decay for mismatches: 0.85^mismatches
        for (int i = 0; i < c.mismatches; i++) conf *= 0.85;
        conf = std::max(0.25, conf);

        std::array<double, 4> probs;
        double rest = (1.0 - conf) / 3.0;
        for (int b = 0; b < 4; b++) {
            probs[b] = (b == predicted) ? conf : rest;
        }
        return probs;
    }

    int numActiveChannels() const {
        int count = 0;
        for (int ch = 0; ch < NUM_CHANNELS; ch++) {
            if (channels_[ch].active) count++;
        }
        return count;
    }

    // Update copy state with actual base, then re-anchor if needed
    void update(uint64_t fullContext, int base) {
        // Advance all active channels
        for (int ch = 0; ch < NUM_CHANNELS; ch++) {
            auto& c = channels_[ch];
            if (!c.active) continue;
            if (c.sourcePos < 0 || c.sourcePos >= static_cast<int64_t>(historyLen_)) {
                c.active = false;
                continue;
            }

            int expected = getBase(static_cast<size_t>(c.sourcePos));
            if (c.v4Transform == Transform::C || c.v4Transform == Transform::CR)
                expected = 3 - expected;

            // Slide the mismatch window
            bool isMismatch = (expected != base);
            c.mismatchWindow = (c.mismatchWindow << 1) | (isMismatch ? 1ULL : 0ULL);
            c.mismatches = __builtin_popcountll(c.mismatchWindow & ((1ULL << WINDOW_SIZE) - 1));

            if (c.mismatches > MAX_MISMATCHES) {
                c.active = false;
            } else {
                c.matchLen++;
                c.sourcePos++;
            }
        }

        // Re-anchor: try to fill inactive channels from hash table
        int activeCount = numActiveChannels();
        if (activeCount < NUM_CHANNELS && historyLen_ >= static_cast<size_t>(anchorOrder_)) {
            // Canonicalize context for hash lookup
            uint64_t mask = (1ULL << (2 * anchorOrder_)) - 1;
            uint64_t ctx = fullContext & mask;
            auto ctxT = canonicalWithTransform(ctx, anchorOrder_);
            bool currentCompl = (ctxT.second == Transform::C || ctxT.second == Transform::CR);
            uint32_t hash = hashCtx(ctxT.first);

            const auto& entry = hashTable_[hash];
            for (int e = 0; e < entry.count && activeCount < NUM_CHANNELS; e++) {
                int32_t storedPos = entry.positions[e];
                if (storedPos < 0 || storedPos >= static_cast<int32_t>(historyLen_)) continue;
                // Don't re-use very recent positions (avoid self-match)
                if (static_cast<int64_t>(historyLen_) - storedPos < anchorOrder_) continue;

                // Check if already tracked
                bool alreadyTracked = false;
                for (int ch = 0; ch < NUM_CHANNELS; ch++) {
                    if (channels_[ch].active && std::abs(channels_[ch].sourcePos - storedPos) < 4) {
                        alreadyTracked = true;
                        break;
                    }
                }
                if (alreadyTracked) continue;

                // Determine V₄ transform between source and current
                bool storedCompl = entry.complemented[e];
                bool compl_ = (currentCompl != storedCompl);
                Transform t = compl_ ? Transform::C : Transform::I;

                // Verify the base matches
                int storedBase = getBase(static_cast<size_t>(storedPos));
                int expected = compl_ ? (3 - storedBase) : storedBase;
                if (expected == base) {
                    // Find first inactive channel
                    for (int ch = 0; ch < NUM_CHANNELS; ch++) {
                        if (!channels_[ch].active) {
                            channels_[ch] = CopyChannel{
                                storedPos + 1,  // advance past matched base
                                1,              // matchLen
                                0,              // mismatches
                                0,              // mismatchWindow
                                t,
                                true
                            };
                            activeCount++;
                            break;
                        }
                    }
                }
            }
        }

        // Short-anchor re-anchoring: try order-8 hash for diverged repeats
        if (activeCount < NUM_CHANNELS && historyLen_ >= SHORT_ANCHOR_K) {
            uint64_t shortMask = (1ULL << (2 * SHORT_ANCHOR_K)) - 1;
            uint64_t shortCtx = fullContext & shortMask;
            auto shortCtxT = canonicalWithTransform(shortCtx, SHORT_ANCHOR_K);
            bool shortCompl = (shortCtxT.second == Transform::C || shortCtxT.second == Transform::CR);
            uint32_t shortHash = static_cast<uint32_t>(
                (shortCtxT.first * FIBONACCI_HASH_CONSTANT) >> (64 - SHORT_HASH_BITS));

            const auto& shortEntry = shortHashTable_[shortHash];
            for (int e = 0; e < shortEntry.count && activeCount < NUM_CHANNELS; e++) {
                int32_t storedPos = shortEntry.positions[e];
                if (storedPos < 0 || storedPos >= static_cast<int32_t>(historyLen_)) continue;
                if (static_cast<int64_t>(historyLen_) - storedPos < SHORT_ANCHOR_K) continue;

                bool alreadyTracked = false;
                for (int ch = 0; ch < NUM_CHANNELS; ch++) {
                    if (channels_[ch].active && std::abs(channels_[ch].sourcePos - storedPos) < 4) {
                        alreadyTracked = true;
                        break;
                    }
                }
                if (alreadyTracked) continue;

                bool storedCompl = shortEntry.complemented[e];
                bool compl_ = (shortCompl != storedCompl);
                Transform t = compl_ ? Transform::C : Transform::I;

                int storedBase = getBase(static_cast<size_t>(storedPos));
                int expected = compl_ ? (3 - storedBase) : storedBase;
                if (expected == base) {
                    for (int ch = 0; ch < NUM_CHANNELS; ch++) {
                        if (!channels_[ch].active) {
                            channels_[ch] = CopyChannel{
                                storedPos + 1, 1, 0, 0, t, true
                            };
                            activeCount++;
                            break;
                        }
                    }
                }
            }
        }

        // Store current position in both hash tables
        if (historyLen_ >= static_cast<size_t>(anchorOrder_)) {
            uint64_t mask = (1ULL << (2 * anchorOrder_)) - 1;
            uint64_t ctx = fullContext & mask;
            auto ctxT = canonicalWithTransform(ctx, anchorOrder_);
            bool currentCompl = (ctxT.second == Transform::C || ctxT.second == Transform::CR);
            uint32_t hash = hashCtx(ctxT.first);
            hashTable_[hash].push(static_cast<int32_t>(historyLen_), currentCompl);
        }
        if (historyLen_ >= SHORT_ANCHOR_K) {
            uint64_t shortMask = (1ULL << (2 * SHORT_ANCHOR_K)) - 1;
            uint64_t shortCtx = fullContext & shortMask;
            auto shortCtxT = canonicalWithTransform(shortCtx, SHORT_ANCHOR_K);
            bool shortCompl = (shortCtxT.second == Transform::C || shortCtxT.second == Transform::CR);
            uint32_t shortHash = static_cast<uint32_t>(
                (shortCtxT.first * FIBONACCI_HASH_CONSTANT) >> (64 - SHORT_HASH_BITS));
            shortHashTable_[shortHash].push(static_cast<int32_t>(historyLen_), shortCompl);
        }

        // Store base in packed history
        pushBase(base);
    }

    void reset() {
        for (auto& e : hashTable_) e = HashEntry{};
        for (auto& e : shortHashTable_) e = HashEntry{};
        history_.clear();
        historyLen_ = 0;
        for (int ch = 0; ch < NUM_CHANNELS; ch++) {
            channels_[ch] = CopyChannel{};
        }
    }

    // Soft reset: keep history + hash table, just deactivate channels
    void resetChannels() {
        for (int ch = 0; ch < NUM_CHANNELS; ch++) {
            channels_[ch] = CopyChannel{};
        }
    }

private:
    int getBase(size_t pos) const {
        return (history_[pos >> 2] >> (2 * (pos & 3))) & 3;
    }

    void pushBase(int base) {
        // Bound history to prevent unbounded memory growth on large files.
        // When limit is reached, clear history and deactivate channels.
        // Hash table entries become stale (rejected by position validity checks).
        // Encoder and decoder hit this at the same position, staying synchronized.
        if (historyLen_ >= MAX_HISTORY_BASES) {
            history_.clear();
            historyLen_ = 0;
            for (int ch = 0; ch < NUM_CHANNELS; ch++) {
                channels_[ch] = CopyChannel{};
            }
        }
        if ((historyLen_ & 3) == 0) history_.push_back(0);
        history_.back() |= static_cast<uint8_t>(base << (2 * (historyLen_ & 3)));
        historyLen_++;
    }

    static uint32_t hashCtx(uint64_t ctx) {
        return static_cast<uint32_t>((ctx * FIBONACCI_HASH_CONSTANT) >> (64 - HASH_BITS));
    }

    int anchorOrder_;
    mutable CopyChannel channels_[NUM_CHANNELS];
    std::vector<HashEntry> hashTable_;       // Long anchor (order-16)
    std::vector<HashEntry> shortHashTable_;  // Short anchor (order-8) for diverged repeats
    std::vector<uint8_t> history_;     // 2-bit packed
    size_t historyLen_ = 0;
};

// ============================================================================
// OnlineContextModel — adaptive PPM* with neural mixer and V₄ symmetry
// ============================================================================

class OnlineContextModel {
public:
    static constexpr int NUM_MIXER_CTX = 64;  // Context-dependent mixer array size (V6+)
    static constexpr int K_MATCH_SLOTS = MatchModel::K_MATCHES;  // Multi-match slots
    static constexpr int K_COPY_CHANNELS = CopyModel::NUM_CHANNELS;  // Long-range copy channels

    // modelVersion 7: adds V₄-canonical long-range copy model
    // modelVersion 6: SSE chain + context mixer array + multi-match
    // modelVersion 5: TwoLayerMixer + APM + MatchModel + V₄-symmetric enhancements
    OnlineContextModel(int maxOrder = 16, int poolingOrder = -1,
                       bool /*unused*/ = true, int modelVersion = 7)
        : maxOrder_(maxOrder),
          poolingOrder_(poolingOrder < 0 ? maxOrder : poolingOrder),
          modelVersion_(std::max(3, modelVersion)),
          subs13_(std::min(13, maxOrder), 1),
          subs17_(maxOrder, 3),
          apm_(16),
          ctx2_(0) {
        fullMask_ = (1ULL << (2 * maxOrder_)) - 1;

        // Mixer inputs: orders 0..maxOrder + subs13 + subs17 [+ match(es)] [+ coset]
        numModels_ = maxOrder_ + 3;
        subs13Idx_ = maxOrder_ + 1;
        subs17Idx_ = maxOrder_ + 2;
        if (modelVersion_ >= 4) {
            matchIdx_ = numModels_++;
            if (modelVersion_ >= 6) {
                // Multi-match: K-1 additional match slots as extra mixer inputs
                for (int m = 1; m < K_MATCH_SLOTS; m++) {
                    extraMatchIdx_[m] = numModels_++;
                }
                matchModel_ = MatchModel(maxOrder_, true);
            } else {
                matchModel_ = MatchModel(modelVersion_ >= 5 ? maxOrder_ : 0, false);
            }
        }
        if (modelVersion_ >= 5) {
            cosetIdx_ = numModels_++;
            subs13_.v4EditClass = true;
            subs17_.v4EditClass = true;
        }
        if (modelVersion_ >= 7) {
            // Long-range copy model: 4 channels as additional mixer inputs
            for (int ch = 0; ch < K_COPY_CHANNELS; ch++) {
                copyModelIdx_[ch] = numModels_++;
            }
            copyModel_ = CopyModel(maxOrder_);
        }

        // Initialize mixer(s)
        if (modelVersion_ >= 6) {
            // Context-dependent mixer array
            for (int i = 0; i < NUM_MIXER_CTX; i++) {
                mixers_[i] = TwoLayerMixer(numModels_);
            }
            // SSE chain: APM₂ (8 ctx), APM₃ (8 ctx)
            apm2_ = AdaptiveProbabilityMap(8);
            apm3_ = AdaptiveProbabilityMap(8);
        } else {
            mixers_[0] = TwoLayerMixer(numModels_);
        }

        // Flat arrays for orders 0-5 (direct-indexed, cache-friendly)
        flatMaxOrder_ = std::min(5, maxOrder_);
        flatCounts_.resize(flatMaxOrder_ + 1);
        for (int o = 0; o <= flatMaxOrder_; o++) {
            size_t sz = 1ULL << (2 * o);
            flatCounts_[o].assign(sz, {1, 1, 1, 1, 4});
        }

        // Precompute orbit sizes for orbit-stratified forgetting (modelVersion 5+)
        if (modelVersion_ >= 5) {
            flatOrbitSizes_.resize(flatMaxOrder_ + 1);
            for (int o = 0; o <= flatMaxOrder_; o++) {
                size_t sz = 1ULL << (2 * o);
                flatOrbitSizes_[o].resize(sz, 4);
                for (uint64_t ctx = 0; ctx < sz; ctx++) {
                    flatOrbitSizes_[o][ctx] = static_cast<uint8_t>(orbitSize(ctx, o));
                }
            }
        }

        // Hash maps for orders 6+
        // All orders are bounded to prevent unbounded memory growth.
        // Orders 6-9: capped at 2× theoretical k-mer space (4^k).
        // Orders 10+: capped at 2M entries (eviction keeps most-used contexts).
        static constexpr size_t MAX_HASH_ENTRIES = 1 << 21;  // 2M entries per order
        if (maxOrder_ > flatMaxOrder_) {
            hashCounts_.resize(maxOrder_ - flatMaxOrder_);
            for (int i = 0; i < static_cast<int>(hashCounts_.size()); i++) {
                hashCounts_[i].reserve(1 << 16);
                int order = flatMaxOrder_ + 1 + i;
                // Graduated caps: lower orders have fewer possible contexts
                size_t cap = MAX_HASH_ENTRIES;
                if (order < 10) {
                    // 4^k possible k-mers, ×2 safety margin
                    cap = std::min(MAX_HASH_ENTRIES, static_cast<size_t>(2) << (2 * order));
                }
                hashCounts_[i].setMaxCapacity(cap);
            }
        }
    }

    void update(uint64_t fullContext, int base, int validOrder) {
        updateV3(fullContext, base, validOrder);
    }

    void update(uint64_t fullContext, int base) {
        update(fullContext, base, maxOrder_);
    }

    std::array<uint32_t, 5> getCDF(uint64_t fullContext, int validOrder) const {
        return getCDFV3(fullContext, validOrder);
    }

    std::array<uint32_t, 5> getCDF(uint64_t fullContext) const {
        return getCDF(fullContext, maxOrder_);
    }

    void resetPerReadState() {
        if (!persistSubs_) {
            subs13_.reset();
            subs17_.reset();
        }
        if (modelVersion_ >= 4) {
            matchModel_.reset();
        }
        if (modelVersion_ >= 5) {
            transformEMA_ = {0.25, 0.25, 0.25, 0.25};
        }
        if (modelVersion_ >= 7) {
            copyModel_.resetChannels();  // Keep history, just deactivate channels
        }
        numModelsActive_ = 0;
        ctx2_ = 0;
        lastMixerCtx_ = 0;
    }

    void setPersistSubs(bool persist) { persistSubs_ = persist; }

    void reset() {
        for (auto& fc : flatCounts_) {
            std::fill(fc.begin(), fc.end(), CompactCounts{1, 1, 1, 1, 4});
        }
        for (auto& hm : hashCounts_) hm.clear();
        if (modelVersion_ >= 6) {
            for (int i = 0; i < NUM_MIXER_CTX; i++) {
                mixers_[i] = TwoLayerMixer(numModels_);
            }
            apm2_ = AdaptiveProbabilityMap(8);
            apm3_ = AdaptiveProbabilityMap(8);
        } else {
            mixers_[0] = TwoLayerMixer(numModels_);
        }
        apm_ = AdaptiveProbabilityMap(16);
        if (modelVersion_ >= 4) {
            matchModel_.reset();
        }
        if (modelVersion_ >= 5) {
            transformEMA_ = {0.25, 0.25, 0.25, 0.25};
        }
        if (modelVersion_ >= 7) {
            copyModel_.reset();
        }
        subs13_.reset();
        subs17_.reset();
        numModelsActive_ = 0;
        ctx2_ = 0;
        lastMixerCtx_ = 0;
    }

    int maxOrder() const { return maxOrder_; }
    int modelVersion() const { return modelVersion_; }

private:
    // ========== Core model: TwoLayerMixer + APM + compact counts ==========

    void updateV3(uint64_t fullContext, int base, int validOrder) {
        if (numModelsActive_ > 0) {
            // Compute top prediction for subs context update
            int topPred = 0;
            double maxProb = lastMixedProbs_[0];
            for (int b = 1; b < 4; b++) {
                if (lastMixedProbs_[b] > maxProb) {
                    maxProb = lastMixedProbs_[b];
                    topPred = b;
                }
            }

            // Update APM chain
            apm_.update(lastPreAPMProbs_, ctx2_, base);
            if (modelVersion_ >= 6) {
                apm2_.update(lastPostAPM1Probs_, lastAPM2Ctx_, base);
                apm3_.update(lastPostAPM2Probs_, lastAPM3Ctx_, base);
            }

            // Gradient step on mixer (context-dependent in V6+)
            if (modelVersion_ >= 6) {
                mixers_[lastMixerCtx_].update(lastOrderProbs_.data(), lastMixedProbs_, base);
            } else {
                mixers_[0].update(lastOrderProbs_.data(), lastMixedProbs_, base);
            }

            // Update substitutional contexts
            subs13_.update(topPred, base);
            subs17_.update(topPred, base);
        }

        // Update match model (V4+)
        if (modelVersion_ >= 4) {
            matchModel_.update(fullContext, base);
        }

        // Update copy model (V7+)
        if (modelVersion_ >= 7) {
            copyModel_.update(fullContext, base);
        }

        // Update count tables at all valid orders
        int maxOrd = std::min(maxOrder_, validOrder);
        for (int order = maxOrd; order >= 0; order--) {
            auto ctxT = getContextWithTransform(fullContext, order);
            int b = needsComplement(ctxT.second) ? complementBase(base) : base;
            auto& counts = getCompactCounts(order, ctxT.first);
            counts[b]++;
            counts[4]++;
            // Orbit-stratified forgetting (V5+): use precomputed orbit sizes for flat orders
            if (modelVersion_ >= 5 && order <= flatMaxOrder_) {
                applyForgetting(counts, flatOrbitSizes_[order][ctxT.first]);
            } else {
                applyForgetting(counts);
            }
        }

        // Update coset (strand orientation) EMA (V5+)
        if (modelVersion_ >= 5) {
            auto ctxT = getContextWithTransform(fullContext, maxOrd);
            int tIdx = static_cast<int>(ctxT.second);
            constexpr double COSET_ALPHA = 0.05;
            for (int t = 0; t < 4; t++) {
                transformEMA_[t] = (1.0 - COSET_ALPHA) * transformEMA_[t]
                                 + (t == tIdx ? COSET_ALPHA : 0.0);
            }
        }

        // Update secondary context for APM
        ctx2_ = ((ctx2_ << 2) | base) & 0xF;  // Last 2 bases
    }

    std::array<uint32_t, 5> getCDFV3(uint64_t fullContext, int validOrder) const {
        static const std::array<double, 4> UNIFORM = {0.25, 0.25, 0.25, 0.25};
        for (int i = 0; i < numModels_; i++) lastOrderProbs_[i] = UNIFORM;
        numModelsActive_ = numModels_;

        int maxOrd = std::min(maxOrder_, validOrder);

        // Prefetch next order's hash entry while processing current
        for (int order = 0; order <= maxOrd; order++) {
            auto ctxT = getContextWithTransform(fullContext, order);

            // Prefetch for next order if it's a hash-based order
            if (order + 1 <= maxOrd && order + 1 > flatMaxOrder_) {
                auto nextCtxT = getContextWithTransform(fullContext, order + 1);
                int hmIdx = order + 1 - flatMaxOrder_ - 1;
                if (hmIdx >= 0 && hmIdx < static_cast<int>(hashCounts_.size())) {
                    hashCounts_[hmIdx].prefetch(nextCtxT.first);
                }
            }

            const CompactCounts* counts = findCompactCounts(order, ctxT.first);
            if (counts) {
                uint32_t total = (*counts)[0] + (*counts)[1] + (*counts)[2] + (*counts)[3];
                if (total > 0) [[likely]] {
                    bool compl_ = needsComplement(ctxT.second);
                    std::array<double, 4> probs;
                    for (int b = 0; b < 4; b++) {
                        int idx = compl_ ? complementBase(b) : b;
                        probs[b] = static_cast<double>((*counts)[idx]) / total;
                    }
                    lastOrderProbs_[order] = probs;
                }
            }
        }

        // Substitutional contexts
        if (validOrder >= subs13_.order && subs13_.isActive()) {
            auto ctxT = canonicalWithTransform(subs13_.editedCtx, subs13_.order);
            const CompactCounts* counts = findCompactCounts(subs13_.order, ctxT.first);
            if (counts) {
                uint32_t total = (*counts)[0] + (*counts)[1] + (*counts)[2] + (*counts)[3];
                if (total > 0) {
                    bool compl_ = needsComplement(ctxT.second);
                    std::array<double, 4> probs;
                    for (int b = 0; b < 4; b++) {
                        int idx = compl_ ? complementBase(b) : b;
                        probs[b] = static_cast<double>((*counts)[idx]) / total;
                    }
                    lastOrderProbs_[subs13Idx_] = probs;
                }
            }
        }

        if (validOrder >= subs17_.order && subs17_.isActive()) {
            auto ctxT = canonicalWithTransform(subs17_.editedCtx, subs17_.order);
            const CompactCounts* counts = findCompactCounts(subs17_.order, ctxT.first);
            if (counts) {
                uint32_t total = (*counts)[0] + (*counts)[1] + (*counts)[2] + (*counts)[3];
                if (total > 0) {
                    bool compl_ = needsComplement(ctxT.second);
                    std::array<double, 4> probs;
                    for (int b = 0; b < 4; b++) {
                        int idx = compl_ ? complementBase(b) : b;
                        probs[b] = static_cast<double>((*counts)[idx]) / total;
                    }
                    lastOrderProbs_[subs17Idx_] = probs;
                }
            }
        }

        // Match model prediction (V4+, V₄-canonical in V5+)
        if (modelVersion_ >= 4) {
            lastOrderProbs_[matchIdx_] = matchModel_.predict(fullContext, 0);
            // Multi-match: additional match slots (V6+)
            if (modelVersion_ >= 6) {
                for (int m = 1; m < K_MATCH_SLOTS; m++) {
                    lastOrderProbs_[extraMatchIdx_[m]] = matchModel_.predict(fullContext, m);
                }
            }
        }

        // Coset (strand orientation) feature (V5+)
        if (modelVersion_ >= 5) {
            lastOrderProbs_[cosetIdx_] = transformEMA_;
        }

        // Long-range copy model predictions (V7+)
        if (modelVersion_ >= 7) {
            for (int ch = 0; ch < K_COPY_CHANNELS; ch++) {
                lastOrderProbs_[copyModelIdx_[ch]] = copyModel_.predict(ch);
            }
        }

        // Context-dependent mixer selection (V6+) or single mixer (V5-)
        if (modelVersion_ >= 6) {
            // Mixer context = (order-2 context bits, match bucket)
            int matchBucket = 0;
            int mLen = matchModel_.matchLen(0);
            if (mLen >= 8) matchBucket = 3;
            else if (mLen >= 4) matchBucket = 2;
            else if (mLen >= 2) matchBucket = 1;

            lastMixerCtx_ = ((static_cast<int>(fullContext) & 0xF) << 2) | matchBucket;
            lastMixerCtx_ = lastMixerCtx_ % NUM_MIXER_CTX;
            lastPreAPMProbs_ = mixers_[lastMixerCtx_].predict(lastOrderProbs_.data(), numModels_);
        } else {
            lastPreAPMProbs_ = mixers_[0].predict(lastOrderProbs_.data(), numModels_);
        }

        // Soft V₄-equivariant blending (V5+)
        if (modelVersion_ >= 5) {
            constexpr double SYM_WEIGHT = 0.05;
            double avg_AT = (lastPreAPMProbs_[0] + lastPreAPMProbs_[3]) * 0.5;
            double avg_CG = (lastPreAPMProbs_[1] + lastPreAPMProbs_[2]) * 0.5;
            lastPreAPMProbs_[0] += SYM_WEIGHT * (avg_AT - lastPreAPMProbs_[0]);
            lastPreAPMProbs_[3] += SYM_WEIGHT * (avg_AT - lastPreAPMProbs_[3]);
            lastPreAPMProbs_[1] += SYM_WEIGHT * (avg_CG - lastPreAPMProbs_[1]);
            lastPreAPMProbs_[2] += SYM_WEIGHT * (avg_CG - lastPreAPMProbs_[2]);
        }

        // APM₁ refinement (last-2-bases context, 16 ctx)
        lastMixedProbs_ = apm_.refine(lastPreAPMProbs_, ctx2_);

        // SSE chain: APM₂ and APM₃ (V6+)
        if (modelVersion_ >= 6) {
            lastPostAPM1Probs_ = lastMixedProbs_;

            // APM₂: order-1 context (4) × match_active (2) = 8 contexts
            int matchActive = (matchModel_.matchLen(0) >= 2) ? 1 : 0;
            lastAPM2Ctx_ = (static_cast<int>(fullContext) & 3) * 2 + matchActive;
            lastMixedProbs_ = apm2_.refine(lastMixedProbs_, lastAPM2Ctx_);
            lastPostAPM2Probs_ = lastMixedProbs_;

            // APM₃: best-order bucket (which orders dominate)
            lastAPM3Ctx_ = computeBestOrderBucket(lastOrderProbs_.data(), maxOrd);
            lastMixedProbs_ = apm3_.refine(lastMixedProbs_, lastAPM3Ctx_);
        }

        return probsToCDF(lastMixedProbs_);
    }

    // Compact count access: flat array for orders 0-5, hash map for 6+
    CompactCounts& getCompactCounts(int order, uint64_t ctx) {
        if (order <= flatMaxOrder_) {
            return flatCounts_[order][ctx];
        }
        int hmIdx = order - flatMaxOrder_ - 1;
        auto* val = hashCounts_[hmIdx].find(ctx);
        if (!val) {
            hashCounts_[hmIdx][ctx] = {1, 1, 1, 1, 4};
            return hashCounts_[hmIdx][ctx];
        }
        return *val;
    }

    const CompactCounts* findCompactCounts(int order, uint64_t ctx) const {
        if (order <= flatMaxOrder_) {
            return &flatCounts_[order][ctx];
        }
        int hmIdx = order - flatMaxOrder_ - 1;
        return hashCounts_[hmIdx].find(ctx);
    }

    // ========== Utilities ==========

    // Compute best-order bucket for APM₃ context (8 buckets)
    // Partitions based on which order range has the sharpest (most confident) prediction
    static int computeBestOrderBucket(const std::array<double, 4>* orderProbs, int maxOrd) {
        // Find the order with the highest max probability (sharpest prediction)
        double bestConf = 0.0;
        int bestOrder = 0;
        int limit = maxOrd + 1;
        for (int o = 0; o < limit; o++) {
            double maxP = *std::max_element(orderProbs[o].begin(), orderProbs[o].end());
            if (maxP > bestConf) {
                bestConf = maxP;
                bestOrder = o;
            }
        }
        // Map best order to 8 buckets: 0-1, 2-3, 4-5, 6-7, 8-9, 10-11, 12-13, 14+
        return std::min(7, bestOrder / 2);
    }

    static std::array<uint32_t, 5> probsToCDF(const std::array<double, 4>& probs) {
        std::array<uint32_t, 5> cdf;
        cdf[0] = 0;
        double cum = 0.0;
        for (int i = 0; i < 3; i++) {
            cum += probs[i];
            uint32_t val = static_cast<uint32_t>(cum * PROB_SCALE + 0.5);
            cdf[i+1] = std::max(cdf[i] + 1, val);
        }
        cdf[4] = PROB_SCALE;
        if (cdf[3] >= cdf[4]) {
            cdf[3] = cdf[4] - 1;
            if (cdf[2] >= cdf[3]) cdf[2] = cdf[3] - 1;
            if (cdf[1] >= cdf[2]) cdf[1] = cdf[2] - 1;
        }
        return cdf;
    }

    std::pair<uint64_t, Transform> getContextWithTransform(
            uint64_t fullContext, int order) const {
        uint64_t mask = (1ULL << (2 * order)) - 1;
        uint64_t ctx = fullContext & mask;
        if (order > 0) {
            return canonicalWithTransform(ctx, order);
        }
        return {ctx, Transform::I};
    }

    static bool needsComplement(Transform t) {
        return t == Transform::C || t == Transform::CR;
    }

    static int complementBase(int base) {
        return 3 - base;
    }

    // Member variables
    int maxOrder_;
    int poolingOrder_;
    int modelVersion_;
    mutable SubstitutionalContext subs13_;
    mutable SubstitutionalContext subs17_;

    // APM chain: APM₁ (16 ctx), APM₂ (8 ctx), APM₃ (8 ctx) — V6+
    mutable AdaptiveProbabilityMap apm_;
    mutable AdaptiveProbabilityMap apm2_{8};
    mutable AdaptiveProbabilityMap apm3_{8};

    mutable int ctx2_ = 0;

    // Derived / initialized in constructor body
    uint64_t fullMask_;
    int numModels_ = 0;
    int subs13Idx_ = 0;
    int subs17Idx_ = 0;
    int matchIdx_ = 0;
    int extraMatchIdx_[K_MATCH_SLOTS] = {};  // Multi-match mixer input indices (V6+)
    int cosetIdx_ = 0;
    bool persistSubs_ = false;

    // V3 data structures (compact counts + flat arrays + hash maps)
    int flatMaxOrder_ = 5;
    mutable std::vector<std::vector<CompactCounts>> flatCounts_;
    mutable std::vector<FlatHashMap<CompactCounts>> hashCounts_;

    // Context-dependent mixer array (V6: 64 mixers; V5-: only mixers_[0] used)
    mutable std::array<TwoLayerMixer, NUM_MIXER_CTX> mixers_;
    mutable int lastMixerCtx_ = 0;
    mutable std::array<double, 4> lastPreAPMProbs_ = {0.25, 0.25, 0.25, 0.25};
    mutable std::array<double, 4> lastPostAPM1Probs_ = {0.25, 0.25, 0.25, 0.25};
    mutable std::array<double, 4> lastPostAPM2Probs_ = {0.25, 0.25, 0.25, 0.25};
    mutable int lastAPM2Ctx_ = 0;
    mutable int lastAPM3Ctx_ = 0;

    // V₄-symmetric enhancements
    std::vector<std::vector<uint8_t>> flatOrbitSizes_;
    mutable std::array<double, 4> transformEMA_ = {0.25, 0.25, 0.25, 0.25};

    // Match model (single or multi-match)
    mutable MatchModel matchModel_;

    // Long-range copy model (V7+)
    int copyModelIdx_[K_COPY_CHANNELS] = {};
    mutable CopyModel copyModel_;

    // Shared state
    static constexpr int MAX_MODELS = 40;  // Must be >= numModels_
    mutable int numModelsActive_ = 0;
    mutable std::array<std::array<double, 4>, MAX_MODELS> lastOrderProbs_;
    mutable std::array<double, 4> lastMixedProbs_ = {0.25, 0.25, 0.25, 0.25};
};

// ============================================================================
// Interleaved rANS encoder/decoder (4-way)
// ============================================================================

class InterleavedRANSEncoder {
public:
    static constexpr int NUM_STATES = 4;

    static std::vector<uint8_t> encodeAll(
        const std::vector<int>& symbols,
        const std::vector<std::array<uint32_t, 5>>& cdfs
    ) {
        if (symbols.empty()) {
            // Write 4 initial states
            std::vector<uint8_t> out(NUM_STATES * 4, 0);
            return out;
        }

        // N separate encoders
        std::array<uint32_t, NUM_STATES> states;
        std::array<std::vector<uint8_t>, NUM_STATES> outputs;
        for (int s = 0; s < NUM_STATES; s++) {
            states[s] = RANS_L;
        }

        // Encode in reverse order, round-robin assignment
        for (int i = static_cast<int>(symbols.size()) - 1; i >= 0; i--) {
            int stateIdx = i % NUM_STATES;
            uint8_t sym = static_cast<uint8_t>(symbols[i]);
            uint32_t start = cdfs[i][sym];
            uint32_t freq = cdfs[i][sym + 1] - start;
            assert(freq > 0 && "CDF must not contain zero-frequency symbols");

            uint32_t limit = ((RANS_L >> PROB_BITS) << 8) * freq;
            while (states[stateIdx] >= limit) {
                outputs[stateIdx].push_back(static_cast<uint8_t>(states[stateIdx] & 0xFF));
                states[stateIdx] >>= 8;
            }
            states[stateIdx] = ((states[stateIdx] / freq) << PROB_BITS) +
                               (states[stateIdx] % freq) + start;
        }

        // Build output: [4 final states][interleaved byte streams]
        std::vector<uint8_t> result;
        // Write states
        for (int s = 0; s < NUM_STATES; s++) {
            for (int b = 3; b >= 0; b--) {
                result.push_back(static_cast<uint8_t>((states[s] >> (b * 8)) & 0xFF));
            }
        }
        // Write stream sizes
        for (int s = 0; s < NUM_STATES; s++) {
            uint32_t sz = static_cast<uint32_t>(outputs[s].size());
            result.push_back(sz & 0xFF);
            result.push_back((sz >> 8) & 0xFF);
            result.push_back((sz >> 16) & 0xFF);
            result.push_back((sz >> 24) & 0xFF);
        }
        // Write streams (reversed, as rANS output is reversed)
        for (int s = 0; s < NUM_STATES; s++) {
            auto& out = outputs[s];
            std::reverse(out.begin(), out.end());
            result.insert(result.end(), out.begin(), out.end());
        }

        return result;
    }
};

class InterleavedRANSDecoder {
public:
    static constexpr int NUM_STATES = 4;

    template<typename CDFGetter, typename Emitter>
    static void decodeIncremental(
        const uint8_t* data, size_t len, size_t count,
        CDFGetter getCDF,
        Emitter emit
    ) {
        if (count == 0) return;

        // Read 4 initial states
        size_t pos = 0;
        std::array<uint32_t, NUM_STATES> states;
        for (int s = 0; s < NUM_STATES; s++) {
            states[s] = 0;
            for (int b = 0; b < 4 && pos < len; b++) {
                states[s] = (states[s] << 8) | data[pos++];
            }
        }

        // Read stream sizes
        std::array<uint32_t, NUM_STATES> streamSizes;
        for (int s = 0; s < NUM_STATES; s++) {
            streamSizes[s] = 0;
            for (int b = 0; b < 4 && pos < len; b++) {
                streamSizes[s] |= static_cast<uint32_t>(data[pos++]) << (b * 8);
            }
        }

        // Compute stream offsets
        std::array<size_t, NUM_STATES> streamPos;
        std::array<size_t, NUM_STATES> streamEnd;
        size_t off = pos;
        for (int s = 0; s < NUM_STATES; s++) {
            streamPos[s] = off;
            off += streamSizes[s];
            streamEnd[s] = off;
        }

        // Decode symbols round-robin
        for (size_t i = 0; i < count; i++) {
            int stateIdx = static_cast<int>(i % NUM_STATES);
            auto cdf = getCDF(i);

            uint32_t slot = states[stateIdx] & (PROB_SCALE - 1);
            uint8_t symbol = static_cast<uint8_t>(searchCDF(cdf, 4, slot));

            uint32_t start = cdf[symbol];
            uint32_t freq = cdf[symbol + 1] - start;
            assert(freq > 0 && "CDF must not contain zero-frequency symbols");

            states[stateIdx] = freq * (states[stateIdx] >> PROB_BITS) + slot - start;

            // Renormalize from this state's stream
            while (states[stateIdx] < RANS_L && streamPos[stateIdx] < streamEnd[stateIdx]) {
                states[stateIdx] = (states[stateIdx] << 8) | data[streamPos[stateIdx]++];
            }

            emit(i, symbol);
        }
    }
};

// ============================================================================
// CompressorV4 — FASTA compression with version-aware format
// ============================================================================

class CompressorV4 {
public:
    explicit CompressorV4(int maxOrder = 16, int poolingOrder = 16)
        : maxOrder_(maxOrder), poolingOrder_(poolingOrder) {}

    std::vector<uint8_t> compress(const std::string& sequence, CompressionStats* stats = nullptr) {
        return compressWithHeader("", sequence, stats);
    }

    std::vector<uint8_t> compressWithHeader(const std::string& header,
                                             const std::string& rawSequence,
                                             CompressionStats* stats = nullptr) {
        if (rawSequence.empty()) {
            return createEmptyOutput(header, stats);
        }

        // Separate DNA from non-DNA characters
        std::string dnaSequence;
        std::vector<NonDNARun> nonDNARuns;
        dnaSequence.reserve(rawSequence.size());

        size_t i = 0;
        while (i < rawSequence.size()) {
            char c = rawSequence[i];
            if (isValidBase(c)) {
                dnaSequence.push_back(toUpper(c));
                i++;
            } else {
                size_t runStart = i;
                char runChar = c;
                while (i < rawSequence.size() && rawSequence[i] == runChar) {
                    dnaSequence.push_back('A');
                    i++;
                }
                nonDNARuns.push_back({runStart, i - runStart, runChar});
            }
        }

        if (dnaSequence.size() <= static_cast<size_t>(maxOrder_)) {
            return compressShortSequence(header, rawSequence, nonDNARuns, stats);
        }

        // V7 model: SSE chain + context mixer + multi-match + copy model
        OnlineContextModel model(maxOrder_, poolingOrder_, true, 7);
        model.setPersistSubs(true);
        uint64_t fullMask = (1ULL << (2 * maxOrder_)) - 1;

        std::vector<int> symbols;
        std::vector<std::array<uint32_t, 5>> cdfs;
        symbols.reserve(dnaSequence.size());
        cdfs.reserve(dnaSequence.size());

        std::array<uint32_t, 5> uniformCDF = {
            0, PROB_SCALE/4, PROB_SCALE/2, 3*PROB_SCALE/4, PROB_SCALE
        };

        uint64_t ctx = 0;
        for (size_t ii = 0; ii < dnaSequence.size(); ii++) {
            int base = charToBase(dnaSequence[ii]);
            symbols.push_back(base);

            if (ii == 0) {
                cdfs.push_back(uniformCDF);
            } else {
                int validOrder = static_cast<int>(std::min(ii, static_cast<size_t>(maxOrder_)));
                cdfs.push_back(model.getCDF(ctx, validOrder));
            }

            if (ii > 0) {
                int updateOrder = static_cast<int>(std::min(ii - 1, static_cast<size_t>(maxOrder_)));
                model.update(ctx, base, updateOrder);
            }

            ctx = ((ctx << 2) | base) & fullMask;
        }

        // Use interleaved rANS for V4
        auto encodedData = InterleavedRANSEncoder::encodeAll(symbols, cdfs);

        // Build output
        BinaryWriter writer;
        writer.writeU32LE(MAGIC_V4);
        writer.writeU8(VERSION_V4_MAJOR);
        writer.writeU8(VERSION_V4_MINOR);
        writer.writeU8(static_cast<uint8_t>(maxOrder_));
        writer.writeU8(static_cast<uint8_t>(poolingOrder_));
        writer.writeU64LE(rawSequence.size());

        writer.writeU32LE(static_cast<uint32_t>(header.size()));
        if (!header.empty()) {
            writer.writeBytes(std::vector<uint8_t>(header.begin(), header.end()));
        }

        writer.writeU32LE(static_cast<uint32_t>(nonDNARuns.size()));
        if (!nonDNARuns.empty()) {
            std::vector<uint8_t> nonDNAData;
            uint64_t prevPos = 0;
            for (const auto& run : nonDNARuns) {
                writeVarint(nonDNAData, run.position - prevPos);
                writeVarint(nonDNAData, run.length);
                nonDNAData.push_back(static_cast<uint8_t>(run.character));
                prevPos = run.position + run.length;
            }
            writer.writeU32LE(static_cast<uint32_t>(nonDNAData.size()));
            writer.writeBytes(nonDNAData);
        }

        writer.writeU32LE(static_cast<uint32_t>(encodedData.size()));
        writer.writeBytes(encodedData);

        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());

        if (stats) {
            stats->originalBytes = rawSequence.size();
            stats->compressedBytes = writer.data().size();
            stats->numKmers = dnaSequence.size();
            stats->numOrbits = nonDNARuns.size();
            stats->bitsPerBase = (writer.data().size() * 8.0) / rawSequence.size();
            stats->compressionRatio = static_cast<double>(rawSequence.size()) / writer.data().size();
        }

        return writer.data();
    }

    std::string decompress(const std::vector<uint8_t>& data) {
        std::string header;
        return decompressWithHeader(data, header);
    }

    std::string decompressWithHeader(const std::vector<uint8_t>& data, std::string& header) {
        BinaryReader reader(data);

        uint32_t magic = reader.readU32LE();
        if (magic != MAGIC_V4) {
            throw std::runtime_error("Invalid magic for V4");
        }

        uint8_t versionMajor = reader.readU8();
        uint8_t versionMinor = reader.readU8();
        int maxOrder = reader.readU8();
        int poolingOrder = reader.readU8();
        uint64_t originalLength = reader.readU64LE();

        if (originalLength == 0) {
            header = "";
            return "";
        }

        // Validate CRC
        size_t dataEnd = data.size() - 4;
        CRC32 crc;
        crc.update(data.data(), dataEnd);
        if (crc.finalize() != BinaryReader(data.data() + dataEnd, 4).readU32LE()) {
            throw std::runtime_error("CRC mismatch");
        }

        (void)versionMajor;
        int fastaModelVersion = (versionMinor >= 8) ? 7 : (versionMinor >= 7) ? 6 : (versionMinor >= 6) ? 5 : 3;

        // Read FASTA header
        uint32_t headerLen = reader.readU32LE();
        if (headerLen > 0) {
            auto headerBytes = reader.readBytes(headerLen);
            header = std::string(headerBytes.begin(), headerBytes.end());
        } else {
            header = "";
        }

        // Read non-DNA runs (N characters and other ambiguous bases)
        std::vector<NonDNARun> nonDNARuns;
        uint32_t nonDNARunCount = reader.readU32LE();
        if (nonDNARunCount > 0) {
            uint32_t nonDNADataSize = reader.readU32LE();
            auto nonDNAData = reader.readBytes(nonDNADataSize);
            const uint8_t* p = nonDNAData.data();
            const uint8_t* end = p + nonDNAData.size();
            uint64_t pos = 0;
            for (uint32_t i = 0; i < nonDNARunCount && p < end; i++) {
                uint64_t delta = readVarint(p);
                pos += delta;
                uint64_t length = readVarint(p);
                char c = static_cast<char>(*p++);
                nonDNARuns.push_back({pos, length, c});
                pos += length;
            }
        }

        if (originalLength <= static_cast<size_t>(maxOrder)) {
            uint32_t rawSize = reader.readU32LE();
            auto rawData = reader.readBytes(rawSize);
            std::string result;
            int bitPos = 0;
            size_t bytePos = 0;
            for (size_t i = 0; i < originalLength; i++) {
                int base = (rawData[bytePos] >> (6 - bitPos)) & 3;
                result.push_back(baseToChar(base));
                bitPos += 2;
                if (bitPos == 8) { bitPos = 0; bytePos++; }
            }
            for (const auto& run : nonDNARuns) {
                for (uint64_t j = 0; j < run.length && run.position + j < result.size(); j++) {
                    result[run.position + j] = run.character;
                }
            }
            return result;
        }

        uint32_t encodedSize = reader.readU32LE();
        auto encodedData = reader.readBytes(encodedSize);

        std::string result;
        result.resize(originalLength);

        OnlineContextModel model(maxOrder, poolingOrder, true, fastaModelVersion);
        model.setPersistSubs(true);

        uint64_t fullMask = (1ULL << (2 * maxOrder)) - 1;
        uint64_t ctx = 0;

        std::array<uint32_t, 5> uniformCDF = {
            0, PROB_SCALE/4, PROB_SCALE/2, 3*PROB_SCALE/4, PROB_SCALE
        };

        auto getCDFLambda = [&](size_t i) -> std::array<uint32_t, 5> {
            if (i == 0) return uniformCDF;
            int validOrder = static_cast<int>(std::min(i, static_cast<size_t>(maxOrder)));
            return model.getCDF(ctx, validOrder);
        };

        auto emitLambda = [&](size_t i, int base) {
            result[i] = baseToChar(base);
            if (i > 0) {
                int updateOrder = static_cast<int>(std::min(i - 1, static_cast<size_t>(maxOrder)));
                model.update(ctx, base, updateOrder);
            }
            ctx = ((ctx << 2) | base) & fullMask;
        };

        InterleavedRANSDecoder::decodeIncremental(
            encodedData.data(), encodedData.size(), originalLength,
            getCDFLambda, emitLambda);

        for (const auto& run : nonDNARuns) {
            for (uint64_t j = 0; j < run.length && run.position + j < result.size(); j++) {
                result[run.position + j] = run.character;
            }
        }

        return result;
    }

    std::vector<uint8_t> compressFile(const std::string& path, CompressionStats* stats = nullptr) {
        FastaReader fastaReader(path);
        auto records = fastaReader.readAll();

        if (records.empty()) {
            return createEmptyOutput("", stats);
        }

        // Always use multi-record format (handles 1 or more records uniformly)
        return compressMultiRecord(records, stats);
    }

    std::vector<uint8_t> compressMultiRecord(const std::vector<FastaRecord>& records,
                                              CompressionStats* stats = nullptr) {
        std::vector<std::string> headers;
        headers.reserve(records.size());
        for (const auto& rec : records) headers.push_back(rec.header);

        IdentifierEncoderV2 idEncoder;
        auto identifierStream = idEncoder.encode(headers);

        std::vector<uint8_t> seqLengthsData;
        int64_t prevLen = 0;
        for (const auto& rec : records) {
            int64_t delta = static_cast<int64_t>(rec.sequence.size()) - prevLen;
            writeSignedVarint(seqLengthsData, delta);
            prevLen = static_cast<int64_t>(rec.sequence.size());
        }

        // RLE ambiguous encoding: merge consecutive same-char positions into runs
        std::vector<uint8_t> ambiguousData;
        std::vector<std::vector<AmbiguousRun>> allAmbiguousRuns;
        allAmbiguousRuns.reserve(records.size());
        for (const auto& rec : records) {
            auto runs = ambiguousToRuns(rec.ambiguous);
            writeVarint(ambiguousData, runs.size());
            uint64_t prevEnd = 0;
            for (const auto& run : runs) {
                writeVarint(ambiguousData, run.position - prevEnd);
                writeVarint(ambiguousData, run.length);
                ambiguousData.push_back(static_cast<uint8_t>(run.character));
                prevEnd = run.position + run.length;
            }
            allAmbiguousRuns.push_back(std::move(runs));
        }

        // Build ACGT-only DNA string, skipping ambiguous positions
        // Also build boundaries: context resets at record starts + after each N-gap
        size_t totalSeqLen = 0;
        for (const auto& rec : records) totalSeqLen += rec.sequence.size();

        std::string acgtDNA;
        acgtDNA.reserve(totalSeqLen);  // upper bound
        std::vector<size_t> allBoundaries;

        size_t acgtPos = 0;
        for (size_t r = 0; r < records.size(); r++) {
            const auto& seq = records[r].sequence;
            const auto& runs = allAmbiguousRuns[r];
            allBoundaries.push_back(acgtPos);  // record start = context reset

            size_t seqPos = 0;
            for (const auto& run : runs) {
                // Copy ACGT bases before this run
                for (size_t j = seqPos; j < run.position && j < seq.size(); j++) {
                    acgtDNA += seq[j];
                    acgtPos++;
                }
                seqPos = run.position + run.length;
                // Context reset after N-gap (if there are ACGT bases after)
                if (seqPos < seq.size()) {
                    allBoundaries.push_back(acgtPos);
                }
            }
            // Copy remaining ACGT bases after last run
            for (size_t j = seqPos; j < seq.size(); j++) {
                acgtDNA += seq[j];
                acgtPos++;
            }
        }

        uint64_t acgtDNALen = acgtDNA.size();
        auto encodedDNA = encodeDNAMultiRecord(acgtDNA, allBoundaries);

        BinaryWriter writer;
        writer.writeU32LE(MAGIC_V4);
        writer.writeU8(VERSION_V4_MAJOR);
        writer.writeU8(VERSION_V4_MINOR);
        writer.writeU8(static_cast<uint8_t>(maxOrder_));
        writer.writeU8(static_cast<uint8_t>(poolingOrder_));
        writer.writeU64LE(totalSeqLen);

        writer.writeU32LE(static_cast<uint32_t>(records.size()));

        writer.writeU32LE(static_cast<uint32_t>(identifierStream.size()));
        writer.writeBytes(identifierStream);
        writer.writeU32LE(static_cast<uint32_t>(seqLengthsData.size()));
        writer.writeBytes(seqLengthsData);
        writer.writeU32LE(static_cast<uint32_t>(ambiguousData.size()));
        writer.writeBytes(ambiguousData);
        writer.writeU64LE(acgtDNALen);  // V5: ACGT-only DNA length
        writer.writeU32LE(static_cast<uint32_t>(encodedDNA.size()));
        writer.writeBytes(encodedDNA);

        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());

        if (stats) {
            stats->originalBytes = totalSeqLen;
            stats->compressedBytes = writer.data().size();
            stats->numKmers = totalSeqLen;
            stats->numOrbits = records.size();
            stats->bitsPerBase = totalSeqLen > 0 ? (writer.data().size() * 8.0) / totalSeqLen : 0;
            stats->compressionRatio = writer.data().size() > 0 ?
                static_cast<double>(totalSeqLen) / writer.data().size() : 0;
        }

        return writer.data();
    }

    std::vector<FastaRecord> decompressToRecords(const std::vector<uint8_t>& data) {
        BinaryReader reader(data);
        uint32_t magic = reader.readU32LE();
        if (magic != MAGIC_V4) throw std::runtime_error("Invalid magic for V4");

        reader.readU8();  // versionMajor
        uint8_t vMinor = reader.readU8();

        int maxOrder = reader.readU8();
        int poolingOrder = reader.readU8();
        uint64_t totalSeqLen = reader.readU64LE();
        uint32_t recordCount = reader.readU32LE();
        int fastaModelVersion = (vMinor >= 8) ? 7 : (vMinor >= 7) ? 6 : (vMinor >= 6) ? 5 : 3;
        bool chunkedFastaDNA = (vMinor >= 9);

        size_t dataEnd = data.size() - 4;
        CRC32 crc;
        crc.update(data.data(), dataEnd);
        if (crc.finalize() != BinaryReader(data.data() + dataEnd, 4).readU32LE()) {
            throw std::runtime_error("CRC mismatch");
        }
        return decompressMultiRecordImpl(reader, maxOrder, poolingOrder,
                                          totalSeqLen, recordCount, fastaModelVersion,
                                          chunkedFastaDNA);
    }

    void decompressFile(const std::string& inputPath, const std::string& outputPath) {
        auto data = BinaryReader::readFile(inputPath);
        auto records = decompressToRecords(data);
        FastaWriter writer(outputPath);
        writer.writeAll(records);
    }

private:
    // Helper to append a little-endian u32 to a byte vector
    static void appendU32LE(std::vector<uint8_t>& out, uint32_t v) {
        out.push_back(v & 0xFF);
        out.push_back((v >> 8) & 0xFF);
        out.push_back((v >> 16) & 0xFF);
        out.push_back((v >> 24) & 0xFF);
    }

    // Helper to read a little-endian u32 from a pointer and advance it
    static uint32_t readU32LEPtr(const uint8_t*& p) {
        uint32_t v = p[0] | (static_cast<uint32_t>(p[1]) << 8)
                          | (static_cast<uint32_t>(p[2]) << 16)
                          | (static_cast<uint32_t>(p[3]) << 24);
        p += 4;
        return v;
    }

    // V9 chunked encoding: processes CHUNK_SIZE bases at a time to limit peak memory.
    // Each chunk is independently rANS-encoded; model persists across chunks.
    // Format: numChunks(u32) + for each chunk: baseCount(u32) + encodedSize(u32) + data
    // Peak memory: ~96 MB per chunk (4M bases × 24 bytes for symbols+CDFs) instead of
    // O(N) for entire genome (which would be 72 GB for a 3B-base human genome).
    std::vector<uint8_t> encodeDNAMultiRecord(const std::string& allDNA,
                                               const std::vector<size_t>& recordBoundaries) {
        if (allDNA.empty()) return {};

        static constexpr size_t FASTA_DNA_CHUNK_SIZE = 4 * 1024 * 1024;  // 4M bases

        OnlineContextModel model(maxOrder_, poolingOrder_, true, 7);
        model.setPersistSubs(true);
        uint64_t fullMask = (1ULL << (2 * maxOrder_)) - 1;

        std::array<uint32_t, 5> uniformCDF = {
            0, PROB_SCALE/4, PROB_SCALE/2, 3*PROB_SCALE/4, PROB_SCALE
        };

        uint32_t numChunks = static_cast<uint32_t>(
            (allDNA.size() + FASTA_DNA_CHUNK_SIZE - 1) / FASTA_DNA_CHUNK_SIZE);
        if (numChunks == 0) numChunks = 1;

        std::vector<uint8_t> output;
        appendU32LE(output, numChunks);

        uint64_t ctx = 0;
        size_t nextBoundary = 0;
        int basesInRecord = 0;

        for (size_t chunkStart = 0; chunkStart < allDNA.size(); chunkStart += FASTA_DNA_CHUNK_SIZE) {
            size_t chunkEnd = std::min(chunkStart + FASTA_DNA_CHUNK_SIZE, allDNA.size());
            size_t chunkLen = chunkEnd - chunkStart;

            std::vector<int> symbols;
            std::vector<std::array<uint32_t, 5>> cdfs;
            symbols.reserve(chunkLen);
            cdfs.reserve(chunkLen);

            for (size_t i = chunkStart; i < chunkEnd; i++) {
                while (nextBoundary < recordBoundaries.size() && i == recordBoundaries[nextBoundary]) {
                    ctx = 0;
                    basesInRecord = 0;
                    model.resetPerReadState();
                    nextBoundary++;
                }

                int base = charToBase(allDNA[i]);
                symbols.push_back(base);

                if (basesInRecord == 0) {
                    cdfs.push_back(uniformCDF);
                } else {
                    int validOrder = std::min(basesInRecord, maxOrder_);
                    cdfs.push_back(model.getCDF(ctx, validOrder));
                }

                if (basesInRecord > 0) {
                    model.update(ctx, base, std::min(basesInRecord - 1, maxOrder_));
                }

                ctx = ((ctx << 2) | base) & fullMask;
                basesInRecord++;
            }

            auto encoded = InterleavedRANSEncoder::encodeAll(symbols, cdfs);
            appendU32LE(output, static_cast<uint32_t>(chunkLen));
            appendU32LE(output, static_cast<uint32_t>(encoded.size()));
            output.insert(output.end(), encoded.begin(), encoded.end());
        }

        return output;
    }

    std::vector<FastaRecord> decompressMultiRecordImpl(BinaryReader& reader,
                                                        int maxOrder, int poolingOrder,
                                                        uint64_t /*totalSeqLen*/,
                                                        uint32_t recordCount,
                                                        int fastaModelVersion = 6,
                                                        bool chunkedDNA = false) {
        uint32_t identifierStreamSize = reader.readU32LE();
        auto identifierData = reader.readBytes(identifierStreamSize);
        uint32_t seqLengthsSize = reader.readU32LE();
        auto seqLengthsData = reader.readBytes(seqLengthsSize);
        uint32_t ambiguousDataSize = reader.readU32LE();
        auto ambiguousDataBytes = reader.readBytes(ambiguousDataSize);

        uint64_t acgtDNALen = reader.readU64LE();

        uint32_t encodedDNASize = reader.readU32LE();
        auto encodedDNA = reader.readBytes(encodedDNASize);

        IdentifierDecoderV2 idDecoder;
        auto headers = idDecoder.decode(identifierData.data(), identifierData.size());

        // Delta-decode sequence lengths
        std::vector<size_t> seqLengths;
        seqLengths.reserve(recordCount);
        {
            const uint8_t* p = seqLengthsData.data();
            int64_t prevLen = 0;
            for (uint32_t i = 0; i < recordCount; i++) {
                int64_t delta = readSignedVarint(p);
                prevLen += delta;
                seqLengths.push_back(static_cast<size_t>(prevLen));
            }
        }

        // Decode RLE ambiguous runs per record
        std::vector<std::vector<AmbiguousRun>> perRecordAmbiguousRuns(recordCount);
        {
            const uint8_t* p = ambiguousDataBytes.data();
            const uint8_t* end = p + ambiguousDataBytes.size();
            for (uint32_t r = 0; r < recordCount && p < end; r++) {
                uint64_t numRuns = readVarint(p);
                uint64_t prevEnd = 0;
                for (uint64_t j = 0; j < numRuns && p < end; j++) {
                    uint64_t startDelta = readVarint(p);
                    uint64_t length = readVarint(p);
                    char c = static_cast<char>(*p++);
                    uint64_t pos = prevEnd + startDelta;
                    perRecordAmbiguousRuns[r].push_back({pos, length, c});
                    prevEnd = pos + length;
                }
            }
        }

        // Compute boundaries: record starts + after each N-gap (in ACGT-space)
        std::vector<size_t> allBoundaries;
        {
            size_t acgtPos = 0;
            for (uint32_t r = 0; r < recordCount; r++) {
                allBoundaries.push_back(acgtPos);
                const auto& runs = perRecordAmbiguousRuns[r];
                size_t seqPos = 0;
                for (const auto& run : runs) {
                    size_t acgtBefore = (run.position > seqPos) ? (run.position - seqPos) : 0;
                    acgtPos += acgtBefore;
                    seqPos = run.position + run.length;
                    if (seqPos < seqLengths[r]) {
                        allBoundaries.push_back(acgtPos);
                    }
                }
                size_t acgtAfter = (seqLengths[r] > seqPos) ? (seqLengths[r] - seqPos) : 0;
                acgtPos += acgtAfter;
            }
        }

        // Decode DNA with adaptive context model
        std::string decodedDNA;
        decodedDNA.resize(acgtDNALen);

        OnlineContextModel model(maxOrder, poolingOrder, true, fastaModelVersion);
        model.setPersistSubs(true);

        uint64_t fullMask = (1ULL << (2 * maxOrder)) - 1;
        uint64_t ctx = 0;

        std::array<uint32_t, 5> uniformCDF = {
            0, PROB_SCALE/4, PROB_SCALE/2, 3*PROB_SCALE/4, PROB_SCALE
        };

        size_t nextBoundary = 0;
        int basesInRecord = 0;

        if (chunkedDNA) {
            // V9+ chunked DNA decoding: parse chunk structure
            const uint8_t* p = encodedDNA.data();
            const uint8_t* end = p + encodedDNA.size();

            uint32_t numChunks = readU32LEPtr(p);
            size_t globalPos = 0;

            for (uint32_t c = 0; c < numChunks && p + 8 <= end; c++) {
                uint32_t chunkBaseCount = readU32LEPtr(p);
                uint32_t chunkEncodedSize = readU32LEPtr(p);
                if (p + chunkEncodedSize > end) {
                    throw std::runtime_error("Chunked DNA: encoded size exceeds data");
                }

                size_t chunkGlobalStart = globalPos;

                InterleavedRANSDecoder::decodeIncremental(
                    p, chunkEncodedSize, chunkBaseCount,
                    [&](size_t localI) -> std::array<uint32_t, 5> {
                        size_t gi = chunkGlobalStart + localI;
                        while (nextBoundary < allBoundaries.size() && gi == allBoundaries[nextBoundary]) {
                            ctx = 0;
                            basesInRecord = 0;
                            model.resetPerReadState();
                            nextBoundary++;
                        }
                        if (basesInRecord == 0) return uniformCDF;
                        int validOrder = std::min(basesInRecord, maxOrder);
                        return model.getCDF(ctx, validOrder);
                    },
                    [&](size_t localI, int base) {
                        size_t gi = chunkGlobalStart + localI;
                        decodedDNA[gi] = baseToChar(base);
                        if (basesInRecord > 0) {
                            model.update(ctx, base, std::min(basesInRecord - 1, maxOrder));
                        }
                        ctx = ((ctx << 2) | base) & fullMask;
                        basesInRecord++;
                    }
                );

                p += chunkEncodedSize;
                globalPos += chunkBaseCount;
            }
        } else {
            // V8 and earlier: single rANS block
            InterleavedRANSDecoder::decodeIncremental(
                encodedDNA.data(), encodedDNA.size(), acgtDNALen,
                [&](size_t i) -> std::array<uint32_t, 5> {
                    while (nextBoundary < allBoundaries.size() && i == allBoundaries[nextBoundary]) {
                        ctx = 0;
                        basesInRecord = 0;
                        model.resetPerReadState();
                        nextBoundary++;
                    }
                    if (basesInRecord == 0) return uniformCDF;
                    int validOrder = std::min(basesInRecord, maxOrder);
                    return model.getCDF(ctx, validOrder);
                },
                [&](size_t i, int base) {
                    decodedDNA[i] = baseToChar(base);
                    if (basesInRecord > 0) {
                        model.update(ctx, base, std::min(basesInRecord - 1, maxOrder));
                    }
                    ctx = ((ctx << 2) | base) & fullMask;
                    basesInRecord++;
                }
            );
        }

        // Reconstruct full sequences by interleaving ACGT stream with ambiguous runs
        std::vector<FastaRecord> result;
        result.reserve(recordCount);

        size_t acgtPos = 0;
        for (uint32_t r = 0; r < recordCount; r++) {
            FastaRecord rec;
            rec.header = r < headers.size() ? headers[r] : "";
            rec.sequence.resize(seqLengths[r]);
            const auto& runs = perRecordAmbiguousRuns[r];
            size_t seqPos = 0;
            for (const auto& run : runs) {
                while (seqPos < run.position && seqPos < seqLengths[r]) {
                    rec.sequence[seqPos++] = decodedDNA[acgtPos++];
                }
                for (uint64_t k = 0; k < run.length && seqPos < seqLengths[r]; k++) {
                    rec.sequence[seqPos++] = run.character;
                }
            }
            while (seqPos < seqLengths[r]) {
                rec.sequence[seqPos++] = decodedDNA[acgtPos++];
            }
            result.push_back(std::move(rec));
        }

        return result;
    }

    static int charToBase(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
        }
        return 0;
    }

    static char baseToChar(int b) { return "ACGT"[b & 3]; }

    static char toUpper(char c) {
        if (c >= 'a' && c <= 'z') return c - 32;
        return c;
    }

    static bool isValidBase(char c) {
        return c == 'A' || c == 'C' || c == 'G' || c == 'T' ||
               c == 'a' || c == 'c' || c == 'g' || c == 't';
    }

    std::vector<uint8_t> createEmptyOutput(const std::string& header, CompressionStats* stats) {
        BinaryWriter writer;
        writer.writeU32LE(MAGIC_V4);
        writer.writeU8(VERSION_V4_MAJOR);
        writer.writeU8(VERSION_V4_MINOR);
        writer.writeU8(maxOrder_);
        writer.writeU8(poolingOrder_);
        writer.writeU64LE(0);
        writer.writeU32LE(static_cast<uint32_t>(header.size()));
        if (!header.empty()) {
            writer.writeBytes(std::vector<uint8_t>(header.begin(), header.end()));
        }
        writer.writeU32LE(0);
        writer.writeU32LE(0);
        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());
        if (stats) *stats = {0, writer.data().size(), 0, 0, 0, 0};
        return writer.data();
    }

    std::vector<uint8_t> compressShortSequence(const std::string& header,
                                                const std::string& rawSeq,
                                                const std::vector<NonDNARun>& nonDNARuns,
                                                CompressionStats* stats) {
        BinaryWriter writer;
        writer.writeU32LE(MAGIC_V4);
        writer.writeU8(VERSION_V4_MAJOR);
        writer.writeU8(VERSION_V4_MINOR);
        writer.writeU8(maxOrder_);
        writer.writeU8(poolingOrder_);
        writer.writeU64LE(rawSeq.size());

        writer.writeU32LE(static_cast<uint32_t>(header.size()));
        if (!header.empty()) {
            writer.writeBytes(std::vector<uint8_t>(header.begin(), header.end()));
        }

        writer.writeU32LE(static_cast<uint32_t>(nonDNARuns.size()));
        if (!nonDNARuns.empty()) {
            std::vector<uint8_t> nonDNAData;
            uint64_t prevPos = 0;
            for (const auto& run : nonDNARuns) {
                writeVarint(nonDNAData, run.position - prevPos);
                writeVarint(nonDNAData, run.length);
                nonDNAData.push_back(static_cast<uint8_t>(run.character));
                prevPos = run.position + run.length;
            }
            writer.writeU32LE(static_cast<uint32_t>(nonDNAData.size()));
            writer.writeBytes(nonDNAData);
        }

        std::vector<uint8_t> raw;
        uint8_t byte = 0;
        int bits = 0;
        for (char c : rawSeq) {
            int base = isValidBase(c) ? charToBase(c) : 0;
            byte = (byte << 2) | base;
            bits += 2;
            if (bits == 8) { raw.push_back(byte); byte = 0; bits = 0; }
        }
        if (bits > 0) { raw.push_back(byte << (8 - bits)); }

        writer.writeU32LE(static_cast<uint32_t>(raw.size()));
        writer.writeBytes(raw);

        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());

        if (stats) {
            stats->originalBytes = rawSeq.size();
            stats->compressedBytes = writer.data().size();
            stats->bitsPerBase = rawSeq.size() > 0 ? (writer.data().size() * 8.0) / rawSeq.size() : 0;
            stats->compressionRatio = rawSeq.size() > 0 ? (double)rawSeq.size() / writer.data().size() : 0;
        }
        return writer.data();
    }

    int maxOrder_;
    int poolingOrder_;
};

} // namespace v4zip

#endif // V4ZIP_COMPRESSOR_V4_HPP
