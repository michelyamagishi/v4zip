/**
 * ArithmeticCoder.hpp - rANS (range Asymmetric Numeral Systems) encoder/decoder
 *
 * Parameters:
 *   PROB_BITS = 14 (probability precision)
 *   PROB_SCALE = 16384
 *   RANS_L = 2^23 (normalization threshold)
 */

#ifndef V4ZIP_ARITHMETICCODER_HPP
#define V4ZIP_ARITHMETICCODER_HPP

#include "ContextModel.hpp"  // For PROB_BITS, PROB_SCALE
#include <vector>
#include <array>
#include <cstdint>
#include <cassert>
#include <algorithm>


namespace v4zip {

// RANS_L defined in ContextModel.hpp (included above)

/**
 * Binary search a CDF array to find symbol where cdf[sym] <= slot < cdf[sym+1].
 * Works with any CDF container supporting operator[] (array, vector, etc.).
 */
template<typename CDF>
inline int searchCDF(const CDF& cdf, int alphabetSize, uint32_t slot) {
    int lo = 0, hi = alphabetSize;
    while (lo < hi) {
        int mid = (lo + hi) / 2;
        if (cdf[mid + 1] <= slot) {
            lo = mid + 1;
        } else {
            hi = mid;
        }
    }
    return lo;
}

class RANSEncoder {
public:
    RANSEncoder() : state_(RANS_L) {}

    // Encode a symbol given its CDF (reverse order for rANS)
    void encode(uint8_t symbol, const std::array<uint32_t, 5>& cdf) {
        uint32_t start = cdf[symbol];
        uint32_t freq = cdf[symbol + 1] - start;
        assert(freq > 0 && "CDF must not contain zero-frequency symbols");

        uint32_t limit = ((RANS_L >> PROB_BITS) << 8) * freq;
        while (state_ >= limit) {
            output_.push_back(static_cast<uint8_t>(state_ & 0xFF));
            state_ >>= 8;
        }

        state_ = ((state_ / freq) << PROB_BITS) + (state_ % freq) + start;
        assert(state_ >= RANS_L && "rANS state underflow after encode");
    }

    // Finalize and return encoded bytes (reversed)
    std::vector<uint8_t> finish() {
        // Flush final state
        for (int i = 0; i < 4; i++) {
            output_.push_back(static_cast<uint8_t>(state_ & 0xFF));
            state_ >>= 8;
        }

        std::reverse(output_.begin(), output_.end());
        return std::move(output_);
    }

    // Encode multiple symbols with their CDFs (in reverse order internally)
    static std::vector<uint8_t> encodeAll(
        const std::vector<int>& symbols,
        const std::vector<std::array<uint32_t, 5>>& cdfs
    ) {
        if (symbols.empty()) {
            // Return valid initial state (RANS_L in big-endian)
            return {static_cast<uint8_t>((RANS_L >> 24) & 0xFF),
                    static_cast<uint8_t>((RANS_L >> 16) & 0xFF),
                    static_cast<uint8_t>((RANS_L >> 8) & 0xFF),
                    static_cast<uint8_t>(RANS_L & 0xFF)};
        }

        RANSEncoder encoder;
        // Encode in reverse order
        for (int i = static_cast<int>(symbols.size()) - 1; i >= 0; i--) {
            encoder.encode(static_cast<uint8_t>(symbols[i]), cdfs[i]);
        }
        return encoder.finish();
    }

private:
    uint32_t state_;
    std::vector<uint8_t> output_;
};

class RANSDecoder {
public:
    RANSDecoder(const uint8_t* data, size_t len)
        : data_(data), len_(len), pos_(0), state_(0) {
        // Initialize state from first 4 bytes
        for (int i = 0; i < 4 && pos_ < len_; i++) {
            state_ = (state_ << 8) | data_[pos_++];
        }
    }

    // Decode a symbol given its CDF
    uint8_t decode(const std::array<uint32_t, 5>& cdf) {
        uint32_t slot = state_ & (PROB_SCALE - 1);

        // Binary search for symbol
        uint8_t symbol = 0;
        while (symbol < 4 && cdf[symbol + 1] <= slot) {
            symbol++;
        }

        uint32_t start = cdf[symbol];
        uint32_t freq = cdf[symbol + 1] - start;
        assert(freq > 0 && "CDF must not contain zero-frequency symbols");

        state_ = freq * (state_ >> PROB_BITS) + slot - start;

        // Renormalize
        while (state_ < RANS_L && pos_ < len_) {
            state_ = (state_ << 8) | data_[pos_++];
        }

        return symbol;
    }

    // Decode with incremental CDF generation (template to avoid std::function overhead)
    template<typename CDFGetter, typename Emitter>
    static void decodeIncremental(
        const uint8_t* data, size_t len, size_t count,
        CDFGetter getCDF,
        Emitter emit
    ) {
        if (count == 0) return;

        RANSDecoder decoder(data, len);
        for (size_t i = 0; i < count; i++) {
            auto cdf = getCDF(i);
            int symbol = decoder.decode(cdf);
            emit(i, symbol);
        }
    }

private:
    const uint8_t* data_;
    size_t len_;
    size_t pos_;
    uint32_t state_;
};

// ============================================================================
// Generic rANS encoder/decoder for arbitrary alphabet sizes
// ============================================================================

/**
 * GenericRANSEncoder - Template-based rANS encoder for any alphabet size
 *
 * @tparam AlphabetSize The number of symbols in the alphabet (e.g., 94 for quality scores)
 */
template<size_t AlphabetSize>
class GenericRANSEncoder {
public:
    GenericRANSEncoder() : state_(RANS_L) {}

    // Encode a symbol given its CDF (CDF has AlphabetSize+1 entries: [0, cumsum..., PROB_SCALE])
    void encode(size_t symbol, const std::array<uint32_t, AlphabetSize + 1>& cdf) {
        if (symbol >= AlphabetSize) {
            symbol = AlphabetSize - 1;  // Clamp to valid range
        }

        uint32_t start = cdf[symbol];
        uint32_t freq = cdf[symbol + 1] - start;
        assert(freq > 0 && "CDF must not contain zero-frequency symbols");

        uint32_t limit = ((RANS_L >> PROB_BITS) << 8) * freq;
        while (state_ >= limit) {
            output_.push_back(static_cast<uint8_t>(state_ & 0xFF));
            state_ >>= 8;
        }

        state_ = ((state_ / freq) << PROB_BITS) + (state_ % freq) + start;
    }

    // Finalize and return encoded bytes (reversed)
    std::vector<uint8_t> finish() {
        // Flush final state
        for (int i = 0; i < 4; i++) {
            output_.push_back(static_cast<uint8_t>(state_ & 0xFF));
            state_ >>= 8;
        }

        std::reverse(output_.begin(), output_.end());
        return std::move(output_);
    }

    // Encode multiple symbols with their CDFs (in reverse order internally)
    static std::vector<uint8_t> encodeAll(
        const std::vector<size_t>& symbols,
        const std::vector<std::array<uint32_t, AlphabetSize + 1>>& cdfs
    ) {
        if (symbols.empty()) {
            // Return valid initial state (RANS_L in big-endian)
            return {static_cast<uint8_t>((RANS_L >> 24) & 0xFF),
                    static_cast<uint8_t>((RANS_L >> 16) & 0xFF),
                    static_cast<uint8_t>((RANS_L >> 8) & 0xFF),
                    static_cast<uint8_t>(RANS_L & 0xFF)};
        }

        GenericRANSEncoder<AlphabetSize> encoder;
        // Encode in reverse order for rANS
        for (int i = static_cast<int>(symbols.size()) - 1; i >= 0; i--) {
            encoder.encode(symbols[i], cdfs[i]);
        }
        return encoder.finish();
    }

private:
    uint32_t state_;
    std::vector<uint8_t> output_;
};

/**
 * GenericRANSDecoder - Template-based rANS decoder for any alphabet size
 *
 * Uses binary search for efficient symbol lookup in large alphabets.
 *
 * @tparam AlphabetSize The number of symbols in the alphabet
 */
template<size_t AlphabetSize>
class GenericRANSDecoder {
public:
    GenericRANSDecoder(const uint8_t* data, size_t len)
        : data_(data), len_(len), pos_(0), state_(0) {
        // Initialize state from first 4 bytes
        for (int i = 0; i < 4 && pos_ < len_; i++) {
            state_ = (state_ << 8) | data_[pos_++];
        }
    }

    // Decode a symbol given its CDF using binary search
    size_t decode(const std::array<uint32_t, AlphabetSize + 1>& cdf) {
        uint32_t slot = state_ & (PROB_SCALE - 1);
        size_t symbol = static_cast<size_t>(searchCDF(cdf, static_cast<int>(AlphabetSize), slot));

        uint32_t start = cdf[symbol];
        uint32_t freq = cdf[symbol + 1] - start;
        assert(freq > 0 && "CDF must not contain zero-frequency symbols");

        state_ = freq * (state_ >> PROB_BITS) + slot - start;

        // Renormalize
        while (state_ < RANS_L && pos_ < len_) {
            state_ = (state_ << 8) | data_[pos_++];
        }

        return symbol;
    }

    // Decode with incremental CDF generation (template to avoid std::function overhead)
    template<typename CDFGetter, typename Emitter>
    static void decodeIncremental(
        const uint8_t* data, size_t len, size_t count,
        CDFGetter getCDF,
        Emitter emit
    ) {
        if (count == 0) return;

        GenericRANSDecoder<AlphabetSize> decoder(data, len);
        for (size_t i = 0; i < count; i++) {
            auto cdf = getCDF(i);
            size_t symbol = decoder.decode(cdf);
            emit(i, symbol);
        }
    }

private:
    const uint8_t* data_;
    size_t len_;
    size_t pos_;
    uint32_t state_;
};

// Convenience aliases for common alphabet sizes
using QualityRANSEncoder = GenericRANSEncoder<94>;   // Phred+33 quality (33-126)
using QualityRANSDecoder = GenericRANSDecoder<94>;

} // namespace v4zip

#endif // V4ZIP_ARITHMETICCODER_HPP
