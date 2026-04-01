/**
 * IdentifierCodec.hpp - Tokenization-based identifier compression
 *
 * FASTQ read identifiers often have regular patterns:
 *   @HWUSI-EAS100R:6:73:941:1973#0/1
 *   @SRR001666.1 071112_SLXA-EAS1_s_7:5:1:817:345 length=72
 *
 * Strategy:
 *   1. Tokenize on delimiters (. : _ / # space)
 *   2. Build dictionary of frequent tokens
 *   3. Delta-encode sequential numeric tokens
 *   4. Byte-level fallback for novel tokens
 *
 * Expected: 0.1-0.5 bits/char vs 5-6 for raw ASCII
 */

#ifndef V4ZIP_IDENTIFIERCODEC_HPP
#define V4ZIP_IDENTIFIERCODEC_HPP

#include "Varint.hpp"
#include "ContextModel.hpp"  // For PROB_BITS, PROB_SCALE
#include <string>
#include <vector>
#include <algorithm>
#include <array>
#include <unordered_map>
#include <cstdint>
#include <cassert>

namespace v4zip {

// Token types
enum class TokenType : uint8_t {
    LITERAL = 0,      // Raw byte sequence (fallback)
    DICT_REF = 1,     // Reference to dictionary entry
    NUMBER = 2,       // Numeric value (varint encoded)
    DELTA_NUM = 3,    // Delta-encoded number (from previous)
    DELIMITER = 4,    // Single delimiter character
    REPEAT_PREV = 5,  // Repeat previous token
};

// Common delimiters in FASTQ identifiers
inline bool isDelimiter(char c) {
    return c == '.' || c == ':' || c == '_' || c == '/' ||
           c == '#' || c == ' ' || c == '\t' || c == '@' ||
           c == '-' || c == '|';
}

// Check if string is a number
inline bool isNumeric(const std::string& s) {
    if (s.empty()) return false;
    for (char c : s) {
        if (c < '0' || c > '9') return false;
    }
    return true;
}

// Maximum digits that safely fit in uint64_t without overflow checks
constexpr int MAX_SAFE_NUMBER_DIGITS = 18;  // 999999999999999999 < 2^64

// Parse number from string (returns 0 if not numeric or overflow)
inline uint64_t parseNumber(const std::string& s) {
    if (s.size() > MAX_SAFE_NUMBER_DIGITS) return 0;  // Overflow: treat as literal
    uint64_t val = 0;
    for (char c : s) {
        if (c < '0' || c > '9') return 0;
        val = val * 10 + (c - '0');
    }
    return val;
}

// Format number as string
inline std::string formatNumber(uint64_t val) {
    if (val == 0) return "0";
    std::string s;
    while (val > 0) {
        s.push_back('0' + (val % 10));
        val /= 10;
    }
    std::reverse(s.begin(), s.end());
    return s;
}

// Format number with specified width (preserving leading zeros)
inline std::string formatNumberWithWidth(uint64_t val, uint8_t width) {
    std::string s = formatNumber(val);
    if (width > 0 && s.size() < width) {
        s = std::string(width - s.size(), '0') + s;
    }
    return s;
}

/**
 * Token representation
 */
struct Token {
    TokenType type;
    std::string value;      // For LITERAL, DICT_REF key
    uint64_t numValue;      // For NUMBER, DELTA_NUM, DICT_REF index
    uint8_t numWidth;       // For NUMBER: original digit count (to preserve leading zeros)
    char delimiter;         // For DELIMITER

    Token() : type(TokenType::LITERAL), numValue(0), numWidth(0), delimiter(0) {}
    Token(TokenType t, const std::string& v) : type(t), value(v), numValue(0), numWidth(0), delimiter(0) {}
    Token(TokenType t, uint64_t n, uint8_t w = 0) : type(t), numValue(n), numWidth(w), delimiter(0) {}
    Token(char d) : type(TokenType::DELIMITER), numValue(0), numWidth(0), delimiter(d) {}
};

/**
 * Tokenize an identifier string
 */
inline std::vector<Token> tokenize(const std::string& identifier) {
    std::vector<Token> tokens;
    std::string current;

    for (size_t i = 0; i < identifier.size(); i++) {
        char c = identifier[i];
        if (isDelimiter(c)) {
            if (!current.empty()) {
                if (isNumeric(current) && current.size() <= MAX_SAFE_NUMBER_DIGITS) {
                    // Store width to preserve leading zeros
                    uint8_t width = static_cast<uint8_t>(std::min(current.size(), size_t(255)));
                    tokens.push_back(Token(TokenType::NUMBER, parseNumber(current), width));
                } else {
                    tokens.push_back(Token(TokenType::LITERAL, current));
                }
                current.clear();
            }
            tokens.push_back(Token(c));
        } else {
            current.push_back(c);
        }
    }

    // Don't forget the last token
    if (!current.empty()) {
        if (isNumeric(current) && current.size() <= MAX_SAFE_NUMBER_DIGITS) {
            uint8_t width = static_cast<uint8_t>(std::min(current.size(), size_t(255)));
            tokens.push_back(Token(TokenType::NUMBER, parseNumber(current), width));
        } else {
            tokens.push_back(Token(TokenType::LITERAL, current));
        }
    }

    return tokens;
}

/**
 * Reconstruct identifier from tokens
 */
inline std::string detokenize(const std::vector<Token>& tokens) {
    std::string result;
    for (const auto& tok : tokens) {
        switch (tok.type) {
            case TokenType::LITERAL:
            case TokenType::DICT_REF:
                result += tok.value;
                break;
            case TokenType::NUMBER:
            case TokenType::DELTA_NUM:
                result += formatNumberWithWidth(tok.numValue, tok.numWidth);
                break;
            case TokenType::DELIMITER:
                result += tok.delimiter;
                break;
            case TokenType::REPEAT_PREV:
                // Handled during decoding
                break;
        }
    }
    return result;
}

// ============================================================================
// Entropy-coded identifier codec using rANS
// ============================================================================

// Number of token types for entropy coding
constexpr int NUM_TOKEN_TYPES = 6;  // LITERAL, DICT_REF, NUMBER, DELTA_NUM, DELIMITER, REPEAT_PREV

// Number of delimiter characters
constexpr int NUM_DELIMITERS = 11;  // . : _ / # space tab @ - |

// Map delimiter character to index (0-10)
inline int delimiterToIndex(char c) {
    switch (c) {
        case '.': return 0;
        case ':': return 1;
        case '_': return 2;
        case '/': return 3;
        case '#': return 4;
        case ' ': return 5;
        case '\t': return 6;
        case '@': return 7;
        case '-': return 8;
        case '|': return 9;
        default: return 10;  // Other (shouldn't happen)
    }
}

// Map index back to delimiter character (index 10 = "other", requires raw char from stream)
inline char indexToDelimiter(int idx) {
    const char delims[] = ".:_/# \t@-|";
    return (idx >= 0 && idx < 10) ? delims[idx] : '?';  // '?' placeholder; callers must handle idx==10
}

/**
 * Adaptive model for entropy coding small alphabets (token types, delimiters, dict refs)
 */
template<int AlphabetSize>
class AdaptiveSymbolModel {
public:
    AdaptiveSymbolModel() : activeSize_(AlphabetSize) {
        reset();
    }

    void reset() {
        for (int i = 0; i < activeSize_; i++) counts_[i] = 1;
        for (int i = activeSize_; i < AlphabetSize; i++) counts_[i] = 0;
        total_ = activeSize_;
    }

    void setActiveSize(int size) {
        activeSize_ = std::min(size, AlphabetSize);
        reset();
    }

    std::vector<uint32_t> getCDF() const {
        return normalizeCDFVec(counts_.data(), activeSize_, PROB_SCALE);
    }

    void update(int symbol) {
        if (symbol >= 0 && symbol < activeSize_) {
            counts_[symbol]++;
            total_++;
            if (total_ > 65535) {  // COUNT_RESCALE_THRESHOLD
                rescale();
            }
        }
    }

    int activeSize() const { return activeSize_; }

private:
    void rescale() {
        rescaleCounts(counts_, activeSize_, total_);
    }

    std::array<uint32_t, AlphabetSize> counts_;
    uint32_t total_;
    int activeSize_;
};

/**
 * rANS encoder for identifier streams
 */
class IdentifierRANSEncoder {
public:
    static std::vector<uint8_t> encode(
        const std::vector<int>& symbols,
        const std::vector<std::vector<uint32_t>>& cdfs,
        int /*alphabetSize*/) {

        if (symbols.empty()) return {};

        uint32_t state = RANS_L;
        std::vector<uint8_t> output;

        for (int i = static_cast<int>(symbols.size()) - 1; i >= 0; i--) {
            int sym = symbols[i];
            const auto& cdf = cdfs[i];

            uint32_t start = cdf[sym];
            uint32_t freq = cdf[sym + 1] - start;
            assert(freq > 0 && "CDF must not contain zero-frequency symbols");

            uint32_t max_state = ((RANS_L / PROB_SCALE) * freq) << 8;
            while (state >= max_state) {
                output.push_back(state & 0xFF);
                state >>= 8;
            }

            state = ((state / freq) << PROB_BITS) + (state % freq) + start;
        }

        // Append final state (LSB first; reverse will put MSB first)
        output.push_back((state >> 0) & 0xFF);
        output.push_back((state >> 8) & 0xFF);
        output.push_back((state >> 16) & 0xFF);
        output.push_back((state >> 24) & 0xFF);

        std::reverse(output.begin(), output.end());
        return output;
    }
};

/**
 * rANS decoder for identifier streams
 */
class IdentifierRANSDecoder {
public:

    template<typename CDFGetter, typename SymbolHandler>
    static void decode(
        const uint8_t* data, size_t len,
        size_t count, int alphabetSize,
        CDFGetter getCDF, SymbolHandler handleSymbol) {

        if (count == 0 || len < 4) return;

        const uint8_t* ptr = data;
        const uint8_t* end = data + len;

        uint32_t state = (static_cast<uint32_t>(ptr[0]) << 24) |
                         (static_cast<uint32_t>(ptr[1]) << 16) |
                         (static_cast<uint32_t>(ptr[2]) << 8) |
                         static_cast<uint32_t>(ptr[3]);
        ptr += 4;

        for (size_t i = 0; i < count; i++) {
            auto cdf = getCDF(i);

            uint32_t slot = state & (PROB_SCALE - 1);
            int sym = searchCDF(cdf, alphabetSize, slot);

            handleSymbol(i, sym);

            uint32_t start = cdf[sym];
            uint32_t freq = cdf[sym + 1] - start;
            assert(freq > 0 && "CDF must not contain zero-frequency symbols");
            state = freq * (state >> PROB_BITS) + slot - start;

            while (state < RANS_L && ptr < end) {
                state = (state << 8) | *ptr++;
            }
        }
    }
};

/**
 * IdentifierEncoderV2 - Entropy-coded identifier encoder using rANS
 *
 * Improvements over V1:
 * - Token types are entropy coded (6 symbols, highly skewed distribution)
 * - Dictionary references are entropy coded (up to 256 symbols)
 * - Delimiter characters are entropy coded (11 possible delimiters)
 *
 * Expected: 10-20% better compression than V1
 */
class IdentifierEncoderV2 {
public:
    IdentifierEncoderV2() : prevNumber_(0), prevWidth_(0) {}

    std::vector<uint8_t> encode(const std::vector<std::string>& identifiers) {
        if (identifiers.empty()) {
            return {};
        }

        // First pass: collect token statistics and build dictionary
        std::unordered_map<std::string, uint32_t> tokenCounts;
        for (const auto& id : identifiers) {
            auto tokens = tokenize(id);
            for (const auto& tok : tokens) {
                if (tok.type == TokenType::LITERAL) {
                    tokenCounts[tok.value]++;
                }
            }
        }

        // Build dictionary from frequent tokens
        std::vector<std::pair<std::string, uint32_t>> sortedTokens(
            tokenCounts.begin(), tokenCounts.end());
        std::sort(sortedTokens.begin(), sortedTokens.end(),
            [](const auto& a, const auto& b) {
                return a.second * a.first.size() > b.second * b.first.size();
            });

        dictionary_.clear();
        dictIndex_.clear();
        for (size_t i = 0; i < sortedTokens.size() && dictionary_.size() < 256; i++) {
            if (sortedTokens[i].second >= 2) {
                dictIndex_[sortedTokens[i].first] = dictionary_.size();
                dictionary_.push_back(sortedTokens[i].first);
            }
        }

        // Second pass: collect all tokens with their types
        std::vector<int> tokenTypeSymbols;
        std::vector<std::vector<uint32_t>> tokenTypeCDFs;
        std::vector<int> delimiterSymbols;
        std::vector<std::vector<uint32_t>> delimiterCDFs;
        std::vector<int> dictRefSymbols;
        std::vector<std::vector<uint32_t>> dictRefCDFs;
        std::vector<uint8_t> rawData;  // Numbers, literals, widths

        AdaptiveSymbolModel<NUM_TOKEN_TYPES> tokenTypeModel;
        AdaptiveSymbolModel<NUM_DELIMITERS> delimiterModel;
        AdaptiveSymbolModel<256> dictRefModel;
        int dictSize = std::max(1, static_cast<int>(dictionary_.size()));
        dictRefModel.setActiveSize(dictSize);

        prevNumber_ = 0;
        prevWidth_ = 0;

        for (const auto& id : identifiers) {
            auto tokens = tokenize(id);
            writeVarint(rawData, tokens.size());  // Token count

            for (const auto& tok : tokens) {
                int tokenType;
                if (tok.type == TokenType::DELIMITER) {
                    tokenType = static_cast<int>(TokenType::DELIMITER);
                    tokenTypeCDFs.push_back(tokenTypeModel.getCDF());
                    tokenTypeSymbols.push_back(tokenType);
                    tokenTypeModel.update(tokenType);

                    // Entropy code the delimiter
                    int delimIdx = delimiterToIndex(tok.delimiter);
                    delimiterCDFs.push_back(delimiterModel.getCDF());
                    delimiterSymbols.push_back(delimIdx);
                    delimiterModel.update(delimIdx);

                    // Store raw char for unknown delimiters (index 10 = "other")
                    if (delimIdx == 10) {
                        rawData.push_back(static_cast<uint8_t>(tok.delimiter));
                    }
                } else if (tok.type == TokenType::NUMBER) {
                    int64_t delta = static_cast<int64_t>(tok.numValue) - static_cast<int64_t>(prevNumber_);
                    if (delta >= -128 && delta <= 127 && prevNumber_ != 0 && tok.numWidth == prevWidth_) {
                        tokenType = static_cast<int>(TokenType::DELTA_NUM);
                        tokenTypeCDFs.push_back(tokenTypeModel.getCDF());
                        tokenTypeSymbols.push_back(tokenType);
                        tokenTypeModel.update(tokenType);
                        rawData.push_back(static_cast<uint8_t>(static_cast<int8_t>(delta)));
                    } else {
                        tokenType = static_cast<int>(TokenType::NUMBER);
                        tokenTypeCDFs.push_back(tokenTypeModel.getCDF());
                        tokenTypeSymbols.push_back(tokenType);
                        tokenTypeModel.update(tokenType);
                        rawData.push_back(tok.numWidth);
                        writeVarint(rawData, tok.numValue);
                    }
                    prevNumber_ = tok.numValue;
                    prevWidth_ = tok.numWidth;
                } else if (tok.type == TokenType::LITERAL) {
                    auto it = dictIndex_.find(tok.value);
                    if (it != dictIndex_.end()) {
                        tokenType = static_cast<int>(TokenType::DICT_REF);
                        tokenTypeCDFs.push_back(tokenTypeModel.getCDF());
                        tokenTypeSymbols.push_back(tokenType);
                        tokenTypeModel.update(tokenType);

                        // Entropy code the dictionary reference
                        int dictIdx = static_cast<int>(it->second);
                        dictRefCDFs.push_back(dictRefModel.getCDF());
                        dictRefSymbols.push_back(dictIdx);
                        dictRefModel.update(dictIdx);
                    } else {
                        tokenType = static_cast<int>(TokenType::LITERAL);
                        tokenTypeCDFs.push_back(tokenTypeModel.getCDF());
                        tokenTypeSymbols.push_back(tokenType);
                        tokenTypeModel.update(tokenType);
                        writeVarint(rawData, tok.value.size());
                        for (char c : tok.value) {
                            rawData.push_back(static_cast<uint8_t>(c));
                        }
                    }
                }
            }
        }

        // Encode the streams
        auto tokenTypeData = IdentifierRANSEncoder::encode(tokenTypeSymbols, tokenTypeCDFs, NUM_TOKEN_TYPES);
        auto delimiterData = IdentifierRANSEncoder::encode(delimiterSymbols, delimiterCDFs, NUM_DELIMITERS);
        auto dictRefData = IdentifierRANSEncoder::encode(dictRefSymbols, dictRefCDFs, dictSize);

        // Entropy-code the raw data stream using order-0 adaptive model
        size_t rawDataOrigSize = rawData.size();
        std::vector<uint8_t> rawDataCompressed;
        if (!rawData.empty()) {
            AdaptiveSymbolModel<256> rawDataModel;
            std::vector<int> rawDataSymbols;
            std::vector<std::vector<uint32_t>> rawDataCDFs;
            rawDataSymbols.reserve(rawData.size());
            rawDataCDFs.reserve(rawData.size());
            for (uint8_t b : rawData) {
                rawDataCDFs.push_back(rawDataModel.getCDF());
                rawDataSymbols.push_back(b);
                rawDataModel.update(b);
            }
            rawDataCompressed = IdentifierRANSEncoder::encode(rawDataSymbols, rawDataCDFs, 256);
        }

        // Build output
        std::vector<uint8_t> output;

        // Write dictionary
        writeVarint(output, dictionary_.size());
        for (const auto& entry : dictionary_) {
            writeVarint(output, entry.size());
            for (char c : entry) {
                output.push_back(static_cast<uint8_t>(c));
            }
        }

        // Write record count
        writeVarint(output, identifiers.size());

        // Write stream sizes
        writeVarint(output, tokenTypeSymbols.size());
        writeVarint(output, tokenTypeData.size());
        writeVarint(output, delimiterSymbols.size());
        writeVarint(output, delimiterData.size());
        writeVarint(output, dictRefSymbols.size());
        writeVarint(output, dictRefData.size());
        writeVarint(output, rawDataOrigSize);           // Uncompressed size
        writeVarint(output, rawDataCompressed.size());  // Compressed size

        // Write streams
        output.insert(output.end(), tokenTypeData.begin(), tokenTypeData.end());
        output.insert(output.end(), delimiterData.begin(), delimiterData.end());
        output.insert(output.end(), dictRefData.begin(), dictRefData.end());
        output.insert(output.end(), rawDataCompressed.begin(), rawDataCompressed.end());

        return output;
    }

    void reset() {
        dictionary_.clear();
        dictIndex_.clear();
        prevNumber_ = 0;
        prevWidth_ = 0;
    }

private:
    std::vector<std::string> dictionary_;
    std::unordered_map<std::string, size_t> dictIndex_;
    uint64_t prevNumber_;
    uint8_t prevWidth_;
};

/**
 * IdentifierDecoderV2 - Entropy-coded identifier decoder
 */
class IdentifierDecoderV2 {
public:
    IdentifierDecoderV2() : prevNumber_(0), prevWidth_(0) {}

    std::vector<std::string> decode(const uint8_t* data, size_t len) {
        if (len == 0) {
            return {};
        }

        const uint8_t* p = data;
        const uint8_t* end = data + len;

        // Read dictionary
        dictionary_.clear();
        size_t dictSize = readVarint(p);
        dictionary_.reserve(dictSize);
        for (size_t i = 0; i < dictSize && p < end; i++) {
            size_t entryLen = readVarint(p);
            std::string entry;
            entry.reserve(entryLen);
            for (size_t j = 0; j < entryLen && p < end; j++) {
                entry.push_back(static_cast<char>(*p++));
            }
            dictionary_.push_back(std::move(entry));
        }

        // Read record count
        size_t recordCount = readVarint(p);

        // Read stream sizes
        size_t tokenTypeCount = readVarint(p);
        size_t tokenTypeDataLen = readVarint(p);
        size_t delimiterCount = readVarint(p);
        size_t delimiterDataLen = readVarint(p);
        size_t dictRefCount = readVarint(p);
        size_t dictRefDataLen = readVarint(p);
        size_t rawDataOrigSize = readVarint(p);
        size_t rawDataCompressedLen = readVarint(p);

        // Get stream pointers
        const uint8_t* tokenTypeData = p;
        p += tokenTypeDataLen;
        const uint8_t* delimiterData = p;
        p += delimiterDataLen;
        const uint8_t* dictRefData = p;
        p += dictRefDataLen;
        const uint8_t* rawDataCompressed = p;

        // Decompress raw data stream
        std::vector<uint8_t> rawDataBuf(rawDataOrigSize);
        if (rawDataOrigSize > 0) {
            AdaptiveSymbolModel<256> rawDataModel;
            IdentifierRANSDecoder::decode(
                rawDataCompressed, rawDataCompressedLen, rawDataOrigSize, 256,
                [&](size_t) { return rawDataModel.getCDF(); },
                [&](size_t i, int sym) { rawDataBuf[i] = static_cast<uint8_t>(sym); rawDataModel.update(sym); }
            );
        }
        const uint8_t* rawData = rawDataBuf.data();
        const uint8_t* rawEnd = rawData + rawDataOrigSize;

        // Decode token types
        std::vector<int> tokenTypes(tokenTypeCount);
        AdaptiveSymbolModel<NUM_TOKEN_TYPES> tokenTypeModel;
        IdentifierRANSDecoder::decode(
            tokenTypeData, tokenTypeDataLen, tokenTypeCount, NUM_TOKEN_TYPES,
            [&](size_t /*i*/) { return tokenTypeModel.getCDF(); },
            [&](size_t i, int sym) {
                tokenTypes[i] = sym;
                tokenTypeModel.update(sym);
            }
        );

        // Decode delimiters
        std::vector<int> delimiters(delimiterCount);
        AdaptiveSymbolModel<NUM_DELIMITERS> delimiterModel;
        IdentifierRANSDecoder::decode(
            delimiterData, delimiterDataLen, delimiterCount, NUM_DELIMITERS,
            [&](size_t /*i*/) { return delimiterModel.getCDF(); },
            [&](size_t i, int sym) {
                delimiters[i] = sym;
                delimiterModel.update(sym);
            }
        );

        // Decode dictionary references (right-sized model)
        std::vector<int> dictRefs(dictRefCount);
        int dictAlphabetSize = std::max(1, static_cast<int>(dictionary_.size()));
        AdaptiveSymbolModel<256> dictRefModel;
        dictRefModel.setActiveSize(dictAlphabetSize);
        IdentifierRANSDecoder::decode(
            dictRefData, dictRefDataLen, dictRefCount, dictAlphabetSize,
            [&](size_t /*i*/) { return dictRefModel.getCDF(); },
            [&](size_t i, int sym) {
                dictRefs[i] = sym;
                dictRefModel.update(sym);
            }
        );

        // Reconstruct identifiers
        std::vector<std::string> results;
        results.reserve(recordCount);

        size_t tokenIdx = 0;
        size_t delimIdx = 0;
        size_t dictRefIdx = 0;
        prevNumber_ = 0;
        prevWidth_ = 0;

        for (size_t rec = 0; rec < recordCount && rawData < rawEnd; rec++) {
            size_t tokenCount = readVarint(rawData);
            std::string result;

            for (size_t t = 0; t < tokenCount && tokenIdx < tokenTypes.size(); t++, tokenIdx++) {
                TokenType type = static_cast<TokenType>(tokenTypes[tokenIdx]);

                switch (type) {
                    case TokenType::DELIMITER:
                        if (delimIdx < delimiters.size()) {
                            int dIdx = delimiters[delimIdx++];
                            if (dIdx == 10 && rawData < rawEnd) {
                                // Unknown delimiter: read raw char from data stream
                                result.push_back(static_cast<char>(*rawData++));
                            } else {
                                result.push_back(indexToDelimiter(dIdx));
                            }
                        }
                        break;

                    case TokenType::NUMBER: {
                        uint8_t width = *rawData++;
                        uint64_t num = readVarint(rawData);
                        result += formatNumberWithWidth(num, width);
                        prevNumber_ = num;
                        prevWidth_ = width;
                        break;
                    }

                    case TokenType::DELTA_NUM: {
                        int8_t delta = static_cast<int8_t>(*rawData++);
                        uint64_t num = static_cast<uint64_t>(
                            static_cast<int64_t>(prevNumber_) + delta);
                        result += formatNumberWithWidth(num, prevWidth_);
                        prevNumber_ = num;
                        break;
                    }

                    case TokenType::DICT_REF:
                        if (dictRefIdx < dictRefs.size()) {
                            size_t idx = dictRefs[dictRefIdx++];
                            if (idx < dictionary_.size()) {
                                result += dictionary_[idx];
                            }
                        }
                        break;

                    case TokenType::LITERAL: {
                        size_t literalLen = readVarint(rawData);
                        for (size_t j = 0; j < literalLen && rawData < rawEnd; j++) {
                            result.push_back(static_cast<char>(*rawData++));
                        }
                        break;
                    }

                    case TokenType::REPEAT_PREV:
                        // Not used
                        break;
                }
            }

            results.push_back(std::move(result));
        }

        return results;
    }

    void reset() {
        dictionary_.clear();
        prevNumber_ = 0;
        prevWidth_ = 0;
    }

private:
    std::vector<std::string> dictionary_;
    uint64_t prevNumber_;
    uint8_t prevWidth_;
};

} // namespace v4zip

#endif // V4ZIP_IDENTIFIERCODEC_HPP
