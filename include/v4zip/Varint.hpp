/**
 * Varint.hpp - Variable-length integer encoding/decoding utilities
 *
 * Uses standard LEB128 (Little Endian Base 128) encoding:
 * - Each byte uses 7 bits for data
 * - High bit (0x80) indicates continuation
 *
 * Format:
 *   Values 0-127:       1 byte  (0x00-0x7F)
 *   Values 128-16383:   2 bytes (0x80-0xFF, 0x00-0x7F)
 *   Values 16384+:      3+ bytes
 */

#ifndef V4ZIP_VARINT_HPP
#define V4ZIP_VARINT_HPP

#include <cstdint>
#include <vector>
#include <stdexcept>

namespace v4zip {

/**
 * Write a variable-length integer to a byte vector
 *
 * @param out Output vector to append to
 * @param v Value to encode (up to 64 bits)
 */
inline void writeVarint(std::vector<uint8_t>& out, uint64_t v) {
    do {
        out.push_back((v & 0x7F) | (v >= 128 ? 0x80 : 0));
        v >>= 7;
    } while (v);
}

/**
 * Read a variable-length integer from a byte stream
 *
 * @param p Pointer to current position (updated after read)
 * @return Decoded value
 */
inline uint64_t readVarint(const uint8_t*& p) {
    uint64_t v = 0;
    int shift = 0;
    uint8_t byte;
    do {
        byte = *p++;
        v |= (uint64_t)(byte & 0x7F) << shift;
        shift += 7;
    } while (byte & 0x80);
    return v;
}

/**
 * Read a variable-length integer with bounds checking
 *
 * @param p Pointer to current position (updated after read)
 * @param end Pointer past the end of the buffer
 * @return Decoded value
 * @throws std::runtime_error on truncated or overlong varint
 */
inline uint64_t readVarint(const uint8_t*& p, const uint8_t* end) {
    uint64_t v = 0;
    int shift = 0;
    while (p < end) {
        uint8_t byte = *p++;
        v |= static_cast<uint64_t>(byte & 0x7F) << shift;
        if (!(byte & 0x80)) return v;
        shift += 7;
        if (shift >= 64) throw std::runtime_error("Varint overflow: too many continuation bytes");
    }
    throw std::runtime_error("Truncated varint: unexpected end of data");
}

/**
 * Calculate the encoded size of a varint
 *
 * @param v Value to measure
 * @return Number of bytes needed
 */
inline size_t varintSize(uint64_t v) {
    size_t size = 1;
    while (v >= 128) {
        v >>= 7;
        size++;
    }
    return size;
}

/**
 * Zigzag encoding for signed integers
 * Maps signed values to unsigned: 0→0, -1→1, 1→2, -2→3, 2→4, ...
 */
inline uint64_t zigzagEncode(int64_t v) {
    return static_cast<uint64_t>((v << 1) ^ (v >> 63));
}

inline int64_t zigzagDecode(uint64_t v) {
    return static_cast<int64_t>((v >> 1) ^ -(static_cast<int64_t>(v) & 1));
}

// Write a zigzag-encoded signed integer
inline void writeSignedVarint(std::vector<uint8_t>& out, int64_t v) {
    writeVarint(out, zigzagEncode(v));
}

// Read a zigzag-encoded signed integer
inline int64_t readSignedVarint(const uint8_t*& p) {
    return zigzagDecode(readVarint(p));
}

// Read a zigzag-encoded signed integer with bounds checking
inline int64_t readSignedVarint(const uint8_t*& p, const uint8_t* end) {
    return zigzagDecode(readVarint(p, end));
}

} // namespace v4zip

#endif // V4ZIP_VARINT_HPP
