#pragma once
#include "v4zip/dna.hpp"
#include <array>
#include <limits>
#include <memory>
#include <stdexcept>
#include <zstd.h>
namespace v4zip::detail {
[[noreturn]] inline void fail(const char* message) { throw std::runtime_error(message); }
inline std::uint32_t crc32(const std::uint8_t* data, std::size_t n) {
    static const auto table = [] {
        std::array<std::uint32_t, 256> values{};
        for (std::uint32_t i = 0; i < 256; ++i) {
            auto value = i;
            for (unsigned j = 0; j < 8; ++j) value = (value >> 1) ^ (0xedb88320U & (0U - (value & 1U)));
            values[i] = value;
        }
        return values;
    }();
    std::uint32_t value = 0xffffffffU;
    for (std::size_t i = 0; i < n; ++i) value = table[(value ^ data[i]) & 255U] ^ (value >> 8);
    return value ^ 0xffffffffU;
}
inline std::uint32_t crc32(const Bytes& b) { return crc32(b.data(), b.size()); }
inline void put(Bytes& data, std::size_t pos, std::uint64_t value, unsigned n) {
    for (unsigned i = 0; i < n; ++i) data[pos + i] = static_cast<std::uint8_t>(value >> (8 * i));
}
inline std::uint64_t get(const Bytes& data, std::size_t pos, unsigned n) {
    std::uint64_t result = 0;
    for (unsigned i = 0; i < n; ++i) result |= std::uint64_t(data[pos + i]) << (8 * i);
    return result;
}
inline void varint(Bytes& data, std::uint64_t value) {
    do {
        const auto low = static_cast<std::uint8_t>(value & 127);
        value >>= 7;
        data.push_back(static_cast<std::uint8_t>(low | (value ? 128 : 0)));
    } while (value);
}
inline std::uint64_t read_uint64(const Bytes& data, std::size_t& pos) {
    std::uint64_t value = 0;
    unsigned shift = 0;
    for (;;) {
        if (pos == data.size() || shift > 63) fail("invalid layout length");
        const auto byte = data[pos++];
        if (shift == 63 && (byte & 127) > 1) fail("layout length overflow");
        value |= std::uint64_t(byte & 127) << shift;
        if (!(byte & 128)) {
            if (shift && !(byte & 127)) fail("noncanonical layout length");
            return value;
        }
        shift += 7;
    }
}
inline std::size_t read_varint(const Bytes& data,std::size_t& pos) {
    const auto value=read_uint64(data,pos);
    if(!value || value>std::numeric_limits<std::size_t>::max())fail("invalid layout length");
    return static_cast<std::size_t>(value);
}
inline Bytes pack(const Bytes& data) {
    Bytes result(ZSTD_compressBound(data.size()));
    const auto n = ZSTD_compress(result.data(), result.size(), data.data(), data.size(), 9);
    if (ZSTD_isError(n)) fail("zstd compression failed");
    result.resize(n);
    return result;
}
inline Bytes unpack(const std::uint8_t* data, std::size_t size, std::size_t expected) {
    const auto content_size = ZSTD_getFrameContentSize(data, size);
    if (content_size == ZSTD_CONTENTSIZE_ERROR || content_size == ZSTD_CONTENTSIZE_UNKNOWN || content_size != expected)
        fail("invalid zstd content size");
    const auto frame_size = ZSTD_findFrameCompressedSize(data, size);
    if (ZSTD_isError(frame_size) || frame_size != size) fail("invalid zstd frame extent");
    Bytes result(expected);
    std::unique_ptr<ZSTD_DCtx, decltype(&ZSTD_freeDCtx)> ctx(ZSTD_createDCtx(), &ZSTD_freeDCtx);
    if (!ctx || ZSTD_isError(ZSTD_DCtx_setParameter(ctx.get(), ZSTD_d_windowLogMax, 27))) fail("zstd initialization failed");
    const auto n = ZSTD_decompressDCtx(ctx.get(), result.data(), result.size(), data, size);
    if (ZSTD_isError(n) || n != expected) fail("zstd decompression failed");
    return result;
}

inline std::size_t read_size(const Bytes& data,std::size_t& pos) {
    if(pos<data.size() && data[pos]==0){++pos;return 0;}
    return read_varint(data,pos);
}
}
