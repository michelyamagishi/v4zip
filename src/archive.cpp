#include "v4zip/archive.hpp"
#include "detail/bytes.hpp"
#include "detail/dna_config.hpp"

#include <algorithm>
#include <array>
#include <cstring>
#include <limits>
#include <memory>
#include <stdexcept>
#include <zstd.h>

namespace v4zip {
namespace {
using namespace detail;
constexpr std::size_t header_size = 96;
constexpr std::array<std::uint8_t, 8> magic{{'V','4','Z','I','P','N','E','W'}};
void check_layout(const Bytes& layout, std::size_t original, std::size_t dna, std::size_t qualities, bool fastq) {
    std::size_t pos = 0, out = 0, bases = 0, quality_count = 0;
    while (pos < layout.size()) {
        const auto type = layout[pos++];
        if (type > 3 || (!fastq && type == 3)) fail("invalid layout token");
        const auto count = read_varint(layout, pos);
        if (count > original - out) fail("layout exceeds original length");
        out += count;
        if (type == 0) {
            if (count > layout.size() - pos) fail("truncated literal");
            pos += count;
        } else if (type == 3) {
            if (count > qualities - quality_count) fail("layout exceeds quality count");
            quality_count += count;
        } else {
            if (count > dna - bases) fail("layout exceeds DNA count");
            bases += count;
        }
    }
    if (out != original || bases != dna || quality_count != qualities) fail("layout counts differ");
}
}

Bytes decode_archive_v1(const Bytes& archive, std::size_t max_output) {
    if (archive.size() < header_size || !std::equal(magic.begin(), magic.end(), archive.begin()))
        fail("not a new V4ZIP archive (legacy archives require the original tool)");
    if (get(archive, 8, 4) != 1 || get(archive, 12, 4) > 1 || get(archive, 92, 4)) fail("unsupported archive version or flags");
    if (get(archive, 88, 4) != crc32(archive.data(), 88)) fail("archive header checksum mismatch");
    const auto original64 = get(archive, 16, 8);
    if (original64 > max_output || original64 > std::numeric_limits<std::size_t>::max()) fail("output exceeds --max-output");
    const auto original = static_cast<std::size_t>(original64);
    const auto dna64 = get(archive, 24, 8), layout64 = get(archive, 32, 8), quality64 = get(archive, 40, 8);
    if (dna64 > original || quality64 > original - dna64 || (!get(archive, 12, 4) && quality64)) fail("invalid stream counts");
    // Each nonempty run uses a tag and at most one varint byte per represented
    // byte, plus literal bytes. Thus layout <= 3 * original, without overflow.
    if (layout64 > std::numeric_limits<std::size_t>::max() || (layout64 / 3 + (layout64 % 3 != 0) > original))
        fail("invalid layout size");
    std::array<std::size_t, 3> offsets{}, sizes{};
    std::size_t pos = header_size;
    for (unsigned i = 0; i < offsets.size(); ++i) {
        const auto n = get(archive, 48 + 8 * i, 8);
        if (n > archive.size() - pos) fail("truncated archive stream");
        offsets[i] = pos;
        sizes[i] = static_cast<std::size_t>(n);
        pos += sizes[i];
    }
    if (pos != archive.size()) fail("trailing archive data");
    for (unsigned i = 0; i < offsets.size(); ++i)
        if (get(archive, 76 + 4 * i, 4) != crc32(archive.data() + offsets[i], sizes[i])) fail("archive stream checksum mismatch");
    const auto layout = unpack(archive.data() + offsets[1], sizes[1], static_cast<std::size_t>(layout64));
    check_layout(layout, original, static_cast<std::size_t>(dna64), static_cast<std::size_t>(quality64), get(archive, 12, 4) != 0);
    const auto qualities = unpack(archive.data() + offsets[2], sizes[2], static_cast<std::size_t>(quality64));
    const Bytes encoded_dna(archive.begin() + offsets[0], archive.begin() + offsets[0] + sizes[0]);
    const auto dna = detail::decode_dna(encoded_dna, static_cast<std::size_t>(dna64), detail::DnaConfig{1,{}});
    if (dna.size() != dna64) fail("decoded DNA count differs");
    Bytes result;
    result.reserve(original);
    std::size_t cursor = 0, base = 0, quality = 0;
    while (cursor < layout.size()) {
        const auto type = layout[cursor++];
        const auto n = read_varint(layout, cursor);
        if (type == 0) {
            result.insert(result.end(), layout.begin() + cursor, layout.begin() + cursor + n);
            cursor += n;
        } else if (type == 3) {
            result.insert(result.end(), qualities.begin() + quality, qualities.begin() + quality + n);
            quality += n;
        } else {
            for (std::size_t i = 0; i < n; ++i) {
                const auto b = dna[base++];
                if (b > 3) fail("invalid decoded DNA symbol");
                result.push_back(static_cast<std::uint8_t>((type == 1 ? "ACGT" : "acgt")[b]));
            }
        }
    }
    if (crc32(result) != get(archive, 72, 4)) fail("original checksum mismatch");
    return result;
}
}
