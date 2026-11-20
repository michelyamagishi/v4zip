#pragma once

#include "bytes.hpp"
#include "v4zip/archive.hpp"
#include <charconv>
#include <cstdio>
#include <filesystem>
#include <sstream>
#include <string>
#include <vector>

namespace v4zip::detail {

inline std::size_t parse_byte_size(const std::string& text, const char* flag) {
    if (text.empty()) throw std::runtime_error(std::string(flag) + " requires a size");
    std::size_t mult = 1;
    char suffix = text.back();
    std::string num_str = text;
    if (suffix == 'K' || suffix == 'k') { mult = 1024ULL; num_str.pop_back(); }
    else if (suffix == 'M' || suffix == 'm') { mult = 1048576ULL; num_str.pop_back(); }
    else if (suffix == 'G' || suffix == 'g') { mult = 1073741824ULL; num_str.pop_back(); }
    if (num_str.empty()) throw std::runtime_error(std::string(flag) + " invalid number");
    std::uint64_t val = 0;
    auto res = std::from_chars(num_str.data(), num_str.data() + num_str.size(), val);
    if (res.ec != std::errc() || res.ptr != num_str.data() + num_str.size())
        throw std::runtime_error(std::string(flag) + " invalid number");
    if (val > std::numeric_limits<std::size_t>::max() / mult)
        throw std::runtime_error(std::string(flag) + " size overflow");
    std::size_t result = static_cast<std::size_t>(val * mult);
    if (result == 0) throw std::runtime_error(std::string(flag) + " must be positive");
    return result;
}

inline std::string split_part_path(const std::string& base, std::size_t index) {
    char buf[32];
    if (index < 1000) std::snprintf(buf, sizeof(buf), ".%03zu", index);
    else std::snprintf(buf, sizeof(buf), ".%zu", index);
    return base + buf;
}

inline bool is_split_first_part(const std::string& path) {
    return path.size() > 4 && path.substr(path.size() - 4) == ".001";
}

inline bool is_manifest_path(const std::string& path) {
    return path.size() > 5 && path.substr(path.size() - 5) == ".v4zm";
}

inline std::string chunk_part_name(const std::string& manifest_path, std::size_t index) {
    char buf[32];
    std::snprintf(buf, sizeof(buf), ".part%04zu.v4z", index);
    return std::filesystem::path(manifest_path).stem().string() + buf;
}

inline std::vector<std::pair<std::size_t, std::size_t>> plan_chunks(
    const Bytes& input, bool fastq, std::size_t chunk_size) {
    if (input.empty()) return {{0, 0}};
    std::vector<std::pair<std::size_t, std::size_t>> chunks;
    std::size_t start = 0, pos = 0, n = input.size();
    if (fastq) {
        std::size_t lines = 0;
        while (pos < n) {
            if (input[pos] == '\n') {
                ++lines;
                if (lines % 4 == 0) {
                    std::size_t end = pos + 1;
                    if (end - start >= chunk_size && end < n) {
                        chunks.emplace_back(start, end - start);
                        start = end;
                    }
                }
            }
            ++pos;
        }
    } else {
        while (pos < n) {
            if (input[pos] == '\n') {
                std::size_t next = pos + 1;
                if (next < n && input[next] == '>' && next - start >= chunk_size) {
                    chunks.emplace_back(start, next - start);
                    start = next;
                }
            }
            ++pos;
        }
    }
    if (start < n) chunks.emplace_back(start, n - start);
    return chunks;
}

struct ManifestEntry {
    std::string filename;
    std::size_t raw_size = 0;
    std::uint32_t raw_crc = 0;
};

struct Manifest {
    bool fastq = false;
    std::size_t total_size = 0;
    std::uint32_t total_crc = 0;
    std::vector<ManifestEntry> parts;
};

inline Bytes serialize_manifest(const Manifest& m) {
    std::ostringstream oss;
    oss << "V4ZM1 " << (m.fastq ? 1 : 0) << ' ' << m.total_size << ' '
        << m.total_crc << ' ' << m.parts.size() << '\n';
    for (const auto& p : m.parts)
        oss << p.filename << ' ' << p.raw_size << ' ' << p.raw_crc << '\n';
    std::string s = oss.str();
    return Bytes(s.begin(), s.end());
}

inline Manifest parse_manifest(const Bytes& data) {
    std::string text(data.begin(), data.end());
    std::istringstream iss(text);
    std::string magic;
    int fq_int = 0;
    std::size_t declared_parts = 0;
    Manifest m;
    if (!(iss >> magic >> fq_int >> m.total_size >> m.total_crc >> declared_parts))
        throw std::runtime_error("malformed manifest header");
    if (magic != "V4ZM1" || (fq_int != 0 && fq_int != 1))
        throw std::runtime_error("invalid manifest format or flag");
    m.fastq = (fq_int != 0);
    if (declared_parts == 0 || declared_parts > 10000)
        throw std::runtime_error("invalid manifest part count");
    m.parts.reserve(declared_parts);
    std::size_t sum_size = 0;
    for (std::size_t i = 0; i < declared_parts; ++i) {
        ManifestEntry entry;
        if (!(iss >> entry.filename >> entry.raw_size >> entry.raw_crc))
            throw std::runtime_error("incomplete manifest parts");
        if (entry.filename.empty() || entry.filename == "." || entry.filename == ".." ||
            entry.filename.find('/') != std::string::npos ||
            entry.filename.find('\\') != std::string::npos)
            throw std::runtime_error("unsafe manifest part filename: " + entry.filename);
        if (std::numeric_limits<std::size_t>::max() - sum_size < entry.raw_size)
            throw std::runtime_error("manifest total size overflow");
        sum_size += entry.raw_size;
        m.parts.push_back(std::move(entry));
    }
    if (sum_size != m.total_size)
        throw std::runtime_error("manifest size mismatch: parts sum != total_size");
    std::string extra;
    if (iss >> extra) throw std::runtime_error("extra data in manifest");
    return m;
}

} // namespace v4zip::detail
