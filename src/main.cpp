#include "v4zip/archive.hpp"
#include "detail/partition.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <future>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <string>
#include <system_error>
#include <vector>

#include <fcntl.h>
#include <sys/stat.h>
#include <unistd.h>

namespace {
using v4zip::Bytes;
constexpr std::size_t default_limit = std::size_t{1} << 30;
void system_failure(const char* action) { throw std::runtime_error(std::string(action) + ": " + std::strerror(errno)); }
class Fd {
public:
    explicit Fd(int fd) : value(fd) {}
    ~Fd() { if (value >= 0) ::close(value); }
    Fd(const Fd&) = delete;
    Fd& operator=(const Fd&) = delete;
    int value;
};
Bytes read_file(const std::string& path) {
    // The descriptor fixes which inode is read, even when its path changes.
    Fd fd(::open(path.c_str(), O_RDONLY | O_NONBLOCK | O_CLOEXEC));
    if (fd.value < 0) system_failure("open input");
    struct stat info{};
    if (::fstat(fd.value, &info) < 0) system_failure("stat input");
    if (!S_ISREG(info.st_mode)) throw std::runtime_error("input must be a regular file");
    if (info.st_size < 0 || static_cast<std::uintmax_t>(info.st_size) > Bytes().max_size())
        throw std::runtime_error("input file is too large for this platform");
    Bytes result;
    result.reserve(static_cast<std::size_t>(info.st_size));
    std::array<std::uint8_t, 65536> buffer{};
    for (;;) {
        const auto n = ::read(fd.value, buffer.data(), buffer.size());
        if (n < 0) { if (errno == EINTR) continue; system_failure("read input"); }
        if (n == 0) break;
        result.insert(result.end(), buffer.begin(), buffer.begin() + n);
    }
    return result;
}
class Destination {
public:
    explicit Destination(const std::string& path) : parent(open_parent(path)), name(std::filesystem::path(path).filename().string()) {
        if (name.empty() || name == "." || name == "..") throw std::runtime_error("invalid output filename");
        struct stat info{};
        if (::fstatat(parent.value, name.c_str(), &info, AT_SYMLINK_NOFOLLOW) == 0)
            throw std::runtime_error("output already exists; refusing to overwrite");
        if (errno != ENOENT) system_failure("inspect output");
    }
    ~Destination() { if (!temporary.empty()) ::unlinkat(parent.value, temporary.c_str(), 0); }
    void publish(const Bytes& data) {
        int file = -1;
        const auto nonce = std::chrono::steady_clock::now().time_since_epoch().count();
        for (unsigned i = 0; i < 100; ++i) {
            temporary = ".v4zip-" + std::to_string(::getpid()) + "-" + std::to_string(nonce) + "-" + std::to_string(i);
            file = ::openat(parent.value, temporary.c_str(), O_WRONLY | O_CREAT | O_EXCL | O_CLOEXEC | O_NOFOLLOW, 0600);
            if (file >= 0) break;
            if (errno != EEXIST) { temporary.clear(); system_failure("create temporary output"); }
            // A colliding file is not ours and must never be cleaned up.
            temporary.clear();
        }
        if (file < 0) throw std::runtime_error("cannot reserve temporary output filename");
        Fd fd(file);
        std::size_t pos = 0;
        while (pos < data.size()) {
            const auto n = ::write(fd.value, data.data() + pos, std::min<std::size_t>(data.size() - pos, 1U << 20));
            if (n < 0) { if (errno == EINTR) continue; system_failure("write output"); }
            if (!n) throw std::runtime_error("zero-length output write");
            pos += static_cast<std::size_t>(n);
        }
        if (::fsync(fd.value) < 0) system_failure("flush output");
        // linkat publishes atomically and fails if *any* directory entry now
        // occupies the destination. No exists/rename overwrite race is possible.
        if (::linkat(parent.value, temporary.c_str(), parent.value, name.c_str(), 0) < 0) system_failure("publish output (overwrite refused)");
        ::unlinkat(parent.value, temporary.c_str(), 0);
        temporary.clear();
        // Publication has succeeded. Some filesystems do not support dir fsync.
        (void)::fsync(parent.value);
    }
private:
    static int open_parent(const std::string& path) {
        auto directory = std::filesystem::path(path).parent_path();
        if (directory.empty()) directory = ".";
        const auto fd = ::open(directory.c_str(), O_RDONLY | O_DIRECTORY | O_CLOEXEC);
        if (fd < 0) system_failure("open output directory");
        return fd;
    }
    Fd parent;
    std::string name, temporary;
};
void help() {
    std::cout << "v4zip " V4ZIP_VERSION " - exact FASTA/FASTQ compression\n"
              << "Usage: v4zip [fastq] compress INPUT OUTPUT [--split SIZE | --chunk-size SIZE] [--verbose]\n"
              << "       v4zip [fastq] decompress INPUT OUTPUT [--max-output BYTES] [--verbose]\n"
              << "       v4zip --help | --version\n"
              << "Default decompression limit: 1073741824 bytes (1 GiB).\n"
              << "Existing outputs are never overwritten. Legacy archives require the old tool.\n";
}
std::size_t parse_limit(const std::string& text) {
    std::uint64_t value = 0;
    const auto result = std::from_chars(text.data(), text.data() + text.size(), value);
    if (text.empty() || result.ec != std::errc() || result.ptr != text.data() + text.size() || value > std::numeric_limits<std::size_t>::max())
        throw std::runtime_error("--max-output requires a nonnegative decimal byte count");
    return static_cast<std::size_t>(value);
}
}
int main(int argc, char** argv) {
    try {
        if (argc == 2 && std::string(argv[1]) == "--help") { help(); return 0; }
        if (argc == 2 && std::string(argv[1]) == "--version") { std::cout << "v4zip " V4ZIP_VERSION "\n"; return 0; }
        int arg = 1;
        bool fastq = false, verbose = false, limit_set = false;
        std::size_t limit = default_limit, split_size = 0, chunk_size = 0;
        if (arg < argc && std::string(argv[arg]) == "fastq") { fastq = true; ++arg; }
        if (arg == argc) throw std::runtime_error("expected compress or decompress; use --help");
        const std::string command = argv[arg++];
        if (command != "compress" && command != "decompress") throw std::runtime_error("unknown command; use --help");
        std::vector<std::string> paths;
        bool options = true;
        for (; arg < argc; ++arg) {
            const std::string value = argv[arg];
            if (options && value == "--") options = false;
            else if (options && value == "--verbose") {
                if (verbose) throw std::runtime_error("duplicate --verbose");
                verbose = true;
            } else if (options && value == "--max-output") {
                if (limit_set || command != "decompress" || ++arg == argc) throw std::runtime_error("invalid --max-output usage");
                limit = parse_limit(argv[arg]);
                limit_set = true;
            } else if (options && value == "--split") {
                if (split_size || command != "compress" || ++arg == argc) throw std::runtime_error("invalid --split usage");
                split_size = v4zip::detail::parse_byte_size(argv[arg], "--split");
            } else if (options && value == "--chunk-size") {
                if (chunk_size || command != "compress" || ++arg == argc) throw std::runtime_error("invalid --chunk-size usage");
                chunk_size = v4zip::detail::parse_byte_size(argv[arg], "--chunk-size");
            } else if (options && !value.empty() && value[0] == '-') throw std::runtime_error("unknown option: " + value);
            else paths.push_back(value);
        }
        if (split_size && chunk_size) throw std::runtime_error("--split and --chunk-size are mutually exclusive");
        if (paths.size() != 2) throw std::runtime_error("expected exactly INPUT and OUTPUT; use --help");
        const auto start = std::chrono::steady_clock::now();
        std::size_t in_bytes = 0, out_bytes = 0;
        if (command == "compress") {
            const auto input = read_file(paths[0]);
            in_bytes = input.size();
            if (split_size) {
                const auto archive = v4zip::compress_archive(input, fastq);
                out_bytes = archive.size();
                std::size_t part_idx = 1;
                for (std::size_t off = 0; off < archive.size(); off += split_size, ++part_idx) {
                    const std::size_t len = std::min(split_size, archive.size() - off);
                    Destination part_dest(v4zip::detail::split_part_path(paths[1], part_idx));
                    part_dest.publish(Bytes(archive.begin() + off, archive.begin() + off + len));
                }
            } else if (chunk_size) {
                if (!v4zip::detail::is_manifest_path(paths[1])) throw std::runtime_error("manifest output must end with .v4zm");
                Destination manifest_dest(paths[1]);
                const auto chunks = v4zip::detail::plan_chunks(input, fastq, chunk_size);
                const auto parent_dir = std::filesystem::path(paths[1]).parent_path();
                v4zip::detail::Manifest manifest{fastq, input.size(), v4zip::detail::crc32(input), {}};
                for (std::size_t i = 0; i < chunks.size(); i += 4) {
                    const std::size_t batch_end = std::min(chunks.size(), i + 4);
                    std::vector<std::future<Bytes>> tasks;
                    for (std::size_t j = i; j < batch_end; ++j) {
                        tasks.push_back(std::async(std::launch::async, [&input, &chunks, j, fastq]() {
                            return v4zip::compress_archive(Bytes(input.begin() + chunks[j].first, input.begin() + chunks[j].first + chunks[j].second), fastq);
                        }));
                    }
                    for (std::size_t j = i; j < batch_end; ++j) {
                        const auto comp = tasks[j - i].get();
                        out_bytes += comp.size();
                        const auto part_name = v4zip::detail::chunk_part_name(paths[1], j + 1);
                        const auto part_path = parent_dir.empty() ? part_name : (parent_dir / part_name).string();
                        Destination part_dest(part_path);
                        part_dest.publish(comp);
                        Bytes slice(input.begin() + chunks[j].first, input.begin() + chunks[j].first + chunks[j].second);
                        manifest.parts.push_back({part_name, chunks[j].second, v4zip::detail::crc32(slice)});
                    }
                }
                const auto man_bytes = v4zip::detail::serialize_manifest(manifest);
                out_bytes += man_bytes.size();
                manifest_dest.publish(man_bytes);
            } else {
                Destination destination(paths[1]);
                const auto output = v4zip::compress_archive(input, fastq);
                out_bytes = output.size();
                destination.publish(output);
            }
        } else {
            Destination destination(paths[1]);
            Bytes output;
            if (v4zip::detail::is_split_first_part(paths[0])) {
                Bytes archive;
                const auto base = paths[0].substr(0, paths[0].size() - 4);
                const std::size_t max_archive = limit <= (std::numeric_limits<std::size_t>::max() - 65536) / 2 ? limit * 2 + 65536 : std::numeric_limits<std::size_t>::max();
                for (std::size_t idx = 1; idx <= 10000; ++idx) {
                    const auto part_path = v4zip::detail::split_part_path(base, idx);
                    std::error_code ec;
                    if (!std::filesystem::is_regular_file(part_path, ec)) {
                        if (idx == 1) throw std::runtime_error("split part 1 not found");
                        break;
                    }
                    const auto part_data = read_file(part_path);
                    in_bytes += part_data.size();
                    if (archive.size() + part_data.size() > max_archive) throw std::runtime_error("archive exceeds size limit");
                    archive.insert(archive.end(), part_data.begin(), part_data.end());
                }
                output = v4zip::decompress_archive(archive, limit);
            } else if (v4zip::detail::is_manifest_path(paths[0])) {
                const auto man_data = read_file(paths[0]);
                in_bytes += man_data.size();
                const auto manifest = v4zip::detail::parse_manifest(man_data);
                if (manifest.total_size > limit) throw std::runtime_error("manifest total size exceeds limit");
                output.reserve(manifest.total_size);
                const auto parent_dir = std::filesystem::path(paths[0]).parent_path();
                for (const auto& entry : manifest.parts) {
                    const auto part_path = parent_dir.empty() ? entry.filename : (parent_dir / entry.filename).string();
                    const auto part_data = read_file(part_path);
                    in_bytes += part_data.size();
                    const auto chunk_raw = v4zip::decompress_archive(part_data, limit - output.size());
                    if (chunk_raw.size() != entry.raw_size) throw std::runtime_error("chunk raw size mismatch");
                    if (v4zip::detail::crc32(chunk_raw) != entry.raw_crc) throw std::runtime_error("chunk crc mismatch");
                    output.insert(output.end(), chunk_raw.begin(), chunk_raw.end());
                }
                if (output.size() != manifest.total_size) throw std::runtime_error("manifest total size mismatch");
                if (v4zip::detail::crc32(output) != manifest.total_crc) throw std::runtime_error("manifest total crc mismatch");
            } else {
                const auto input = read_file(paths[0]);
                in_bytes = input.size();
                output = v4zip::decompress_archive(input, limit);
            }
            out_bytes = output.size();
            destination.publish(output);
        }
        if (verbose) {
            const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            std::cerr << in_bytes << " input bytes -> " << out_bytes << " output bytes in " << elapsed << " seconds\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v4zip: " << error.what() << '\n';
        return 1;
    }
}
