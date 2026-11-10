#include "v4zip/archive.hpp"

#include <algorithm>
#include <array>
#include <cerrno>
#include <charconv>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <iostream>
#include <limits>
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
              << "Usage: v4zip [fastq] compress INPUT OUTPUT [--verbose]\n"
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
        std::size_t limit = default_limit;
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
            } else if (options && !value.empty() && value[0] == '-') throw std::runtime_error("unknown option: " + value);
            else paths.push_back(value);
        }
        if (paths.size() != 2) throw std::runtime_error("expected exactly INPUT and OUTPUT; use --help");
        Destination destination(paths[1]);
        const auto start = std::chrono::steady_clock::now();
        const auto input = read_file(paths[0]);
        const auto output = command == "compress" ? v4zip::compress_archive(input, fastq) : v4zip::decompress_archive(input, limit);
        destination.publish(output);
        if (verbose) {
            const auto elapsed = std::chrono::duration<double>(std::chrono::steady_clock::now() - start).count();
            std::cerr << input.size() << " input bytes -> " << output.size() << " output bytes in " << elapsed << " seconds\n";
        }
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "v4zip: " << error.what() << '\n';
        return 1;
    }
}
