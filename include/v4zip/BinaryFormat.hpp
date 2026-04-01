/**
 * BinaryFormat.hpp - .v4z file format specification and I/O
 *
 * V4ZIP FASTA compressed format (.v4z):
 * +----------------------+--------+----------------------------------+
 * | Field                | Size   | Description                      |
 * +----------------------+--------+----------------------------------+
 * | Magic "VGRM"         | 4 bytes| 0x56 0x47 0x52 0x4D (LE: 0x4D524756)|
 * | Version major        | 1 byte | Format version (currently 1)     |
 * | Version minor        | 1 byte | Minor version (currently 0)      |
 * | k (context order)    | 1 byte | Context model order (2-12)       |
 * | Flags                | 1 byte | Reserved (currently 0)           |
 * | Original length      | 8 bytes| Uncompressed size (LE uint64)    |
 * | Context model size   | 4 bytes| Context model stream (LE uint32) |
 * | Context model        | var    | Serialized context model         |
 * | Prob model size      | 4 bytes| Probability model (LE uint32)    |
 * | Prob model           | var    | Serialized prob model            |
 * | Grammar stream size  | 4 bytes| V4 grammar data (LE uint32)      |
 * | Grammar stream       | var    | rANS-encoded grammar symbols     |
 * | Symmetry stream size | 4 bytes| Symmetry bits (LE uint32)        |
 * | Symmetry stream      | var    | rANS-encoded symmetry flags      |
 * | CRC32                | 4 bytes| IEEE CRC32 of all preceding data |
 * +----------------------+--------+----------------------------------+
 *
 * Byte ordering: All multi-byte integers are little-endian (LE).
 */

#ifndef V4ZIP_BINARYFORMAT_HPP
#define V4ZIP_BINARYFORMAT_HPP

#include <cstdint>
#include <vector>
#include <string>
#include <fstream>
#include <stdexcept>
#include <cstring>

namespace v4zip {

constexpr uint32_t MAGIC = 0x4D524756;  // "VGRM" in little-endian
constexpr uint8_t VERSION_MAJOR = 1;
constexpr uint8_t VERSION_MINOR = 0;

// Flags
constexpr uint8_t FLAG_NONE = 0x00;

struct FileHeader {
    uint32_t magic;
    uint8_t versionMajor;
    uint8_t versionMinor;
    uint8_t k;
    uint8_t flags;
    uint64_t originalLength;
};

// CRC32 implementation (IEEE 802.3 polynomial)
static constexpr uint32_t CRC32_POLYNOMIAL = 0xEDB88320;  // Reversed bit-order

class CRC32 {
public:
    CRC32() : crc_(0xFFFFFFFF) {
        initTable();
    }

    void update(const uint8_t* data, size_t len) {
        for (size_t i = 0; i < len; i++) {
            crc_ = table_[(crc_ ^ data[i]) & 0xFF] ^ (crc_ >> 8);
        }
    }

    void update(const std::vector<uint8_t>& data) {
        update(data.data(), data.size());
    }

    uint32_t finalize() const {
        return crc_ ^ 0xFFFFFFFF;
    }

    void reset() {
        crc_ = 0xFFFFFFFF;
    }

private:
    void initTable() {
        for (uint32_t i = 0; i < 256; i++) {
            uint32_t c = i;
            for (int j = 0; j < 8; j++) {
                c = (c & 1) ? (CRC32_POLYNOMIAL ^ (c >> 1)) : (c >> 1);
            }
            table_[i] = c;
        }
    }

    uint32_t crc_;
    uint32_t table_[256];
};

class BinaryWriter {
public:
    void writeU8(uint8_t v) {
        data_.push_back(v);
    }

    void writeU32LE(uint32_t v) {
        data_.push_back(v & 0xFF);
        data_.push_back((v >> 8) & 0xFF);
        data_.push_back((v >> 16) & 0xFF);
        data_.push_back((v >> 24) & 0xFF);
    }

    void writeU64LE(uint64_t v) {
        for (int i = 0; i < 8; i++) {
            data_.push_back((v >> (i * 8)) & 0xFF);
        }
    }

    void writeBytes(const std::vector<uint8_t>& bytes) {
        data_.insert(data_.end(), bytes.begin(), bytes.end());
    }

    void writeBytes(const uint8_t* bytes, size_t len) {
        data_.insert(data_.end(), bytes, bytes + len);
    }

    const std::vector<uint8_t>& data() const { return data_; }
    std::vector<uint8_t>& data() { return data_; }

    void writeToFile(const std::string& path) const {
        std::ofstream f(path, std::ios::binary);
        if (!f) {
            throw std::runtime_error("Cannot create file: " + path);
        }
        f.write(reinterpret_cast<const char*>(data_.data()), data_.size());
    }

private:
    std::vector<uint8_t> data_;
};

class BinaryReader {
public:
    BinaryReader(const uint8_t* data, size_t len)
        : data_(data), len_(len), pos_(0) {}

    explicit BinaryReader(const std::vector<uint8_t>& data)
        : data_(data.data()), len_(data.size()), pos_(0) {}

    uint8_t readU8() {
        if (pos_ >= len_) throw std::runtime_error("Unexpected end of file");
        return data_[pos_++];
    }

    uint32_t readU32LE() {
        if (pos_ + 4 > len_) throw std::runtime_error("Unexpected end of file");
        uint32_t v = data_[pos_] |
                     (static_cast<uint32_t>(data_[pos_ + 1]) << 8) |
                     (static_cast<uint32_t>(data_[pos_ + 2]) << 16) |
                     (static_cast<uint32_t>(data_[pos_ + 3]) << 24);
        pos_ += 4;
        return v;
    }

    uint64_t readU64LE() {
        if (pos_ + 8 > len_) throw std::runtime_error("Unexpected end of file");
        uint64_t v = 0;
        for (int i = 0; i < 8; i++) {
            v |= static_cast<uint64_t>(data_[pos_ + i]) << (i * 8);
        }
        pos_ += 8;
        return v;
    }

    std::vector<uint8_t> readBytes(size_t count) {
        if (pos_ + count > len_) throw std::runtime_error("Unexpected end of file");
        std::vector<uint8_t> result(data_ + pos_, data_ + pos_ + count);
        pos_ += count;
        return result;
    }

    const uint8_t* currentPtr() const { return data_ + pos_; }
    void skip(size_t count) { pos_ += count; }
    size_t remaining() const { return len_ - pos_; }
    size_t position() const { return pos_; }

    static std::vector<uint8_t> readFile(const std::string& path) {
        std::ifstream f(path, std::ios::binary | std::ios::ate);
        if (!f) {
            throw std::runtime_error("Cannot open file: " + path);
        }
        size_t size = f.tellg();
        f.seekg(0);
        std::vector<uint8_t> data(size);
        f.read(reinterpret_cast<char*>(data.data()), size);
        return data;
    }

private:
    const uint8_t* data_;
    size_t len_;
    size_t pos_;
};

} // namespace v4zip

#endif // V4ZIP_BINARYFORMAT_HPP
