/**
 * FastqIO.hpp - FASTQ file parsing and writing
 *
 * FASTQ format: 4-line records
 *   Line 1: @identifier (header)
 *   Line 2: DNA sequence
 *   Line 3: + (optionally followed by identifier)
 *   Line 4: Quality string (Phred+33 ASCII)
 */

#ifndef V4ZIP_FASTQIO_HPP
#define V4ZIP_FASTQIO_HPP

#include <string>
#include <vector>
#include <fstream>
#include <stdexcept>
#include <algorithm>

namespace v4zip {

struct FastqRecord {
    std::string identifier;    // Header line (without '@')
    std::string sequence;      // DNA sequence (ACGTN)
    std::string plusLine;      // Plus line content (may be empty, without '+')
    std::string quality;       // Quality string (Phred+33 ASCII)

    bool isValid() const {
        return !sequence.empty() && sequence.size() == quality.size();
    }

    void clear() {
        identifier.clear();
        sequence.clear();
        plusLine.clear();
        quality.clear();
    }
};

class FastqReader {
public:
    explicit FastqReader(const std::string& path) : file_(path), lineNum_(0) {
        if (!file_) {
            throw std::runtime_error("Cannot open FASTQ file: " + path);
        }
        path_ = path;
    }

    // Read next record, returns false at EOF
    bool next(FastqRecord& record) {
        record.clear();
        std::string line;

        // Line 1: @identifier
        if (!readLine(line)) return false;
        if (line.empty() || line[0] != '@') {
            throw std::runtime_error("FASTQ format error at line " +
                std::to_string(lineNum_) + ": expected '@' header");
        }
        record.identifier = line.substr(1);

        // Line 2: sequence
        if (!readLine(line)) {
            throw std::runtime_error("FASTQ format error: unexpected EOF after header at line " +
                std::to_string(lineNum_));
        }
        record.sequence = std::move(line);

        // Line 3: + (optionally followed by identifier)
        if (!readLine(line)) {
            throw std::runtime_error("FASTQ format error: unexpected EOF after sequence at line " +
                std::to_string(lineNum_));
        }
        if (line.empty() || line[0] != '+') {
            throw std::runtime_error("FASTQ format error at line " +
                std::to_string(lineNum_) + ": expected '+' separator");
        }
        record.plusLine = line.size() > 1 ? line.substr(1) : "";

        // Line 4: quality
        if (!readLine(line)) {
            throw std::runtime_error("FASTQ format error: unexpected EOF after '+' at line " +
                std::to_string(lineNum_));
        }
        record.quality = std::move(line);

        // Validate lengths match
        if (record.sequence.size() != record.quality.size()) {
            throw std::runtime_error("FASTQ format error at line " +
                std::to_string(lineNum_) + ": sequence length (" +
                std::to_string(record.sequence.size()) + ") != quality length (" +
                std::to_string(record.quality.size()) + ")");
        }

        return true;
    }

    // Read all records
    std::vector<FastqRecord> readAll() {
        std::vector<FastqRecord> records;
        FastqRecord record;
        while (next(record)) {
            records.push_back(std::move(record));
            record.clear();
        }
        return records;
    }

    // Get current line number (for error reporting)
    size_t lineNumber() const { return lineNum_; }

    // Read a chunk of records (up to maxRecords or maxBases, whichever comes first)
    // Returns true if any records were read, false at EOF
    // The records vector is cleared before reading
    bool readChunk(std::vector<FastqRecord>& records, size_t maxRecords, size_t maxBases) {
        records.clear();
        size_t totalBases = 0;

        FastqRecord record;
        while (records.size() < maxRecords && totalBases < maxBases) {
            if (!next(record)) {
                break;  // EOF
            }
            totalBases += record.sequence.size();
            records.push_back(std::move(record));
            record.clear();
        }

        return !records.empty();
    }

private:
    bool readLine(std::string& line) {
        if (!std::getline(file_, line)) {
            return false;
        }
        lineNum_++;
        // Strip carriage return if present (Windows line endings)
        if (!line.empty() && line.back() == '\r') {
            line.pop_back();
        }
        return true;
    }

    std::ifstream file_;
    std::string path_;
    size_t lineNum_;
};

class FastqWriter {
public:
    explicit FastqWriter(const std::string& path) : file_(path) {
        if (!file_) {
            throw std::runtime_error("Cannot create FASTQ file: " + path);
        }
    }

    void write(const FastqRecord& record) {
        file_ << '@' << record.identifier << '\n';
        file_ << record.sequence << '\n';
        file_ << '+' << record.plusLine << '\n';
        file_ << record.quality << '\n';
    }

    void writeAll(const std::vector<FastqRecord>& records) {
        for (const auto& record : records) {
            write(record);
        }
    }

    void flush() {
        file_.flush();
    }

private:
    std::ofstream file_;
};

// Utility: count records without loading all into memory
inline size_t countFastqRecords(const std::string& path) {
    FastqReader reader(path);
    FastqRecord record;
    size_t count = 0;
    while (reader.next(record)) {
        count++;
        record.clear();
    }
    return count;
}

// Utility: validate FASTQ file format
inline bool validateFastqFile(const std::string& path, std::string& error) {
    try {
        FastqReader reader(path);
        FastqRecord record;
        while (reader.next(record)) {
            // Validation happens in next()
        }
        return true;
    } catch (const std::exception& e) {
        error = e.what();
        return false;
    }
}

} // namespace v4zip

#endif // V4ZIP_FASTQIO_HPP
