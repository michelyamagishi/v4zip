/**
 * FastaIO.hpp - FASTA file parsing and writing
 */

#ifndef V4ZIP_FASTAIO_HPP
#define V4ZIP_FASTAIO_HPP

#include "Alphabet.hpp"
#include <string>
#include <vector>
#include <fstream>
#include <sstream>
#include <stdexcept>

namespace v4zip {

struct FastaRecord {
    std::string header;     // Header line (without '>')
    std::string sequence;   // Pure ACGT sequence
    std::vector<std::pair<size_t, char>> ambiguous;  // (position, original char)
};

class FastaReader {
public:
    explicit FastaReader(const std::string& path) : file_(path) {
        if (!file_) {
            throw std::runtime_error("Cannot open FASTA file: " + path);
        }
    }

    std::vector<FastaRecord> readAll() {
        std::vector<FastaRecord> records;
        FastaRecord current;
        std::string line;

        while (std::getline(file_, line)) {
            if (line.empty()) continue;

            if (line[0] == '>') {
                if (!current.sequence.empty()) {
                    records.push_back(std::move(current));
                    current = FastaRecord();
                }
                current.header = line.substr(1);
            } else {
                // Process sequence line
                for (char c : line) {
                    if (c == '\r' || c == '\n') continue;

                    if (isValidBase(c)) {
                        current.sequence.push_back(toChar(fromChar(c)));
                    } else {
                        // Store ambiguous base position and replace with 'A'
                        current.ambiguous.emplace_back(current.sequence.size(), c);
                        current.sequence.push_back('A');
                    }
                }
            }
        }

        if (!current.sequence.empty()) {
            records.push_back(std::move(current));
        }

        return records;
    }

    // Read only the pure DNA sequence (replace ambiguous with A)
    static std::string readSequenceOnly(const std::string& path) {
        std::ifstream file(path);
        if (!file) {
            throw std::runtime_error("Cannot open FASTA file: " + path);
        }

        std::string sequence;
        std::string line;

        while (std::getline(file, line)) {
            if (line.empty() || line[0] == '>') continue;

            for (char c : line) {
                if (c == '\r' || c == '\n') continue;
                if (isValidBase(c)) {
                    sequence.push_back(toChar(fromChar(c)));
                } else {
                    sequence.push_back('A');
                }
            }
        }

        return sequence;
    }

private:
    std::ifstream file_;
};

class FastaWriter {
public:
    explicit FastaWriter(const std::string& path) : file_(path) {
        if (!file_) {
            throw std::runtime_error("Cannot create FASTA file: " + path);
        }
    }

    void write(const FastaRecord& record, int lineWidth = 80) {
        file_ << '>' << record.header << '\n';

        // Restore ambiguous bases
        std::string seq = record.sequence;
        for (const auto& [pos, c] : record.ambiguous) {
            if (pos < seq.size()) {
                seq[pos] = c;
            }
        }

        // Write sequence with line wrapping
        for (size_t i = 0; i < seq.size(); i += lineWidth) {
            file_ << seq.substr(i, lineWidth) << '\n';
        }
    }

    void write(const std::string& header, const std::string& sequence, int lineWidth = 80) {
        file_ << '>' << header << '\n';
        for (size_t i = 0; i < sequence.size(); i += lineWidth) {
            file_ << sequence.substr(i, lineWidth) << '\n';
        }
    }

    void writeAll(const std::vector<FastaRecord>& records, int lineWidth = 80) {
        for (const auto& record : records) {
            write(record, lineWidth);
        }
    }

private:
    std::ofstream file_;
};

} // namespace v4zip

#endif // V4ZIP_FASTAIO_HPP
