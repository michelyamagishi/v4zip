/**
 * FastqCompressor.hpp - Multi-stream FASTQ compressor/decompressor
 *
 * Architecture:
 *   Stream 1: Identifiers  — rANS-coded token/delimiter/dict-ref (IdentifierCodecV2)
 *   Stream 2: DNA sequences — PPM* + TwoLayerMixer + V₄ orbit pooling (CompressorV4)
 *   Stream 3: Quality scores — position-bucketed PPM, sparse order-3, per-qual RLE (V5)
 *   Stream 4: Metadata      — delta-varint lengths, plus-line bitmask, non-DNA runs
 *
 * Binary format: .v4zq (V4ZIP-FASTQ), chunked for memory-efficient large files
 */

#ifndef V4ZIP_FASTQCOMPRESSOR_HPP
#define V4ZIP_FASTQCOMPRESSOR_HPP

#include "FastqIO.hpp"
#include "QualityModelV4.hpp"
#include "IdentifierCodec.hpp"
#include "CompressorV4.hpp"
#include "BinaryFormat.hpp"
#include "Alphabet.hpp"
#include "Varint.hpp"

#include <string>
#include <vector>
#include <cstdint>
#include <stdexcept>
#include <future>
#include <numeric>
#include <algorithm>

namespace v4zip {

// Non-DNA run record for lossless N character handling
struct NonDNARunFQ {
    uint64_t position;
    uint64_t length;
    char character;
};

// VGFQ file format constants
constexpr uint32_t MAGIC_VGFQ = 0x51464756;  // "VGFQ" in little-endian
constexpr uint8_t VERSION_VGFQ_MAJOR = 2;     // Bumped for chunked format
constexpr uint8_t VERSION_VGFQ_MINOR = 1;     // Added maxOrder/poolingOrder fields

// Feature flags (bitmask in file header)
constexpr uint16_t FLAG_VGFQ_HAS_PLUSLINE = 0x0001;    // Plus lines contain data beyond "+"
constexpr uint16_t FLAG_VGFQ_CHUNKED = 0x0002;         // Block-based chunked format
constexpr uint16_t FLAG_VGFQ_REORDERED = 0x0004;       // Reads sorted by minimizer within blocks
constexpr uint16_t FLAG_VGFQ_IDENTIFIER_V2 = 0x0010;   // rANS-coded identifier tokens
constexpr uint16_t FLAG_VGFQ_QUALITY_V6 = 0x0020;      // DNA-base-conditioned quality (32 bins)
constexpr uint16_t FLAG_VGFQ_QUALITY_V5 = 0x0040;      // Position-bucketed PPM + per-qual RLE
constexpr uint16_t FLAG_VGFQ_DNA_V2 = 0x0080;          // TwoLayerMixer + substitutional contexts
constexpr uint16_t FLAG_VGFQ_DNA_PERSIST = 0x0100;     // DNA model persists across chunks
constexpr uint16_t FLAG_VGFQ_DNA_CROSS_READ = 0x0200;  // DNA context carries across read boundaries
constexpr uint16_t FLAG_VGFQ_DNA_MATCH = 0x0400;       // Match model for repetitive sequences
constexpr uint16_t FLAG_VGFQ_DNA_V4SYM = 0x0800;       // V₄-symmetric enhancements
constexpr uint16_t FLAG_VGFQ_DNA_V6MODEL = 0x1000;    // SSE chain + context mixer + multi-match
constexpr uint16_t FLAG_VGFQ_DNA_COPYMODEL = 0x2000;  // Long-range V₄-canonical copy model

// Chunk configuration for memory-efficient processing
struct ChunkConfig {
    size_t maxRecordsPerBlock = 10000;   // ~1.5M bases typical
    size_t maxBasesPerBlock = 2000000;   // 2M base hard limit
};

// Block header for chunked format
struct BlockHeader {
    uint32_t recordCount;      // Number of records in block
    uint32_t totalBases;       // Total bases in block
    uint32_t identifierSize;   // Size of identifier stream
    uint32_t dnaSize;          // Size of DNA stream
    uint32_t qualitySize;      // Size of quality stream
    uint32_t metadataSize;     // Size of metadata stream
};

// Encoded block data
struct EncodedBlockData {
    std::vector<uint8_t> identifierData;
    std::vector<uint8_t> dnaData;
    std::vector<uint8_t> qualityData;
    std::vector<uint8_t> metadataData;
};

struct FastqCompressionStats {
    size_t recordCount;
    size_t totalBases;
    size_t totalQualityChars;
    size_t identifierBytes;
    size_t dnaBytes;
    size_t qualityBytes;
    size_t metadataBytes;
    size_t totalCompressedBytes;
    double bitsPerBase;
    double bitsPerQuality;
    double bitsPerIdentifierChar;
    double compressionRatio;
};

// ============================================================================
// Read reordering: minimizer-based block-local sorting for improved DNA compression
// ============================================================================

// Compute RC-canonical minimizer: the lexicographically smallest canonical k-mer
// in a DNA sequence. Canonical = min(kmer, reverse_complement(kmer)).
inline uint64_t rcCanonicalMinimizer(const std::string& seq, int k) {
    if (static_cast<int>(seq.size()) < k) return UINT64_MAX;

    uint64_t mask = (1ULL << (2 * k)) - 1;
    uint64_t fwd = 0, rev = 0;
    uint64_t best = UINT64_MAX;
    int valid = 0;

    for (size_t i = 0; i < seq.size(); i++) {
        int base = -1;
        switch (seq[i]) {
            case 'A': case 'a': base = 0; break;
            case 'C': case 'c': base = 1; break;
            case 'G': case 'g': base = 2; break;
            case 'T': case 't': base = 3; break;
        }

        if (base < 0) { valid = 0; fwd = rev = 0; continue; }

        fwd = ((fwd << 2) | base) & mask;
        rev = (rev >> 2) | (static_cast<uint64_t>(3 - base) << (2 * (k - 1)));
        valid++;

        if (valid >= k) {
            uint64_t canonical = std::min(fwd, rev);
            if (canonical < best) best = canonical;
        }
    }
    return best;
}

// Sort records by minimizer. Returns permutation[sortedPos] = originalIndex.
// Records are reordered in-place.
inline std::vector<uint32_t> sortByMinimizer(std::vector<FastqRecord>& records) {
    size_t n = records.size();
    if (n <= 1) return {};

    // Compute sort keys: primary k=15, secondary k=7
    struct SortKey { uint64_t primary; uint64_t secondary; };
    std::vector<SortKey> keys(n);
    for (size_t i = 0; i < n; i++) {
        keys[i].primary = rcCanonicalMinimizer(records[i].sequence, 15);
        keys[i].secondary = rcCanonicalMinimizer(records[i].sequence, 7);
    }

    // Sort indices by minimizer
    std::vector<uint32_t> indices(n);
    std::iota(indices.begin(), indices.end(), 0);
    std::sort(indices.begin(), indices.end(), [&](uint32_t a, uint32_t b) {
        if (keys[a].primary != keys[b].primary) return keys[a].primary < keys[b].primary;
        return keys[a].secondary < keys[b].secondary;
    });

    // Reorder records in-place using the permutation
    std::vector<FastqRecord> sorted(n);
    for (size_t i = 0; i < n; i++) {
        sorted[i] = std::move(records[indices[i]]);
    }
    records = std::move(sorted);

    return indices;  // permutation[sortedPos] = originalIndex
}

// Restore original order. permutation[sortedPos] = originalIndex.
inline void unsortRecords(std::vector<FastqRecord>& records,
                          const std::vector<uint32_t>& permutation) {
    size_t n = records.size();
    if (permutation.size() != n || n == 0) return;

    std::vector<FastqRecord> restored(n);
    for (size_t i = 0; i < n; i++) {
        restored[permutation[i]] = std::move(records[i]);
    }
    records = std::move(restored);
}

/**
 * FASTQ Compressor — multi-stream chunked architecture
 *
 * Splits FASTQ records into four independently-coded streams (identifiers,
 * DNA, quality, metadata), processes them in blocks for bounded memory usage,
 * and writes a .v4zq container with CRC32 integrity checking.
 */
class FastqCompressor {
public:
    explicit FastqCompressor(int maxOrder = 16, int poolingOrder = 16,
                             const ChunkConfig& config = ChunkConfig())
        : maxOrder_(maxOrder), poolingOrder_(poolingOrder), chunkConfig_(config) {}

    void setChunkConfig(const ChunkConfig& config) { chunkConfig_ = config; }
    const ChunkConfig& chunkConfig() const { return chunkConfig_; }

    // Compress from file using chunked format (memory-efficient)
    std::vector<uint8_t> compressFile(const std::string& path,
                                       FastqCompressionStats* stats = nullptr) {
        return compressFileChunked(path, stats);
    }

    // Compress from file using chunked format (memory-efficient)
    std::vector<uint8_t> compressFileChunked(const std::string& path,
                                              FastqCompressionStats* stats = nullptr) {
        FastqReader reader(path);

        BinaryWriter writer;

        // Write header (will update some fields later)
        writer.writeU32LE(MAGIC_VGFQ);
        writer.writeU8(VERSION_VGFQ_MAJOR);
        writer.writeU8(VERSION_VGFQ_MINOR);

        // Flags placeholder - will update after scanning for plus lines
        size_t flagsPos = writer.data().size();
        writer.writeU8(0);  // Flags low
        writer.writeU8(0);  // Flags high

        // DNA context model parameters (v2.1+)
        writer.writeU8(static_cast<uint8_t>(maxOrder_));
        writer.writeU8(static_cast<uint8_t>(poolingOrder_));

        // Record count and total bases placeholders
        size_t recordCountPos = writer.data().size();
        writer.writeU64LE(0);  // Record count placeholder
        writer.writeU64LE(0);  // Total bases placeholder

        // Block count placeholder
        size_t blockCountPos = writer.data().size();
        writer.writeU32LE(0);  // Block count placeholder

        // Process chunks
        std::vector<FastqRecord> chunk;
        std::vector<uint64_t> blockOffsets;
        size_t totalRecords = 0;
        size_t totalBases = 0;
        size_t totalIdentifierChars = 0;
        bool hasPlusLineContent = false;

        // Accumulate stats across all blocks
        size_t totalIdentifierBytes = 0;
        size_t totalDnaBytes = 0;
        size_t totalQualityBytes = 0;
        size_t totalMetadataBytes = 0;

        // Persistent DNA model across chunks (eliminates cold start)
        OnlineContextModel persistentDNAModel(maxOrder_, poolingOrder_);

        while (reader.readChunk(chunk, chunkConfig_.maxRecordsPerBlock, chunkConfig_.maxBasesPerBlock)) {
            // Record block start offset
            blockOffsets.push_back(writer.data().size());

            // Check for plus line content and count stats
            size_t blockBases = 0;
            bool blockHasPlusLine = false;
            for (const auto& record : chunk) {
                blockBases += record.sequence.size();
                totalIdentifierChars += record.identifier.size();
                if (!record.plusLine.empty()) {
                    blockHasPlusLine = true;
                    hasPlusLineContent = true;
                }
            }

            // Sort reads by minimizer for better DNA compression
            std::vector<uint32_t> permutation;
            if (chunk.size() > 1) {
                permutation = sortByMinimizer(chunk);
            }

            // Encode block with persistent DNA model
            auto blockData = encodeBlock(chunk, blockHasPlusLine, &persistentDNAModel,
                                          /*crossRead=*/true, permutation);

            // Write block header
            BlockHeader header;
            header.recordCount = static_cast<uint32_t>(chunk.size());
            header.totalBases = static_cast<uint32_t>(blockBases);
            header.identifierSize = static_cast<uint32_t>(blockData.identifierData.size());
            header.dnaSize = static_cast<uint32_t>(blockData.dnaData.size());
            header.qualitySize = static_cast<uint32_t>(blockData.qualityData.size());
            header.metadataSize = static_cast<uint32_t>(blockData.metadataData.size());

            writer.writeU32LE(header.recordCount);
            writer.writeU32LE(header.totalBases);
            writer.writeU32LE(header.identifierSize);
            writer.writeU32LE(header.dnaSize);
            writer.writeU32LE(header.qualitySize);
            writer.writeU32LE(header.metadataSize);

            // Write block data
            writer.writeBytes(blockData.identifierData);
            writer.writeBytes(blockData.dnaData);
            writer.writeBytes(blockData.qualityData);
            writer.writeBytes(blockData.metadataData);

            // Accumulate stats
            totalRecords += chunk.size();
            totalBases += blockBases;
            totalIdentifierBytes += blockData.identifierData.size();
            totalDnaBytes += blockData.dnaData.size();
            totalQualityBytes += blockData.qualityData.size();
            totalMetadataBytes += blockData.metadataData.size();

            // Free memory
            chunk.clear();
            chunk.shrink_to_fit();
        }

        // Update header fields
        auto& data = writer.data();

        // All current features always enabled
        uint16_t flags = FLAG_VGFQ_CHUNKED | FLAG_VGFQ_REORDERED
                       | FLAG_VGFQ_IDENTIFIER_V2 | FLAG_VGFQ_QUALITY_V5 | FLAG_VGFQ_QUALITY_V6
                       | FLAG_VGFQ_DNA_V2 | FLAG_VGFQ_DNA_PERSIST | FLAG_VGFQ_DNA_CROSS_READ
                       | FLAG_VGFQ_DNA_MATCH | FLAG_VGFQ_DNA_V4SYM | FLAG_VGFQ_DNA_V6MODEL
                       | FLAG_VGFQ_DNA_COPYMODEL;
        if (hasPlusLineContent) flags |= FLAG_VGFQ_HAS_PLUSLINE;
        data[flagsPos] = flags & 0xFF;
        data[flagsPos + 1] = (flags >> 8) & 0xFF;

        // Update record count
        for (int i = 0; i < 8; i++) {
            data[recordCountPos + i] = (totalRecords >> (i * 8)) & 0xFF;
        }

        // Update total bases
        for (int i = 0; i < 8; i++) {
            data[recordCountPos + 8 + i] = (totalBases >> (i * 8)) & 0xFF;
        }

        // Update block count
        uint32_t blockCount = static_cast<uint32_t>(blockOffsets.size());
        data[blockCountPos] = blockCount & 0xFF;
        data[blockCountPos + 1] = (blockCount >> 8) & 0xFF;
        data[blockCountPos + 2] = (blockCount >> 16) & 0xFF;
        data[blockCountPos + 3] = (blockCount >> 24) & 0xFF;

        // Write block offset index (for random access)
        for (uint64_t offset : blockOffsets) {
            writer.writeU64LE(offset);
        }

        // CRC32
        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());

        // Statistics
        if (stats) {
            stats->recordCount = totalRecords;
            stats->totalBases = totalBases;
            stats->totalQualityChars = totalBases;
            stats->identifierBytes = totalIdentifierBytes;
            stats->dnaBytes = totalDnaBytes;
            stats->qualityBytes = totalQualityBytes;
            stats->metadataBytes = totalMetadataBytes;
            stats->totalCompressedBytes = writer.data().size();

            // Estimate original size
            size_t originalSize = totalBases * 2  // sequence + quality
                                + totalIdentifierChars + totalRecords * 4;  // @\n+\n per record

            stats->bitsPerBase = totalBases > 0 ?
                (totalDnaBytes * 8.0) / totalBases : 0;
            stats->bitsPerQuality = totalBases > 0 ?
                (totalQualityBytes * 8.0) / totalBases : 0;
            stats->bitsPerIdentifierChar = totalIdentifierChars > 0 ?
                (totalIdentifierBytes * 8.0) / totalIdentifierChars : 0;
            stats->compressionRatio = writer.data().size() > 0 ?
                static_cast<double>(originalSize) / writer.data().size() : 0;
        }

        return writer.data();
    }

    // Compress records
    std::vector<uint8_t> compress(const std::vector<FastqRecord>& records,
                                   FastqCompressionStats* stats = nullptr) {
        if (records.empty()) {
            return createEmptyOutput(stats);
        }

        // Collect data for each stream
        std::vector<std::string> identifiers;
        std::vector<std::string> sequences;
        std::vector<std::string> plusLines;
        std::vector<std::string> qualities;
        std::vector<size_t> seqLengths;

        identifiers.reserve(records.size());
        sequences.reserve(records.size());
        plusLines.reserve(records.size());
        qualities.reserve(records.size());
        seqLengths.reserve(records.size());

        size_t totalBases = 0;
        size_t totalIdentifierChars = 0;
        bool hasPlusLineContent = false;

        for (const auto& record : records) {
            identifiers.push_back(record.identifier);
            sequences.push_back(record.sequence);
            plusLines.push_back(record.plusLine);
            qualities.push_back(record.quality);
            seqLengths.push_back(record.sequence.size());
            totalBases += record.sequence.size();
            totalIdentifierChars += record.identifier.size();
            if (!record.plusLine.empty()) {
                hasPlusLineContent = true;
            }
        }

        // Encode each stream
        std::vector<NonDNARunFQ> nonDNARuns;
        auto identifierData = encodeIdentifiers(identifiers);
        auto dnaData = encodeDNA(sequences, nonDNARuns);
        auto qualityData = encodeQualities(qualities, seqLengths);
        auto metadataData = encodeMetadata(seqLengths, plusLines, hasPlusLineContent, nonDNARuns);

        // Build output
        BinaryWriter writer;

        // Header
        writer.writeU32LE(MAGIC_VGFQ);
        writer.writeU8(VERSION_VGFQ_MAJOR);
        writer.writeU8(VERSION_VGFQ_MINOR);

        uint16_t flags = FLAG_VGFQ_IDENTIFIER_V2 | FLAG_VGFQ_QUALITY_V5
                       | FLAG_VGFQ_DNA_V2 | FLAG_VGFQ_DNA_MATCH | FLAG_VGFQ_DNA_V4SYM
                       | FLAG_VGFQ_DNA_V6MODEL | FLAG_VGFQ_DNA_COPYMODEL;
        if (hasPlusLineContent) flags |= FLAG_VGFQ_HAS_PLUSLINE;
        writer.writeU8(flags & 0xFF);
        writer.writeU8((flags >> 8) & 0xFF);

        // DNA context model parameters (v2.1+)
        writer.writeU8(static_cast<uint8_t>(maxOrder_));
        writer.writeU8(static_cast<uint8_t>(poolingOrder_));

        writer.writeU64LE(records.size());  // Record count
        writer.writeU64LE(totalBases);       // Total bases

        // Stream 1: Identifiers
        writer.writeU32LE(static_cast<uint32_t>(identifierData.size()));
        writer.writeBytes(identifierData);

        // Stream 2: DNA sequences
        writer.writeU32LE(static_cast<uint32_t>(dnaData.size()));
        writer.writeBytes(dnaData);

        // Stream 3: Quality scores
        writer.writeU32LE(static_cast<uint32_t>(qualityData.size()));
        writer.writeBytes(qualityData);

        // Stream 4: Metadata
        writer.writeU32LE(static_cast<uint32_t>(metadataData.size()));
        writer.writeBytes(metadataData);

        // CRC32
        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());

        // Statistics
        if (stats) {
            stats->recordCount = records.size();
            stats->totalBases = totalBases;
            stats->totalQualityChars = totalBases;  // Same as bases
            stats->identifierBytes = identifierData.size();
            stats->dnaBytes = dnaData.size();
            stats->qualityBytes = qualityData.size();
            stats->metadataBytes = metadataData.size();
            stats->totalCompressedBytes = writer.data().size();

            size_t originalSize = 0;
            for (const auto& r : records) {
                originalSize += 1 + r.identifier.size() + 1;  // @id\n
                originalSize += r.sequence.size() + 1;         // seq\n
                originalSize += 1 + r.plusLine.size() + 1;     // +plus\n
                originalSize += r.quality.size() + 1;          // qual\n
            }

            stats->bitsPerBase = totalBases > 0 ?
                (dnaData.size() * 8.0) / totalBases : 0;
            stats->bitsPerQuality = totalBases > 0 ?
                (qualityData.size() * 8.0) / totalBases : 0;
            stats->bitsPerIdentifierChar = totalIdentifierChars > 0 ?
                (identifierData.size() * 8.0) / totalIdentifierChars : 0;
            stats->compressionRatio = writer.data().size() > 0 ?
                static_cast<double>(originalSize) / writer.data().size() : 0;
        }

        return writer.data();
    }

    // Decompress to file (memory-efficient for chunked format)
    void decompressFile(const std::string& inputPath, const std::string& outputPath) {
        auto data = BinaryReader::readFile(inputPath);
        BinaryReader reader(data);

        // Read and validate header
        uint32_t magic = reader.readU32LE();
        if (magic != MAGIC_VGFQ) {
            throw std::runtime_error("Invalid VGFQ magic number");
        }

        reader.readU8();  // versionMajor
        reader.readU8();  // versionMinor

        uint16_t flags = reader.readU8() | (static_cast<uint16_t>(reader.readU8()) << 8);
        bool hasPlusLineContent = (flags & FLAG_VGFQ_HAS_PLUSLINE) != 0;
        bool isChunked = (flags & FLAG_VGFQ_CHUNKED) != 0;
        bool isReordered = (flags & FLAG_VGFQ_REORDERED) != 0;
        bool useDNAPersist = (flags & FLAG_VGFQ_DNA_PERSIST) != 0;
        bool useDNACrossRead = (flags & FLAG_VGFQ_DNA_CROSS_READ) != 0;
        bool useDNAMatch = (flags & FLAG_VGFQ_DNA_MATCH) != 0;
        bool useDNAV4Sym = (flags & FLAG_VGFQ_DNA_V4SYM) != 0;
        bool useDNAV6Model = (flags & FLAG_VGFQ_DNA_V6MODEL) != 0;
        bool useDNACopyModel = (flags & FLAG_VGFQ_DNA_COPYMODEL) != 0;
        int dnaModelVersion = useDNACopyModel ? 7 : (useDNAV6Model ? 6 : (useDNAV4Sym ? 5 : (useDNAMatch ? 4 : 3)));

        // Read DNA model parameters
        maxOrder_ = reader.readU8();
        poolingOrder_ = reader.readU8();

        uint64_t recordCount = reader.readU64LE();
        reader.readU64LE();  // totalBases (unused)

        // Validate CRC
        size_t dataEnd = data.size() - 4;
        CRC32 crc;
        crc.update(data.data(), dataEnd);
        uint32_t expectedCRC = BinaryReader(data.data() + dataEnd, 4).readU32LE();
        if (crc.finalize() != expectedCRC) {
            throw std::runtime_error("CRC mismatch");
        }

        FastqWriter writer(outputPath);

        if (isChunked) {
            uint32_t blockCount = reader.readU32LE();

            OnlineContextModel persistentDNAModel(maxOrder_, poolingOrder_, true, dnaModelVersion);
            OnlineContextModel* modelPtr = useDNAPersist ? &persistentDNAModel : nullptr;

            for (uint32_t b = 0; b < blockCount; b++) {
                auto blockRecords = decodeBlock(reader, hasPlusLineContent, modelPtr, useDNACrossRead, dnaModelVersion, isReordered);
                for (const auto& record : blockRecords) {
                    writer.write(record);
                }
            }
        } else {
            auto records = decompressV1(reader, recordCount, hasPlusLineContent, dnaModelVersion);
            writer.writeAll(records);
        }

        writer.flush();
    }

    // Decompress to records (in-memory)
    std::vector<FastqRecord> decompress(const std::vector<uint8_t>& data) {
        BinaryReader reader(data);

        uint32_t magic = reader.readU32LE();
        if (magic != MAGIC_VGFQ) {
            throw std::runtime_error("Invalid VGFQ magic number");
        }

        reader.readU8();  // versionMajor
        reader.readU8();  // versionMinor

        uint16_t flags = reader.readU8() | (static_cast<uint16_t>(reader.readU8()) << 8);
        bool hasPlusLineContent = (flags & FLAG_VGFQ_HAS_PLUSLINE) != 0;
        bool isChunked = (flags & FLAG_VGFQ_CHUNKED) != 0;
        bool isReordered = (flags & FLAG_VGFQ_REORDERED) != 0;
        bool useDNAPersist = (flags & FLAG_VGFQ_DNA_PERSIST) != 0;
        bool useDNACrossRead = (flags & FLAG_VGFQ_DNA_CROSS_READ) != 0;
        bool useDNAMatch = (flags & FLAG_VGFQ_DNA_MATCH) != 0;
        bool useDNAV4Sym = (flags & FLAG_VGFQ_DNA_V4SYM) != 0;
        bool useDNAV6Model = (flags & FLAG_VGFQ_DNA_V6MODEL) != 0;
        bool useDNACopyModel = (flags & FLAG_VGFQ_DNA_COPYMODEL) != 0;
        int dnaModelVersion = useDNACopyModel ? 7 : (useDNAV6Model ? 6 : (useDNAV4Sym ? 5 : (useDNAMatch ? 4 : 3)));

        maxOrder_ = reader.readU8();
        poolingOrder_ = reader.readU8();

        uint64_t recordCount = reader.readU64LE();
        reader.readU64LE();  // totalBases

        if (recordCount == 0) {
            return {};
        }

        // Validate CRC
        size_t dataEnd = data.size() - 4;
        CRC32 crc;
        crc.update(data.data(), dataEnd);
        uint32_t expectedCRC = BinaryReader(data.data() + dataEnd, 4).readU32LE();
        if (crc.finalize() != expectedCRC) {
            throw std::runtime_error("CRC mismatch");
        }

        if (isChunked) {
            return decompressChunked(reader, recordCount, hasPlusLineContent, useDNAPersist, useDNACrossRead, dnaModelVersion, isReordered);
        } else {
            return decompressV1(reader, recordCount, hasPlusLineContent, dnaModelVersion);
        }
    }

private:
    // Decompress single-block (non-chunked) format
    std::vector<FastqRecord> decompressV1(BinaryReader& reader, uint64_t recordCount,
                                           bool hasPlusLineContent, int dnaModelVersion) {
        uint32_t identifierSize = reader.readU32LE();
        if (identifierSize > reader.remaining()) {
            throw std::runtime_error("Identifier stream size exceeds remaining file data");
        }
        auto identifierData = reader.readBytes(identifierSize);

        uint32_t dnaSize = reader.readU32LE();
        if (dnaSize > reader.remaining()) {
            throw std::runtime_error("DNA stream size exceeds remaining file data");
        }
        auto dnaData = reader.readBytes(dnaSize);

        uint32_t qualitySize = reader.readU32LE();
        if (qualitySize > reader.remaining()) {
            throw std::runtime_error("Quality stream size exceeds remaining file data");
        }
        auto qualityData = reader.readBytes(qualitySize);

        uint32_t metadataSize = reader.readU32LE();
        if (metadataSize > reader.remaining()) {
            throw std::runtime_error("Metadata stream size exceeds remaining file data");
        }
        auto metadataData = reader.readBytes(metadataSize);

        std::vector<size_t> seqLengths;
        std::vector<std::string> plusLines;
        std::vector<NonDNARunFQ> nonDNARuns;
        std::vector<uint32_t> permutation;
        decodeMetadata(metadataData.data(), metadataData.size(),
                       recordCount, hasPlusLineContent, seqLengths, plusLines, nonDNARuns,
                       permutation, false);

        auto identifiers = decodeIdentifiers(identifierData.data(), identifierData.size());
        auto sequences = decodeDNA(dnaData.data(), dnaData.size(), seqLengths, nonDNARuns, nullptr, false, dnaModelVersion);
        auto qualities = decodeQualities(qualityData.data(), qualityData.size(), seqLengths);

        std::vector<FastqRecord> records;
        records.reserve(recordCount);
        for (size_t i = 0; i < recordCount; i++) {
            FastqRecord record;
            record.identifier = i < identifiers.size() ? identifiers[i] : "";
            record.sequence = i < sequences.size() ? sequences[i] : "";
            record.plusLine = i < plusLines.size() ? plusLines[i] : "";
            record.quality = i < qualities.size() ? qualities[i] : "";
            records.push_back(std::move(record));
        }
        return records;
    }

    // Decompress chunked format to vector of records
    std::vector<FastqRecord> decompressChunked(BinaryReader& reader, uint64_t recordCount,
                                                bool hasPlusLineContent,
                                                bool useDNAPersist, bool useDNACrossRead,
                                                int dnaModelVersion,
                                                bool isReordered = false) {
        uint32_t blockCount = reader.readU32LE();

        OnlineContextModel persistentDNAModel(maxOrder_, poolingOrder_, true, dnaModelVersion);
        OnlineContextModel* modelPtr = useDNAPersist ? &persistentDNAModel : nullptr;

        std::vector<FastqRecord> allRecords;
        allRecords.reserve(recordCount);

        for (uint32_t b = 0; b < blockCount; b++) {
            auto blockRecords = decodeBlock(reader, hasPlusLineContent, modelPtr, useDNACrossRead, dnaModelVersion, isReordered);
            for (auto& record : blockRecords) {
                allRecords.push_back(std::move(record));
            }
        }
        return allRecords;
    }

    // Encode a block of records
    EncodedBlockData encodeBlock(const std::vector<FastqRecord>& records, bool hasPlusLineContent,
                                  OnlineContextModel* persistentDNAModel = nullptr,
                                  bool crossRead = true,
                                  const std::vector<uint32_t>& permutation = {}) {
        EncodedBlockData result;

        if (records.empty()) {
            return result;
        }

        // Collect data for each stream
        std::vector<std::string> identifiers;
        std::vector<std::string> sequences;
        std::vector<std::string> plusLines;
        std::vector<std::string> qualities;
        std::vector<size_t> seqLengths;

        identifiers.reserve(records.size());
        sequences.reserve(records.size());
        plusLines.reserve(records.size());
        qualities.reserve(records.size());
        seqLengths.reserve(records.size());

        for (const auto& record : records) {
            identifiers.push_back(record.identifier);
            sequences.push_back(record.sequence);
            plusLines.push_back(record.plusLine);
            qualities.push_back(record.quality);
            seqLengths.push_back(record.sequence.size());
        }

        // Encode streams in parallel (identifiers + quality are independent of DNA)
        std::vector<NonDNARunFQ> nonDNARuns;

        // Launch identifiers and quality on async threads
        auto idFuture = std::async(std::launch::async,
            [&]{ return encodeIdentifiers(identifiers); });
        auto qualFuture = std::async(std::launch::async,
            [&]{ return encodeQualities(qualities, seqLengths, sequences); });

        // DNA runs on current thread (needed for nonDNARuns before metadata)
        result.dnaData = encodeDNA(sequences, nonDNARuns, persistentDNAModel, crossRead);

        // Metadata depends on nonDNARuns from DNA encoding
        result.metadataData = encodeMetadata(seqLengths, plusLines, hasPlusLineContent, nonDNARuns, permutation);

        // Collect async results
        result.identifierData = idFuture.get();
        result.qualityData = qualFuture.get();

        return result;
    }

    // Decode a single block
    std::vector<FastqRecord> decodeBlock(BinaryReader& reader, bool hasPlusLineContent,
                                          OnlineContextModel* persistentDNAModel = nullptr,
                                          bool crossRead = false, int dnaModelVersion = 5,
                                          bool isReordered = false) {
        uint32_t recordCount = reader.readU32LE();
        reader.readU32LE();  // totalBases
        uint32_t identifierSize = reader.readU32LE();
        uint32_t dnaSize = reader.readU32LE();
        uint32_t qualitySize = reader.readU32LE();
        uint32_t metadataSize = reader.readU32LE();

        // Validate stream sizes against remaining data
        uint64_t totalStreamSize = static_cast<uint64_t>(identifierSize) + dnaSize + qualitySize + metadataSize;
        if (totalStreamSize > reader.remaining()) {
            throw std::runtime_error("Block stream sizes exceed remaining file data");
        }

        auto identifierData = reader.readBytes(identifierSize);
        auto dnaData = reader.readBytes(dnaSize);
        auto qualityData = reader.readBytes(qualitySize);
        auto metadataData = reader.readBytes(metadataSize);

        std::vector<size_t> seqLengths;
        std::vector<std::string> plusLines;
        std::vector<NonDNARunFQ> nonDNARuns;
        std::vector<uint32_t> permutation;
        decodeMetadata(metadataData.data(), metadataData.size(),
                       recordCount, hasPlusLineContent, seqLengths, plusLines, nonDNARuns,
                       permutation, isReordered);

        auto identifiers = decodeIdentifiers(identifierData.data(), identifierData.size());
        auto sequences = decodeDNA(dnaData.data(), dnaData.size(), seqLengths, nonDNARuns, persistentDNAModel, crossRead, dnaModelVersion);
        auto qualities = decodeQualities(qualityData.data(), qualityData.size(), seqLengths, sequences);

        std::vector<FastqRecord> records;
        records.reserve(recordCount);
        for (size_t i = 0; i < recordCount; i++) {
            FastqRecord record;
            record.identifier = i < identifiers.size() ? identifiers[i] : "";
            record.sequence = i < sequences.size() ? sequences[i] : "";
            record.plusLine = i < plusLines.size() ? plusLines[i] : "";
            record.quality = i < qualities.size() ? qualities[i] : "";
            records.push_back(std::move(record));
        }

        // Restore original read order if reordered
        if (isReordered && !permutation.empty()) {
            unsortRecords(records, permutation);
        }

        return records;
    }

    // Encode identifiers using entropy-coded V2 codec
    std::vector<uint8_t> encodeIdentifiers(const std::vector<std::string>& identifiers) {
        IdentifierEncoderV2 encoder;
        return encoder.encode(identifiers);
    }

    // Decode identifiers using V2 entropy-coded codec
    std::vector<std::string> decodeIdentifiers(const uint8_t* data, size_t len) {
        IdentifierDecoderV2 decoder;
        return decoder.decode(data, len);
    }

    // Encode DNA sequences using adaptive context model
    // Returns encoded data and populates nonDNARuns with positions of N and other non-ACGT characters
    // If persistentModel is provided, uses it instead of creating a local model
    std::vector<uint8_t> encodeDNA(const std::vector<std::string>& sequences,
                                    std::vector<NonDNARunFQ>& nonDNARuns,
                                    OnlineContextModel* persistentModel = nullptr,
                                    bool crossRead = false) {
        nonDNARuns.clear();

        // Concatenate all sequences, tracking non-DNA positions
        std::string allDNA;
        size_t totalLen = 0;
        for (const auto& seq : sequences) {
            totalLen += seq.size();
        }
        allDNA.reserve(totalLen);

        size_t pos = 0;
        for (const auto& seq : sequences) {
            size_t i = 0;
            while (i < seq.size()) {
                char c = seq[i];
                if (isValidBase(c)) {
                    allDNA.push_back(toUpper(c));
                    pos++;
                    i++;
                } else {
                    // Start of non-DNA run
                    size_t runStart = pos;
                    char runChar = c;
                    while (i < seq.size() && seq[i] == runChar) {
                        allDNA.push_back('A');  // Placeholder
                        pos++;
                        i++;
                    }
                    nonDNARuns.push_back({runStart, pos - runStart, runChar});
                }
            }
        }

        if (allDNA.empty()) {
            return {};
        }

        // Build read boundaries (cumulative positions where each read starts)
        std::vector<size_t> readBoundaries;
        {
            size_t cumLen = 0;
            for (const auto& seq : sequences) {
                readBoundaries.push_back(cumLen);
                cumLen += seq.size();
            }
        }

        // Use adaptive context model with orbit pooling
        // If a persistent model is provided, use it; otherwise create a local one
        OnlineContextModel localModel(maxOrder_, poolingOrder_);
        OnlineContextModel& model = persistentModel ? *persistentModel : localModel;
        uint64_t fullMask = (1ULL << (2 * maxOrder_)) - 1;

        std::vector<int> symbols;
        std::vector<std::array<uint32_t, 5>> cdfs;
        symbols.reserve(allDNA.size());
        cdfs.reserve(allDNA.size());

        std::array<uint32_t, 5> uniformCDF = {
            0, PROB_SCALE/4, PROB_SCALE/2, 3*PROB_SCALE/4, PROB_SCALE
        };

        uint64_t ctx = 0;
        size_t nextBoundary = 0;
        int basesInRead = 0;

        for (size_t i = 0; i < allDNA.size(); i++) {
            // Check for read boundary (handle empty reads with while)
            while (nextBoundary < readBoundaries.size() && i == readBoundaries[nextBoundary]) {
                if (!crossRead) {
                    ctx = 0;
                    basesInRead = 0;
                    model.resetPerReadState();
                }
                nextBoundary++;
            }

            int base = charToBase(allDNA[i]);
            symbols.push_back(base);

            if (basesInRead == 0) {
                cdfs.push_back(uniformCDF);
            } else {
                int validOrder = std::min(basesInRead, maxOrder_);
                cdfs.push_back(model.getCDF(ctx, validOrder));
            }

            if (basesInRead > 0) {
                model.update(ctx, base, std::min(basesInRead - 1, maxOrder_));
            }

            ctx = ((ctx << 2) | base) & fullMask;
            basesInRead++;
        }

        return InterleavedRANSEncoder::encodeAll(symbols, cdfs);
    }

    // Helper to check if character is valid DNA base
    static bool isValidBase(char c) {
        switch (c) {
            case 'A': case 'a':
            case 'C': case 'c':
            case 'G': case 'g':
            case 'T': case 't':
                return true;
            default:
                return false;
        }
    }

    // Helper to convert to uppercase
    static char toUpper(char c) {
        if (c >= 'a' && c <= 'z') return c - 32;
        return c;
    }

    // Decode DNA sequences, restoring non-DNA characters from runs
    std::vector<std::string> decodeDNA(const uint8_t* data, size_t len,
                                        const std::vector<size_t>& lengths,
                                        const std::vector<NonDNARunFQ>& nonDNARuns,
                                        OnlineContextModel* persistentModel = nullptr,
                                        bool crossRead = false,
                                        int dnaModelVersion = 5) {
        size_t totalLen = 0;
        for (size_t l : lengths) {
            totalLen += l;
        }

        if (totalLen == 0 || len == 0) {
            std::vector<std::string> result(lengths.size());
            return result;
        }

        // Build read boundaries from cumulative lengths
        std::vector<size_t> readBoundaries;
        {
            size_t cumLen = 0;
            for (size_t l : lengths) {
                readBoundaries.push_back(cumLen);
                cumLen += l;
            }
        }

        // Decode all DNA as one stream
        std::string allDNA;
        allDNA.resize(totalLen);

        // Use persistent model if provided, otherwise create local
        OnlineContextModel localModel(maxOrder_, poolingOrder_, true, dnaModelVersion);
        OnlineContextModel& model = persistentModel ? *persistentModel : localModel;
        uint64_t fullMask = (1ULL << (2 * maxOrder_)) - 1;
        uint64_t ctx = 0;

        std::array<uint32_t, 5> uniformCDF = {
            0, PROB_SCALE/4, PROB_SCALE/2, 3*PROB_SCALE/4, PROB_SCALE
        };

        size_t nextBoundary = 0;
        int basesInRead = 0;

        InterleavedRANSDecoder::decodeIncremental(
            data, len, totalLen,
            [&](size_t i) -> std::array<uint32_t, 5> {
                // Check for read boundary (handle empty reads with while)
                while (nextBoundary < readBoundaries.size() && i == readBoundaries[nextBoundary]) {
                    if (!crossRead) {
                        ctx = 0;
                        basesInRead = 0;
                        model.resetPerReadState();
                    }
                    nextBoundary++;
                }

                if (basesInRead == 0) {
                    return uniformCDF;
                }
                int validOrder = std::min(basesInRead, maxOrder_);
                return model.getCDF(ctx, validOrder);
            },
            [&](size_t i, int base) {
                allDNA[i] = baseToChar(base);

                if (basesInRead > 0) {
                    model.update(ctx, base, std::min(basesInRead - 1, maxOrder_));
                }

                ctx = ((ctx << 2) | base) & fullMask;
                basesInRead++;
            }
        );

        // Restore non-DNA characters
        for (const auto& run : nonDNARuns) {
            for (uint64_t j = 0; j < run.length && run.position + j < allDNA.size(); j++) {
                allDNA[run.position + j] = run.character;
            }
        }

        // Split into separate sequences
        std::vector<std::string> results;
        results.reserve(lengths.size());
        size_t pos = 0;
        for (size_t l : lengths) {
            results.push_back(allDNA.substr(pos, l));
            pos += l;
        }

        return results;
    }

    // Encode quality scores using V5 model
    std::vector<uint8_t> encodeQualities(const std::vector<std::string>& qualities,
                                          const std::vector<size_t>& /*lengths*/,
                                          const std::vector<std::string>& sequences = {}) {
        QualityEncoderV5 encoder;
        if (!sequences.empty()) {
            return encoder.encodeAll(qualities, sequences);  // V6: DNA-conditioned
        }
        return encoder.encodeAll(qualities);  // V5: position-only
    }

    // Decode quality scores (dispatches V5 vs V6 based on format byte)
    std::vector<std::string> decodeQualities(const uint8_t* data, size_t len,
                                              const std::vector<size_t>& lengths,
                                              const std::vector<std::string>& sequences = {}) {
        QualityDecoderV5 decoder;
        if (len > 0 && data[0] == 0x04 && !sequences.empty()) {
            return decoder.decodeAll(data, len, lengths, sequences);  // V6
        }
        return decoder.decodeAll(data, len, lengths);  // V5
    }

    // Encode metadata (lengths + plus lines + non-DNA runs)
    std::vector<uint8_t> encodeMetadata(const std::vector<size_t>& lengths,
                                         const std::vector<std::string>& plusLines,
                                         bool hasPlusLineContent,
                                         const std::vector<NonDNARunFQ>& nonDNARuns,
                                         const std::vector<uint32_t>& permutation = {}) {
        std::vector<uint8_t> output;

        // Delta-encode lengths using zigzag varint (most reads have similar lengths)
        int64_t prevLen = 0;
        for (size_t len : lengths) {
            int64_t delta = static_cast<int64_t>(len) - prevLen;
            writeSignedVarint(output, delta);
            prevLen = static_cast<int64_t>(len);
        }

        // Plus lines: bitmask for presence + only non-empty content stored
        if (hasPlusLineContent) {
            // Write bitmask: 1 bit per record, packed into bytes
            size_t numBytes = (plusLines.size() + 7) / 8;
            std::vector<uint8_t> bitmask(numBytes, 0);
            for (size_t i = 0; i < plusLines.size(); i++) {
                if (!plusLines[i].empty()) {
                    bitmask[i / 8] |= (1 << (i % 8));
                }
            }
            output.insert(output.end(), bitmask.begin(), bitmask.end());

            // Write only non-empty plus lines
            for (const auto& plus : plusLines) {
                if (!plus.empty()) {
                    writeVarint(output, plus.size());
                    for (char c : plus) {
                        output.push_back(static_cast<uint8_t>(c));
                    }
                }
            }
        }

        // Encode non-DNA runs with delta encoding
        writeVarint(output, nonDNARuns.size());
        if (!nonDNARuns.empty()) {
            uint64_t prevPos = 0;
            for (const auto& run : nonDNARuns) {
                writeVarint(output, run.position - prevPos);
                writeVarint(output, run.length);
                output.push_back(static_cast<uint8_t>(run.character));
                prevPos = run.position + run.length;
            }
        }

        // Append permutation for read reordering (varint-encoded original indices)
        if (!permutation.empty()) {
            writeVarint(output, permutation.size());
            for (uint32_t idx : permutation) {
                writeVarint(output, idx);
            }
        }

        return output;
    }

    // Decode metadata
    void decodeMetadata(const uint8_t* data, size_t len,
                        size_t recordCount, bool hasPlusLineContent,
                        std::vector<size_t>& lengths,
                        std::vector<std::string>& plusLines,
                        std::vector<NonDNARunFQ>& nonDNARuns,
                        std::vector<uint32_t>& permutation,
                        bool isReordered = false) {
        const uint8_t* p = data;
        const uint8_t* end = data + len;

        // Read delta-encoded lengths (bounded varint for safety)
        lengths.clear();
        lengths.reserve(recordCount);
        int64_t prevLen = 0;
        for (size_t i = 0; i < recordCount && p < end; i++) {
            int64_t delta = readSignedVarint(p, end);
            prevLen += delta;
            lengths.push_back(static_cast<size_t>(prevLen));
        }

        // Read plus lines with bitmask
        plusLines.clear();
        plusLines.reserve(recordCount);
        if (hasPlusLineContent) {
            // Read bitmask
            size_t numBytes = (recordCount + 7) / 8;
            std::vector<uint8_t> bitmask(numBytes, 0);
            for (size_t i = 0; i < numBytes && p < end; i++) {
                bitmask[i] = *p++;
            }

            // Read non-empty plus lines
            for (size_t i = 0; i < recordCount; i++) {
                bool hasContent = (bitmask[i / 8] >> (i % 8)) & 1;
                if (hasContent && p < end) {
                    size_t plusLen = readVarint(p, end);
                    std::string plus;
                    plus.reserve(plusLen);
                    for (size_t j = 0; j < plusLen && p < end; j++) {
                        plus.push_back(static_cast<char>(*p++));
                    }
                    plusLines.push_back(std::move(plus));
                } else {
                    plusLines.push_back("");
                }
            }
        } else {
            for (size_t i = 0; i < recordCount; i++) {
                plusLines.push_back("");
            }
        }

        // Read non-DNA runs if present
        nonDNARuns.clear();
        if (p < end) {
            uint64_t runCount = readVarint(p, end);
            nonDNARuns.reserve(std::min(runCount, static_cast<uint64_t>(end - p)));
            uint64_t pos = 0;
            for (uint64_t i = 0; i < runCount && p < end; i++) {
                uint64_t delta = readVarint(p, end);
                pos += delta;
                uint64_t length = readVarint(p, end);
                char c = static_cast<char>(*p++);
                nonDNARuns.push_back({pos, length, c});
                pos += length;
            }
        }

        // Read permutation for read reordering
        permutation.clear();
        if (isReordered && p < end) {
            uint64_t permSize = readVarint(p, end);
            permutation.reserve(std::min(permSize, static_cast<uint64_t>(end - p)));
            for (uint64_t i = 0; i < permSize && p < end; i++) {
                permutation.push_back(static_cast<uint32_t>(readVarint(p, end)));
            }
        }
    }

    std::vector<uint8_t> createEmptyOutput(FastqCompressionStats* stats) {
        BinaryWriter writer;
        writer.writeU32LE(MAGIC_VGFQ);
        writer.writeU8(VERSION_VGFQ_MAJOR);
        writer.writeU8(VERSION_VGFQ_MINOR);
        writer.writeU8(0);  // Flags low
        writer.writeU8(0);  // Flags high
        writer.writeU64LE(0);  // Record count
        writer.writeU64LE(0);  // Total bases
        writer.writeU32LE(0);  // Identifier stream size
        writer.writeU32LE(0);  // DNA stream size
        writer.writeU32LE(0);  // Quality stream size
        writer.writeU32LE(0);  // Metadata stream size

        CRC32 crc;
        crc.update(writer.data());
        writer.writeU32LE(crc.finalize());

        if (stats) {
            *stats = {};
        }

        return writer.data();
    }

    static int charToBase(char c) {
        switch (c) {
            case 'A': case 'a': return 0;
            case 'C': case 'c': return 1;
            case 'G': case 'g': return 2;
            case 'T': case 't': return 3;
            default: return 0;  // N and others map to A
        }
    }

    static char baseToChar(int b) {
        return "ACGT"[b & 3];
    }

    int maxOrder_;
    int poolingOrder_;
    ChunkConfig chunkConfig_;
};

} // namespace v4zip

#endif // V4ZIP_FASTQCOMPRESSOR_HPP
