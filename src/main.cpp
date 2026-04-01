/**
 * V4ZIP CLI — V₄ Group-Based DNA/FASTQ Compressor
 *
 * Usage:
 *   v4zip compress [-o N] [--verbose] input.fasta output.v4z
 *   v4zip decompress input.v4z output.fasta
 *   v4zip fastq compress [-o N] [--block-size N] [--verbose] [--stats] input.fastq output.v4zq
 *   v4zip fastq decompress input.v4zq output.fastq
 */

#include <v4zip/CompressorV4.hpp>
#include <v4zip/BinaryFormat.hpp>
#include <v4zip/FastaIO.hpp>
#include <v4zip/FastqCompressor.hpp>

#include <iostream>
#include <fstream>
#include <chrono>
#include <cstring>
#include <iomanip>

void printUsage(const char* prog) {
    std::cerr << "V4ZIP - V4 Group-Based DNA Compressor\n\n";
    std::cerr << "Usage:\n";
    std::cerr << "  " << prog << " compress [options] input.fasta output.v4z\n";
    std::cerr << "  " << prog << " decompress input.v4z output.fasta\n";
    std::cerr << "  " << prog << " fastq compress [options] input.fastq output.v4zq\n";
    std::cerr << "  " << prog << " fastq decompress input.v4zq output.fastq\n";
    std::cerr << "\nOptions:\n";
    std::cerr << "  -o N            Context order (default: 16)\n";
    std::cerr << "  --block-size N  Records per block for chunked compression (default: 10000)\n";
    std::cerr << "  --verbose       Print compression statistics\n";
    std::cerr << "  --stats         Show detailed stream statistics (fastq only)\n";
    std::cerr << "  -h, --help      Show this help message\n";
}

int doCompress(int argc, char** argv) {
    int order = 16;
    bool verbose = false;
    std::string inputPath;
    std::string outputPath;

    for (int i = 2; i < argc; i++) {
        if ((strcmp(argv[i], "-o") == 0 || strcmp(argv[i], "-k") == 0) && i + 1 < argc) {
            order = std::atoi(argv[++i]);
            if (order < 2 || order > 20) {
                std::cerr << "Error: order must be between 2 and 20\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
            verbose = true;
        } else if (inputPath.empty()) {
            inputPath = argv[i];
        } else {
            outputPath = argv[i];
        }
    }

    if (inputPath.empty() || outputPath.empty()) {
        std::cerr << "Error: Missing input or output file\n";
        printUsage(argv[0]);
        return 1;
    }

    try {
        auto startTime = std::chrono::high_resolution_clock::now();

        v4zip::CompressionStats stats;
        v4zip::CompressorV4 compressor(order, order);
        auto compressed = compressor.compressFile(inputPath, &stats);

        std::ofstream ofs(outputPath, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("Cannot open output file: " + outputPath);
        }
        ofs.write(reinterpret_cast<const char*>(compressed.data()), compressed.size());

        auto endTime = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(endTime - startTime).count();

        if (verbose) {
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "Input:       " << inputPath << "\n";
            std::cout << "Output:      " << outputPath << "\n";
            std::cout << "Order:       " << order << "\n";
            std::cout << "Original:    " << stats.originalBytes << " bytes\n";
            std::cout << "Compressed:  " << stats.compressedBytes << " bytes\n";
            std::cout << "Bits/base:   " << stats.bitsPerBase << "\n";
            std::cout << "Ratio:       " << stats.compressionRatio << "x\n";
            std::cout << "Time:        " << elapsed << " s\n";
        } else {
            std::cout << stats.originalBytes << " -> " << stats.compressedBytes
                      << " bytes (" << stats.bitsPerBase << " bpb, "
                      << std::setprecision(2) << stats.compressionRatio << "x) in "
                      << std::setprecision(3) << elapsed << "s\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

int doDecompress(int argc, char** argv) {
    std::string inputPath;
    std::string outputPath;

    for (int i = 2; i < argc; i++) {
        if (inputPath.empty()) {
            inputPath = argv[i];
        } else {
            outputPath = argv[i];
        }
    }

    if (inputPath.empty() || outputPath.empty()) {
        std::cerr << "Error: Missing input or output file\n";
        printUsage(argv[0]);
        return 1;
    }

    try {
        auto startTime = std::chrono::high_resolution_clock::now();

        auto inputData = v4zip::BinaryReader::readFile(inputPath);

        v4zip::CompressorV4 compressor;
        auto records = compressor.decompressToRecords(inputData);

        v4zip::FastaWriter writer(outputPath);
        writer.writeAll(records);

        size_t totalLen = 0;
        for (const auto& rec : records) {
            totalLen += rec.sequence.size();
        }

        auto endTime = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(endTime - startTime).count();

        std::cout << inputData.size() << " bytes -> " << totalLen
                  << " bases (" << records.size() << " records) in " << elapsed << "s\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

int doFastqCompress(int argc, char** argv) {
    int order = 16;
    int blockSize = 10000;
    bool verbose = false;
    bool stats = false;
    std::string inputPath;
    std::string outputPath;

    for (int i = 3; i < argc; i++) {
        if ((strcmp(argv[i], "-o") == 0) && i + 1 < argc) {
            order = std::atoi(argv[++i]);
            if (order < 2 || order > 20) {
                std::cerr << "Error: order must be between 2 and 20\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--block-size") == 0 && i + 1 < argc) {
            blockSize = std::atoi(argv[++i]);
            if (blockSize < 100 || blockSize > 1000000) {
                std::cerr << "Error: block-size must be between 100 and 1000000\n";
                return 1;
            }
        } else if (strcmp(argv[i], "--verbose") == 0 || strcmp(argv[i], "-v") == 0) {
            verbose = true;
        } else if (strcmp(argv[i], "--stats") == 0) {
            stats = true;
        } else if (inputPath.empty()) {
            inputPath = argv[i];
        } else {
            outputPath = argv[i];
        }
    }

    if (inputPath.empty() || outputPath.empty()) {
        std::cerr << "Error: Missing input or output file\n";
        printUsage(argv[0]);
        return 1;
    }

    try {
        auto startTime = std::chrono::high_resolution_clock::now();

        v4zip::ChunkConfig config;
        config.maxRecordsPerBlock = blockSize;
        config.maxBasesPerBlock = blockSize * 200;

        v4zip::FastqCompressionStats compStats;
        v4zip::FastqCompressor compressor(order, order, config);
        auto compressed = compressor.compressFile(inputPath, &compStats);

        std::ofstream ofs(outputPath, std::ios::binary);
        if (!ofs) {
            throw std::runtime_error("Cannot open output file: " + outputPath);
        }
        ofs.write(reinterpret_cast<const char*>(compressed.data()), compressed.size());

        auto endTime = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(endTime - startTime).count();

        if (verbose || stats) {
            std::cout << std::fixed << std::setprecision(4);
            std::cout << "Input:         " << inputPath << "\n";
            std::cout << "Output:        " << outputPath << "\n";
            std::cout << "Order:         " << order << "\n";
            std::cout << "Records:       " << compStats.recordCount << "\n";
            std::cout << "Total bases:   " << compStats.totalBases << "\n";
            std::cout << "Compressed:    " << compStats.totalCompressedBytes << " bytes\n";
            std::cout << "Ratio:         " << compStats.compressionRatio << "x\n";
            std::cout << "Time:          " << elapsed << " s\n";

            if (stats) {
                std::cout << "\nStream breakdown:\n";
                std::cout << "  Identifiers: " << compStats.identifierBytes << " bytes ("
                          << compStats.bitsPerIdentifierChar << " bits/char)\n";
                std::cout << "  DNA:         " << compStats.dnaBytes << " bytes ("
                          << compStats.bitsPerBase << " bits/base)\n";
                std::cout << "  Quality:     " << compStats.qualityBytes << " bytes ("
                          << compStats.bitsPerQuality << " bits/qual)\n";
                std::cout << "  Metadata:    " << compStats.metadataBytes << " bytes\n";
            }
        } else {
            std::cout << compStats.recordCount << " records, "
                      << compStats.totalBases << " bases -> "
                      << compStats.totalCompressedBytes << " bytes ("
                      << std::setprecision(2) << compStats.compressionRatio << "x) in "
                      << std::setprecision(3) << elapsed << "s\n";
        }

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

int doFastqDecompress(int argc, char** argv) {
    std::string inputPath;
    std::string outputPath;

    for (int i = 3; i < argc; i++) {
        if (inputPath.empty()) {
            inputPath = argv[i];
        } else {
            outputPath = argv[i];
        }
    }

    if (inputPath.empty() || outputPath.empty()) {
        std::cerr << "Error: Missing input or output file\n";
        printUsage(argv[0]);
        return 1;
    }

    try {
        auto startTime = std::chrono::high_resolution_clock::now();

        v4zip::FastqCompressor compressor;
        compressor.decompressFile(inputPath, outputPath);

        auto inputData = v4zip::BinaryReader::readFile(inputPath);

        auto endTime = std::chrono::high_resolution_clock::now();
        double elapsed = std::chrono::duration<double>(endTime - startTime).count();

        std::cout << inputData.size() << " bytes -> FASTQ in " << elapsed << "s\n";

        return 0;
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }
}

int doFastq(int argc, char** argv) {
    if (argc < 3) {
        std::cerr << "Error: Missing fastq subcommand\n";
        printUsage(argv[0]);
        return 1;
    }

    const char* subcmd = argv[2];

    if (strcmp(subcmd, "compress") == 0 || strcmp(subcmd, "c") == 0) {
        return doFastqCompress(argc, argv);
    }

    if (strcmp(subcmd, "decompress") == 0 || strcmp(subcmd, "d") == 0) {
        return doFastqDecompress(argc, argv);
    }

    std::cerr << "Unknown fastq subcommand: " << subcmd << "\n";
    printUsage(argv[0]);
    return 1;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        printUsage(argv[0]);
        return 1;
    }

    const char* cmd = argv[1];

    if (strcmp(cmd, "-h") == 0 || strcmp(cmd, "--help") == 0) {
        printUsage(argv[0]);
        return 0;
    }

    if (strcmp(cmd, "compress") == 0 || strcmp(cmd, "c") == 0) {
        return doCompress(argc, argv);
    }

    if (strcmp(cmd, "decompress") == 0 || strcmp(cmd, "d") == 0) {
        return doDecompress(argc, argv);
    }

    if (strcmp(cmd, "fastq") == 0 || strcmp(cmd, "fq") == 0) {
        return doFastq(argc, argv);
    }

    std::cerr << "Unknown command: " << cmd << "\n";
    printUsage(argv[0]);
    return 1;
}
