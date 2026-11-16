#pragma once
#include "hash.hpp"
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <vector>

namespace v4zip::detail {

enum class DnaStrategy : std::uint8_t {
    coding_small = 0,   // Compact bacterial / viral genomes (< 1.5M bases)
    coding_large = 1,   // Larger coding genomes (>= 1.5M bases)
    short_reads  = 2,   // Fragmented reads (FASTQ / short contigs)
    repeat_heavy = 3,   // High-repeat / multi-strain / pan-genome (>= 20% 16-mer repeats)
    low_complex  = 4    // Homopolymer / tandem-repeat dominated sequences
};

struct GenomeProfile {
    std::size_t length = 0;
    double gc_ratio = 0.5;
    double run_rate = 0.0;
    double codon_bias = 0.0;
    double repeat_16_rate = 0.0;
    double avg_record_len = 0.0;
    DnaStrategy strategy = DnaStrategy::coding_large;
};

inline GenomeProfile analyze_genome(const Bytes& dna, const DnaConfig& config) {
    GenomeProfile p;
    p.length = dna.size();
    if (dna.empty()) return p;

    const std::size_t n = dna.size();
    const std::size_t num_rec = std::max<std::size_t>(1, config.record_boundaries.size());
    p.avg_record_len = static_cast<double>(n) / static_cast<double>(num_rec);

    std::array<std::size_t, 4> counts{};
    std::array<std::array<std::size_t, 4>, 3> frame_counts{};
    std::size_t runs = 0;
    for (std::size_t i = 0; i < n; ++i) {
        const unsigned b = dna[i] & 3u;
        ++counts[b];
        ++frame_counts[i % 3][b];
        if (i >= 2 && dna[i] == dna[i - 1] && dna[i - 1] == dna[i - 2]) ++runs;
    }
    p.gc_ratio = static_cast<double>(counts[1] + counts[2]) / static_cast<double>(n);
    p.run_rate = static_cast<double>(runs) / static_cast<double>(n);

    // Measure 3-base codon periodicity via L1 divergence across reading frames
    double frame_div = 0.0;
    for (unsigned b = 0; b < 4; ++b) {
        const double mean_f = static_cast<double>(counts[b]) / static_cast<double>(n);
        for (unsigned f = 0; f < 3; ++f) {
            const double fn = static_cast<double>(frame_counts[f][0] + frame_counts[f][1] + frame_counts[f][2] + frame_counts[f][3]);
            if (fn > 0) frame_div += std::abs(static_cast<double>(frame_counts[f][b]) / fn - mean_f);
        }
    }
    p.codon_bias = frame_div;

    // Sampled 16-mer canonical repeat rate (stride = 4 for fast O(N/4) profiling)
    if (n >= 32) {
        constexpr std::size_t BLOOM_BITS = 1u << 23; // 8M bits = 1 MB
        std::vector<std::uint64_t> bloom(BLOOM_BITS >> 6, 0);
        std::size_t sampled = 0, hits = 0;
        const std::uint64_t mask16 = word_mask(16);
        std::uint64_t w = 0;
        for (std::size_t i = 0; i < n; ++i) {
            w = ((w << 2) | (dna[i] & 3u)) & mask16;
            if (i >= 15 && (i & 3u) == 0) {
                const auto canon = canonical_word(w, 16).word;
                const auto h = hash_word(canon) & (BLOOM_BITS - 1);
                const std::uint64_t bit = std::uint64_t{1} << (h & 63u);
                const std::size_t idx = h >> 6;
                if (bloom[idx] & bit) ++hits;
                else bloom[idx] |= bit;
                ++sampled;
            }
        }
        p.repeat_16_rate = sampled ? static_cast<double>(hits) / static_cast<double>(sampled) : 0.0;
    }

    if (p.run_rate >= 0.25) {
        p.strategy = DnaStrategy::low_complex;
    } else if (p.repeat_16_rate >= 0.30) {
        p.strategy = DnaStrategy::repeat_heavy;
    } else if (!config.record_boundaries.empty() && p.avg_record_len <= 250.0) {
        p.strategy = DnaStrategy::short_reads;
    } else if (n < 1500000) {
        p.strategy = DnaStrategy::coding_small;
    } else {
        p.strategy = DnaStrategy::coding_large;
    }
    return p;
}

} // namespace v4zip::detail
