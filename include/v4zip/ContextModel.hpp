/**
 * ContextModel.hpp - P(base|context) with orbit pooling
 *
 * Exploits Chargaff's second parity rule: contexts in the same
 * V₄ orbit have equivalent base frequency statistics.
 */

#ifndef V4ZIP_CONTEXTMODEL_HPP
#define V4ZIP_CONTEXTMODEL_HPP

#include "V4Group.hpp"
#include "Varint.hpp"
#include <vector>
#include <array>
#include <cstdint>

namespace v4zip {

// rANS parameters (shared across all codec components)
constexpr uint32_t PROB_BITS = 16;  // Increased from 14 for better precision on skewed distributions
constexpr uint32_t PROB_SCALE = 1u << PROB_BITS;
constexpr uint32_t RANS_L = 1u << 23;  // rANS normalization threshold

/**
 * Rescale adaptive counts: halve all counts (keeping minimum of 1) and recompute total.
 * Used by adaptive models to prevent count overflow while preserving relative frequencies.
 */
template<typename CountArray>
inline void rescaleCounts(CountArray& counts, int size, uint32_t& total) {
    total = 0;
    for (int i = 0; i < size; i++) {
        counts[i] = (counts[i] >> 1) | 1;
        total += counts[i];
    }
}

/**
 * Normalize observation counts to a CDF using the largest-remainder method
 * (Hamilton's method). Produces the optimal integer approximation minimizing
 * the maximum per-symbol rounding error.
 *
 * @param counts    Observation counts per symbol
 * @param N         Alphabet size (max 256)
 * @param targetSum Desired CDF total (e.g., PROB_SCALE)
 * @param cdf       Output array of N+1 elements: cdf[0]=0, cdf[N]=targetSum
 */
inline void normalizeCDF(const uint32_t* counts, int N, uint32_t targetSum, uint32_t* cdf) {
    uint64_t total = 0;
    for (int i = 0; i < N; i++) total += counts[i];

    if (total == 0) {
        cdf[0] = 0;
        for (int i = 1; i <= N; i++) {
            cdf[i] = static_cast<uint32_t>((static_cast<uint64_t>(i) * targetSum) / N);
        }
        return;
    }

    uint32_t probs[256];
    uint32_t frac[256];  // fractional remainder numerators (denominator = total)
    uint32_t sum = 0;

    for (int i = 0; i < N; i++) {
        if (counts[i] == 0) {
            probs[i] = 0;
            frac[i] = 0;
        } else {
            uint64_t scaled = static_cast<uint64_t>(counts[i]) * targetSum;
            uint32_t floored = static_cast<uint32_t>(scaled / total);
            uint32_t remainder = static_cast<uint32_t>(scaled % total);
            if (floored == 0) {
                probs[i] = 1;
                frac[i] = 0;  // already rounded up
            } else {
                probs[i] = floored;
                frac[i] = remainder;
            }
        }
        sum += probs[i];
    }

    // Distribute remaining units to symbols with largest fractional parts
    int32_t remaining = static_cast<int32_t>(targetSum) - static_cast<int32_t>(sum);
    while (remaining > 0) {
        uint32_t bestFrac = 0;
        int bestIdx = -1;
        for (int i = 0; i < N; i++) {
            if (frac[i] > bestFrac) {
                bestFrac = frac[i];
                bestIdx = i;
            }
        }
        if (bestIdx < 0) break;
        probs[bestIdx]++;
        frac[bestIdx] = 0;
        remaining--;
    }
    while (remaining < 0) {
        uint32_t bestFrac = UINT32_MAX;
        int bestIdx = -1;
        for (int i = 0; i < N; i++) {
            if (probs[i] > 1 && frac[i] < bestFrac) {
                bestFrac = frac[i];
                bestIdx = i;
            }
        }
        if (bestIdx < 0) break;
        probs[bestIdx]--;
        frac[bestIdx] = UINT32_MAX;
        remaining++;
    }

    cdf[0] = 0;
    for (int i = 0; i < N; i++) {
        cdf[i + 1] = cdf[i] + probs[i];
    }
    cdf[N] = targetSum;
}

// Convenience wrapper returning vector CDF
inline std::vector<uint32_t> normalizeCDFVec(const uint32_t* counts, int N, uint32_t targetSum) {
    std::vector<uint32_t> cdf(N + 1);
    normalizeCDF(counts, N, targetSum, cdf.data());
    return cdf;
}

class ContextModel {
public:
    explicit ContextModel(int k) : k_(k), ctxBits_(2 * (k - 1)) {
        uint32_t numContexts = 1u << ctxBits_;
        counts_.resize(numContexts);
        canonicalMap_.resize(numContexts);

        // Initialize with Laplace smoothing (alpha=1)
        // Apply orbit pooling for lower orders (Chargaff's rule):
        // DNA sequences exhibit reverse-complement symmetry. Contexts in the same
        // V4 orbit share equivalent base distributions. Pooling statistics improves
        // compression for lower orders with fewer observations.
        for (uint32_t ctx = 0; ctx < numContexts; ctx++) {
            canonicalMap_[ctx] = static_cast<uint32_t>(canonical(ctx, k - 1));
            counts_[ctx] = {1, 1, 1, 1};
        }
    }

    int k() const noexcept { return k_; }
    int contextBits() const noexcept { return ctxBits_; }

    void add(uint32_t ctx, int base) {
        if (base >= 0 && base < 4) {
            counts_[canonicalMap_[ctx]][base]++;
        }
    }

    // Get CDF for rANS encoding (5 elements: 0, cumsum[0..3], PROB_SCALE)
    std::array<uint32_t, 5> getCDF(uint32_t ctx) const {
        const auto& c = counts_[canonicalMap_[ctx]];
        std::array<uint32_t, 5> cdf;
        normalizeCDF(c.data(), 4, PROB_SCALE, cdf.data());
        return cdf;
    }

    // Serialize model (only canonical contexts)
    std::vector<uint8_t> serialize() const {
        std::vector<uint8_t> out;
        for (size_t ctx = 0; ctx < counts_.size(); ctx++) {
            if (canonicalMap_[ctx] == ctx) {
                for (int b = 0; b < 4; b++) {
                    writeVarint(out, counts_[ctx][b]);
                }
            }
        }
        return out;
    }

    // Deserialize model
    void deserialize(const uint8_t*& p) {
        for (size_t ctx = 0; ctx < counts_.size(); ctx++) {
            if (canonicalMap_[ctx] == ctx) {
                for (int b = 0; b < 4; b++) {
                    counts_[ctx][b] = readVarint(p);
                }
            }
        }
        // Copy to non-canonical contexts
        for (size_t ctx = 0; ctx < counts_.size(); ctx++) {
            if (canonicalMap_[ctx] != ctx) {
                counts_[ctx] = counts_[canonicalMap_[ctx]];
            }
        }
    }

    // Reset to initial state (Laplace prior)
    void reset() {
        for (auto& c : counts_) {
            c = {1, 1, 1, 1};
        }
    }

private:
    int k_;
    int ctxBits_;
    std::vector<std::array<uint32_t, 4>> counts_;
    std::vector<uint32_t> canonicalMap_;
};

} // namespace v4zip

#endif // V4ZIP_CONTEXTMODEL_HPP
