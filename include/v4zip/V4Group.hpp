/**
 * V4Group.hpp - Klein four-group (V₄) operations on k-mers
 *
 * V₄ = {I, R, C, CR} where:
 *   I:  Identity
 *   R:  Reverse
 *   C:  Complement (XOR with mask)
 *   CR: Reverse-complement
 *
 * All elements are self-inverse: g(g(w)) = w
 * Group is abelian: g*h = h*g
 */

#ifndef V4ZIP_V4GROUP_HPP
#define V4ZIP_V4GROUP_HPP

#include <cstdint>
#include <utility>
#include <algorithm>

namespace v4zip {

enum class Transform : uint8_t {
    I  = 0,  // Identity
    R  = 1,  // Reverse
    C  = 2,  // Complement
    CR = 3   // Reverse-complement
};

// V₄ multiplication table (all elements are self-inverse)
// result[g][h] = g * h
constexpr Transform V4_MULTIPLY[4][4] = {
    // I    R    C    CR
    {Transform::I,  Transform::R,  Transform::C,  Transform::CR},  // I
    {Transform::R,  Transform::I,  Transform::CR, Transform::C },  // R
    {Transform::C,  Transform::CR, Transform::I,  Transform::R },  // C
    {Transform::CR, Transform::C,  Transform::R,  Transform::I }   // CR
};

constexpr Transform multiply(Transform g, Transform h) noexcept {
    return V4_MULTIPLY[static_cast<int>(g)][static_cast<int>(h)];
}

constexpr Transform inverse(Transform g) noexcept {
    return g; // All V₄ elements are self-inverse
}

// Reverse a k-mer (bit pairs reversed)
inline uint64_t reverse(uint64_t kmer, int k) noexcept {
    uint64_t result = 0;
    for (int i = 0; i < k; i++) {
        result = (result << 2) | (kmer & 3);
        kmer >>= 2;
    }
    return result;
}

// Complement a k-mer (XOR with all-1s mask of appropriate width)
inline uint64_t complement(uint64_t kmer, int k) noexcept {
    uint64_t mask = (1ULL << (2 * k)) - 1;
    return kmer ^ mask;
}

// Apply a transform to a k-mer
inline uint64_t apply(uint64_t kmer, int k, Transform h) noexcept {
    switch (h) {
        case Transform::I:  return kmer;
        case Transform::R:  return reverse(kmer, k);
        case Transform::C:  return complement(kmer, k);
        case Transform::CR: return reverse(complement(kmer, k), k);
    }
    return kmer;
}

// Get canonical representative (lex-min) and transform that maps canonical -> original
// Since all V₄ elements are self-inverse: if w = h(g), then g = h(w)
inline std::pair<uint64_t, Transform> canonicalWithTransform(uint64_t kmer, int k) noexcept {
    uint64_t variants[4] = {
        kmer,
        reverse(kmer, k),
        complement(kmer, k),
        reverse(complement(kmer, k), k)
    };

    int minIdx = 0;
    for (int i = 1; i < 4; i++) {
        if (variants[i] < variants[minIdx]) {
            minIdx = i;
        }
    }

    return {variants[minIdx], static_cast<Transform>(minIdx)};
}

inline uint64_t canonical(uint64_t kmer, int k) noexcept {
    return canonicalWithTransform(kmer, k).first;
}

// Check if k-mer is a palindrome (fixed under CR)
inline bool isPalindrome(uint64_t kmer, int k) noexcept {
    return kmer == reverse(complement(kmer, k), k);
}

// Check if k-mer is fixed under R (standard palindrome)
inline bool isReversePalindrome(uint64_t kmer, int k) noexcept {
    return kmer == reverse(kmer, k);
}

// Get orbit size (1, 2, or 4)
inline int orbitSize(uint64_t kmer, int k) noexcept {
    uint64_t r = reverse(kmer, k);
    uint64_t c = complement(kmer, k);
    uint64_t rc = reverse(c, k);

    int size = 1;
    if (r != kmer) size++;
    if (c != kmer && c != r) size++;
    if (rc != kmer && rc != r && rc != c) size++;
    return size;
}

// Get valid transforms for an orbit (depends on stabilizer)
// For size-4 orbits: all 4 transforms
// For size-2 orbits: only 2 transforms produce distinct results
inline int getValidTransforms(uint64_t canonical, int k, Transform* out) noexcept {
    int size = orbitSize(canonical, k);
    if (size == 4) {
        out[0] = Transform::I;
        out[1] = Transform::R;
        out[2] = Transform::C;
        out[3] = Transform::CR;
        return 4;
    } else if (size == 2) {
        // Find which transforms give distinct results
        out[0] = Transform::I;
        int count = 1;
        if (apply(canonical, k, Transform::R) != canonical) {
            out[count++] = Transform::R;
        } else if (apply(canonical, k, Transform::C) != canonical) {
            out[count++] = Transform::C;
        } else {
            out[count++] = Transform::CR;
        }
        return count;
    }
    // Size 1 - only identity (theoretically impossible for DNA k-mers, k>0)
    out[0] = Transform::I;
    return 1;
}

} // namespace v4zip

#endif // V4ZIP_V4GROUP_HPP
