/**
 * Alphabet.hpp - DNA base encoding and complement operations
 *
 * Defines 2-bit encoding: A=0, C=1, G=2, T=3
 * Complement: b XOR 3 (A<->T, C<->G)
 */

#ifndef V4ZIP_ALPHABET_HPP
#define V4ZIP_ALPHABET_HPP

#include <cstdint>
#include <stdexcept>

namespace v4zip {

enum class Base : uint8_t {
    A = 0,
    C = 1,
    G = 2,
    T = 3
};

constexpr Base complement(Base b) noexcept {
    return static_cast<Base>(static_cast<uint8_t>(b) ^ 3);
}

constexpr char toChar(Base b) noexcept {
    constexpr char table[] = {'A', 'C', 'G', 'T'};
    return table[static_cast<uint8_t>(b)];
}

constexpr Base fromChar(char c) {
    switch (c) {
        case 'A': case 'a': return Base::A;
        case 'C': case 'c': return Base::C;
        case 'G': case 'g': return Base::G;
        case 'T': case 't': return Base::T;
        default:
            return Base::A; // Default for ambiguous bases
    }
}

constexpr bool isValidBase(char c) noexcept {
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

constexpr int baseToInt(Base b) noexcept {
    return static_cast<int>(b);
}

constexpr Base intToBase(int i) noexcept {
    return static_cast<Base>(i & 3);
}

/**
 * Convert character to base index (0-3)
 * Returns 0 for invalid characters (maps to A)
 */
constexpr int charToBase(char c) noexcept {
    switch (c) {
        case 'A': case 'a': return 0;
        case 'C': case 'c': return 1;
        case 'G': case 'g': return 2;
        case 'T': case 't': return 3;
        default: return 0;
    }
}

/**
 * Convert base index (0-3) to character
 */
constexpr char baseToChar(int b) noexcept {
    return "ACGT"[b & 3];
}

/**
 * Convert character to uppercase
 */
constexpr char toUpper(char c) noexcept {
    return (c >= 'a' && c <= 'z') ? (c - 32) : c;
}

} // namespace v4zip

#endif // V4ZIP_ALPHABET_HPP
