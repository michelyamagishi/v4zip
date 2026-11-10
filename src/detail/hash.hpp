#pragma once
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

namespace v4zip::detail {
inline std::uint64_t hash_word(std::uint64_t word) {
    word^=word>>30;word*=UINT64_C(0xbf58476d1ce4e5b9);
    word^=word>>27;word*=UINT64_C(0x94d049bb133111eb);
    return word^(word>>31);
}

}
