#pragma once
#include <array>
#include <cstddef>
namespace v4zip::detail {
// Format-2 constants. Changes after format freeze require a new archive version.
constexpr std::array<unsigned,3> context_table_bits{26,25,25};
constexpr std::array<unsigned,2> repeat_table_bits{24,23};
constexpr unsigned context_ways=4;
inline unsigned table_bits(unsigned word_length) {
    return word_length<=9?2*word_length:context_table_bits[word_length==13?0:(word_length==17?1:2)];
}
}
