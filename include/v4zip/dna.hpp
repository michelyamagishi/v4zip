#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

namespace v4zip {
using Bytes=std::vector<std::uint8_t>;
// Input symbols use A=0,C=1,G=2,T=3. Raw streams use the current DNA format.
// The archive API selects compatibility decoding and supplies record boundaries.
Bytes encode_dna(const Bytes& bases);
Bytes decode_dna(const Bytes& encoded,std::size_t count);
} // namespace v4zip
