#pragma once
#include "v4zip/dna.hpp"
namespace v4zip::detail {
Bytes encode_identifiers(const Bytes&, const std::vector<std::size_t>&);
Bytes decode_identifiers(const Bytes&, const std::vector<std::size_t>&, std::size_t);
}
