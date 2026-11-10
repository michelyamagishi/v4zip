#pragma once
#include "v4zip/dna.hpp"
namespace v4zip::detail {
struct DnaConfig {
    unsigned format_version=2;
    // Offsets in extracted ACGT, including duplicates for empty/ambiguous reads.
    std::vector<std::size_t> record_boundaries;
};
Bytes encode_dna(const Bytes&,const DnaConfig&);
Bytes decode_dna(const Bytes&,std::size_t,const DnaConfig&);
}
