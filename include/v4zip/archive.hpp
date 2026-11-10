#pragma once

#include "v4zip/dna.hpp"
#include <cstddef>

namespace v4zip {
// FASTA (false) or standard four-line FASTQ (true), preserving every input byte.
Bytes compress_archive(const Bytes& input, bool fastq);
// Rejects malformed archives and declared outputs exceeding max_output.
Bytes decompress_archive(const Bytes& archive, std::size_t max_output);
}
