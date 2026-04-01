/**
 * Types.hpp - Common types used across vgram
 */

#ifndef V4ZIP_TYPES_HPP
#define V4ZIP_TYPES_HPP

#include <cstddef>

namespace v4zip {

struct CompressionStats {
    size_t originalBytes;
    size_t compressedBytes;
    size_t numKmers;
    size_t numOrbits;
    double bitsPerBase;
    double compressionRatio;
};

} // namespace v4zip

#endif // V4ZIP_TYPES_HPP
