#pragma once

#include <cstdint>
#include <stdexcept>

namespace v4zip {

// I=0, R=1, C=2, RC=3. Composition is XOR; A=0,C=1,G=2,T=3.
inline std::uint64_t word_mask(unsigned length) {
    if (length>32) throw std::invalid_argument("V4 word exceeds 32 bases");
    return length==32 ? UINT64_MAX : (std::uint64_t{1}<<(2*length))-1;
}

inline std::uint64_t reverse_word(std::uint64_t word, unsigned length) {
    if (length>32) throw std::invalid_argument("V4 word exceeds 32 bases");
    if (length==0) return 0;
    word=((word>>2)&UINT64_C(0x3333333333333333))|((word&UINT64_C(0x3333333333333333))<<2);
    word=((word>>4)&UINT64_C(0x0f0f0f0f0f0f0f0f))|((word&UINT64_C(0x0f0f0f0f0f0f0f0f))<<4);
    word=((word>>8)&UINT64_C(0x00ff00ff00ff00ff))|((word&UINT64_C(0x00ff00ff00ff00ff))<<8);
    word=((word>>16)&UINT64_C(0x0000ffff0000ffff))|((word&UINT64_C(0x0000ffff0000ffff))<<16);
    word=(word>>32)|(word<<32);
    return word>>(64-2*length);
}

inline std::uint64_t transform_word(std::uint64_t word,unsigned length,unsigned transform) {
    const auto mask=word_mask(length);
    word&=mask;
    if (transform&1) word=reverse_word(word,length);
    if (transform&2) word^=mask;
    return word;
}

struct CanonicalWord { std::uint64_t word; unsigned transform; };

inline CanonicalWord canonical_word(std::uint64_t word,unsigned length) {
    const auto mask=word_mask(length);
    word&=mask;
    const auto reversed=reverse_word(word,length);
    CanonicalWord result{word,0};
    if(reversed<result.word) result={reversed,1};
    if((word^mask)<result.word) result={word^mask,2};
    if((reversed^mask)<result.word) result={reversed^mask,3};
    return result;
}

} // namespace v4zip
