#pragma once
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

namespace v4zip::detail {
// Integer-only log/exp tables are part of the format. No host libm or floating
// contraction changes an adaptive prediction between encoder and decoder.
struct WeightTables {
    std::array<std::int16_t,65536> log_probability{};
    std::array<std::uint32_t,8193> weight{};
    WeightTables() {
        log_probability[0]=-4096;
        for(std::uint32_t p=1;p<65536;++p) {
            unsigned whole=0;
            for(auto q=p;q>1;q>>=1) ++whole;
            std::uint64_t x=std::uint64_t(p)<<(31-whole);
            unsigned fraction=0;
            for(unsigned bit=0;bit<8;++bit) {
                x=(x*x)>>31;
                fraction<<=1;
                if(x>=(UINT64_C(1)<<32)) {x>>=1;fraction|=1;}
            }
            log_probability[p]=static_cast<std::int16_t>((static_cast<int>(whole)-16)*256+static_cast<int>(fraction));
        }
        std::uint64_t value=UINT32_MAX;
        for(auto& w:weight) {
            w=static_cast<std::uint32_t>(value>>8);
            // round(2^(-1/256) * 2^32)
            value=(value*UINT64_C(4283353945))>>32;
        }
    }
};

inline const WeightTables& weight_tables() {static const WeightTables tables;return tables;}

}
