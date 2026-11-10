#pragma once
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

namespace v4zip::detail {
constexpr std::uint32_t scale=65536;
using Distribution=std::array<std::uint32_t,4>;
using Cdf=std::array<std::uint32_t,5>;

class RangeEncoder {
public:
    void put(unsigned symbol,const Cdf& cdf) {
        const auto unit=range_>>16;
        low_+=std::uint64_t(unit)*cdf[symbol];
        range_=unit*(cdf[symbol+1]-cdf[symbol]);
        while(range_<(1u<<24)) {range_<<=8; shift();}
    }
    Bytes finish() {for(unsigned i=0;i<5;++i) shift(); return std::move(data_);}
private:
    void shift() {
        const auto low32=static_cast<std::uint32_t>(low_);
        const auto carry=static_cast<std::uint8_t>(low_>>32);
        if(low32<0xff000000u || carry!=0) {
            auto byte=cache_;
            do {data_.push_back(static_cast<std::uint8_t>(byte+carry));byte=255;} while(--pending_!=0);
            cache_=static_cast<std::uint8_t>(low32>>24);
        }
        ++pending_;
        low_=static_cast<std::uint32_t>(low32<<8);
    }
    Bytes data_;
    std::uint64_t low_=0;
    std::uint32_t range_=UINT32_MAX;
    std::size_t pending_=1;
    std::uint8_t cache_=0;
};

class RangeDecoder {
public:
    explicit RangeDecoder(const Bytes& data):data_(data) {
        if(data.size()<5 || data[0]!=0) throw std::runtime_error("Invalid DNA range stream");
        for(unsigned i=0;i<5;++i) code_=(code_<<8)|read();
    }
    unsigned get(const Cdf& cdf) {
        const auto unit=range_>>16;
        const auto value=code_/unit;
        if(value>=scale) throw std::runtime_error("Invalid DNA range code");
        unsigned symbol=(value>=cdf[1])+(value>=cdf[2])+(value>=cdf[3]);
        code_-=unit*cdf[symbol];
        range_=unit*(cdf[symbol+1]-cdf[symbol]);
        while(range_<(1u<<24)) {range_<<=8;code_=(code_<<8)|read();}
        return symbol;
    }
    void finish() const {
        if(position_!=data_.size()) throw std::runtime_error("Trailing DNA range bytes");
    }
private:
    std::uint8_t read() {
        if(position_==data_.size()) throw std::runtime_error("Truncated DNA range stream");
        return data_[position_++];
    }
    const Bytes& data_;
    std::size_t position_=0;
    std::uint32_t range_=UINT32_MAX,code_=0;
};

}
