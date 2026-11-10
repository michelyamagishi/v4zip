#pragma once
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

#include "hash.hpp"
namespace v4zip::detail {
struct Row {
    std::uint64_t tag=UINT64_MAX;
    std::array<std::uint32_t,4> counts{};
};

// One row per canonical extended word. Separate orientation counts retain
// directed evidence; reverse-complement and full-orbit experts pool *words*,
// not mismatched left-context/right-extension conditionals.
class WordTable {
public:
#ifdef V4ZIP_DIAGNOSTICS
    std::uint64_t replacements=0;
    std::size_t heap_bytes() const {return rows_.capacity()*sizeof(Row);}
#endif
    explicit WordTable(unsigned length,unsigned ways=1,unsigned bits=20):
        length_(length),ways_(length<=9?1:ways),mask_(word_mask(length)),
        rows_(std::size_t{1}<<std::min(2*length,bits)),index_mask_(rows_.size()-1),dense_(2*length<=bits) {}

    std::array<std::uint64_t,3> counts(std::uint64_t word) {
        word&=mask_;
        const auto canonical=canonical_word(word,length_);
        const auto* found=find(canonical.word,false);
        if(!found)return {};
        const auto& row=*found;
        const unsigned t=canonical.transform;
        const auto r=reverse_word(word,length_);
        const bool reverse_fixed=(r==word);
        const bool rc_fixed=((r^mask_)==word);
        const auto raw=std::uint64_t(row.counts[t]);
        const auto rc=rc_fixed ? raw : raw+row.counts[t^3];
        std::uint64_t pooled=0;
        for(auto c:row.counts) pooled+=c;
        if(reverse_fixed || rc_fixed) pooled/=2;
        return {raw,rc,pooled};
    }

    void observe(std::uint64_t word) {
        word&=mask_;
        const auto canonical=canonical_word(word,length_);
        auto& row=*find(canonical.word,true);
        const auto t=canonical.transform;
        // A conditional compares several word rows. Aging rows independently
        // destroys their relative frequencies; retain wide observation counts.
        // Saturation is synchronized across a word's equivalent orientations.
        if(row.counts[t]!=UINT32_MAX) ++row.counts[t];
        const auto r=reverse_word(word,length_);
        // A stabilizer gives multiple orientation labels to one observation.
        // Replicate its count; the orbit sum divides out that multiplicity.
        if(r==word && row.counts[t^1]!=UINT32_MAX) ++row.counts[t^1];
        else if((r^mask_)==word && row.counts[t^3]!=UINT32_MAX) ++row.counts[t^3];
    }
private:
    Row* find(std::uint64_t key,bool insert) {
        const auto first=index(key)&~std::size_t(ways_-1);
        unsigned found=ways_;
        for(unsigned i=0;i<ways_;++i)if(rows_[first+i].tag==key){found=i;break;}
        if(found==ways_) {
            if(!insert)return nullptr;
            found=ways_-1;
#ifdef V4ZIP_DIAGNOSTICS
            replacements+=rows_[first+found].tag!=UINT64_MAX;
#endif
            rows_[first+found]={key,{}};
        }
        const auto row=rows_[first+found];
        for(unsigned i=found;i>0;--i)rows_[first+i]=rows_[first+i-1];
        rows_[first]=row;
        return &rows_[first];
    }
    std::size_t index(std::uint64_t word) const {
        return (dense_ ? word : hash_word(word))&index_mask_;
    }
    unsigned length_,ways_;
    std::uint64_t mask_;
    std::vector<Row> rows_;
    std::size_t index_mask_;
    bool dense_;
};

}
