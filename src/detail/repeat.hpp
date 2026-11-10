#pragma once
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

#include "range.hpp"
#include "hash.hpp"
namespace v4zip::detail {
struct MatchRow {std::uint64_t tag=UINT64_MAX;std::uint64_t position_transform=0;};
class RepeatPredictor {
public:
    RepeatPredictor():matches_(1u<<20) {}
    Distribution predict(const Bytes& history,std::size_t position,std::uint64_t context) {
        Distribution prediction{16384,16384,16384,16384};
        predicted_base_=4;
        if(position>=match_order) {
            const auto canonical=canonical_word(context&word_mask(match_order),match_order);
            auto& match=matches_[hash_word(canonical.word)&(matches_.size()-1)];
            if(match.tag==canonical.word) {
                const auto previous=match.position_transform>>2;
                const auto transform=canonical.transform^unsigned(match.position_transform&3);
                const bool reversed=(transform&1)!=0;
                if(!reversed || previous>match_order) {
                    const auto source=reversed ? previous-match_order-1 : previous;
                    if(source<position) {
                        predicted_base_=history[static_cast<std::size_t>(source)]^((transform&2)?3:0);
                        const std::uint32_t error=std::max(8u,4096u/(1+match_streak_));
                        prediction.fill(error);
                        prediction[predicted_base_]=scale-3*error;
                    }
                }
            }
            match.tag=canonical.word;
            match.position_transform=(std::uint64_t(position)<<2)|canonical.transform;
        }
        return prediction;
    }
    void observe(unsigned base) {
        match_streak_=predicted_base_==base ? std::min(1023u,match_streak_+1) : 0;
    }
private:
    static constexpr unsigned match_order=12;
    std::vector<MatchRow> matches_;
    unsigned predicted_base_=4,match_streak_=0;
};
}
