#pragma once
#include "v4zip/dna.hpp"
#include "v4zip/v4.hpp"
#include <algorithm>
#include <array>
#include <limits>
#include <stdexcept>
#include <utility>

#include "range.hpp"
#include "weights.hpp"
#include "orbit.hpp"
#include "repeat.hpp"
namespace v4zip::detail {
class Predictor {
public:
    Predictor() {
        tables_.reserve(orders.size());
        for(auto order:orders) tables_.emplace_back(order+1);
    }
    Cdf predict(const Bytes& history,std::size_t position) {
        const auto& lookup=weight_tables();
        unsigned expert=0;
        for(std::size_t i=0;i<orders.size();++i) {
            std::array<std::array<std::uint64_t,4>,3> evidence{};
            if(position>=orders[i]) {
                for(unsigned base=0;base<4;++base) {
                    const auto c=tables_[i].counts((context_<<2)|base);
                    for(unsigned kind=0;kind<3;++kind) evidence[kind][base]=c[kind];
                }
            }
            for(unsigned kind=0;kind<3;++kind) {
                const auto& count=evidence[kind];
                const std::uint64_t total=count[0]+count[1]+count[2]+count[3];
                const auto reciprocal=(UINT64_C(1)<<48)/(std::uint64_t(total)*20+4);
                for(unsigned b=0;b<4;++b)
                    predictions_[expert][b]=static_cast<std::uint32_t>(((std::uint64_t(count[b])*20+1)*reciprocal)>>32);
                ++expert;
            }
        }
        predictions_[expert]=repeat_.predict(history,position,context_);
        std::array<std::uint64_t,4> mixed{};
        for(unsigned e=0;e<expert_count;++e) {
            const auto weight=lookup.weight[static_cast<std::size_t>(-weights_[e])];
            for(unsigned b=0;b<4;++b) mixed[b]+=std::uint64_t(weight)*predictions_[e][b];
        }
        // Keep a one-unit floor per symbol; use exact cumulative rounding.
        Cdf cdf{0,0,0,0,scale};
        const auto sum=mixed[0]+mixed[1]+mixed[2]+mixed[3];
        std::uint64_t cumulative=0;
        for(unsigned b=0;b<3;++b) {
            cumulative+=mixed[b];
            cdf[b+1]=b+1+static_cast<std::uint32_t>(cumulative*(scale-4)/sum);
        }
        return cdf;
    }

    void observe(unsigned base,std::size_t position) {
        const auto& lookup=weight_tables();
        std::int32_t maximum=std::numeric_limits<std::int32_t>::min();
        for(unsigned e=0;e<expert_count;++e) {
            const auto p=std::max(1u,std::min(65535u,predictions_[e][base]));
            // Gamma=65208/65536, rounded toward zero (no signed shifts).
            weights_[e]=weights_[e]*65208/65536+lookup.log_probability[p];
            maximum=std::max(maximum,weights_[e]);
        }
        for(auto& w:weights_) w=std::max(-8192,w-maximum);
        repeat_.observe(base);
        const auto word=(context_<<2)|base;
        for(std::size_t i=0;i<orders.size();++i)
            if(position>=orders[i]) tables_[i].observe(word);
        context_=word&word_mask(20);
    }
private:
    static constexpr std::array<unsigned,7> orders{2,4,6,8,12,16,20};
    static constexpr unsigned expert_count=3*orders.size()+1;
    std::vector<WordTable> tables_;
    RepeatPredictor repeat_;
    std::array<Distribution,expert_count> predictions_{};
    std::array<std::int32_t,expert_count> weights_{};
    std::uint64_t context_=0;
};
}
