#pragma once
#include "orbit.hpp"
#include "repeat_v2.hpp"
#include "weights.hpp"
#include "model_settings.hpp"

namespace v4zip::detail {
class PredictorV2 {
    static constexpr std::array<unsigned,7> orders{2,4,6,8,12,16,20};
    using Weights=std::array<std::int32_t,22>;
    std::vector<WordTable> tables_;
    V4Repeat repeat_;
    std::array<Distribution,22> predictions_{};
    Weights global_{};
    std::array<Weights,56> local_{};
    unsigned bucket_=0;
    std::uint64_t context_=0;
    std::size_t since_reset_=0;
    void update(Weights& weights,unsigned base) {
        const auto& lookup=weight_tables();
        auto maximum=std::numeric_limits<std::int32_t>::min();
        for(unsigned e=0;e<22;++e) {
            const auto p=std::max(1u,std::min(65535u,predictions_[e][base]));
            weights[e]=weights[e]*65208/65536+lookup.log_probability[p];
            maximum=std::max(maximum,weights[e]);
        }
        for(auto& w:weights)w=std::max(-8192,w-maximum);
    }
public:
#ifdef V4ZIP_DIAGNOSTICS
    struct Diagnostics {
        std::uint64_t loss_256bits=0,repeat_positions=0;
        std::array<double,22> expert_contribution{};
        Cdf last{};
    } diagnostics;
    std::uint64_t replacements() const {
        std::uint64_t sum=0;for(const auto& t:tables_)sum+=t.replacements;return sum;
    }
    std::size_t model_storage_bytes() const {
        auto total=sizeof(*this)+sizeof(WeightTables)+tables_.capacity()*sizeof(WordTable)+repeat_.heap_bytes();
        for(const auto& t:tables_)total+=t.heap_bytes();
        return total;
    }
#endif
    PredictorV2(){tables_.reserve(orders.size());for(auto k:orders)tables_.emplace_back(k+1,context_ways,table_bits(k+1));}
    void reset(){context_=0;since_reset_=0;repeat_.reset();}
    Cdf predict(const Bytes& history,std::size_t position) {
        unsigned chosen_order=0,chosen_bin=0;
        for(unsigned i=0;i<orders.size();++i) {
            std::array<std::array<std::uint64_t,4>,3> evidence{};
            if(since_reset_>=orders[i])for(unsigned b=0;b<4;++b) {
                const auto counts=tables_[i].counts((context_<<2)|b);
                for(unsigned kind=0;kind<3;++kind)evidence[kind][b]=counts[kind];
            }
            const auto observations=evidence[0][0]+evidence[0][1]+evidence[0][2]+evidence[0][3];
            if(observations){chosen_order=i;chosen_bin=observations<4?0:(observations<16?1:(observations<64?2:3));}
            for(unsigned kind=0;kind<3;++kind) {
                const auto& c=evidence[kind];const auto total=c[0]+c[1]+c[2]+c[3];
                const auto reciprocal=(UINT64_C(1)<<48)/(total*20+4);
                for(unsigned b=0;b<4;++b)predictions_[3*i+kind][b]=static_cast<std::uint32_t>(((c[b]*20+1)*reciprocal)>>32);
            }
        }
        predictions_[21]=repeat_.predict(history,position,context_,since_reset_);
        // Uniform candidate votes carry no repeat information for this bank.
        const bool repeat_available=predictions_[21]!=Distribution{16384,16384,16384,16384};
        bucket_=(chosen_order*4+chosen_bin)*2+unsigned(repeat_available);
        std::array<std::uint64_t,4> global{},local{},mixed{};
        const auto& lookup=weight_tables();
        for(unsigned e=0;e<22;++e) {
            const auto g=lookup.weight[static_cast<std::size_t>(-global_[e])];
            const auto l=lookup.weight[static_cast<std::size_t>(-local_[bucket_][e])];
            for(unsigned b=0;b<4;++b){global[b]+=std::uint64_t(g)*predictions_[e][b];local[b]+=std::uint64_t(l)*predictions_[e][b];}
        }
        const auto gsum=global[0]+global[1]+global[2]+global[3],lsum=local[0]+local[1]+local[2]+local[3];
#ifdef V4ZIP_DIAGNOSTICS
        diagnostics.repeat_positions+=repeat_available;
        for(unsigned e=0;e<22;++e) {
            const auto total=predictions_[e][0]+predictions_[e][1]+predictions_[e][2]+predictions_[e][3];
            const auto g=lookup.weight[static_cast<std::size_t>(-global_[e])];
            const auto l=lookup.weight[static_cast<std::size_t>(-local_[bucket_][e])];
            diagnostics.expert_contribution[e]+=0.5*total*(double(g)/double(gsum)+double(l)/double(lsum));
        }
#endif
        for(unsigned b=0;b<4;++b)mixed[b]=global[b]*scale/gsum+local[b]*scale/lsum;
        const auto sum=mixed[0]+mixed[1]+mixed[2]+mixed[3];
        Cdf cdf{0,0,0,0,scale};std::uint64_t cumulative=0;
        for(unsigned b=0;b<3;++b){cumulative+=mixed[b];cdf[b+1]=b+1+static_cast<std::uint32_t>(cumulative*(scale-4)/sum);}
#ifdef V4ZIP_DIAGNOSTICS
        diagnostics.last=cdf;
#endif
        return cdf;
    }
    void observe(unsigned base,std::size_t) {
#ifdef V4ZIP_DIAGNOSTICS
        const auto frequency=diagnostics.last[base+1]-diagnostics.last[base];
        diagnostics.loss_256bits+=static_cast<std::uint64_t>(-weight_tables().log_probability[std::min(65535u,frequency)]);
#endif
        update(global_,base);update(local_[bucket_],base);repeat_.observe(base);
        const auto word=(context_<<2)|base;
        for(unsigned i=0;i<orders.size();++i)if(since_reset_>=orders[i])tables_[i].observe(word);
        context_=word&word_mask(20);++since_reset_;
    }
};
}
