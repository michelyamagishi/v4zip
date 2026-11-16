#include "v4zip/dna.hpp"
#include "detail/dna_config.hpp"
#include "detail/genome_profiler.hpp"
#include "detail/predictor.hpp"
#include "detail/predictor_v2.hpp"
#include "detail/predictor_v3.hpp"
#include "detail/repeat_codec.hpp"

namespace v4zip::detail {
namespace {
void validate(const DnaConfig& config,std::size_t count) {
    if(config.format_version<1 || config.format_version>3)throw std::runtime_error("unsupported DNA configuration");
    if(!std::is_sorted(config.record_boundaries.begin(),config.record_boundaries.end()) ||
       (!config.record_boundaries.empty() && config.record_boundaries.back()>count) ||
       (config.format_version==1 && !config.record_boundaries.empty()))throw std::runtime_error("invalid DNA boundaries");
}
}
Bytes encode_dna(const Bytes& bases,const DnaConfig& config) {
    validate(config,bases.size());
    if(config.format_version!=2 && config.format_version!=3)throw std::runtime_error("encoding requires format 2 or 3");
    for(auto b:bases)if(b>3)throw std::runtime_error("DNA symbol is not in ACGT");
    if(!bases.empty() && config.format_version==3) {
        const auto prof=analyze_genome(bases,config);
        Bytes result{static_cast<std::uint8_t>(prof.strategy)};
        if(prof.strategy==DnaStrategy::repeat_heavy) {
            const auto rep=RepeatDnaCodec::encode(bases,config);
            result.insert(result.end(),rep.begin(),rep.end());
            return result;
        }
        RangeEncoder encoder;
        PredictorV3 model(prof.strategy);std::size_t boundary=0;
        for(std::size_t i=0;i<bases.size();++i) {
            while(boundary<config.record_boundaries.size() && config.record_boundaries[boundary]==i){model.reset();++boundary;}
            encoder.put(bases[i],model.predict(bases,i));model.observe(bases[i],i);
        }
        const auto body=encoder.finish();
        result.insert(result.end(),body.begin(),body.end());
        return result;
    }
    RangeEncoder encoder;
    if(!bases.empty()) {
        PredictorV2 model;std::size_t boundary=0;
        for(std::size_t i=0;i<bases.size();++i) {
            while(boundary<config.record_boundaries.size() && config.record_boundaries[boundary]==i){model.reset();++boundary;}
            encoder.put(bases[i],model.predict(bases,i));model.observe(bases[i],i);
        }
    }
    return encoder.finish();
}
Bytes decode_dna(const Bytes& encoded,std::size_t count,const DnaConfig& config) {
    validate(config,count);
    if(count && config.format_version==3) {
        if(encoded.empty())fail("missing v3 DNA strategy header");
        if(encoded[0]>4)fail("unsupported v3 DNA strategy");
        const auto strategy=static_cast<DnaStrategy>(encoded[0]);
        const Bytes payload(encoded.begin()+1,encoded.end());
        if(strategy==DnaStrategy::repeat_heavy) return RepeatDnaCodec::decode(payload,count);
        RangeDecoder decoder(payload);Bytes bases;bases.reserve(count);
        PredictorV3 model(strategy);std::size_t boundary=0;
        for(std::size_t i=0;i<count;++i) {
            while(boundary<config.record_boundaries.size() && config.record_boundaries[boundary]==i){model.reset();++boundary;}
            const auto b=decoder.get(model.predict(bases,i));bases.push_back(static_cast<std::uint8_t>(b));model.observe(b,i);
        }
        decoder.finish();return bases;
    }
    RangeDecoder decoder(encoded);Bytes bases;bases.reserve(count);
    if(count && config.format_version==1) {
        Predictor model;
        for(std::size_t i=0;i<count;++i){const auto b=decoder.get(model.predict(bases,i));bases.push_back(static_cast<std::uint8_t>(b));model.observe(b,i);}
    } else if(count) {
        PredictorV2 model;std::size_t boundary=0;
        for(std::size_t i=0;i<count;++i) {
            while(boundary<config.record_boundaries.size() && config.record_boundaries[boundary]==i){model.reset();++boundary;}
            const auto b=decoder.get(model.predict(bases,i));bases.push_back(static_cast<std::uint8_t>(b));model.observe(b,i);
        }
    }
    decoder.finish();return bases;
}
}
namespace v4zip {
Bytes encode_dna(const Bytes& bases){return detail::encode_dna(bases,{});}
Bytes decode_dna(const Bytes& encoded,std::size_t count){return detail::decode_dna(encoded,count,{});}
}
