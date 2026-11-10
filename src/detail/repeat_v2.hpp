#pragma once
#include "range.hpp"
#include "hash.hpp"
#include "model_settings.hpp"
#include <cstddef>
namespace v4zip::detail {
class V4Repeat {
 struct Seed {std::uint64_t tag=UINT64_MAX;std::array<std::uint64_t,4> previous{};};
 struct Copy {std::size_t next=0,recent=0;int direction=1;unsigned complement=0,length=0,errors=0,base=0;std::uint32_t window=0;};
 std::array<std::vector<Seed>,2> seeds_;
 std::array<Copy,4> active_{};
 unsigned active_count_=0;
 static unsigned pop(std::uint32_t w){unsigned n=0;while(w){w&=w-1;++n;}return n;}
public:
#ifdef V4ZIP_DIAGNOSTICS
 std::size_t heap_bytes() const {return (seeds_[0].capacity()+seeds_[1].capacity())*sizeof(Seed);}
#endif
 V4Repeat(){for(unsigned i=0;i<2;++i)seeds_[i].resize(std::size_t{1}<<repeat_table_bits[i]);}
 void reset(){active_count_=0;}
 Distribution predict(const Bytes& history,std::size_t position,std::uint64_t context,std::size_t since_reset) {
  std::array<Copy,36> choices{};unsigned count=0;
  auto add=[&](Copy c){
   if(c.next>=position)return;
   c.base=history[c.next]^c.complement;
   for(unsigned j=0;j<count;++j)if(choices[j].next==c.next && choices[j].direction==c.direction && choices[j].complement==c.complement){if(c.length>choices[j].length)choices[j]=c;return;}
   if(count<choices.size())choices[count++]=c;
  };
  bool retained_top_four=active_count_==4;
  for(unsigned i=0;i<active_count_;++i){
   add(active_[i]);
   if(int(active_[i].length)-4*int(active_[i].errors)<=64)retained_top_four=false;
  }
  for(unsigned which=0;which<2;++which){
   const unsigned k=which?20:12;if(since_reset<k)continue;
   const auto word=context&word_mask(k);const auto canonical=canonical_word(word,k);
   auto& row=seeds_[which][hash_word(canonical.word)&(seeds_[which].size()-1)];
   if(row.tag==canonical.word){
    for(auto packed:row.previous){
     // Retrieved histories are capped at 64. None can displace four
     // continuations whose verified scores are strictly greater than 64.
     if(retained_top_four)break;
     if(!packed)continue;
     const auto previous=static_cast<std::size_t>(packed>>2);
     const auto g=canonical.transform^static_cast<unsigned>(packed&3);
     const auto reversed=reverse_word(word,k);
     const unsigned stabilizer=reversed==word?1:((reversed^word_mask(k))==word?3:0);
     for(unsigned alternative=0;alternative<(stabilizer?2u:1u);++alternative){
      const auto transform=g^(alternative?stabilizer:0);
      const bool reverse=(transform&1)!=0;
      if(reverse && previous<=k)continue;
      Copy c;c.next=reverse?previous-k-1:previous;c.direction=reverse?-1:1;c.complement=(transform&2)?3:0;c.recent=previous;
      if(c.next>=position)continue;
      bool retained=false;
      for(unsigned j=0;j<active_count_;++j)if(active_[j].next==c.next && active_[j].direction==c.direction && active_[j].complement==c.complement && active_[j].length>=64){retained=true;break;}
      if(retained)continue;
      // Full canonical words (not hash fingerprints) and stored orientations
      // authenticate all k seed bases. Only the older history needs scanning.
      if(transform_word(row.tag,k,canonical.transform)!=word)continue;
      c.length=k;
      for(unsigned distance=k+1;distance<=64 && distance<=since_reset;++distance){
       std::size_t source;
       if(reverse){if(distance>=position-c.next)break;source=c.next+distance;}
       else {if(distance>c.next)break;source=c.next-distance;}
       const bool mismatch=(history[source]^c.complement)!=history[position-distance];
       if(mismatch)break;
       ++c.length;
      }
      add(c);
     }
    }
   }else{row.tag=canonical.word;row.previous={};}
   for(unsigned j=3;j>0;--j)row.previous[j]=row.previous[j-1];
   row.previous[0]=(std::uint64_t(position)<<2)|canonical.transform;
  }
  auto better=[](const Copy& a,const Copy& b){
   const auto sa=int(a.length)-4*int(a.errors),sb=int(b.length)-4*int(b.errors);
   if(sa!=sb)return sa>sb;
   if(a.recent!=b.recent)return a.recent>b.recent;
   if(a.next!=b.next)return a.next>b.next;
   if(a.direction!=b.direction)return a.direction>b.direction;
   return a.complement<b.complement;
  };
  // Bounded insertion selection avoids implementation-dependent sort tie choices.
  active_count_=std::min(4u,count);
  for(unsigned i=0;i<active_count_;++i){unsigned best=i;for(unsigned j=i+1;j<count;++j)if(better(choices[j],choices[best]))best=j;std::swap(choices[i],choices[best]);active_[i]=choices[i];}
  if(!active_count_)return {16384,16384,16384,16384};
  std::array<std::uint64_t,4> votes{};std::uint64_t total=0;
  for(unsigned i=0;i<active_count_;++i){const auto& c=active_[i];const auto strength=std::max(1u,c.length-std::min(c.length,11u));const auto weight=std::max(1u,strength*strength/(1+4*c.errors));const auto error=std::max(8u,4096/(1+strength));for(unsigned b=0;b<4;++b)votes[b]+=std::uint64_t(weight)*(b==c.base?scale-3*error:error);total+=weight;}
  Distribution result{};for(unsigned b=0;b<4;++b)result[b]=static_cast<std::uint32_t>(votes[b]/total);return result;
 }
 void observe(unsigned base){
  unsigned kept=0;
  for(unsigned i=0;i<active_count_;++i){auto c=active_[i];const bool mismatch=base!=c.base;c.window=(c.window<<1)|unsigned(mismatch);c.errors=pop(c.window);if(c.errors>2u)continue;if(c.direction<0 && c.next==0)continue;c.next=c.direction>0?c.next+1:c.next-1;c.length=std::min(1024u,c.length+1);active_[kept++]=c;}
  active_count_=kept;
 }
};
}
