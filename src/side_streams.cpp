#include "detail/side_streams.hpp"
#include "detail/bytes.hpp"
#include <string>
#include <unordered_map>

namespace v4zip::detail {
namespace {
void check_lengths(const std::vector<std::size_t>& lengths,std::size_t expected) {
    std::size_t sum=0;
    for(auto n:lengths){if(n>expected-sum)fail("side-stream lengths overflow");sum+=n;}
    if(sum!=expected)fail("side-stream lengths differ");
}
unsigned category(std::uint8_t c) {
    return c>='0' && c<='9'?1:((c>='A' && c<='Z') || (c>='a' && c<='z')?2:0);
}
bool number(const std::string& s,std::int64_t& value) {
    if(s.empty() || s.size()>18)return false;
    value=0;
    for(unsigned char c:s){if(c<'0' || c>'9')return false;value=value*10+(c-'0');}
    return true;
}
struct Dictionary {
    std::array<std::string,4096> words;
    std::unordered_map<std::string,std::size_t> index;
    std::size_t next=0;
    void insert(const std::string& word) {
        const auto old=index.find(words[next]);
        if(old!=index.end() && old->second==next)index.erase(old);
        words[next]=word;index[word]=next;next=(next+1)%words.size();
    }
};
std::size_t varint_size(std::uint64_t n){std::size_t size=1;while(n>=128){n>>=7;++size;}return size;}


}

Bytes encode_identifiers(const Bytes& input,const std::vector<std::size_t>& lengths) {
    check_lengths(lengths,input.size());Bytes result;Dictionary dictionary;
    std::array<std::int64_t,256> previous{};std::size_t cursor=0;
    for(auto length:lengths) {
        const auto end=cursor+length;unsigned column=0;
        while(cursor<end) {
            const auto start=cursor++;const auto kind=category(input[start]);
            while(cursor<end && cursor-start<256 && category(input[cursor])==kind)++cursor;
            const std::string token(input.begin()+static_cast<std::ptrdiff_t>(start),input.begin()+static_cast<std::ptrdiff_t>(cursor));
            const auto found=dictionary.index.find(token);
            const auto literal_cost=1+varint_size(token.size())+token.size();
            const auto reference_cost=found==dictionary.index.end()?literal_cost+1:1+varint_size(found->second);
            std::int64_t value=0;const bool numeric=column<previous.size() && number(token,value);
            const auto delta=numeric?value-previous[column]:0;
            const auto zigzag=delta<0?std::uint64_t(-delta)*2-1:std::uint64_t(delta)*2;
            const auto numeric_cost=numeric?2+varint_size(zigzag):literal_cost+1;
            if(reference_cost<literal_cost && reference_cost<=numeric_cost){result.push_back(1);varint(result,found->second);}
            else if(numeric_cost<literal_cost){result.push_back(2);result.push_back(static_cast<std::uint8_t>(token.size()));varint(result,zigzag);}
            else {result.push_back(0);varint(result,token.size());result.insert(result.end(),token.begin(),token.end());dictionary.insert(token);}
            if(numeric)previous[column]=value;
            ++column;
        }
    }
    return result;
}

Bytes decode_identifiers(const Bytes& encoded,const std::vector<std::size_t>& lengths,std::size_t expected) {
    check_lengths(lengths,expected);Bytes result;result.reserve(expected);Dictionary dictionary;
    std::array<std::int64_t,256> previous{};std::size_t cursor=0;
    for(auto length:lengths) {
        const auto end=result.size()+length;unsigned column=0;
        while(result.size()<end) {
            if(cursor==encoded.size())fail("truncated identifier tokens");
            const auto type=encoded[cursor++];std::string token;
            if(type==0) {
                const auto n=read_varint(encoded,cursor);
                if(n>256 || n>encoded.size()-cursor || n>end-result.size())fail("invalid identifier literal");
                token.assign(encoded.begin()+static_cast<std::ptrdiff_t>(cursor),encoded.begin()+static_cast<std::ptrdiff_t>(cursor+n));cursor+=n;
                dictionary.insert(token);
            } else if(type==1) {
                const auto id=read_size(encoded,cursor);
                if(id>=dictionary.words.size() || dictionary.words[id].empty())fail("invalid identifier reference");
                token=dictionary.words[id];
            } else if(type==2) {
                if(cursor==encoded.size() || column>=previous.size())fail("invalid identifier delta");
                const unsigned width=encoded[cursor++];const auto zz=read_uint64(encoded,cursor);
                if(width==0 || width>18 || zz>UINT64_C(1999999999999999998))fail("identifier delta overflow");
                const auto delta=(zz&1)?-static_cast<std::int64_t>((zz>>1)+1):static_cast<std::int64_t>(zz>>1);
                const auto value=previous[column]+delta;
                if(value<0 || value>INT64_C(999999999999999999))fail("identifier number overflow");
                token=std::to_string(value);if(token.size()>width)fail("identifier width overflow");
                token.insert(0,width-token.size(),'0');
            } else fail("unknown identifier token");
            if(token.size()>end-result.size())fail("identifier exceeds record");
            std::int64_t value=0;if(column<previous.size() && number(token,value))previous[column]=value;
            result.insert(result.end(),token.begin(),token.end());++column;
        }
    }
    if(cursor!=encoded.size())fail("trailing identifier tokens");
    return result;
}

}
