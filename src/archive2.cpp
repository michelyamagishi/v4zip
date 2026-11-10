#include "v4zip/archive.hpp"
#include "detail/bytes.hpp"
#include "detail/dna_config.hpp"
#include "detail/side_streams.hpp"
#include <algorithm>
#include <array>

namespace v4zip {
Bytes decode_archive_v1(const Bytes&,std::size_t);
namespace {
using namespace detail;
constexpr std::size_t stream_count=8,header_size=48+32*stream_count;
constexpr std::array<std::uint8_t,8> magic{'V','4','Z','I','P','N','E','W'};
enum Stream:std::size_t { dna_stream,structure_stream,identifier_stream,case_stream,
    ambiguity_stream,terminator_stream,separator_stream,quality_stream };
struct Line {std::size_t start,content_end,end;};
Line read_line(const Bytes& input,std::size_t& pos) {
    const auto start=pos;while(pos<input.size() && input[pos]!='\n')++pos;
    auto content_end=pos;
    if(pos<input.size() && pos>start && input[pos-1]=='\r')--content_end;
    if(pos<input.size())++pos;
    return {start,content_end,pos};
}
int symbol(std::uint8_t c) {
    switch(c){case 'A':case 'a':return 0;case 'C':case 'c':return 1;
        case 'G':case 'g':return 2;case 'T':case 't':return 3;default:return -1;}
}
struct RawStreams {
    std::array<Bytes,stream_count> data;
    std::vector<std::size_t> identifiers,separators;
    DnaConfig config;
};
RawStreams decompose(const Bytes& input,bool fastq) {
    RawStreams s;
    auto emit=[&](const Line& line,unsigned type) {
        const auto length=line.content_end-line.start;
        auto& structure=s.data[structure_stream];structure.push_back(static_cast<std::uint8_t>(type));varint(structure,length);
        s.data[terminator_stream].push_back(static_cast<std::uint8_t>(line.end-line.content_end));
        if(type!=1) {
            const auto stream=type==0?identifier_stream:(type==2?separator_stream:quality_stream);
            auto& target=s.data[stream];target.insert(target.end(),input.begin()+static_cast<std::ptrdiff_t>(line.start),input.begin()+static_cast<std::ptrdiff_t>(line.content_end));
            if(type==0)s.identifiers.push_back(length);
            else if(type==2)s.separators.push_back(length);
            return;
        }
        auto& dna=s.data[dna_stream];const auto before=dna.size();
        if(fastq)s.config.record_boundaries.push_back(before);
        std::size_t previous=line.start;
        for(auto i=line.start;i<line.content_end;) {
            const auto base=symbol(input[i]);
            if(base>=0) {
                auto& cases=s.data[case_stream];if(dna.size()%8==0)cases.push_back(0);
                if(input[i]>='a' && input[i]<='z')cases.back()|=static_cast<std::uint8_t>(1u<<(dna.size()%8));
                dna.push_back(static_cast<std::uint8_t>(base));++i;
            } else {
                const auto start=i++;while(i<line.content_end && input[i]==input[start])++i;
                auto& ambiguity=s.data[ambiguity_stream];varint(ambiguity,start-previous);varint(ambiguity,i-start);ambiguity.push_back(input[start]);previous=i;
            }
        }
        varint(structure,dna.size()-before);
    };
    std::size_t pos=0;
    if(!fastq) {
        if(!input.empty() && input[0]!='>')fail("FASTA must start with a header");
        while(pos<input.size()){const auto line=read_line(input,pos);emit(line,input[line.start]=='>'?0:1);}
    } else while(pos<input.size()) {
        const auto header=read_line(input,pos);
        if(header.content_end==header.start || input[header.start]!='@' || pos==input.size())fail("invalid FASTQ header");
        const auto sequence=read_line(input,pos);if(pos==input.size())fail("truncated FASTQ sequence");
        const auto separator=read_line(input,pos);
        if(separator.content_end==separator.start || input[separator.start]!='+' || (pos==input.size() && (sequence.content_end!=sequence.start || separator.end==separator.content_end)))fail("invalid FASTQ separator");
        const auto quality=read_line(input,pos);
        if(sequence.content_end-sequence.start!=quality.content_end-quality.start)fail("FASTQ sequence and quality lengths differ");
        emit(header,0);emit(sequence,1);emit(separator,2);emit(quality,3);
    }
    return s;
}
std::uint32_t header_crc(const Bytes& bytes) {
    Bytes header(bytes.begin(),bytes.begin()+static_cast<std::ptrdiff_t>(header_size));
    put(header,40,0,4);return crc32(header);
}
bool bounded(std::uint64_t size,std::size_t original,unsigned factor) {
    return size<=std::numeric_limits<std::size_t>::max() &&
        (size<=16 || (size-16)/factor+((size-16)%factor!=0)<=original);
}
struct Structure {
    std::vector<std::size_t> identifiers,separators;
    DnaConfig config;
};
Structure check_structure(const std::array<Bytes,stream_count>& data,const std::array<std::size_t,stream_count>& sizes,std::size_t original,bool fastq) {
    Structure info;const auto& structure=data[structure_stream];const auto& ambiguity=data[ambiguity_stream];
    const auto& endings=data[terminator_stream];std::size_t pos=0,lines=0,out=0,bases=0,amb=0,sequence_length=0;
    std::array<std::size_t,stream_count> totals{};
    while(pos<structure.size()) {
        const unsigned type=structure[pos++];
        if(type>3 || (fastq?type!=lines%4:(type>1 || (!lines && type!=0))))fail("invalid structural line type");
        const auto length=read_size(structure,pos);
        if(length>original-out || lines>=endings.size() || endings[lines]>2)fail("invalid structural length");
        if(lines && endings[lines-1]==0)fail("line follows unterminated line");
        out+=length;if(endings[lines]>original-out)fail("terminator exceeds original length");out+=endings[lines];
        if((type==0 || type==2) && !length)fail("empty record identifier");
        if(type==1) {
            sequence_length=length;const auto count=read_size(structure,pos);
            if(count>length || count>sizes[dna_stream]-bases)fail("invalid sequence base count");
            if(fastq)info.config.record_boundaries.push_back(bases);
            bases+=count;auto remaining=length-count;std::size_t cursor=0;
            while(remaining) {
                const auto gap=read_size(ambiguity,amb),run=read_varint(ambiguity,amb);
                if(gap>length-cursor || run>length-cursor-gap || run>remaining || amb==ambiguity.size() || symbol(ambiguity[amb])>=0)fail("invalid ambiguity run");
                ++amb;cursor+=gap+run;remaining-=run;
            }
        } else {
            const auto stream=type==0?identifier_stream:(type==2?separator_stream:quality_stream);
            if(length>sizes[stream]-totals[stream])fail("structural stream count differs");
            totals[stream]+=length;
            if(type==0)info.identifiers.push_back(length);
            else if(type==2)info.separators.push_back(length);
            if(type==3 && length!=sequence_length)fail("quality length differs from sequence");
        }
        ++lines;
    }
    if(out!=original || lines!=endings.size() || (fastq && lines%4) || bases!=sizes[dna_stream] || amb!=ambiguity.size())fail("structural totals differ");
    for(auto stream:{identifier_stream,separator_stream,quality_stream})if(totals[stream]!=sizes[stream])fail("unused side-stream data");
    if(sizes[case_stream]!=bases/8+(bases%8!=0))fail("invalid case mask size");
    if(bases%8 && !data[case_stream].empty() && (data[case_stream].back()>>(bases%8)))fail("nonzero case padding");
    return info;
}
Bytes restore(const std::array<Bytes,stream_count>& data,std::size_t original,bool fastq) {
    Bytes result;result.reserve(original);std::size_t pos=0,line=0,base=0,amb=0;
    std::array<std::size_t,stream_count> cursor{};
    const auto& structure=data[structure_stream];
    auto emit_bases=[&](std::size_t n){for(std::size_t j=0;j<n;++j){const auto b=data[dna_stream][base];const bool lower=(data[case_stream][base/8]>>(base%8))&1;result.push_back(static_cast<std::uint8_t>((lower?"acgt":"ACGT")[b]));++base;}};
    while(pos<structure.size()) {
        const auto type=structure[pos++];const auto length=read_size(structure,pos);
        if(type==1) {
            const auto count=read_size(structure,pos);auto remaining=length-count;std::size_t emitted=0;
            while(remaining){const auto gap=read_size(data[ambiguity_stream],amb),run=read_varint(data[ambiguity_stream],amb);emit_bases(gap);result.insert(result.end(),run,data[ambiguity_stream][amb++]);emitted+=gap+run;remaining-=run;}
            emit_bases(length-emitted);
        } else {
            const auto stream=type==0?identifier_stream:(type==2?separator_stream:quality_stream);const auto& bytes=data[stream];
            if((type==0 && bytes[cursor[stream]]!=(fastq?'@':'>')) || (type==2 && bytes[cursor[stream]]!='+'))fail("invalid restored identifier");
            const auto first=bytes.begin()+static_cast<std::ptrdiff_t>(cursor[stream]);result.insert(result.end(),first,first+static_cast<std::ptrdiff_t>(length));cursor[stream]+=length;
        }
        const auto ending=data[terminator_stream][line++];if(ending==2)result.push_back('\r');if(ending)result.push_back('\n');
    }
    return result;
}
}

Bytes compress_archive(const Bytes& input,bool fastq) {
    auto raw=decompose(input,fastq);std::array<Bytes,stream_count> encoded;
    std::array<unsigned,stream_count> methods{};
    encoded[0]=detail::encode_dna(raw.data[0],raw.config);
    for(std::size_t i=1;i<stream_count;++i){encoded[i]=pack(raw.data[i]);methods[i]=1;}
    for(auto i:{identifier_stream,separator_stream}) {
        const auto tokens=encode_identifiers(raw.data[i],i==identifier_stream?raw.identifiers:raw.separators);
        const auto packed=pack(tokens);Bytes alternative(8);put(alternative,0,tokens.size(),8);alternative.insert(alternative.end(),packed.begin(),packed.end());
        if(alternative.size()<encoded[i].size()){encoded[i]=std::move(alternative);methods[i]=2;}
    }
    Bytes result(header_size);std::copy(magic.begin(),magic.end(),result.begin());put(result,8,2,4);put(result,12,fastq?1:0,4);
    put(result,16,input.size(),8);put(result,24,raw.data[0].size(),8);put(result,32,crc32(input),4);put(result,36,stream_count,4);
    for(std::size_t i=0;i<stream_count;++i){const auto at=48+32*i;put(result,at,i+1,4);put(result,at+4,methods[i],4);put(result,at+8,raw.data[i].size(),8);put(result,at+16,encoded[i].size(),8);put(result,at+24,crc32(encoded[i]),4);}
    put(result,40,header_crc(result),4);
    for(const auto& stream:encoded)result.insert(result.end(),stream.begin(),stream.end());
    return result;
}

Bytes decompress_archive(const Bytes& archive,std::size_t max_output) {
    if(archive.size()<12 || !std::equal(magic.begin(),magic.end(),archive.begin()))fail("not a new V4ZIP archive (legacy archives require the original tool)");
    const auto version=get(archive,8,4);if(version==1)return decode_archive_v1(archive,max_output);
    if(version!=2 || archive.size()<header_size || get(archive,12,4)>1 || get(archive,36,4)!=stream_count || get(archive,44,4))fail("unsupported archive version or flags");
    if(get(archive,40,4)!=header_crc(archive))fail("archive header checksum mismatch");
    const auto original64=get(archive,16,8),dna64=get(archive,24,8);
    if(original64>max_output || original64>std::numeric_limits<std::size_t>::max())fail("output exceeds --max-output");
    const auto original=static_cast<std::size_t>(original64);if(dna64>original)fail("invalid DNA count");const bool fastq=get(archive,12,4)!=0;
    std::array<Bytes,stream_count> data;std::array<std::size_t,stream_count> sizes{},offsets{},lengths{};std::array<unsigned,stream_count> methods{};
    std::size_t cursor=header_size;
    for(std::size_t i=0;i<stream_count;++i) {
        const auto at=48+32*i,decoded=get(archive,at+8,8),encoded=get(archive,at+16,8);const auto method=get(archive,at+4,4);
        if(get(archive,at,4)!=i+1 || get(archive,at+28,4) || !bounded(decoded,original,4))fail("invalid stream descriptor");
        if(i==0?method!=0:(method!=1 && !((i==identifier_stream || i==separator_stream) && method==2)))fail("unsupported stream method");
        if(encoded>archive.size()-cursor)fail("truncated archive stream");
        sizes[i]=static_cast<std::size_t>(decoded);offsets[i]=cursor;lengths[i]=static_cast<std::size_t>(encoded);methods[i]=static_cast<unsigned>(method);
        if(get(archive,at+24,4)!=crc32(archive.data()+cursor,lengths[i]))fail("archive stream checksum mismatch");
        cursor+=lengths[i];
    }
    if(cursor!=archive.size())fail("trailing archive data");
    if(sizes[0]!=dna64 || sizes[quality_stream]>original-static_cast<std::size_t>(dna64))fail("invalid stream counts");
    for(auto i:{structure_stream,case_stream,ambiguity_stream,terminator_stream})data[i]=unpack(archive.data()+offsets[i],lengths[i],sizes[i]);
    const auto structure=check_structure(data,sizes,original,fastq);
    for(auto i:{identifier_stream,separator_stream,quality_stream}) {
        const Bytes encoded(archive.begin()+static_cast<std::ptrdiff_t>(offsets[i]),archive.begin()+static_cast<std::ptrdiff_t>(offsets[i]+lengths[i]));
        if(methods[i]==1)data[i]=unpack(encoded.data(),encoded.size(),sizes[i]);
        else if(methods[i]==2){
            const auto& records=i==identifier_stream?structure.identifiers:structure.separators;
            if(encoded.size()<8)fail("truncated identifier frame");
            const auto n=get(encoded,0,8);if(!bounded(n,sizes[i],5))fail("identifier token size exceeds bound");
            const auto tokens=unpack(encoded.data()+8,encoded.size()-8,static_cast<std::size_t>(n));
            data[i]=decode_identifiers(tokens,records,sizes[i]);
        }
    }
    const Bytes dna(archive.begin()+static_cast<std::ptrdiff_t>(offsets[0]),archive.begin()+static_cast<std::ptrdiff_t>(offsets[0]+lengths[0]));
    data[0]=detail::decode_dna(dna,sizes[0],structure.config);
    auto result=restore(data,original,fastq);
    if(result.size()!=original || crc32(result)!=get(archive,32,4))fail("original checksum mismatch");
    return result;
}
}
