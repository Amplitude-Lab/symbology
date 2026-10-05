#pragma once
#include "symrep.hpp"
#include "numeric_parse.hpp"
#include <iomanip>

namespace symrep {
namespace fs=std::filesystem;
inline Tensor tensor(const Mat& m,const std::vector<size_t>& dims,thread_pool* pool=nullptr) {
    require(dims.size()>=2&&dims[0]==m.nrow,"matrix/tensor row mismatch");size_t n=1;
    for(size_t i=1;i<dims.size();++i)n*=dims[i];require(n==m.ncol,"matrix/tensor columns mismatch");
    Tensor result(dims,std::max(size_t(1),m.nnz()));
    for(size_t i=0;i<m.nrow;++i)result.data.rowptr[i+1]=result.rowptr()[i]+m[i].nnz();
    auto fill=[&](size_t i){size_t entry=result.rowptr()[i];
        for(auto [j,x]:m[i]){auto idx=result.index(entry);size_t pos=j;
            for(size_t d=dims.size()-1;d>0;--d){idx[d-1]=index(pos%dims[d]);pos/=dims[d];}
            result.val(entry++)=x;}
    };
    if(pool){pool->detach_loop(size_t(0),m.nrow,fill);pool->wait();}else for(size_t i=0;i<m.nrow;++i)fill(i);
    if(!result.check_sorted())result.sort_indices(pool);return result;
}
inline Mat flatten(const Tensor& t) {
    size_t nc=1;for(size_t j=1;j<t.rank();++j)nc*=t.dim(j);Mat m(t.dim(0),nc);
    for(size_t i=0;i<m.nrow;++i){m[i].reserve(t.rowptr()[i+1]-t.rowptr()[i]);
        for(size_t p=t.rowptr()[i];p<t.rowptr()[i+1];++p){size_t col=0;auto idx=t.index(p);
            for(size_t j=1;j<t.rank();++j)col=col*t.dim(j)+idx[j-1];m[i].push_back(index(col),t.val(p));}
        normalize(m[i]);}
    return m;
}
// A read-only view of a matrix with tensor-shaped columns. Its sequential
// cursors avoid a second allocation of all rational values at final output.
class MatrixTensorView {
    const Mat& matrix;std::vector<size_t> dimensions,offsets;mutable size_t row=0;
    mutable std::vector<I> coordinates;
    void locate(size_t p)const{if(p>=offsets.back())require(false,"matrix tensor view index out of range");
        if(p<offsets[row])row=0;while(p>=offsets[row+1])++row;}
public:
    MatrixTensorView(const Mat& m,std::vector<size_t> dims):matrix(m),dimensions(std::move(dims)),offsets(m.nrow+1){
        require(dimensions.size()>=2&&dimensions[0]==m.nrow,"matrix tensor view shape mismatch");coordinates.resize(dimensions.size()-1);size_t n=1;
        for(size_t j=1;j<dimensions.size();++j){require(!dimensions[j]||n<=SIZE_MAX/dimensions[j],"matrix tensor view shape overflow");n*=dimensions[j];}
        require(n==m.ncol,"matrix tensor view column mismatch");
        for(size_t i=0;i<m.nrow;++i)offsets[i+1]=offsets[i]+m[i].nnz();}
    const auto& dims()const{return dimensions;}size_t dim(size_t i)const{return dimensions[i];}
    size_t rank()const{return dimensions.size();}size_t nnz()const{return offsets.back();}
    const size_t* rowptr()const{return offsets.data();}
    const Q& val(size_t p)const{locate(p);return matrix[row][p-offsets[row]];}
    const I* index(size_t p)const{locate(p);size_t column=matrix[row](p-offsets[row]);
        for(size_t j=dimensions.size()-1;j>0;--j){coordinates[j-1]=static_cast<I>(column%dimensions[j]);column/=dimensions[j];}return coordinates.data();}
};
inline std::vector<uint8_t> encode_tensor_view(const Tensor& t){return sparse_tensor_write_wxf(t);}
inline std::vector<uint8_t> encode_tensor_view(const MatrixTensorView& t){
    using namespace WXF_PARSER;std::unordered_map<std::string,std::function<void(Encoder&)>> fields;
    std::vector<int64_t> dims(t.dims().begin(),t.dims().end());size_t nnz=t.nnz(),rank=t.rank();
    fields["#dims"]=[&](Encoder& e){e.push_packed_array({rank},dims);};
    fields["#rowptr"]=[&](Encoder& e){e.push_array_info({t.dim(0)+1},WXF_HEAD::array,3);e.push_ustr(t.rowptr(),t.dim(0)+1);};
    fields["#colindex"]=[&](Encoder& e){auto width=minimal_pos_signed_bits(1+*std::max_element(dims.begin()+1,dims.end()));
        e.push_array_info({nnz,rank-1},WXF_HEAD::array,width);
        auto append=[&]<class T>(){std::vector<T> entries(nnz*(rank-1));for(size_t i=0;i<nnz;++i){auto idx=t.index(i);
            for(size_t j=0;j+1<rank;++j)entries[i*(rank-1)+j]=idx[j]+1;}e.push_ustr(entries.data(),entries.size());};
        switch(width){case 0:append.template operator()<int8_t>();break;case 1:append.template operator()<int16_t>();break;
            case 2:append.template operator()<int32_t>();break;case 3:append.template operator()<int64_t>();break;default:require(false,"tensor output dimension too large");}
    };
    fields["#vals"]=[&](Encoder& e){auto integer=[&](const int_t& x){if(x.fits_si())e.push_integer(x.to_si());else e.push_bigint(x.get_str());};
        e.push_function("List",nnz);for(size_t i=0;i<nnz;++i){const auto& x=t.val(i);if(x.is_integer())integer(x.num());
            else{e.push_function("Rational",2);integer(x.num());integer(x.den());}}};
    return fullform_to_wxf("SparseArray[Automatic,#dims,0,{1,{#rowptr,#colindex},#vals}]",fields,true).buffer;
}
// Verify the exact serialized tensor without keeping a token object or a
// second rational tensor for every nonzero. Parser retains all syntax checks.
template<class TensorLike> inline void verify_written_tensor(const fs::path& path,const TensorLike& t) {
    using namespace WXF_PARSER;using namespace SparseRREF::wxf_checked;
    auto bytes=file_to_ustr(path);Parser parser(bytes);size_t position=0,value=0;int rational=0;int_t numerator;
    auto check=[&](bool yes){if(!yes)require(false,"WXF roundtrip mismatch: "+path.string());};
    auto integer_value=[&](const Token& token){exact_integer(token);
        return integer(token)?int_t(token.get_integer()):int_t(std::string(token.get_string_view()).c_str());};
    auto accept=[&](const Q& x){check(value<t.nnz()&&x==t.val(value));++value;};
    parser.parse_each([&](Token&& token){size_t p=position++;
        auto symbol=[&](const char* name){check(token.type==WXF_HEAD::symbol&&token.get_string_view()==name);};
        auto function=[&](size_t n){check(token.type==WXF_HEAD::func&&token.length==n);};
        switch(p){
        case 0:function(4);return;case 1:symbol("SparseArray");return;case 2:symbol("Automatic");return;
        case 3:integer_array(token,1);check(token.dimensions[1]==t.rank());for(size_t i=0;i<t.rank();++i)check(natural(token,i)==t.dim(i));return;
        case 4:check(integer(token)&&token.get_integer()==0);return;
        case 5:function(3);return;case 6:symbol("List");return;case 7:check(integer(token)&&token.get_integer()==1);return;
        case 8:function(2);return;case 9:symbol("List");return;
        case 10:integer_array(token,1);check(token.dimensions[1]==t.dim(0)+1);for(size_t i=0;i<=t.dim(0);++i)check(natural(token,i)==t.rowptr()[i]);return;
        case 11:integer_array(token,2);check(token.dim(0)==t.nnz()&&token.dim(1)==t.rank()-1);
            for(size_t i=0;i<t.nnz();++i)for(size_t j=0;j+1<t.rank();++j)check(natural(token,i*(t.rank()-1)+j)==size_t(t.index(i)[j])+1);return;
        case 12:function(t.nnz());return;case 13:symbol("List");return;
        default:break;
        }
        if(rational==1){symbol("Rational");rational=2;}
        else if(rational==2){numerator=integer_value(token);rational=3;}
        else if(rational==3){auto denominator=integer_value(token);check(denominator!=0);accept(Q(numerator,denominator));rational=0;}
        else if(token.type==WXF_HEAD::func){function(2);rational=1;}
        else accept(Q(integer_value(token)));
    });
    check(position>=14&&rational==0&&value==t.nnz());
}
template<class TensorLike,class Sink> inline void stream_tensor_encoding(const TensorLike& t,Sink&& sink) {
    // Encode directly in bounded chunks. The previous encoder held the whole
    // WXF plus a second complete packed-index array alongside the kernel.
    using namespace WXF_PARSER;Encoder e;
    auto flush=[&](){sink(e.buffer);e.clear();};
    auto maybe_flush=[&](){if(e.buffer.size()>=1024*1024)flush();};
    std::vector<int64_t> dims(t.dims().begin(),t.dims().end());size_t nnz=t.nnz(),rank=t.rank();
    e.push_ustr("8:");e.push_function("SparseArray",4);e.push_symbol("Automatic");e.push_packed_array({rank},dims);e.push_integer(0);
    e.push_function("List",3);e.push_integer(1);e.push_function("List",2);e.push_array_info({t.dim(0)+1},WXF_HEAD::array,3);
    for(size_t i=0;i<=t.dim(0);++i){int64_t p=t.rowptr()[i];e.push_ustr(&p,1);maybe_flush();}
    auto width=nnz?minimal_pos_signed_bits(1+*std::max_element(dims.begin()+1,dims.end())):0;
    e.push_array_info({nnz,rank-1},WXF_HEAD::array,width);
    auto indices=[&]<class T>(){for(size_t i=0;i<nnz;++i){auto idx=t.index(i);for(size_t j=0;j+1<rank;++j){T c=idx[j]+1;e.push_ustr(&c,1);}maybe_flush();}};
    switch(width){case 0:indices.template operator()<int8_t>();break;case 1:indices.template operator()<int16_t>();break;
        case 2:indices.template operator()<int32_t>();break;case 3:indices.template operator()<int64_t>();break;default:require(false,"tensor output dimension too large");}
    e.push_function("List",nnz);auto integer=[&](const int_t& x){if(x.fits_si())e.push_integer(x.to_si());else e.push_bigint(x.get_str());};
    for(size_t i=0;i<nnz;++i){const auto& x=t.val(i);if(x.is_integer())integer(x.num());else{e.push_function("Rational",2);integer(x.num());integer(x.den());}maybe_flush();}
    flush();
}
template<class TensorLike> inline void verify_streamed_encoding(const fs::path& path,const TensorLike& t){
    std::ifstream in(path,std::ios::binary);require(bool(in),"cannot read "+path.string());std::vector<uint8_t> actual;
    stream_tensor_encoding(t,[&](const auto& expected){actual.resize(expected.size());in.read(reinterpret_cast<char*>(actual.data()),actual.size());
        require(size_t(in.gcount())==actual.size()&&actual==expected,"WXF exact byte readback mismatch: "+path.string());});
    require(in.peek()==std::char_traits<char>::eof()&&!in.bad(),"WXF trailing bytes or read error: "+path.string());
}
template<class TensorLike> inline void write_tensor_view(const fs::path& path,const TensorLike& t) {
    require(!fs::exists(path),"refusing to overwrite "+path.string());std::ofstream out(path,std::ios::binary);
    stream_tensor_encoding(t,[&](const auto& bytes){out.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());require(bool(out),"cannot write "+path.string());});
    out.close();
    require(bool(out),"cannot write "+path.string());
    // Re-encode against the file in bounded chunks. This is an exact byte
    // comparison, not a probabilistic hash; the encoding is independently
    // tested against SparseRREF's encoder and parser for all scalar widths.
    verify_streamed_encoding(path,t);
}
inline void write(const fs::path& path,const Tensor& t){write_tensor_view(path,t);}
inline void write(const fs::path& path,const Mat& m){write(path,tensor(m,{m.nrow,m.ncol}));}
inline Tensor read(const fs::path& p){return sparse_tensor_read_wxf<Q,I>(p,QQ);}
inline Mat read_matrix(const fs::path& p){auto t=read(p);require(t.rank()==2,"expected matrix: "+p.string());return flatten(t);}
inline void fresh(const fs::path& root){require(!fs::exists(root),"output directory already exists: "+root.string());fs::create_directories(root);}
inline void save_layout(const fs::path& p,const Layout& layout,const std::vector<Model>& ms) {
    std::ofstream out(p);out<<"irrep\tcopies\tdimension\tdivision_dimension\toffset\n";auto off=layout.offsets(ms);
    for(size_t i=0;i<ms.size();++i)out<<i<<'\t'<<layout.copies[i]<<'\t'<<ms[i].dimension<<'\t'<<ms[i].division_dimension<<'\t'<<off[i]<<'\n';
    out.close();require(bool(out),"cannot write layout");
}
inline Layout load_layout(const fs::path& p,const std::vector<Model>& ms) {
    std::ifstream in(p);require(bool(in),"missing layout "+p.string());std::string header;std::getline(in,header);
    require(header=="irrep\tcopies\tdimension\tdivision_dimension\toffset","unknown layout schema");
    Layout l;l.copies.resize(ms.size());size_t offset=0;
    for(size_t i=0;i<ms.size();++i){size_t id,m,d,e,o;in>>id>>m>>d>>e>>o;
        require(bool(in)&&id==i&&d==ms[i].dimension&&e==ms[i].division_dimension&&o==offset,"invalid layout record");
        l.copies[i]=m;offset+=m*d;index(offset);}
    std::string extra;require(!(in>>extra),"trailing layout records");return l;
}
inline void seal(const fs::path& root) {
    require(!fs::exists(root/"checksums.tsv"),"bundle is already complete");
    std::vector<fs::path> files;for(const auto& entry:fs::directory_iterator(root))if(entry.is_regular_file()&&entry.path().filename()!="checksums.tsv")files.push_back(entry.path());
    std::sort(files.begin(),files.end());auto pending=root/".checksums.pending";
    require(!fs::exists(pending),"incomplete bundle seal already exists");std::ofstream out(pending);
    for(const auto& p:files)out<<p.filename().string()<<'\t'<<native_cache::digest(p)<<'\n';
    out.close();require(bool(out),"cannot seal output");fs::rename(pending,root/"checksums.tsv");
}
inline void verify_seal(const fs::path& root) {
    std::ifstream in(root/"checksums.tsv");require(bool(in),"missing bundle completion checksums");std::string file,hash;
    std::set<std::string> recorded,actual;
    while(in>>file){require(bool(in>>hash),"truncated completion checksums");require(fs::path(file).filename()==fs::path(file)&&file!="checksums.tsv","invalid checksum path");
        require(recorded.insert(file).second,"duplicate checksum record");require(native_cache::digest(root/file)==hash,"changed bundle file: "+file);}
    require(in.eof()&&!recorded.empty(),"malformed completion checksums");
    for(const auto& entry:fs::directory_iterator(root))if(entry.is_regular_file()&&entry.path().filename()!="checksums.tsv")actual.insert(entry.path().filename().string());
    require(actual==recorded,"bundle checksum inventory is incomplete");
}
inline void save_group(const fs::path& root,const Group& g,const std::vector<Model>& ms) {
    std::ofstream info(root/"group.tsv");info<<"symbology-symrep-v1 QQ "<<g.order()<<' '<<g.generators.size()<<'\n';
    for(size_t i=0;i<g.generators.size();++i){info<<g.names[i]<<'\n';write(root/("generator_"+std::to_string(i)+".wxf"),g.generators[i]);}
    info.close();require(bool(info),"cannot write group metadata");
    std::ofstream report(root/"irreps.tsv");report<<"irrep\tQQ_dimension\tcenter_degree\tdivision_dimension\tcharacter_in_group_word_order\n";
    Mat central(ms.size(),g.order()),primitive(ms.size(),g.order());
    for(size_t s=0;s<ms.size();++s){report<<s<<'\t'<<ms[s].dimension<<'\t'<<ms[s].degree<<'\t'<<ms[s].division_dimension<<'\t';
        for(size_t h=0;h<g.order();++h){if(h)report<<',';report<<ms[s].character[h];if(ms[s].central[h]!=0)central[s].push_back(index(h),ms[s].central[h]);if(ms[s].primitive[h]!=0)primitive[s].push_back(index(h),ms[s].primitive[h]);}
        report<<'\n';for(size_t j=0;j<g.names.size();++j)write(root/("irrep_"+std::to_string(s)+"_g"+std::to_string(j)+".wxf"),ms[s].actions[g.generator_ids[j]]);
    }
    report.close();require(bool(report),"cannot write irrep metadata");write(root/"central_idempotents.wxf",central);write(root/"primitive_idempotents.wxf",primitive);
    std::ofstream words(root/"group_words.tsv");words<<"element\tparent\tgenerator\n";for(size_t i=0;i<g.order();++i)words<<i<<'\t'<<g.words[i].first<<'\t'<<g.words[i].second<<'\n';
}
inline Group load_group(const fs::path& root,size_t limit=256) {
    verify_seal(root);std::ifstream in(root/"group.tsv");std::string schema,field;size_t order,n;in>>schema>>field>>order>>n;
    require(bool(in)&&schema=="symbology-symrep-v1"&&field=="QQ"&&n>0&&n<1024,"invalid group schema");
    std::vector<std::string> names(n);std::vector<Mat> gs;for(size_t i=0;i<n;++i){in>>names[i];gs.push_back(read_matrix(root/("generator_"+std::to_string(i)+".wxf")));}
    Group g(names,std::move(gs),limit);require(g.order()==order,"group order mismatch");return g;
}
inline void verify_models(const fs::path& root,const Group& g,const std::vector<Model>& ms) {
    auto centers=read_matrix(root/"central_idempotents.wxf"),primitive=read_matrix(root/"primitive_idempotents.wxf");
    require(centers.nrow==ms.size()&&centers.ncol==g.order()&&primitive.nrow==ms.size()&&primitive.ncol==g.order(),"saved model shape mismatch");
    for(size_t s=0;s<ms.size();++s){
        Mat c(1,g.order()),p(1,g.order());for(size_t h=0;h<g.order();++h){if(ms[s].central[h]!=0)c[0].push_back(index(h),ms[s].central[h]);if(ms[s].primitive[h]!=0)p[0].push_back(index(h),ms[s].primitive[h]);}
        Mat sc(1,g.order()),sp(1,g.order());sc[0]=centers[s];sp[0]=primitive[s];
        require(equal(c,sc)&&equal(p,sp),"saved irrep convention differs from the current decomposition");
        for(size_t j=0;j<g.names.size();++j)require(equal(ms[s].actions[g.generator_ids[j]],read_matrix(root/("irrep_"+std::to_string(s)+"_g"+std::to_string(j)+".wxf"))),"saved irrep action mismatch");
    }
}
// Reuse the prepared primitive projectors, while rechecking the rational
// central decomposition and each irreducibility certificate. Searching the
// subgroup lattice again cannot improve an already certified saved model.
inline std::vector<Model> load_models(const fs::path& root,const Group& g,rref_option_t opt){
    std::vector<size_t> degrees;auto centers=central_idempotents(g,degrees);
    auto primitive=read_matrix(root/"primitive_idempotents.wxf");
    require(primitive.nrow==centers.size()&&primitive.ncol==g.order(),"saved primitive projector shape mismatch");
    std::vector<Model> result;size_t total=0;
    for(size_t s=0;s<centers.size();++s){std::vector<Q> p(g.order());for(auto [j,x]:primitive[s])p[j]=x;
        require(g.product(p,p)==p&&g.product(centers[s],p)==p,"saved primitive projector fails algebra certificate");
        auto m=model_from_projector(g,p,degrees[s],opt);require(m.division_dimension>0,"saved model fails irreducibility certificate");
        m.central=centers[s];total+=m.dimension*m.dimension/m.division_dimension;result.push_back(std::move(m));}
    require(total==g.order(),"saved models do not account for the regular representation");
    verify_models(root,g,result);return result;
}
inline void save_space(const fs::path& root,const std::string& name,const Adapted& a,const Group& g,const std::vector<Model>& ms) {
    write(root/(name+"_basis.wxf"),a.basis);save_layout(root/(name+"_copies.tsv"),a.layout,ms);
    for(size_t j=0;j<g.names.size();++j)write(root/(name+"_g"+std::to_string(j)+".wxf"),action(a.layout,ms,g.generator_ids[j]));
}
} // namespace symrep
