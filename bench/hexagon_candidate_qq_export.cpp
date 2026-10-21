#include "irr_stage.hpp"
#include <regex>

// A format-only exporter. It does not solve, reconstruct, or certify physics.
// Its caller supplies independently checked canonical rationals in the fixed
// six-block column order. Output must be a new, separately provenanced root.
int main(int argc,char**argv)try {
    if(argc==3&&std::string(argv[1])=="--seal") {
        irr_stage::no_links(argv[2]);irr_stage::lock lock(std::filesystem::path(argv[2])/".active.lock");
        irr_stage::complete(argv[2]);std::cout<<"QQ_FORMAT_SEAL_PASS physics_certified=0"<<std::endl;return 0;
    }
    irr::require(argc==7,"candidate_qq_export H K FEC LEC VECTOR_TSV NEW_ROOT");
    Flint::set_memory_functions();
    const size_t h=irr_io::number(argv[1]),k=irr_io::number(argv[2]);
    const std::filesystem::path fec=argv[3],lec=argv[4],vector=argv[5],root=argv[6];
    irr_stage::no_links(root);irr_stage::no_links(vector);
    irr::require(!std::filesystem::exists(root),"Output root must not exist");
    irr::sewing_object_t obj;obj.h=h;obj.k=k;obj.kind="unique";obj.directions=0;
    obj.first=irr::chain_copies(fec,h,true);obj.last=irr::chain_copies(lec,k,false);
    size_t columns=0;for(size_t b=0;b<6;++b)columns+=obj.first[b]*obj.last[b];
    irr::require(columns<=size_t(INT32_MAX),"Column index overflow");
    sparse_mat<rat_t,int32_t> packed(1,columns);
    std::ifstream in(vector);irr::require(bool(in),"Cannot read vector");
    std::string line;size_t i=0;const std::regex numeral("0|-?[1-9][0-9]*"),positive("[1-9][0-9]*");
    while(std::getline(in,line)) {
        if(!line.empty()&&line[0]=='#'){irr::require(i==0,"Comment after data");continue;}
        const auto a=line.find('\t'),b=line.find('\t',a==std::string::npos?a:a+1);
        irr::require(a!=std::string::npos&&b!=std::string::npos&&line.find('\t',b+1)==std::string::npos,"Expected three fields");
        irr::require(line.substr(0,a)==std::to_string(i)&&i<columns,"Noncanonical/out-of-order column index");
        auto n=line.substr(a+1,b-a-1),d=line.substr(b+1);
        irr::require(std::regex_match(n,numeral)&&std::regex_match(d,positive),"Invalid integer or nonpositive denominator");
        rat_t value;irr::require(value.set_str(n+"/"+d,10,true)==0,"Invalid rational");
        irr::require(value.num().get_str()==n&&value.den().get_str()==d,"Rational must be in canonical lowest terms");
        if(value!=0)packed[0].push_back(int32_t(i),std::move(value));++i;
    }
    irr::require(in.eof()&&i==columns,"Truncated vector/read failure");
    auto nnz=packed[0].nnz();std::filesystem::create_directories(root);
    irr_stage::lock lock(root/".active.lock");
    rref_option_t opt;irr_run::threads(opt);auto* pool=&opt->pool;
    irr_io::write_reduced_packed(root,obj,std::move(packed),pool);
    // Independently reassemble the written block matrices into the input order.
    auto reload=irr_io::read_object(root,h,k,fec,lec,pool);
    in.clear();in.seekg(0);size_t b=0,off=0,pos=0,seen=0;
    while(std::getline(in,line)){
        if(!line.empty()&&line[0]=='#')continue;
        auto a=line.find('\t'),t=line.find('\t',a+1);
        rat_t expected;irr::require(expected.set_str(line.substr(a+1,t-a-1)+"/"+line.substr(t+1),10,true)==0,"Reload parse failure");
        while(b<6&&seen>=off+obj.first[b]*obj.last[b]){off+=obj.first[b]*obj.last[b];++b;pos=0;}
        irr::require(b<6,"Reload column overflow");const auto& row=reload.S[b][0];rat_t actual;
        if(pos<row.nnz()&&size_t(row(pos))==seen-off)actual=row[pos++];
        irr::require(expected==actual,"Written object differs from input vector");++seen;
    }
    irr::require(seen==columns,"Reload incomplete");
    std::cout<<"QQ_VECTOR_EXPORT_PASS columns="<<columns<<" nnz="<<nnz<<" h="<<h<<" k="<<k<<" exact_readback=1 physics_certified=0"<<std::endl;
    return 0;
}catch(const std::exception& e){std::cerr<<"ERROR "<<e.what()<<std::endl;return 1;}
