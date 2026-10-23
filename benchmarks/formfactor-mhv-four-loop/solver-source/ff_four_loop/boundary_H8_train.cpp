#include "boundary_tensor_train.hpp"
using namespace boundary_tt;
static double elapsed(std::chrono::steady_clock::time_point x){return std::chrono::duration<double>(std::chrono::steady_clock::now()-x).count();}
int main(int argc,char**argv){try{
 require(argc==4,"usage boundary_H8_train factors.txt controls.txt OUT");fs::path out=argv[3];fresh(out);auto start=std::chrono::steady_clock::now();rref_option_t opt;opt->pool.reset(1);
 std::ifstream in(argv[1]);size_t q;in>>q;std::vector<Train>f;
 for(size_t z=0;z<3;++z){size_t w,nnz;in>>w>>nnz;std::map<Word,Q>words;
  for(size_t k=0;k<nnz;++k){Word v(w);for(auto&i:v){in>>i;require(i>0&&i<=q,"one-based factor index");--i;}std::string s;in>>s;require(words.emplace(v,Q(s)).second,"duplicate factor word");}
  require(bool(in),"factor parse");auto t=from_words(q,words,opt);for(auto&[v,x]:words)require(coefficient(t,v)==x,"full factor word replay");save(out/("factor_"+std::to_string(z)),t);f.push_back(std::move(t));std::cout<<"FACTOR "<<z<<" ranks=";for(auto r:f.back().rank)std::cout<<r<<',';std::cout<<" seconds="<<elapsed(start)<<std::endl;
 }
 auto e2=shuffle(f[0],f[0],opt);std::cout<<"SHUFFLE e2 seconds="<<elapsed(start)<<std::endl;
 auto er3=shuffle(f[0],f[2],opt);std::cout<<"SHUFFLE er3 seconds="<<elapsed(start)<<std::endl;
 auto e2r2=shuffle(e2,f[1],opt);std::cout<<"SHUFFLE e2r2 seconds="<<elapsed(start)<<std::endl;
 auto e4=shuffle(e2,e2,opt);std::cout<<"SHUFFLE e4 seconds="<<elapsed(start)<<std::endl;
 auto r22=shuffle(f[1],f[1],opt);std::cout<<"SHUFFLE r22 seconds="<<elapsed(start)<<std::endl;
 auto h=sum({{std::move(er3),Q(1)},{std::move(e2r2),Q(1,2)},{std::move(e4),Q(1,24)},{std::move(r22),Q(1,2)}},opt);
 std::ifstream ctl(argv[2]);size_t count;ctl>>count;for(size_t i=0;i<count;++i){Word v(8);for(auto&j:v){ctl>>j;--j;}std::string x;ctl>>x;require(coefficient(h,v)==Q(x),"independent H8 coefficient oracle mismatch");}require(bool(ctl),"controls parse");
 save(out/"H8",h);std::ofstream receipt(out/"receipt.json");receipt<<"{\"all_checks_pass\":true,\"boundary_size\":"<<q<<",\"factor_words_checked\":true,\"every_core_factorization_exact\":true,\"independent_oracle_controls\":"<<count<<",\"finite_R2_squared_term_retained\":true,\"F4_target_read\":false,\"ranks\":[";for(size_t k=0;k<h.rank.size();++k){if(k)receipt<<',';receipt<<h.rank[k];}receipt<<"],\"seconds\":"<<elapsed(start)<<"}\n";
 std::cout<<"H8_COMPLETE ranks=";for(auto r:h.rank)std::cout<<r<<',';std::cout<<" seconds="<<elapsed(start)<<std::endl;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
