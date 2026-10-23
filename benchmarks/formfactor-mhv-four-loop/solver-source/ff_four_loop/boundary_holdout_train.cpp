// Convert an independently sourced finite descendant into a compact exact
// tensor train for validation only. No candidate form factor is read.
#include "boundary_tensor_train.hpp"
#include <cstdint>
using namespace boundary_tt;
static uint64_t big(const unsigned char*b,size_t n){uint64_t v=0;for(size_t i=0;i<n;++i)v=(v<<8)|b[i];return v;}
int main(int argc,char**argv){try{
 require(argc==5,"usage boundary_holdout_train table.bin encoding_alphabet boundary_alphabet OUT");size_t base=std::stoul(argv[2]),q=std::stoul(argv[3]),w=8;require(base<=q,"finite alphabet embedding");fs::path out=argv[4];fresh(out);auto start=std::chrono::steady_clock::now();rref_option_t opt;opt->pool.reset(1);
 std::ifstream f(argv[1],std::ios::binary);std::vector<uint64_t>code;Mat child(0,1);unsigned char buf[16];while(f.read(reinterpret_cast<char*>(buf),16)){uint64_t word=big(buf,4);int64_t num=static_cast<int64_t>(big(buf+4,8));uint64_t den=big(buf+12,4);require(den>0&&num!=0,"nonzero rational record");require(code.empty()||code.back()<word,"sorted unique holdout words");code.push_back(word);Vec row;row.push_back(0,Q(int_t(num),int_t(den)));child.rows.push_back(std::move(row));}
 require(f.eof()&&f.gcount()==0&&!code.empty(),"complete holdout records");child.nrow=code.size();const size_t words=code.size();Train t{q,std::vector<size_t>(w+1),std::vector<Mat>(w)};t.rank[w]=1;
 for(size_t k=w;k>0;--k){std::vector<uint64_t>parent;for(auto c:code)if(parent.empty()||parent.back()!=c/base)parent.push_back(c/base);Mat m(parent.size(),q*t.rank[k]);size_t row=0;
  for(size_t i=0;i<code.size();++i){while(parent[row]!=code[i]/base)++row;size_t l=code[i]%base;for(auto[j,x]:child[i])m[row].push_back(index(l*t.rank[k]+j),x);}for(auto&v:m.rows)normalize(v);child.clear();code.clear();auto[c,r]=factor(m,opt);m.clear();t.rank[k-1]=r.nrow;t.core[k-1]=std::move(r);child=std::move(c);code=std::move(parent);std::cout<<"HOLDOUT_FACTOR cut="<<k-1<<" rank="<<t.rank[k-1]<<" elapsed="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
 }
 require(code.size()==1&&code[0]==0&&child.nrow==1&&child.ncol==1&&t.rank[0]==1,"valid weight-eight word range and scalar left boundary");t.core[0]=mul(child,t.core[0]);save(out/"R4",t);std::ofstream receipt(out/"receipt.json");receipt<<"{\"all_exact_factorizations_pass\":true,\"word_records\":"<<words<<",\"encoding_alphabet\":"<<base<<",\"boundary_alphabet\":"<<q<<",\"word_encoding\":\"same sorted big-endian >IqI table as validated holdout\",\"used_for_solve\":false,\"candidate_form_factor_read\":false,\"ranks\":[";for(size_t i=0;i<t.rank.size();++i){if(i)receipt<<',';receipt<<t.rank[i];}receipt<<"],\"seconds\":"<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<"}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
