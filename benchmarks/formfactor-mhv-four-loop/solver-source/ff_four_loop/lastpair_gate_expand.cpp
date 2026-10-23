#include "symrep_bootstrap.hpp"
#include <unordered_map>
#include <fstream>
using namespace symrep;using Word=uint64_t;using Symbol=std::unordered_map<Word,Q>;
static void add(Symbol&s,Word w,const Q&x){auto[it,_]=s.try_emplace(w,Q(0));it->second+=x;}
static void cancel(Symbol&s){std::erase_if(s,[](const auto&x){return x.second==0;});}
static std::vector<size_t> decode(Word c,size_t n){std::vector<size_t>w(n);for(size_t j=n;j-->0;){w[j]=c%88;c/=88;}require(c==0,"word overflow");return w;}
int main(int argc,char**argv){try{require(argc==3,"usage: lastpair_gate_expand OUT reconstructed_F3.wxf");fs::path out=argv[1];auto start=std::chrono::steady_clock::now();auto t=read(argv[2]);require(t.rank()==6&&t.dim(0)==88,"F3 must be six physical-letter axes");std::vector<Mat> maps;for(int s=1;s<=2;++s)maps.push_back(read_matrix(out/("lastpair_gate_M"+std::to_string(s)+".wxf")));std::vector<Symbol> contracted(2);
 for(size_t first=0;first<t.dim(0);++first)for(size_t p=t.rowptr()[first];p<t.rowptr()[first+1];++p){auto idx=t.index(p);Word prefix=first;for(size_t j=0;j<4;++j)prefix=88*prefix+idx[j];size_t last=idx[4];for(size_t s=0;s<2;++s)for(auto[j,x]:maps[s][last])add(contracted[s],88*prefix+j,t.val(p)*x);}
 size_t source_nnz=t.nnz();t.clear();for(auto&s:contracted)cancel(s);std::cout<<"CONTRACT source_words="<<source_nnz<<" terms="<<contracted[0].size()<<","<<contracted[1].size()<<" seconds="<<seconds(start)<<std::endl;
 std::ofstream report(out/"lastpair_gate_expanded.json");report<<"{\"source_words\":"<<source_nnz<<",\"gates\":[";bool all=true;
 for(size_t s=0;s<2;++s){Symbol result;for(auto&[c,x]:contracted[s]){size_t l=c%88;auto w=decode(c/88,5);for(size_t pos=0;pos<=5;++pos){Word z=0;for(size_t j=0;j<=5;++j){size_t a=j==pos?l:w[j<pos?j:j-1];z=88*z+a;}add(result,z,x);}}size_t pre=result.size();cancel(result);bool pass=result.empty();all&=pass;if(s)report<<',';report<<"{\"gate\":"<<s+1<<",\"contracted_derivative_log_words\":"<<contracted[s].size()<<",\"shuffle_accumulator_words_before_zero_removal\":"<<pre<<",\"nonzero_residual_words\":"<<result.size()<<",\"pass\":"<<(pass?"true":"false");
  if(!pass){auto it=std::min_element(result.begin(),result.end(),[](const auto&a,const auto&b){return a.first<b.first;});auto w=decode(it->first,6);report<<",\"witness\":{\"word\":[";for(size_t j=0;j<w.size();++j){if(j)report<<',';report<<w[j]+1;}report<<"],\"coefficient\":\""<<it->second<<"\"}";}report<<'}';std::cout<<"GATE "<<s+1<<" residual_words="<<result.size()<<" elapsed="<<seconds(start)<<std::endl;}
 report<<"],\"all_gates_zero\":"<<(all?"true":"false")<<",\"seconds\":"<<seconds(start)<<",\"coefficient_arithmetic\":\"exact Q, all insertion positions including multiplicity\",\"scope\":\"Necessary last-pair gate only; independently reconstructed F3 read, no F4 candidate or published target read\"}\n";
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
