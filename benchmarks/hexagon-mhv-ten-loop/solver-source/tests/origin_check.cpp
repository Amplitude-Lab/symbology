#include "irr_object.hpp"
#include "irr_runtime.hpp"
#include "irr_sewing.hpp"
using Q=rat_t;using Mat=sparse_mat<Q,int32_t>;namespace fs=std::filesystem;
using irr::require;
// At the simultaneous positive hexagon origin, y_i ~ u_i and
// m_i ~ 1/u_i. These are coefficients of dlog(u), dlog(v), dlog(w).
const std::array<std::array<int,3>,9> pullback{{{1,-1,-1},{-1,1,-1},{-1,-1,1},{-1,0,0},{0,-1,0},{0,0,-1},{1,0,0},{0,1,0},{0,0,1}}};
struct Axis{irr::copies m,off;size_t n=0;explicit Axis(const irr::copies& c):m(c){for(size_t b=0;b<6;++b){off[b]=n;n+=m[b]*irr::specs()[b].irrep_dimension;}}};
// Cache rational recurrences and shared word prefixes. Contract the compact
// multiplicity recurrence through its fixed CG vectors; never expand full
// physical tensors or the high-weight symbol.
struct Chain {
 fs::path root;bool first;size_t maximum;const irr::channels_t& cg;
 std::vector<irr::copies> counts;std::vector<irr::matrices> recurrence;
 std::map<std::string,std::vector<Q>> memo;
 Chain(fs::path r,bool f,size_t w,const irr::channels_t& c,thread_pool* pool):root(r),first(f),maximum(w),cg(c),counts(w+1),recurrence(w+1){
  counts[0]={1,0,0,0,0,0};memo.emplace("",std::vector<Q>{Q(1)});
  for(size_t i=1;i<=w;++i){counts[i]=irr::chain_copies(root,i,first);recurrence[i]=irr::read_recurrence(root,i,first,cg,pool);}
 }
 const std::vector<Q>& evaluate(const std::string& sequence){
  if(auto it=memo.find(sequence);it!=memo.end())return it->second;
  size_t w=sequence.size();require(w&&w<=maximum,"Origin prefix length");
  const auto& previous=evaluate(sequence.substr(0,w-1));Axis a(counts[w-1]),b(counts[w]);
  int letter=sequence.back()-'0';std::vector<Q> current(b.n);
  for(size_t block=0;block<6;++block){const auto& K=recurrence[w][block];
   auto decoded=irr::decode_channels(cg,block,counts[w-1]);require(decoded.size()==K.ncol&&K.nrow==counts[w][block],"Origin recurrence layout");
   for(size_t mu=0;mu<irr::specs()[block].irrep_dimension;++mu){
    std::vector<Q> transfer(K.ncol);
    for(size_t j=0;j<decoded.size();++j){auto d=decoded[j];
     for(const auto& term:cg.cg[block].vectors[d.channel].component[mu]){
      int weight=pullback[term.coordinate%9][letter];if(!weight)continue;
      size_t old=a.off[d.block]+(term.coordinate/9-cg.axis[d.block].offset)*a.m[d.block]+d.copy;
      if(previous[old]!=0)transfer[j]+=previous[old]*Q(term.coefficient)*Q(weight);
     }
    }
    for(size_t row=0;row<K.nrow;++row)for(auto [j,x]:K[row])if(transfer[j]!=0)
     current[b.off[block]+mu*b.m[block]+row]+=x*transfer[j];
   }
  }
  return memo.emplace(sequence,std::move(current)).first->second;
 }
};
int main(int argc,char**argv)try{
 require(argc==7,"origin_symbol_check H K SEW FEC LEC CASES_TSV");size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]);Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(2);auto* pool=&opt->pool;irr::channels_t cg(pool);
 auto object=irr_io::read_object(argv[3],h,k,argv[4],argv[5],pool);require(object.kind=="unique"&&object.directions==0,"Need a unique physical symbol");Axis a(object.first),b(object.last);Chain head(argv[4],true,h,cg,pool),tail(argv[5],false,k,cg,pool);std::ifstream in(argv[6]);size_t loop,count=0;std::string word,expected;
 while(in>>loop>>word>>expected){if(loop!=(h+k)/2)continue;require(word.size()==h+k&&word.find_first_not_of("012")==std::string::npos,"Invalid origin word");
  auto suffix=word.substr(h);std::reverse(suffix.begin(),suffix.end());const auto& f=head.evaluate(word.substr(0,h));const auto& l=tail.evaluate(suffix);Q result=0;
  for(size_t block=0;block<6;++block)for(auto [j,x]:object.S[block][0])for(size_t mu=0;mu<irr::specs()[block].irrep_dimension;++mu){auto fi=a.off[block]+mu*a.m[block]+size_t(j)/b.m[block],li=b.off[block]+mu*b.m[block]+size_t(j)%b.m[block];result+=x*f[fi]*l[li]*Q(mu?3:1);}
  std::cout<<"ORIGIN_WORD loop="<<loop<<" word="<<word<<" expected="<<expected<<" actual="<<result<<std::endl;require(result==Q(expected),"All-order origin prediction mismatch");++count;
 }
 require(in.eof()&&count,"No validation cases");std::cout<<"ORIGIN_CHECK PASS words="<<count<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
