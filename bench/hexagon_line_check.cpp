#include "irr_object.hpp"
#include "irr_runtime.hpp"
using Q=rat_t;using Mat=sparse_mat<Q,int32_t>;namespace fs=std::filesystem;
using irr::require;
// Along v=w=1+epsilon, then epsilon -> 0:
// y_u ~ (1-u)^2/(u epsilon^2), y_v=y_w ~ -epsilon/(1-u).
// Dropping dlog epsilon gives the two-letter symbol pullback below.
const std::array<std::array<int,2>,9> pullback{{{1,0},{-1,0},{-1,0},{-1,1},{0,0},{0,0},{-1,2},{0,-1},{0,-1}}};
struct Axis{irr::copies m,off;size_t n=0;explicit Axis(const irr::copies& c):m(c){for(size_t b=0;b<6;++b){off[b]=n;n+=m[b]*irr::specs()[b].irrep_dimension;}}};
std::vector<Q> evaluate(const fs::path& root,size_t maximum,bool first,const std::string& word,const irr::channels_t& cg,thread_pool* pool){
 irr::copies older{1,0,0,0,0,0};std::vector<Q> previous(1,Q(1));
 for(size_t w=1;w<=maximum;++w){auto next=irr::chain_copies(root,w,first);Axis a(older),b(next);std::vector<Q> current(b.n);int letter=word[first?w-1:word.size()-w]-'0';
  for(size_t block=0;block<6;++block){auto K=joint_eet_flatten_tensor(irr::read_tensor(root/joint_eet_weight_name(w)/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[block].name+".wxf"),pool),pool);auto expanded=irr::expand(cg,block,older,K);
   for(size_t row=0;row<expanded.nrow;++row)for(auto [j,x]:expanded[row]){int weight=pullback[size_t(j)%9][letter];if(weight&&previous[size_t(j)/9]!=0)current[b.off[block]+row]+=x*previous[size_t(j)/9]*Q(weight);}}
  older=next;previous=std::move(current);
 }return previous;
}
int main(int argc,char**argv)try{
 require(argc==7,"line_symbol_check H K SEW FEC LEC CASES_TSV");size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]);Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(2);auto* pool=&opt->pool;irr::channels_t cg(pool);
 auto object=irr_io::read_object(argv[3],h,k,argv[4],argv[5],pool);require(object.kind=="unique"&&object.directions==0,"Need a unique physical symbol");Axis a(object.first),b(object.last);std::ifstream in(argv[6]);size_t loop,count=0;std::string word,expected;
 while(in>>loop>>word>>expected){if(loop!=(h+k)/2)continue;require(word.size()==h+k&&word.find_first_not_of("01")==std::string::npos,"Invalid line word");
  auto f=evaluate(argv[4],h,true,word.substr(0,h),cg,pool),l=evaluate(argv[5],k,false,word.substr(h),cg,pool);Q result=0;
  for(size_t block=0;block<6;++block)for(auto [j,x]:object.S[block][0])for(size_t mu=0;mu<irr::specs()[block].irrep_dimension;++mu){auto fi=a.off[block]+mu*a.m[block]+size_t(j)/b.m[block],li=b.off[block]+mu*b.m[block]+size_t(j)%b.m[block];result+=x*f[fi]*l[li]*Q(mu?3:1);}
  std::cout<<"LINE_WORD loop="<<loop<<" word="<<word<<" expected="<<expected<<" actual="<<result<<std::endl;require(result==Q(expected),"Independent published coefficient mismatch");++count;
 }
 require(in.eof()&&count,"No validation cases");std::cout<<"PUBLISHED_LINE_CHECK PASS words="<<count<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
