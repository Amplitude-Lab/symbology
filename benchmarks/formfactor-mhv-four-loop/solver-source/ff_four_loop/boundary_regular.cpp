// Complete finite-letter subspace of recursively compressed collinear symbols.
// The quotient retains every mixed word containing a spurious letter.
#include "symrep_io.hpp"
using namespace symrep;
static double seconds(std::chrono::steady_clock::time_point start) {
 return std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();
}
int main(int argc,char** argv) {try {
 require(argc==5,"usage: boundary_regular BOUNDARY_DIR alphabet finite_letters max_weight");
 fs::path out=argv[1];size_t q=std::stoul(argv[2]),f=std::stoul(argv[3]),maxw=std::stoul(argv[4]);
 require(f<q,"finite alphabet must be proper");rref_option_t opt;opt->pool.reset(1);
 Mat good(f,q);for(size_t i=0;i<f;++i)good[i].push_back(index(i),Q(1));
 if(fs::exists(out/"regular_prefix_w1.wxf"))require(equal(good,read_matrix(out/"regular_prefix_w1.wxf")),"cached finite-letter seed");
 else write(out/"regular_prefix_w1.wxf",good);
 for(size_t w=2;w<=maxw;++w) {
  auto start=std::chrono::steady_clock::now();
  auto target=out/("regular_prefix_w"+std::to_string(w)+".wxf");
  if(fs::exists(target)){good=read_matrix(target);continue;}
  auto t=read(out/"collinear"/("first_w"+std::to_string(w)+"_basis.wxf"));
  require(t.rank()==3&&t.dim(2)==q&&t.dim(1)==good.ncol,"recursive boundary basis shape");
  const size_t n=t.dim(0),d=t.dim(1),g=good.nrow;auto b=flatten(t);t.clear();
  Chart chart(good,opt);auto normal=mul(chart.inverse,good);
  std::vector<I> free(d,-1),pivot(d,-1);for(size_t k=0;k<g;++k)pivot[chart.pivots[k]]=index(k);
  size_t nf=0;for(size_t p=0;p<d;++p)if(pivot[p]<0)free[p]=index(nf++);
  require(nf==d-g,"regular quotient dimension");
  Mat quotient(d,nf);
  for(size_t p=0;p<d;++p)if(free[p]>=0)quotient[p].push_back(free[p],Q(1));
  else for(auto [j,x]:normal[pivot[p]])if(free[j]>=0)quotient[p].push_back(free[j],-x);
  require(mul(good,quotient).nnz()==0,"finite-prefix quotient annihilates every regular row");
  Mat c(n,nf*f+d*(q-f));
  for(size_t i=0;i<n;++i){
   for(auto [j,x]:b[i]){size_t p=j/q,l=j%q;
    if(l>=f)c[i].push_back(index(nf*f+p*(q-f)+l-f),x);
    else for(auto [k,y]:quotient[p])c[i].push_back(index(size_t(k)*f+l),x*y);
   }normalize(c[i]);
  }b.clear();normal.clear();quotient.clear();
  const auto cnnz=c.nnz();
  std::cout<<"FINITE_CONSTRAINTS weight="<<w<<" image="<<n<<" previous_regular="<<g<<" dimensions="<<c.nrow<<"x"<<c.ncol<<" nnz="<<cnnz<<" seconds="<<seconds(start)<<std::endl;
  write(out/("regular_constraints_w"+std::to_string(w)+".wxf"),c);
  good=kernel(transpose(c),opt);
  require(mul(good,c).nnz()==0,"complete exact regularity residual");
  write(target,good);
  std::ofstream receipt(out/("regular_prefix_w"+std::to_string(w)+".json"));
  receipt<<"{\"weight\":"<<w<<",\"image_dimension\":"<<n<<",\"regular_dimension\":"<<good.nrow<<",\"constraint_nonzeros\":"<<cnnz<<",\"exact_full_residual_zero\":true,\"all_spurious_words_retained\":true,\"finite_target_read\":false,\"method\":\"recursive finite-prefix quotient for finite last letters, raw prefix for every spurious last letter; complete modular-reconstructed QQ kernel\",\"seconds\":"<<seconds(start)<<"}\n";
  std::cout<<"FINITE_COMPLETE weight="<<w<<" dimension="<<good.nrow<<" seconds="<<seconds(start)<<std::endl;
 }
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
