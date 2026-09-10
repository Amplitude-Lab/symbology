// Independent small exact systems and WXF fixture/dump helper for the audit.
#include "incremental_solve.hpp"
#include <random>
#include <sstream>
using Q=rat_t;
using I=int32_t;
using COO=sparse_tensor<Q,I,SPARSE_COO>;
using CSR=sparse_tensor<Q,I,SPARSE_CSR>;

int main(int argc,char**argv) {
  field_t F(FIELD_QQ); rref_option_t opt; opt->pool.reset(2); opt->verbose=false;
  auto* pool=&opt->pool;
  if(argc>2 && std::string(argv[1])=="dump") {
    COO t(sparse_tensor_read_wxf<Q,I>(std::filesystem::path(argv[2]),F,pool));
    std::cout<<"DIMS";for(auto d:t.dims())std::cout<<" "<<d;std::cout<<"\n";
    for(size_t k=0;k<t.nnz();k++) {std::cout<<"ENTRY";for(size_t d=0;d<t.rank();d++)std::cout<<" "<<t.index(k)[d];std::cout<<" = "<<t.val(k)<<"\n";}
    return 0;
  }
  if(argc>2 && std::string(argv[1])=="gen") {
    std::filesystem::path dir=argv[2];std::filesystem::create_directories(dir);
    auto write=[&](std::string name,COO&& t){t.canonicalize();auto bytes=sparse_tensor_write_wxf(CSR(std::move(t),pool));std::ofstream f(dir/name,std::ios::binary);f.write((char*)bytes.data(),bytes.size());};
    COO swap(std::vector<size_t>{2,2});swap.push_back({0,1},Q(1));swap.push_back({1,0},Q(1));write("swap.wxf",std::move(swap));
    COO v(std::vector<size_t>{1,2});v.push_back({0,0},Q(1));v.push_back({0,1},Q(2));write("vector.wxf",std::move(v));
    COO t(std::vector<size_t>{2,2,2});t.push_back({0,0,1},Q(2));t.push_back({1,1,0},Q(3));write("tensor.wxf",std::move(t));
    COO bad(std::vector<size_t>{2,3});bad.push_back({0,0},Q(1));write("rect.wxf",std::move(bad));
    COO cond(std::vector<size_t>{2,3});cond.push_back({0,0},Q(1));cond.push_back({0,2},Q(2));cond.push_back({1,2},Q(1));write("inconsistent.wxf",std::move(cond));
    COO seed(std::vector<size_t>{2,1,2});seed.push_back({0,0,0},Q(1));write("pair-seed.wxf",std::move(seed));
    COO rhs(std::vector<size_t>{1,2});rhs.push_back({0,0},Q(2));rhs.push_back({0,1},Q(1));write("pair-rhs.wxf",std::move(rhs));
    COO good(std::vector<size_t>{2,3});good.push_back({0,0},Q(1));good.push_back({0,2},Q(2));write("consistent.wxf",std::move(good));
    return 0;
  }
  std::mt19937 rng(42);int failures=0,total=0;
  for(int mode=0;mode<2;mode++)for(int trial=0;trial<45;trial++) {
    int n=1+trial%5,m=n+5;
    std::vector<std::vector<Q>> M(m,std::vector<Q>(n,Q(0)));
    std::vector<Q> x(n),b(m,Q(0));
    for(int j=0;j<n;j++)x[j]=Q(j-2,3);
    bool inconsistent=trial%3==2;
    int expected_rank=n;
    for(int i=0;i<m-1;i++)for(int j=0;j<n;j++)M[i][j]=i<n?Q(i==j):Q((int)(rng()%7)-3);
    if(trial%3==1 && n>1){for(int i=0;i<m;i++)M[i][n-1]=Q(2)*M[i][0];expected_rank=n-1;}
    for(int i=0;i<m;i++)for(int j=0;j<n;j++)b[i]+=M[i][j]*x[j];
    if(inconsistent)b[m-1]=Q(1); // explicit zero-coefficient/nonzero-RHS row
    COO a(std::vector<size_t>{(size_t)n,(size_t)m}),rhs(std::vector<size_t>{(size_t)m});
    for(int i=0;i<m;i++){for(int j=0;j<n;j++)if(M[i][j]!=Q(0))a.push_back({j,i},M[i][j]);if(b[i]!=Q(0))rhs.push_back({i},b[i]);}
    a.canonicalize();rhs.canonicalize();
    std::ostringstream silent;auto* old=std::cout.rdbuf(silent.rdbuf());
    bool ok=true;std::string why;
    try {
      auto result=mode==0?solve_linear_system_incremental<Q,I>(CSR(std::move(a),pool),CSR(std::move(rhs),pool),F,opt):solve_linear_system<Q,I>(CSR(std::move(a),pool),CSR(std::move(rhs),pool),F,opt);
      if(result.consistent==inconsistent){ok=false;why="wrong consistency verdict";}
      if(result.consistent) {
        for(int i=0;i<m;i++){Q residual=-b[i];for(size_t k=0;k<result.solution.nnz();k++)residual+=M[i][result.solution(k)]*result.solution[k];if(residual!=Q(0)){ok=false;why="nonzero exact residual";}}
        if(result.unique!=(expected_rank==n)){ok=false;why="wrong uniqueness verdict";}
        if(!result.unique) {
          // sparse_mat_rref_kernel(...).transpose() stores basis vectors as rows.
          if(result.null_space.nrow!=n-expected_rank || result.null_space.ncol!=n){ok=false;why="wrong null-space dimensions "+std::to_string(result.null_space.nrow)+"x"+std::to_string(result.null_space.ncol);}
          else for(int k=0;k<n-expected_rank;k++)for(int i=0;i<m;i++){Q v=Q(0);for(int j=0;j<n;j++){auto* c=result.null_space[k].find(j);if(c)v+=M[i][j]*(*c);}if(v!=Q(0)){ok=false;why="invalid null-space vector";}}
        }
      }
    }catch(std::exception&e){ok=false;why=e.what();}
    std::cout.rdbuf(old);total++;failures+=!ok;
    std::cout<<(ok?"PASS":"FAIL")<<" solver="<<(mode==0?"incremental":"sampled")<<" trial="<<trial<<" unknowns="<<n<<" inconsistent="<<inconsistent<<" "<<why<<"\n";
  }
  std::cout<<"SUMMARY "<<total-failures<<" passed "<<failures<<" failed\n";
  return failures?1:0;
}
