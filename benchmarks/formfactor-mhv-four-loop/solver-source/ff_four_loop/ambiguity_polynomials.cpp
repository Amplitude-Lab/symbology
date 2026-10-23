#include "symrep_bootstrap.hpp"
using namespace symrep;
using W=uint64_t;using Dict=std::map<W,Q>;
int main(int argc,char**argv){try{
 require(argc==2,"usage ambiguity_polynomials OUT");fs::path out=argv[1];auto start=std::chrono::steady_clock::now();
 auto n=flatten(read(out/"sewn_affine_w8_directions_TR92.wxf")),s=read_matrix(out/"suffix_TR92_physical.wxf");require(n.nrow==3&&n.ncol==35598*92,"direction shape");
 std::array<std::set<size_t>,7>needed;std::array<std::map<size_t,Vec>,7>rec;
 for(const auto&r:n.rows)for(auto[c,x]:r)needed[6].insert(c/92);
 for(size_t w=6;w>0;--w){auto t=read(out/("FEC_"+std::to_string(w)+".wxf"));require(t.rank()==3&&t.dim(2)==88,"prefixshape");for(auto r:needed[w]){Vec v;for(size_t p=t.rowptr()[r];p<t.rowptr()[r+1];++p){auto idx=t.index(p);v.push_back(index(88*idx[0]+idx[1]),t.val(p));needed[w-1].insert(idx[0]);}rec[w][r]=std::move(v);}std::cout<<"SELECT weight="<<w<<" rows="<<needed[w].size()<<" parent_rows="<<needed[w-1].size()<<" seconds="<<seconds(start)<<std::endl;}
 std::array<std::map<size_t,Dict>,7>exp;exp[0][0][0]=Q(1);
 for(size_t w=1;w<=6;++w)for(auto&[r,v]:rec[w]){auto&d=exp[w][r];for(auto[c,x]:v)for(auto[q,y]:exp[w-1].at(c/88))d[88*q+c%88]+=x*y;for(auto it=d.begin();it!=d.end();)if(it->second==0)it=d.erase(it);else++it;require(d.size()<=1000000,"bounded expansion cap");}
 std::ofstream js(out/"ambiguity_polynomials.json");js<<"{\"weight\":8,\"affine_direction_dimension\":3,\"alphabet\":\"W1,W2,W3,W4=u1,u2,u3,u4\",\"basis_normalization\":\"Exactly the saved TR92 direction rows; coefficients multiply ordinary monomials in logarithms.\",\"directions\":[";
 for(size_t r=0;r<n.nrow;++r){Dict d;for(auto[c,x]:n[r])for(auto[q,y]:exp[6].at(c/92))for(auto[t,z]:s[c%92])d[88*88*q+t]+=x*y*z;for(auto it=d.begin();it!=d.end();)if(it->second==0)it=d.erase(it);else++it;require(d.size()<=1000000,"bounded final expansion cap");std::map<std::array<unsigned,4>,Q>pol;std::map<std::array<unsigned,4>,size_t>counts;
 for(auto[q,x]:d){auto u=q;std::array<unsigned,4>a{};for(size_t k=0;k<8;++k){require(u%88<4,"ambiguity not contained in u-log sector");++a[u%88];u/=88;}auto it=pol.find(a);if(it==pol.end())pol.emplace(a,x);else require(it->second==x,"word coefficient not fully symmetric");++counts[a];}
 auto fact=[](unsigned k){long z=1;for(unsigned j=2;j<=k;++j)z*=j;return z;};
 if(r)js<<',';js<<"{\"row\":"<<r+1<<",\"word_count\":"<<d.size()<<",\"polynomial_terms\":[";size_t k=0;for(auto[a,x]:pol){long div=1;for(auto z:a)div*=fact(z);require(counts[a]==size_t(fact(8)/div),"incomplete permutation orbit");if(k++)js<<',';js<<"[["<<a[0]<<','<<a[1]<<','<<a[2]<<','<<a[3]<<"],\""<<x/Q(div)<<"\"]";}js<<"]}";std::cout<<"POLYNOMIAL row="<<r+1<<" words="<<d.size()<<" terms="<<pol.size()<<" exact_shuffle_reconstruction=pass seconds="<<seconds(start)<<std::endl;}
 js<<"],\"all_checks_pass\":true,\"all_words_reconstructed_from_log_polynomials\":true,\"no_target_read\":true,\"seconds\":"<<seconds(start)<<"}\n";
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
