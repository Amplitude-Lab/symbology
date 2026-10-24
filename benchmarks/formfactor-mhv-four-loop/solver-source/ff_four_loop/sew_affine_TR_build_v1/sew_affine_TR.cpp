#include "symrep_bootstrap.hpp"
#include "product_invariance.hpp"
using namespace symrep;
// Optional exact short-row scheduling only. The complete source owns all
// modular discovery rows, residual blocks and integral-height certificates.
template<class Full> struct GaloisPresolve {
 const Full& full;const staged_kernel::ProductInvarianceRows& short_rows;
 size_t ncol()const{return full.ncol();}
 auto masks(bool remember)const{return full.masks(remember);}
 template<class Emit>void operator()(structured_kernel::Relations&r,size_t limit,Emit&&emit)const{
  if(limit==2)short_rows(r,limit,emit);else full(r,limit,emit);
 }
 template<class Emit>void discover(structured_kernel::Relations&r,const std::vector<I>&columns,const field_t&fp,const staged_kernel::Selection&select,std::vector<std::vector<unsigned char>>&seen,Emit&&emit,thread_pool*pool)const{full.discover(r,columns,fp,select,seen,std::forward<Emit>(emit),pool);}
 size_t residual_blocks()const{return full.residual_blocks();}
 template<class Emit>void residual(const staged_kernel::MM&k,const field_t&fp,size_t block,Emit&&emit,thread_pool*pool)const{full.residual(k,fp,block,std::forward<Emit>(emit),pool);}
 int_t integral_residual_bound(const Mat&k)const{return full.integral_residual_bound(k);}
};

struct ZeroColumn {
 const staged_kernel::ProductInvarianceRows& source;
 explicit ZeroColumn(const staged_kernel::ProductInvarianceRows&s):source(s){}
 size_t ncol()const{return source.ncol()+1;}size_t nrow()const{return source.nrow();}
 template<class Emit>void operator()(structured_kernel::Relations&r,size_t l,Emit&&e)const{source(r,l,std::forward<Emit>(e));}
 template<class Emit>void modular(structured_kernel::Relations&r,const std::vector<I>&c,const field_t&f,Emit&&e,thread_pool*p=nullptr,size_t mx=SIZE_MAX,size_t mn=0,const streaming_echelon::Substitution*s=nullptr,std::vector<unsigned char>*seen=nullptr,unsigned shift=0)const{source.modular(r,c,f,std::forward<Emit>(e),p,mx,mn,s,seen,shift);}
 size_t residual_blocks()const{return source.residual_blocks();}
 template<class Emit>void residual(const staged_kernel::MM&k,const field_t&f,size_t b,Emit&&e,thread_pool*p)const{auto q=k;--q.ncol;for(auto&r:q.rows)if(r.nnz()&&size_t(r(r.nnz()-1))==q.ncol)r.resize(r.nnz()-1);source.residual(q,f,b,std::forward<Emit>(e),p);}
 int_t integral_residual_bound(const Mat&k)const{return source.integral_residual_bound(k);}
};
int main(int argc,char**argv){try{
 require(argc>=5,"usage: sew_selective_presolve OUT prefix_weight full|galois output.wxf [reference.wxf]");
 fs::path out=argv[1],dest=argv[4];size_t w=std::stoul(argv[2]);std::string mode=argv[3];require(mode=="full"||mode=="galois"||mode=="masked","bad mode");size_t threads=1;std::string referencePath;for(int i=5;i<argc;++i){if(std::string(argv[i])=="--threads"){require(i+1<argc,"missing threads");threads=std::stoul(argv[++i]);}else referencePath=argv[i];}rref_option_t opt;opt->pool.reset(threads);auto started=std::chrono::steady_clock::now();
 auto f=read(out/("FEC_"+std::to_string(w)+".wxf"));size_t nf=f.dim(0),nl=92;
 std::vector<Mat>left,right,short_left,short_right;
 for(size_t k=1;k<=5;++k){left.push_back(read_matrix(out/("galois"+std::to_string(k)+"_w"+std::to_string(w)+".wxf")));right.push_back(read_matrix(out/"actions"/("suffix_TR_"+std::array<std::string,5>{"galois1a","galois1b","galois2a","galois2b","galois3"}[k-1]+"_w2.wxf")));}
 short_left=left;short_right=right;
 for(std::string name:{"cycmat","flipmat"}){left.push_back(read_matrix(out/"actions"/("prefix_"+name+"_w"+std::to_string(w)+".wxf")));right.push_back(read_matrix(out/"actions"/("suffix_TR_"+name+"_w2.wxf")));}
 if(mode=="masked"){
 auto grades=[](const std::vector<Mat>&a){std::vector<unsigned>g(a[0].nrow);for(size_t k=0;k<5;++k)for(size_t i=0;i<g.size();++i){require(a[k][i].nnz()==1&&a[k][i](0)==index(i)&&(a[k][i][0]==1||a[k][i][0]==-1),"non-diagonal Galois action");if(a[k][i][0]==-1)g[i]|=1u<<k;}return g;};
 auto lg=grades(left),rg=grades(right);std::vector<Mat>old,letters;std::array<std::string,5>gn{"galois1a","galois1b","galois2a","galois2b","galois3"};for(size_t k=1;k<=5;++k){old.push_back(read_matrix(out/("galois"+std::to_string(k)+"_w"+std::to_string(w-1)+".wxf")));letters.push_back(read_matrix(out/"data"/(gn[k-1]+".wxf")));}
 tensor_kernel::ff_old_grades=grades(old);auto ag=grades(letters);tensor_kernel::ff_local_column_grades.resize(88*nl);for(size_t a=0;a<88;++a)for(size_t b=0;b<nl;++b)tensor_kernel::ff_local_column_grades[a*nl+b]=ag[a]^rg[b];
 for(size_t i=0;i<nf;++i)for(size_t p=f.rowptr()[i];p<f.rowptr()[i+1];++p){auto c=f.index(p);require(lg[i]==(tensor_kernel::ff_old_grades[c[0]]^ag[c[1]]),"prefix recurrence grade mismatch");}
 for(size_t h=5;h<7;++h){std::array<int,32>trans;trans.fill(-1);size_t terms=0;auto check=[&](const Mat&a,const std::vector<unsigned>&g){for(size_t i=0;i<a.nrow;++i)for(auto[j,x]:a[i]){auto&v=trans[g[i]];if(v<0)v=g[j];require(v==int(g[j]),"D4 action not grade homogeneous");++terms;}};check(left[h],lg);check(right[h],rg);std::set<int>images;for(int v:trans)if(v>=0)require(images.insert(v).second,"D4 grade map not injective");std::cout<<"GRADE_NORMALITY generator="<<h<<" exact_support_terms="<<terms<<" injective=1"<<std::endl;}
 staged_kernel::ff_active_output.resize(nf*nl);size_t count=0;for(size_t a=0;a<nf;++a)for(size_t b=0;b<nl;++b)count+=(staged_kernel::ff_active_output[a*nl+b]=lg[a]==rg[b]);std::cout<<"GRADE_MASK active_product="<<count<<" of="<<nf*nl<<std::endl;
 }
 staged_kernel::ProductInvarianceRows invariants(std::move(left),std::move(right)),galois(std::move(short_left),std::move(short_right));
 auto local=read_matrix(out/"boundary_local_TR_rref.wxf");auto rhs=read_matrix(out/("affine_TR_rhs_w"+std::to_string(w+2)+"_residual.wxf"));require(local.ncol==88*nl,"local cache dimension mismatch");
 tensor_kernel::AffineConstraintRows seam(f,std::move(local),nl,std::move(rhs));f.clear();ZeroColumn inv(invariants);auto full=staged_kernel::stack(inv,seam);
 std::cout<<"READY mode="<<mode<<" prefix_weight="<<w<<" columns="<<full.ncol()<<" generators=7 full_residual_blocks="<<full.residual_blocks()<<" setup_seconds="<<seconds(started)<<std::endl;
 structured_kernel::Statistics stats;Mat k;
 require(mode=="masked","affine TR requires masked mode");k=staged_kernel::solve(full,opt,{},&stats);
 for(auto&row:k.rows)vec_cancel_divisor(row);k.sort_rows_by_nnz();write(dest,k);
 std::cout<<"COMPLETE mode="<<mode<<" dimension="<<k.nrow<<" nnz="<<k.nnz()<<" presolve_seconds="<<stats.presolve_s<<" solve_seconds="<<stats.solve_s<<" total_seconds="<<seconds(started)<<std::endl;
 if(!referencePath.empty()){auto reference=read_matrix(referencePath);require(reference.nrow==k.nrow&&reference.ncol==k.ncol,"reference dimensions differ");PrivateCarrierChart chart(reference);require(chart.complete,"reference private chart incomplete");std::vector<I>which(k.ncol,-1);for(size_t i=0;i<chart.columns.size();++i)which[chart.columns[i]]=index(i);Mat coefficients(k.nrow,k.nrow);
  for(size_t r=0;r<k.nrow;++r)for(auto[c,x]:k[r])if(which[c]>=0)coefficients[r].push_back(which[c],x/chart.values[which[c]]);
  for(auto&row:coefficients.rows)normalize(row);require(equal(mul(coefficients,reference,&opt->pool),k),"exact reference row-space residual");write(fs::path(dest.string()+".reference_coefficients.wxf"),coefficients);std::cout<<"EXACT_REFERENCE_SPACE_EQUALITY pass=1 dimension="<<k.nrow<<std::endl;}
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
