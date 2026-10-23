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
int main(int argc,char**argv){try{
 require(argc>=5&&argc<=6,"usage: sew_selective_presolve OUT prefix_weight full|galois output.wxf [reference.wxf]");
 fs::path out=argv[1],dest=argv[4];size_t w=std::stoul(argv[2]);std::string mode=argv[3];require(mode=="full"||mode=="galois","bad mode");rref_option_t opt;opt->pool.reset(1);auto started=std::chrono::steady_clock::now();
 auto f=read(out/("FEC_"+std::to_string(w)+".wxf"));size_t nf=f.dim(0),nl=175;
 std::vector<Mat>left,right,short_left,short_right;
 for(size_t k=1;k<=5;++k){left.push_back(read_matrix(out/("galois"+std::to_string(k)+"_w"+std::to_string(w)+".wxf")));right.push_back(read_matrix(out/("suffix_galois"+std::to_string(k)+"_w2.wxf")));}
 short_left=left;short_right=right;
 for(std::string name:{"cycmat","flipmat"}){left.push_back(read_matrix(out/"actions"/("prefix_"+name+"_w"+std::to_string(w)+".wxf")));right.push_back(read_matrix(out/"actions"/("suffix_"+name+"_w2.wxf")));}
 staged_kernel::ProductInvarianceRows invariants(std::move(left),std::move(right)),galois(std::move(short_left),std::move(short_right));
 auto local=read_matrix(out/"boundary_local_rref.wxf");require(local.ncol==88*nl,"local cache dimension mismatch");
 tensor_kernel::ConstraintRows seam(f,std::move(local),nl);f.clear();auto full=staged_kernel::stack(invariants,seam);
 std::cout<<"READY mode="<<mode<<" prefix_weight="<<w<<" columns="<<full.ncol()<<" generators=7 full_residual_blocks="<<full.residual_blocks()<<" setup_seconds="<<seconds(started)<<std::endl;
 structured_kernel::Statistics stats;Mat k;
 if(mode=="full")k=staged_kernel::solve(full,opt,{},&stats);
 else{GaloisPresolve source{full,galois};k=staged_kernel::solve(source,opt,{},&stats);}
 for(auto&row:k.rows)vec_cancel_divisor(row);k.sort_rows_by_nnz();write(dest,k);
 std::cout<<"COMPLETE mode="<<mode<<" dimension="<<k.nrow<<" nnz="<<k.nnz()<<" presolve_seconds="<<stats.presolve_s<<" solve_seconds="<<stats.solve_s<<" total_seconds="<<seconds(started)<<std::endl;
 if(argc==6){auto reference=read_matrix(argv[5]);require(reference.nrow==k.nrow&&reference.ncol==k.ncol,"reference dimensions differ");PrivateCarrierChart chart(reference);require(chart.complete,"reference private chart incomplete");std::vector<I>which(k.ncol,-1);for(size_t i=0;i<chart.columns.size();++i)which[chart.columns[i]]=index(i);Mat coefficients(k.nrow,k.nrow);
  for(size_t r=0;r<k.nrow;++r)for(auto[c,x]:k[r])if(which[c]>=0)coefficients[r].push_back(which[c],x/chart.values[which[c]]);
  for(auto&row:coefficients.rows)normalize(row);require(equal(mul(coefficients,reference,&opt->pool),k),"exact reference row-space residual");write(fs::path(dest.string()+".reference_coefficients.wxf"),coefficients);std::cout<<"EXACT_REFERENCE_SPACE_EQUALITY pass=1 dimension="<<k.nrow<<std::endl;}
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
