#include "irr_bulk.hpp"
#include "shared/irrep_seam_rows.hpp"
using namespace irr_bulk;
int main(int argc,char**argv)try{
 require(argc==8,"seam_prime_benchmark left|right|factored H K FEC LEC DATA THREADS");std::string mode=argv[1];require(mode=="left"||mode=="right"||mode=="factored","Unknown prime benchmark method");size_t h=std::stoul(argv[2]),k=std::stoul(argv[3]);fs::path fec=argv[4],lec=argv[5],data=argv[6];
 Flint::set_memory_functions();rref_option_t opt;opt->pool.reset(std::stoul(argv[7]));opt->method=1;opt->is_back_sub=true;opt->verbose=false;auto* pool=&opt->pool;irr::channels_t cg(pool);check_lec_convention(lec,cg,pool);std::ofstream discard("/dev/null");
 auto mf=read_copies(fec,h,true),ml=read_copies(lec,k,false),of=read_copies(fec,h-1,true),ol=read_copies(lec,k-1,false);auto first=read_chain(fec,h,true,mf,discard,pool),last=read_chain(lec,k,false,ml,discard,pool);auto W=irr::read_tensor(data/"dlogmat_full_ca_pd3.wxf",pool);auto cb=nmhv_sew_read_irrep_axis(data/"dlogmat_full_ca_pd3_blocks.tsv");auto t=irr::make_seam_template(cg,W,cb,opt);auto layout=irr::seam_layout(t,of,ol,mf,ml);
 field_t F(FIELD_Fp,n_nextprime(1ULL<<60,0));auto start=std::chrono::steady_clock::now();size_t nr=0,nnz=0;sparse_mat<ulong,int32_t> K;
 if(mode=="factored"){IrrepSeamRows rows(cg,t,layout,of,ol,first,last,F,opt);nr=rows.nrow();K=factored_staged_kernel(rows,F,opt);}
 else{auto A=irr::partial_seam(mode=="left",cg,t,layout,of,ol,first,last,F,opt);nr=A.nrow;nnz=A.nnz();A.sort_rows_by_nnz();nmhv_sew_pivot_stage_t<int32_t> stage;K=nmhv_sew_modular_kernel_rows(std::move(A),stage,F,opt);}
 auto duration=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();std::cout<<"PRIME_BENCH_RESULT {\"method\":\""<<mode<<"\",\"prime\":"<<F.mod.n<<",\"columns\":"<<layout.columns<<",\"rows\":"<<nr<<",\"stored_constraint_nnz\":"<<nnz<<",\"nullity\":"<<K.nrow<<",\"kernel_nnz\":"<<K.nnz()<<",\"solve_seconds\":"<<duration<<"}"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
