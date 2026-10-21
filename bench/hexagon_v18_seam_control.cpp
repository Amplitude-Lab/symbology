// One matched binary contains the frozen V17 finishing function and V18.
#define main frozen_v17_main
#include "joint_physics_v17.cpp"
#undef main
#include "shared/hexagon_final_seam_v18.hpp"

int main(int argc,char**argv)try{
    require(argc==11,"v18_seam_control H K FEC LEC DATA CHECKPOINT REFERENCE PRIME v17|v18|audit PLAN_DIRECTORY");
    Flint::set_memory_functions();size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]);
    fs::path fec=argv[3],lec=argv[4],data=argv[5],checkpoint=argv[6],reference=argv[7],plans=argv[10];
    ulong p=std::stoull(argv[8]);std::string mode=argv[9];
    require(mode=="v17"||mode=="v18"||mode=="audit","Invalid finishing mode");
    irr_stage::lock fec_guard(fec/".active.lock",true),lec_guard(lec/".active.lock",true);
    field_t F(FIELD_Fp,p);rref_option_t opt;irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;
    auto* pool=&opt->pool;irr::channels_t cg(pool);irr_bulk::check_lec_convention(lec,cg,pool);std::ofstream discard("/dev/null");
    auto mf=irr_bulk::read_copies(fec,h,true),ml=irr_bulk::read_copies(lec,k,false);
    auto of=irr_bulk::read_copies(fec,h-1,true),ol=irr_bulk::read_copies(lec,k-1,false);
    std::unique_ptr<IrrepSeamRows> storage;
    {
        auto first=irr_bulk::read_chain(fec,h,true,mf,discard,pool),last=irr_bulk::read_chain(lec,k,false,ml,discard,pool);
        auto W=irr::read_tensor(data/"dlogmat_full_ca_pd3.wxf",pool);
        auto t=irr::make_seam_template(cg,W,nmhv_sew_read_irrep_axis(data/"dlogmat_full_ca_pd3_blocks.tsv"),opt);
        auto layout=irr::seam_layout(t,of,ol,mf,ml);
        storage=std::make_unique<IrrepSeamRows>(cg,t,layout,of,ol,first,last,F,opt,false);
    }
    hexagon_native_memory::release_unused_integer_caches(pool);
    auto basis=irr_stage::basis_identity(fec,lec,h,k);irr_stage::verify_basis(checkpoint,basis,h,k);
    irr_stage::verify_basis(reference,basis,h,k);auto current=hexagon_joint_checkpoint::load(checkpoint,p,h,k);
    auto& rows=*storage;size_t initial_rows=current.nrow;
    if(mode=="audit"){
        auto kt=current.transpose();size_t tested=0;
        for(size_t b=0;b<rows.residual_blocks();b+=std::max<size_t>(1,rows.residual_blocks()/32)){
            size_t id=b*128;
            rows.residual_transpose(kt,current.nrow,b,[&](auto&& old){
                auto direct=rows.residual_direct(id++,kt,current.nrow);++tested;
                require(old.nnz()==direct.nnz(),"V18 independent contraction nnz");
                for(size_t j=0;j<old.nnz();++j)require(old(j)==direct(j)&&old[j]==direct[j],"V18 independent contraction coefficient");
            },pool);
        }
        std::cout<<"V18_CONTRACTION_AUDIT rows="<<tested<<" PASS"<<std::endl;return 0;
    }
    auto identity=hex_sha256_file(checkpoint/"input_identity.txt")+hex_sha256_file(data/"dlogmat_full_ca_pd3.wxf");
    auto started=std::chrono::steady_clock::now();
    auto result=mode=="v17"?finish_seam(rows,std::move(current),F,opt):
        hexagon_final_seam_v18::finish(rows,std::move(current),F,opt,plans,identity);
    double seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-started).count();
    auto expected=read_cache(reference,p,h,k);require(result.nrow==expected.nrow&&result.ncol==expected.ncol,"V18 reference shape");
    for(size_t r=0;r<result.nrow;++r){std::vector<ulong> a(result.ncol),b(result.ncol);
        for(auto[j,x]:result[r])a[j]=x;for(auto[j,x]:expected[r])b[j]=x;require(a==b,"V18 every final coefficient equality");}
    // Independent original implementation, after the timed region.
    auto kt=result.transpose();size_t checked=0;
    for(size_t block=0;block<rows.residual_blocks();++block)
        rows.residual_transpose(kt,result.nrow,block,[&](auto&& row){++checked;require(!row.nnz(),"Independent original sewing row failed");},pool);
    require(checked==rows.nrow(),"Independent seam audit incomplete");
    std::cout<<"V18_CONTROL_PASS {\"mode\":\""<<mode<<"\",\"phase_seconds\":"<<seconds
        <<",\"initial_rows\":"<<initial_rows<<",\"rows\":"<<result.nrow<<",\"columns\":"<<result.ncol
        <<",\"nnz\":"<<result.nnz()<<",\"all_coefficients_equal\":true,\"original_rows_verified\":"<<checked<<"}"<<std::endl;
}catch(const std::exception& e){std::cerr<<"V18_CONTROL_ERROR "<<e.what()<<std::endl;return 1;}
