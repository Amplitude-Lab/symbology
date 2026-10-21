// Replay the frozen solver's actual final-seam function. No equation changes.
#define main archived_v17_main
#include "joint_physics_v17.cpp"
#undef main
int main(int argc,char**argv)try{
    require(argc==9,"seam_replay H K FEC LEC DATA CHECKPOINT_ROOT FINAL_NATIVE_ROOT PRIME");
    Flint::set_memory_functions();size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]);fs::path fec=argv[3],lec=argv[4],data=argv[5],checkpoint=argv[6],reference=argv[7];ulong p=std::stoull(argv[8]);
    field_t F(FIELD_Fp,p);rref_option_t opt;irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;auto* pool=&opt->pool;irr::channels_t cg(pool);
    irr_stage::lock fec_guard(fec/".active.lock",true),lec_guard(lec/".active.lock",true);
    irr_bulk::check_lec_convention(lec,cg,pool);std::ofstream discard("/dev/null");
    auto mf=irr_bulk::read_copies(fec,h,true),ml=irr_bulk::read_copies(lec,k,false),of=irr_bulk::read_copies(fec,h-1,true),ol=irr_bulk::read_copies(lec,k-1,false);
    std::unique_ptr<IrrepSeamRows> rows;
    {
        auto first=irr_bulk::read_chain(fec,h,true,mf,discard,pool),last=irr_bulk::read_chain(lec,k,false,ml,discard,pool);
        auto W=irr::read_tensor(data/"dlogmat_full_ca_pd3.wxf",pool);auto t=irr::make_seam_template(cg,W,nmhv_sew_read_irrep_axis(data/"dlogmat_full_ca_pd3_blocks.tsv"),opt);auto layout=irr::seam_layout(t,of,ol,mf,ml);
        rows=std::make_unique<IrrepSeamRows>(cg,t,layout,of,ol,first,last,F,opt,false);
    }
    hexagon_native_memory::release_unused_integer_caches(pool);
    auto identity=irr_stage::basis_identity(fec,lec,h,k);
    irr_stage::verify_basis(checkpoint,identity,h,k);irr_stage::verify_basis(reference,identity,h,k);
    auto current=hexagon_joint_checkpoint::load(checkpoint,p,h,k);
    std::cout<<"REPLAY_INPUT rows="<<current.nrow<<" columns="<<current.ncol<<" nnz="<<current.nnz()<<std::endl;
    auto start=std::chrono::steady_clock::now();auto result=finish_seam(*rows,std::move(current),F,opt);
    auto seconds=std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count();auto expected=read_cache(reference,p,h,k);
    require(result.nrow==expected.nrow&&result.ncol==expected.ncol,"Reference shape");
    for(size_t i=0;i<result.nrow;++i){std::vector<ulong>a(result.ncol),b(result.ncol);for(auto[j,x]:result[i])a[j]=x;for(auto[j,x]:expected[i])b[j]=x;require(a==b,"Full coefficient equality");}
    std::cout<<"REPLAY_PASS {\"columns\":"<<result.ncol<<",\"rows\":"<<result.nrow<<",\"nnz\":"<<result.nnz()<<",\"phase_seconds\":"<<seconds<<",\"all_coefficients_equal\":true,\"equation_changes\":false}"<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}
