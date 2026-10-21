#include "hexagon_mhv_double_bonus.hpp"
using irr::require;
int main(int argc,char** argv)try {
    require(argc>=3,"mhv_double_bonus_audit LEC FRESH_OUTPUT [--native-check H K FEC CACHE_ROOT PRIME]...");
    fs::path lec=fs::absolute(argv[1]),output=fs::absolute(argv[2]);irr_stage::no_links(lec);irr_stage::no_links(output);
    require(!fs::exists(output),"Bonus audit output must be fresh; production inputs are read-only");
    struct Control{size_t h,k;fs::path fec,native;ulong prime;};std::vector<Control> controls;
    for(int i=3;i<argc;){require(std::string(argv[i])=="--native-check"&&i+5<argc,"Invalid bonus native-control arguments");
        controls.push_back({std::stoul(argv[i+1]),std::stoul(argv[i+2]),fs::absolute(argv[i+3]),fs::absolute(argv[i+4]),std::stoull(argv[i+5])});i+=6;}
    for(const auto& control:controls)require(control.k>=2&&control.k<=4,"Bonus controls support LEC2--4");
    Flint::set_memory_functions();rref_option_t opt;irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;
    irr::channels_t cg(&opt->pool);fs::create_directories(output);irr_stage::lock guard(output/".active.lock");
    irr_stage::write(output/"status.tsv",{{"schema","mhv-double-bonus-conditional-audit-v1"},{"complete","false"},{"empirical_assumption","true"},
        {"author_source","https://arxiv.org/html/2308.08199#S7.SS2"},{"equations","76,77,78,79"},{"assumption_scope","double-final relations extrapolated to L10"}});
    auto basis=irr_stage::basis_identity(lec,lec,0,4);irr_stage::write(output/"lec_basis_identity.tsv",basis);
    auto B=mhv_double_bonus::alphabet_matrix(),base=mhv_double_bonus::published_bonus_relations();
    auto orbit=mhv_double_bonus::full_d3_orbit(base),constraints=mhv_double_bonus::coefficient_dual_transform(orbit,B);
    auto words=mhv_double_bonus::expand_lec_words(lec,4,cg,opt);
    mhv_double_bonus::calibrate_single_final_entry(words[1],B);
    auto ordinary=mhv_double_bonus::coefficient_dual_transform(mhv_double_bonus::published_qbar_calibration(),B);
    auto calibration=mhv_double_bonus::constraint_images(words[2],irr::chain_copies(lec,2,false),2,ordinary);
    for(const auto& block:calibration)require(!block.nnz(),"Published ordinary Qbar pair equations fail generic LEC2 calibration");
    std::cout<<"BONUS_ALPHABET_CALIBRATION single_final_entry=PASS ordinary_qbar_pairs=PASS coefficient_dual=B_R_Bt chronological_suffix=PASS"<<std::endl;
    mhv_double_bonus::export_matrix(output/"alphabet_new_to_old_dlog.wxf",B,opt);
    mhv_double_bonus::export_matrix(output/"published_double_bonus_old.wxf",base,opt);
    mhv_double_bonus::export_matrix(output/"double_bonus_D3_current.wxf",constraints,opt);
    std::vector<mhv_double_bonus::Restriction> restrictions;
    std::ofstream dimensions(output/"dimensions.tsv");dimensions<<"weight\told_dimension\treduced_dimension\tblock\told_copies\tnew_copies\n";
    for(size_t w=2;w<=4;++w){auto restriction=mhv_double_bonus::restrict_tail(words[w],irr::chain_copies(lec,w,false),w,constraints,opt);
        size_t olddim=mhv_double_bonus::dimension(restriction.original),newdim=mhv_double_bonus::dimension(restriction.reduced);
        std::cout<<"BONUS_TAIL_DIMENSION weight="<<w<<" original="<<olddim<<" conditional="<<newdim<<std::endl;
        if(w==2)require(olddim==25&&newdim==21,"Double-bonus LEC2 dimension should be 25 to 21");
        if(w==3)require(olddim==83&&newdim==65,"Propagated double-bonus LEC3 dimension should be 83 to 65");
        auto dir=output/joint_eet_weight_name(w);fs::create_directories(dir);
        for(size_t b=0;b<6;++b){auto name=irr::specs()[b].name;
            dimensions<<w<<'\t'<<olddim<<'\t'<<newdim<<'\t'<<name<<'\t'<<restriction.original[b]<<'\t'<<restriction.reduced[b]<<'\n';
            mhv_double_bonus::export_matrix(dir/(name+"_constraint_images.wxf"),restriction.images[b],opt);
            mhv_double_bonus::export_matrix(dir/(name+"_new_to_original_copy_lift.wxf"),restriction.lift[b],opt);
            mhv_double_bonus::export_matrix(dir/(name+"_original_to_new_copy_coordinates.wxf"),restriction.coordinates[b],opt);
        }restrictions.push_back(std::move(restriction));
    }dimensions.close();
    for(size_t i=0;i<controls.size();++i){auto& control=controls[i];auto receipt=mhv_double_bonus::native_check(restrictions[control.k-2],control.h,control.fec,lec,control.native,control.prime,opt);
        irr_stage::write(output/("native_control_"+std::to_string(i+1)+".tsv"),receipt);}
    irr_stage::fields receipt={{"schema","mhv-double-bonus-conditional-audit-v1"},{"complete","true"},{"empirical_assumption","true"},
        {"author_source","https://arxiv.org/html/2308.08199#S7.SS2"},{"equations","76,77,78,79"},{"full_D3_orbit","true"},
        {"shared_E_copy_map","true"},{"exact_QQ_copy_kernel_verified","true"},{"native_controls",std::to_string(controls.size())},
        {"assumption_scope","double-final bonus relations extrapolated to L10"},
        {"audit_executable_sha256",hex_sha256_file(fs::canonical("/proc/self/exe"))},
        {"all_orders_proof","false"},{"ten_loop_result","false"},{"production_inputs_modified","false"}};
    for(const auto& entry:fs::recursive_directory_iterator(output))if(entry.is_regular_file()&&entry.path().filename()!=".active.lock"&&entry.path().filename()!="status.tsv")
        receipt[fs::relative(entry.path(),output).generic_string()+"_sha256"]=hex_sha256_file(entry.path());
    irr_stage::write(output/"audit_complete.tsv",receipt);irr_stage::write(output/"status.tsv",{{"schema","mhv-double-bonus-conditional-audit-v1"},{"complete","true"},{"empirical_assumption","true"},
        {"audit_receipt_sha256",hex_sha256_file(output/"audit_complete.tsv")},{"native_controls",std::to_string(controls.size())}});
    std::cout<<"BONUS_CONDITIONAL_AUDIT_COMPLETE native_controls="<<controls.size()<<" source=author_equations empirical_assumption=1 no_ten_loop_result=1"<<std::endl;
}catch(const std::exception& e){std::cerr<<"BONUS_AUDIT_ERROR "<<e.what()<<std::endl;return 1;}
