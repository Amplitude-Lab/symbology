#include "irr_physics_pipeline.hpp"

int main(int argc,char** argv)try {
    require(argc==2,"stage_options NEW_TEST_DIRECTORY");fs::path root=fs::absolute(argv[1]);
    require(!fs::exists(root),"Use a new test directory");fs::create_directories(root);
    irr_stage::fixed_text(root/"first.tsv","projection\ttail-exact\r\ndata\tfirst_data\r\n");
    irr_stage::fixed_text(root/"second.tsv","data\tsecond_data\ncoverage\tMHV-small-S-T1T2-all\n");
    irr_stage::fixed_text(root/"plan.tsv","one collinear first.tsv\ntwo ope_small_s second.tsv\n");
    settings s;s.helicity="MHV";s.h=10;s.k=2;s.data=root/"default_data";s.plan=root/"plan.tsv";
    auto ops=operations(s);
    require(ops.size()==2&&ops[0].projection=="tail-exact"&&ops[1].projection=="scalar","Stage projection not independent");
    require(ops[0].data==root/"first_data"&&ops[1].data==root/"second_data","Data path not relative to configuration");
    require(ops[1].coverage=="MHV-small-S-T1T2-all","Coverage lost");
    s.k=4;bool failed=false;try{operations(s);}catch(const std::exception&){failed=true;}
    require(failed,"Unsupported exact-tail weight accepted");
    std::cout<<"PASS independent stage options, relative data, CRLF, unsupported projection\n";
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
