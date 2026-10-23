#include "symrep_bootstrap.hpp"
using namespace symrep;
int main(int argc,char**argv){try{
 require(argc==2,"usage ES_suffix_audit RESEARCH_ROOT");fs::path r=argv[1],out=r/"results/ff_four_loop/ES_comparison";
 rref_option_t opt;opt->pool.reset(1);auto start=std::chrono::steady_clock::now();
 auto tr=flatten(read(r/"results/ff_endpoint_literature_LEC_R.wxf"));
 auto es=flatten(read(out/"LEC_TR_2.wxf")),ca=flatten(read(r/"results/ff_four_loop/LEC_TR_2.wxf"));
 require(tr.nrow==12&&tr.ncol==88&&es.nrow==120&&es.ncol==1056&&ca.nrow==92&&ca.ncol==1056,"suffix basis dimensions");
 auto expand=kron(identity(88),tr),physical=mul(es,expand),caPhysical=mul(ca,expand);
 Chart chart(es,opt);auto inclusion=chart.coordinates(ca);require(equal(mul(inclusion,es),ca),"exact CA92 to ES120 inclusion");
 auto conditions=condition_rows(read(r/"results/ff_bootstrap/data/interior.wxf"));
 require(mul(physical,transpose(conditions)).nnz()==0,"all physical ES and integrability equations");
 write(out/"CA92_to_ES120.wxf",inclusion);write(out/"ES120_physical_pair.wxf",physical);
 std::ofstream receipt(out/"suffix_independent_audit.json");receipt<<"{\"all_checks_pass\":true,\"ES_suffix_dimension\":120,\"proposed_suffix_dimension\":92,\"inclusion_nonzeros\":"<<inclusion.nnz()<<",\"full_7744_pair_coefficients_verified\":true,\"all_ES_integrability_equations_verified\":true,\"boundaries\":[";
 for(size_t b=0;b<2;++b){size_t q=b?29:22,f=b?9:6;auto map=read_matrix(r/(b?"results/ff_boot_triple_map.wxf":"results/ff_boot_boundary_map.wxf"));require(map.nrow==88&&map.ncol==q,"boundary convention");auto projection=mul(physical,kron(map,map));write(out/(b?"suffix_TR120_TC_raw2.wxf":"suffix_TR120_DC_raw2.wxf"),projection);auto basis=row_basis(projection);Mat spur(basis.nrow,q*q);for(size_t i=0;i<basis.nrow;++i)for(auto[j,x]:basis[i])if(size_t(j)/q>=f||size_t(j)%q>=f)spur[i].push_back(j,x);auto bad=row_basis(spur);if(b)receipt<<',';receipt<<"{\"name\":\""<<(b?"TC":"DC")<<"\",\"raw_letters\":"<<q<<",\"finite_letters\":"<<f<<",\"image_rank\":"<<basis.nrow<<",\"regular_image_rank\":"<<basis.nrow-bad.nrow<<",\"zero_image_kernel_dimension\":"<<120-basis.nrow<<",\"projected_nonzeros\":"<<projection.nnz()<<"}";}
 receipt<<"],\"scope\":\"Exact suffix audit only; no weaker ES prefix6 or four-loop solve. Native suffix completeness certificate is separate in suffix_TR_w2.log.\",\"seconds\":"<<seconds(start)<<"}\n";
 std::cout<<"PASS ES120 contains CA92; all equations checked; seconds="<<seconds(start)<<std::endl;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
