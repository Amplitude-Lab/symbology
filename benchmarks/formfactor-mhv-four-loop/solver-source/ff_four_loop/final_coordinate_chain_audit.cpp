#include "symrep_bootstrap.hpp"
using namespace symrep;
int main(int argc,char**argv){try{
 require(argc==2,"usage final_coordinate_chain_audit OUTPUT");fs::path out=argv[1];auto begin=std::chrono::steady_clock::now();
 auto origin=flatten(read(out/"sewn_affine_w8_origin.wxf")),dirs=flatten(read(out/"sewn_affine_w8_directions_TR92.wxf"));auto j=read_matrix(out/"suffix_TR92_to_TF175.wxf");
 auto change=read_matrix(out/"boundary/DC/affine_w8_mixed/complete_parameter_map.wxf");require(change.nrow==1&&change.ncol==4&&dirs.nrow==3,"DC parameter shape");auto alpha=change[0].find(0);require(alpha&&*alpha!=0,"nonempty affine DC section");
 Vec correction;for(auto[c,x]:change[0])if(c>0)for(auto[k,y]:dirs[c-1])correction.push_back(k,x*y/(*alpha));normalize(correction);
 Mat expected=origin;for(auto[k,x]:correction)for(auto[t,y]:j[k%92])expected[0].push_back(index((k/92)*175+t),x*y);normalize(expected[0]);
 auto dc=flatten(read(out/"boundary/DC/affine_w8_mixed/origin.wxf")),tc=flatten(read(out/"boundary/TC/affine_w8_mixed/origin.wxf"));require(equal(expected,dc),"all final DC coordinates reconstructed from complete parameter map");require(equal(dc,tc),"TC preserves complete DC origin");
 auto cm=read_matrix(out/"boundary/TC/affine_w8_mixed/complete_parameter_map.wxf");require(cm.nrow==1&&cm.ncol==1&&cm[0].find(0)&&*cm[0].find(0)!=0,"TC affine section");
 auto dd=read(out/"boundary/DC/affine_w8_mixed/directions_TR92.wxf"),td=read(out/"boundary/TC/affine_w8_mixed/directions_TR92.wxf");require(dd.rank()==3&&td.rank()==3&&dd.dim(0)==0&&td.dim(0)==0,"both final homogeneous kernels empty");
 std::ofstream f(out/"final_coordinate_chain_independent_audit.json");f<<"{\"all_checks_pass\":true,\"DC_parameter_map_dimensions\":[1,4],\"TC_parameter_map_dimensions\":[1,1],\"DC_origin_reconstructed_from_original_sewn_family\":true,\"TC_origin_exactly_equal_DC_origin\":true,\"DC_and_TC_affine_sections_nonempty\":true,\"both_final_homogeneous_dimensions\":0,\"final_origin_nonzeros\":"<<dc.nnz()<<",\"all_origin_coefficients_checked\":true,\"scope\":\"Exact coordinate transport audit, separate from boundary block and physical holdout certificates. No full word expansion.\",\"seconds\":"<<seconds(begin)<<"}\n";std::cout<<"PASS final origin and empty directions; seconds="<<seconds(begin)<<std::endl;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
