// Exact affine regularity, one suffix-pair block at a time. The source and
// all previously imposed equations are transported by the same row maps.
#include "symrep_io.hpp"
using namespace symrep;
static double elapsed(std::chrono::steady_clock::time_point s){return std::chrono::duration<double>(std::chrono::steady_clock::now()-s).count();}
int main(int argc,char**argv){try{
 require(argc==8,"usage boundary_affine_solve EXPERIMENT BOUNDARY prefix_weight finite_letters origin.wxf directions.wxf OUTPUT");fs::path exp=argv[1],dir=exp/"boundary"/argv[2],out=argv[7];fresh(out);size_t pw=std::stoul(argv[3]),finite=std::stoul(argv[4]);auto start=std::chrono::steady_clock::now();rref_option_t opt;opt->pool.reset(1);
 auto a=read(argv[5]),b=read(argv[6]);require(a.rank()==3&&a.dim(0)==1&&b.rank()==3&&a.dim(1)==b.dim(1)&&a.dim(2)==175&&b.dim(2)==175,"candidate prefix/suffix shapes");size_t original=b.dim(0)+1,nf=a.dim(1);Mat z(original,nf*175);z[0]=flatten(a)[0];auto bb=flatten(b);for(size_t i=0;i<bb.nrow;++i)z[i+1]=std::move(bb[i]);a.clear();b.clear();bb.clear();Mat change=identity(original);
 auto p=read_matrix(dir/"collinear"/("first_w"+std::to_string(pw)+".wxf"));require(p.nrow==nf,"prefix map input dimension");auto suffix=read_matrix(dir/"suffix_TF175_raw2.wxf"),quotient=read_matrix(dir/("prefix"+std::to_string(pw)+"_mod_finite.wxf"));auto h0=read(dir/("H"+std::to_string(pw+2)+"_mod_finite_prefix"+std::to_string(pw)+"_raw2.wxf"));size_t q=h0.dim(1);require(h0.rank()==3&&h0.dim(0)==p.ncol&&h0.dim(2)==q&&suffix.ncol==q*q&&suffix.nrow==175&&quotient.nrow==p.ncol,"boundary preparation shapes");auto h=transpose(flatten(h0)),s=transpose(suffix);h0.clear();suffix.clear();
 std::vector<size_t>pairs;size_t skipped=0;for(size_t c=0;c<q*q;++c)if(s[c].nnz()||h[c].nnz())pairs.push_back(c);else++skipped;
 auto bad=[&](size_t c){return c/q>=finite||c%q>=finite;};std::stable_sort(pairs.begin(),pairs.end(),[&](size_t i,size_t j){if(bad(i)!=bad(j))return bad(i)>bad(j);return s[i].nnz()+h[i].nnz()>s[j].nnz()+h[j].nnz();});
 std::ofstream report(out/"blocks.jsonl");size_t tested=0,totalcomponents=0,changes=0;bool consistent=true;
 for(size_t pair:pairs){auto tick=std::chrono::steady_clock::now();std::vector<Q>weights(175);for(auto[j,x]:s[pair])weights[j]=x;Mat raw(z.nrow,nf);
  for(size_t r=0;r<z.nrow;++r){for(auto[j,x]:z[r])if(weights[j%175]!=0)raw[r].push_back(index(j/175),x*weights[j%175]);normalize(raw[r]);}
  auto condition=mul(raw,p);raw.clear();for(size_t r=0;r<condition.nrow;++r)if(auto alpha=change[r].find(0)){for(auto[j,x]:h[pair])condition[r].push_back(j,-(*alpha)*x);normalize(condition[r]);}
  if(!bad(pair))condition=mul(condition,quotient);++tested;totalcomponents+=condition.ncol;size_t before=z.nrow;
  if(condition.nnz()){
   auto k=kernel(transpose(condition),opt);require(mul(k,condition).nnz()==0,"complete exact affine boundary block residual");
   change=mul(k,change);z=mul(k,z);++changes;write(out/("change_after_block"+std::to_string(tested)+".wxf"),change);
   bool affine=false;for(auto&row:change.rows)if(auto c=row.find(0);c&&*c!=0)affine=true;
   if(!affine){consistent=false;write(out/"inconsistent_block.wxf",condition);write(out/"inconsistent_section_map.wxf",change);}
  }
  report<<"{\"block\":"<<tested<<",\"pair\":["<<pair/q+1<<','<<pair%q+1<<"],\"spurious_suffix\":"<<(bad(pair)?"true":"false")<<",\"components\":"<<condition.ncol<<",\"before\":"<<before<<",\"after\":"<<z.nrow<<",\"exact_kernel_residual_zero\":true,\"seconds\":"<<elapsed(tick)<<"}\n";report.flush();
  if(before!=z.nrow||tested%20==0)std::cout<<"BOUNDARY_BLOCK checked="<<tested<<'/'<<pairs.size()<<" augmented_dimension="<<z.nrow<<" elapsed="<<elapsed(start)<<std::endl;
  if(!consistent)break;
 }
 if(consistent){size_t point=z.nrow;for(size_t r=0;r<change.nrow;++r)if(auto c=change[r].find(0);c&&*c!=0){point=r;break;}require(point<z.nrow,"affine section exists");Q scale=*change[point].find(0);Mat origin(1,z.ncol),directions(z.nrow-1,z.ncol);origin[0]=z[point];for(size_t j=0;j<origin[0].nnz();++j)origin[0][j]/=scale;size_t d=0;
  for(size_t r=0;r<z.nrow;++r)if(r!=point){directions[d]=z[r];if(auto alpha=change[r].find(0))for(auto[j,x]:origin[0])directions[d].push_back(j,-(*alpha)*x);normalize(directions[d++]);}
  write_tensor_view(out/"origin.wxf",MatrixTensorView(origin,{1,nf,175}));write_tensor_view(out/"directions.wxf",MatrixTensorView(directions,{directions.nrow,nf,175}));write(out/"complete_parameter_map.wxf",change);
 }
 std::ofstream receipt(out/"receipt.json");receipt<<"{\"consistent\":"<<(consistent?"true":"false")<<",\"prefix_weight\":"<<pw<<",\"initial_affine_dimension\":"<<original-1<<",\"final_augmented_dimension\":"<<z.nrow<<",\"all_suffix_pairs_checked\":"<<(tested==pairs.size()?"true":"false")<<",\"nonzero_input_pair_blocks\":"<<pairs.size()<<",\"checked_pair_blocks\":"<<tested<<",\"exact_zero_input_pair_blocks\":"<<skipped<<",\"full_coordinate_components_checked\":"<<totalcomponents<<",\"coordinate_changes\":"<<changes<<",\"all_spurious_words_retained\":true,\"all_block_kernels_exactly_verified\":true,\"finite_targets_read\":false,\"affine_RHS_transported_with_each_change\":true,\"seconds\":"<<elapsed(start)<<"}\n";
 std::cout<<"BOUNDARY_COMPLETE consistent="<<consistent<<" augmented_dimension="<<z.nrow<<" seconds="<<elapsed(start)<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
