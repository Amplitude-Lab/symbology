// Remove only certified wholly finite corrections from the affine source.
#include "boundary_tensor_train.hpp"
using namespace boundary_tt;
int main(int argc,char**argv){try{
 require(argc==4||argc==6,"usage boundary_affine_prepare EXPERIMENT BOUNDARY finite_count [prefix_weight SOURCE]");size_t weight=argc==6?std::stoul(argv[4]):6;std::string source=argc==6?argv[5]:"H8";fs::path experiment=argv[1],dir=experiment/"boundary"/argv[2],joint=dir/(source=="H8"?"joint":"joint_H6");size_t f=std::stoul(argv[3]);rref_option_t opt;opt->pool.reset(1);
 auto pt=read(dir/"collinear"/("first_w"+std::to_string(weight)+".wxf"));require(pt.rank()==2,"prefix projection rank");size_t np=pt.dim(1);pt.clear();
 auto ht=read(joint/(source+"_joint"+std::to_string(weight)+"_raw2.wxf"));size_t nj=ht.dim(0),q=ht.dim(1),ne=nj-np;require(ht.rank()==3&&ht.dim(2)==q&&nj>=np,"joint H8 shape");auto h=flatten(ht);ht.clear();
 auto gj=read_matrix(joint/("regular_prefix_w"+std::to_string(weight)+".wxf")),gp=read_matrix(dir/("regular_prefix_w"+std::to_string(weight)+".wxf"));require(gj.ncol==nj&&gp.ncol==np,"finite prefix shape");
 Mat extra(gj.nrow,ne);for(size_t i=0;i<gj.nrow;++i)for(auto[j,x]:gj[i])if(size_t(j)>=np)extra[i].push_back(index(j-np),x);
 Mat hx(ne,q*q);for(size_t i=0;i<ne;++i){hx[i]=h[np+i];for(auto[j,x]:hx[i])require(size_t(j)/q<f&&size_t(j)%q<f,"source outside candidate prefix has a spurious suffix: affine boundary obstruction");}
 Mat finite_lift(ne,nj),correction(nj,q*q);
 if(ne){Chart chart(transpose(extra),opt);require(chart.pivots.size()==ne,"all added source directions have finite representatives");Mat square(ne,ne),selected(ne,nj);for(size_t i=0;i<ne;++i){square[i]=extra[chart.pivots[i]];selected[i]=gj[chart.pivots[i]];}
  finite_lift=mul(chart_inverse(square,opt),selected);Mat check(ne,ne);for(size_t i=0;i<ne;++i)for(auto[j,x]:finite_lift[i])if(size_t(j)>=np)check[i].push_back(index(j-np),x);require(equal(check,identity(ne)),"finite lift auxiliary block identity");correction=mul(transpose(finite_lift),hx);
 }
 for(size_t i=0;i<nj;++i){for(auto[j,x]:correction[i])h[i].push_back(j,-x);normalize(h[i]);if(i>=np)require(h[i].nnz()==0,"finite correction removes full auxiliary source");}
 h.rows.resize(np);h.nrow=np;
 write(joint/"finite_correction_prefix_lift.wxf",finite_lift);write(joint/"finite_correction_suffix_coefficients.wxf",hx);write(dir/(source+"_mod_finite_prefix"+std::to_string(weight)+"_raw2.wxf"),tensor(h,{np,q,q}));
 Chart gc(gp,opt);auto normal=mul(gc.inverse,gp);std::vector<I>piv(np,-1),free(np,-1);for(size_t k=0;k<gc.pivots.size();++k)piv[gc.pivots[k]]=index(k);size_t nf=0;for(size_t p=0;p<np;++p)if(piv[p]<0)free[p]=index(nf++);Mat quotient(np,nf);
 for(size_t p=0;p<np;++p)if(free[p]>=0)quotient[p].push_back(free[p],Q(1));else for(auto[j,x]:normal[piv[p]])if(free[j]>=0)quotient[p].push_back(free[j],-x);require(mul(gp,quotient).nnz()==0,"finite prefix quotient");write(dir/("prefix"+std::to_string(weight)+"_mod_finite.wxf"),quotient);
 auto map=read_matrix(dir/"data/colmat42.wxf"),physical=read_matrix(experiment/"suffix_TF175_physical.wxf");require(map.nrow==88&&map.ncol==q&&physical.nrow==175&&physical.ncol==88*88,"suffix physical map shapes");Mat suffix(175,q*q);
 for(size_t r=0;r<175;++r){for(auto[j,x]:physical[r])for(auto[a,y]:map[j/88])for(auto[b,z]:map[j%88])suffix[r].push_back(index(size_t(a)*q+b),x*y*z);normalize(suffix[r]);}if(!fs::exists(dir/"suffix_TF175_raw2.wxf"))write(dir/"suffix_TF175_raw2.wxf",suffix);else require(equal(suffix,read_matrix(dir/"suffix_TF175_raw2.wxf")),"cached boundary suffix map");
 std::ofstream manifest(dir/("boundary_affine_preparation_w"+std::to_string(weight+2)+".json"));manifest<<"{\"all_checks_pass\":true,\"boundary_alphabet\":"<<q<<",\"finite_letters\":"<<f<<",\"prefix_dimension\":"<<np<<",\"regular_prefix_dimension\":"<<gp.nrow<<",\"additional_source_directions\":"<<ne<<",\"finite_correction_suffix_nonzeros\":"<<hx.nnz()<<",\"corrected_source_nonzeros\":"<<h.nnz()<<",\"source_auxiliary_spurious_suffix_zero\":true,\"source_correction_wholly_finite\":true,\"full_finite_correction_saved\":true,\"F4_target_read\":false,\"constraint_recipe\":\"Every suffix pair containing a spurious letter keeps all prefix coordinates; every finite suffix pair is projected by prefix6_mod_finite.\"}\n";
 std::cout<<"BOUNDARY_AFFINE_READY prefix="<<np<<" regular="<<gp.nrow<<" finite_correction_rank="<<ne<<" source_nnz="<<h.nnz()<<std::endl;
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
