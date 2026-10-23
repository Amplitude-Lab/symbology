// Joint prefix coordinates for the candidate boundary image and the complete
// affine normalization source. All new directions are retained, including
// finite corrections and any source direction outside the candidate image.
#include "boundary_tensor_train.hpp"
using namespace boundary_tt;
static double elapsed(std::chrono::steady_clock::time_point s){return std::chrono::duration<double>(std::chrono::steady_clock::now()-s).count();}
static void save_or_check(const fs::path&p,const Tensor&t){if(fs::exists(p))require(equal(flatten(read(p)),flatten(t)),"cached joint tensor mismatch");else write(p,t);}
int main(int argc,char**argv){try{
 require(argc==4||argc==5,"usage boundary_joint_prefix BOUNDARY_DIR maxweight OUTPUT_DIR [H6|H8]");std::string source=argc==5?argv[4]:"H8";require(source=="H6"||source=="H8","normalization source");fs::path input=argv[1],out=argv[3];size_t maxw=std::stoul(argv[2]);fs::create_directories(out/"collinear");rref_option_t opt;opt->pool.reset(1);auto start=std::chrono::steady_clock::now();
 std::ifstream cfg(input/"train"/source/"shape.txt");size_t q,w;cfg>>q>>w;std::vector<size_t>rank(w+1);for(auto&r:rank)cfg>>r;require(bool(cfg)&&(w==8||w==6)&&maxw<=w-2,"H8 train shape");
 auto core=read(input/"train"/source/"core_1.wxf");auto c=flatten(core);Mat u(rank[1],q);for(auto[j,x]:c[0])u[j%rank[1]].push_back(index(j/rank[1]),x);for(auto&row:u.rows)normalize(row);
 size_t old=q;std::ofstream report(out/"joint_prefix_progress.txt",std::ios::app);
 for(size_t k=2;k<=maxw;++k){
  auto t=read(input/"collinear"/("first_w"+std::to_string(k)+"_basis.wxf"));require(t.dim(2)==q&&t.dim(1)<=old,"candidate boundary recurrence");size_t np=t.dim(0),prev=t.dim(1);Mat b=flatten(t);t.clear();b.ncol=old*q;
  // A complete private-coordinate chart of the existing independent rows.
  std::vector<unsigned char>count(prev*q,0);for(auto&row:b.rows)for(auto[j,x]:row)if(count[j]<2)++count[j];
  std::vector<I>piv(np,-1),col(old*q,-1);std::vector<Q>scale(np);
  for(size_t i=0;i<np;++i)for(auto[j,x]:b[i])if(count[j]==1){piv[i]=j;scale[i]=x;col[j]=index(i);break;}
  require(std::all_of(piv.begin(),piv.end(),[](I i){return i>=0;}),"candidate recurrence has complete private chart");
  core=read(input/"train"/source/("core_"+std::to_string(k)+".wxf"));c=flatten(core);require(c.nrow==u.nrow&&c.ncol==q*rank[k],"source train recurrence");
  Mat next(rank[k],old*q);
  for(size_t a=0;a<c.nrow;++a)for(auto[z,x]:c[a]){size_t l=z/rank[k],s=z%rank[k];for(auto[j,y]:u[a])next[s].push_back(index(size_t(j)*q+l),x*y);}
  for(auto&row:next.rows)normalize(row);
  Mat coeff(next.nrow,np);for(size_t s=0;s<next.nrow;++s)for(auto[j,x]:next[s])if(col[j]>=0)coeff[s].push_back(col[j],x/scale[col[j]]);for(auto&row:coeff.rows)normalize(row);
  auto reconstructed=mul(coeff,b);Mat residual=next;
  for(size_t i=0;i<residual.nrow;++i){for(auto[j,x]:reconstructed[i])residual[i].push_back(j,-x);normalize(residual[i]);}reconstructed.clear();next.clear();
  auto [extra_coeff,extra_basis]=factor(residual,opt);const auto ne=extra_basis.nrow;residual.clear();
  Mat joint(np+ne,old*q);for(size_t i=0;i<np;++i)joint[i]=std::move(b[i]);for(size_t i=0;i<ne;++i)joint[np+i]=std::move(extra_basis[i]);b.clear();extra_basis.clear();
  u=Mat(rank[k],np+ne);for(size_t s=0;s<u.nrow;++s){u[s]=std::move(coeff[s]);for(auto[j,x]:extra_coeff[s])u[s].push_back(index(np+j),x);normalize(u[s]);}
  auto path=out/"collinear"/("first_w"+std::to_string(k)+"_basis.wxf");save_or_check(path,tensor(joint,{np+ne,old,q}));
  auto upath=out/("source_prefix_w"+std::to_string(k)+".wxf");if(fs::exists(upath))require(equal(u,read_matrix(upath)),"cached source prefix map");else write(upath,u);
  std::ofstream js(out/("joint_prefix_w"+std::to_string(k)+".json"));js<<"{\"weight\":"<<k<<",\"candidate_dimension\":"<<np<<",\"additional_source_dimension\":"<<ne<<",\"joint_dimension\":"<<np+ne<<",\"source_prefix_state_count\":"<<rank[k]<<",\"source_coordinate_nonzeros\":"<<u.nnz()<<",\"exact_full_source_reconstruction\":true,\"all_source_directions_retained\":true,\"finite_target_read\":false,\"seconds_since_start\":"<<elapsed(start)<<"}\n";
  report<<k<<' '<<np<<' '<<ne<<' '<<u.nnz()<<' '<<elapsed(start)<<std::endl;std::cout<<"JOINT weight="<<k<<" candidate="<<np<<" extra="<<ne<<" source_nnz="<<u.nnz()<<" seconds="<<elapsed(start)<<std::endl;old=np+ne;
 }
 if(maxw==w-2){auto c7=flatten(read(input/"train"/source/("core_"+std::to_string(w-1)+".wxf"))),c8=flatten(read(input/"train"/source/("core_"+std::to_string(w)+".wxf")));Mat suffix(rank[w-2],q*q);
  for(size_t a=0;a<rank[w-2];++a){for(auto[j,x]:c7[a])for(auto[l,y]:c8[j%rank[w-1]])suffix[a].push_back(index((j/rank[w-1])*q+l),x*y);normalize(suffix[a]);}
  auto h=mul(transpose(u),suffix);save_or_check(out/(source+"_joint"+std::to_string(maxw)+"_raw2.wxf"),tensor(h,{old,q,q}));
  std::ofstream receipt(out/(source+"_joint"+std::to_string(maxw)+"_raw2.json"));receipt<<"{\"all_exact_prefix_reconstructions_pass\":true,\"joint_dimension\":"<<old<<",\"raw_suffix_dimension\":"<<q*q<<",\"full_source_nonzeros\":"<<h.nnz()<<",\"finite_R2_squared_term_retained\":"<<(source=="H8"?"true":"null")<<",\"no_finite_correction_discarded\":true,\"F4_target_read\":false,\"seconds\":"<<elapsed(start)<<"}\n";
 }
}catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
