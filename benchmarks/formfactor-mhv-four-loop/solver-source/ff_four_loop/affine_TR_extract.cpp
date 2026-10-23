#include "symrep_bootstrap.hpp"
using namespace symrep;
int main(int argc,char**argv){try{
 require(argc==5||argc==6,"usage affine_TR_extract OUT prefix_weight kernel.wxf originTF.wxf [tf|tr|profile]");
 fs::path out=argv[1];size_t w=std::stoul(argv[2]);std::string mode=argc==6?argv[5]:"tf";require(mode=="tf"||mode=="tr"||mode=="profile","invalid extraction mode");
 auto k=read_matrix(argv[3]),o=read_matrix(argv[4]),j=read_matrix(out/"suffix_TR92_to_TF175.wxf");size_t nf=o.nrow,nx=nf*92;require(k.ncol==nx+1&&o.ncol==175,"affine extraction dimensions");
 size_t point=k.nrow,alphaRows=0,rawTR=0,joinedTerms=0;std::vector<size_t>rowTerms(k.nrow);
 for(size_t r=0;r<k.nrow;++r){if(k[r].find(index(nx))){++alphaRows;if(point==k.nrow)point=r;}for(auto[c,x]:k[r])if(size_t(c)<nx){++rawTR;rowTerms[r]+=j[c%92].nnz();}joinedTerms+=rowTerms[r];}
 std::string stem="sewn_affine_w"+std::to_string(w+2);std::ofstream profile(out/(stem+"_extraction_profile.json"));profile<<"{\"augmented_kernel_dimension\":"<<k.nrow<<",\"kernel_nonzeros\":"<<k.nnz()<<",\"alpha_supported_rows\":"<<alphaRows<<",\"raw_TR_coefficients\":"<<rawTR<<",\"TR_to_TF_raw_product_terms\":"<<joinedTerms;
 if(point<k.nrow){size_t pterms=rowTerms[point],pn=k[point].nnz()-1;profile<<",\"origin_TF_product_term_upper_bound\":"<<o.nnz()+pterms<<",\"directions_TF_product_term_upper_bound\":"<<joinedTerms-pterms+(alphaRows-1)*pterms<<",\"directions_TR_term_upper_bound\":"<<rawTR-pn+(alphaRows-1)*pn;}
 profile<<",\"scope\":\"Support-cost preflight before embedding; upper bounds precede exact cancellation. No physical target read.\"}\n";profile.close();if(mode=="profile")return 0;
 std::ofstream receipt(out/(stem+(mode=="tr"?"_TR.json":".json")));
 if(point==k.nrow){receipt<<"{\"consistent\":false,\"augmented_kernel_dimension\":"<<k.nrow<<",\"complete_kernel_file\":\""<<argv[3]<<"\",\"certificate\":\"The complete certified augmented kernel has identically zero affine coordinate, so its t=1 section is empty.\"}\n";std::cout<<"INCONSISTENT complete kernel has zero affine coordinate\n";return 0;}
 auto embed=[&](const Vec&v){Vec z;for(auto[c,x]:v)if(size_t(c)<nx)for(auto[t,y]:j[c%92])z.push_back(index((c/92)*175+t),x*y);else for(size_t b=0;b<nf;++b)for(auto[t,y]:o[b])z.push_back(index(b*175+t),x*y);normalize(z);return z;};
 Q scale=*k[point].find(index(nx));Vec p=std::move(k[point]);for(size_t i=0;i<p.nnz();++i)p[i]/=scale;Mat P(1,nf*175);P[0]=embed(p);size_t tail=mode=="tr"?92:175;Mat N(k.nrow-1,nf*tail);size_t n=0;
 for(size_t r=0;r<k.nrow;++r)if(r!=point){auto v=std::move(k[r]);if(auto a=v.find(index(nx))){Q x=*a;for(auto[c,y]:p)v.push_back(c,-x*y);normalize(v);}require(!v.find(index(nx)),"homogeneous affine coordinate");N[n++]=mode=="tr"?std::move(v):embed(v);}
 MatrixTensorView pt(P,{1,nf,175});write_tensor_view(out/(stem+"_origin.wxf"),pt);MatrixTensorView nt(N,{N.nrow,nf,tail});write_tensor_view(out/(stem+(mode=="tr"?"_directions_TR92.wxf":"_directions.wxf")),nt);
 receipt<<"{\"consistent\":true,\"augmented_kernel_dimension\":"<<k.nrow<<",\"affine_dimension\":"<<N.nrow<<",\"prefix_dimension\":"<<nf<<",\"origin_suffix_TF_dimension\":175,\"direction_suffix_dimension\":"<<tail<<",\"origin_nonzeros\":"<<P.nnz()<<",\"directions_nonzeros\":"<<N.nnz()<<",\"coordinate_map\":\"S=O_invariant*t+X.J; t=1, directions t=0\",\"direction_embedding\":\""<<(mode=="tr"?"Native TR92; multiply suffix_TR92_to_TF175 before adding to TF origin":"Already embedded in TF175")<<"\",\"input_rows_released_during_extraction\":true,\"scope\":\"Complete symmetry, proposed adjacency, integrability and continued remainder-tail equations. DC/TC boundary regularity pending.\"}\n";
 std::cout<<"COMPLETE affine_dimension="<<N.nrow<<" direction_suffix="<<tail<<" origin_nnz="<<P.nnz()<<" directions_nnz="<<N.nnz()<<std::endl;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
