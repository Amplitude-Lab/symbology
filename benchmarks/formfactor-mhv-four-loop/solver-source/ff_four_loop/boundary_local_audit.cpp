#include "symrep_bootstrap.hpp"
#include <fstream>
using namespace symrep;
int main(int argc,char**argv){try{
 require(argc==2,"usage: boundary_local_audit OUT");fs::path out=argv[1];rref_option_t opt;opt->pool.reset(1);
 auto start=std::chrono::steady_clock::now();auto d=read(out/"data/interior.wxf"),l=read(out/"LEC_2.wxf");
 auto raw=tensor_kernel::right_boundary_conditions(d,l,tensor_kernel::BoundaryReduction::raw);d.clear();l.clear();
 write(out/"boundary_local_raw.wxf",raw);std::cout<<"RAW saved rows="<<raw.nrow<<" columns="<<raw.ncol<<" nnz="<<raw.nnz()<<" seconds="<<seconds(start)<<std::endl;
 auto begin=std::chrono::steady_clock::now();auto rr=raw;rr.sort_rows_by_nnz();auto piv=sparse_mat_rref_reconstruct(rr,opt);
 std::map<I,size_t> where;std::vector<I> oldrows;for(auto&batch:piv)for(auto [r,c]:batch){where[c]=oldrows.size();oldrows.push_back(r);}
 Mat basis(oldrows.size(),raw.ncol);for(size_t r=0;r<oldrows.size();++r)basis[r]=std::move(rr[oldrows[r]]);rr.clear();
 for(size_t r=0;r<basis.nrow;++r)for(auto[c,x]:basis[r])if(where.contains(c))require(where[c]==r&&x==1,"RREF pivot block is not the identity");
 std::cout<<"RREF rank="<<basis.nrow<<" nnz="<<basis.nnz()<<" seconds="<<seconds(begin)<<std::endl;write(out/"boundary_local_rref.wxf",basis);
 // Exact full-row certificate: raw = raw[:, pivot columns] * basis.
 // The independently computed modular rank in rref_reconstruct equals the
 // number of rows, so this containment and the identity pivot block certify
 // equality of the two row spaces, not merely a sampled residual check.
 begin=std::chrono::steady_clock::now();size_t verified=0;Mat coefficients(raw.nrow,basis.nrow);
 for(size_t r=0;r<raw.nrow;++r){Vec image;for(auto[c,x]:raw[r])if(where.contains(c)){auto k=where[c];coefficients[r].push_back(index(k),x);for(auto[j,y]:basis[k])image.push_back(j,x*y);}normalize(image);normalize(coefficients[r]);require(image.nnz()==raw[r].nnz(),"row-space residual support");for(size_t j=0;j<image.nnz();++j)require(image(j)==raw[r](j)&&image[j]==raw[r][j],"row-space residual coefficient");++verified;}
 write(out/"boundary_local_coefficients.wxf",coefficients);std::cout<<"EXACT raw_equals_coefficients_times_rref rows="<<verified<<" pass=1 seconds="<<seconds(begin)<<" total_seconds="<<seconds(start)<<std::endl;
 }catch(const std::exception&e){std::cerr<<e.what()<<std::endl;return 1;}return 0;}
