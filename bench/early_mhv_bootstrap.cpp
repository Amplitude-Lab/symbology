// Identify coefficients with necessary small exact collinear images. Then
// independently check EVERY word containing any divergent letter against the
// original shuffle boundary. Discovery projections alone are not a certificate.
#include "compute_rhs.hpp"
#include "packed_shuffle.hpp"
using namespace symrep;
using COO=sparse_tensor<Q,I,SPARSE_COO>;
size_t power(size_t n,size_t w){size_t r=1;for(size_t i=0;i<w;++i){require(r<=SIZE_MAX/n,"word index overflow");r*=n;}return r;}
Mat project_basis(const Tensor& s,const Tensor& last2,const Tensor& last1,const Mat& letter,const fs::path& data,const fs::path& chain,size_t weight,thread_pool* pool){
    size_t alphabet=letter.ncol,suffix=alphabet*alphabet,n=s.dim(0);auto right=mul(flatten(last2),kron(letter,mul(flatten(last1),letter,pool)),pool);
    Mat coefficients(n*suffix,s.dim(1));for(size_t b=0;b<n;++b)for(size_t p=s.rowptr()[b];p<s.rowptr()[b+1];++p){auto i=s.index(p);for(auto [j,x]:right[i[1]])coefficients[b*suffix+j].push_back(i[0],s.val(p)*x);}
    for(auto& row:coefficients.rows)normalize(row);
    for(size_t w=weight;w>0;--w){auto f=read(w==1?data/"FEC_1.wxf":chain/("FEC_"+std::to_string(w)+".wxf"));size_t previous=f.dim(1),letters=f.dim(2);
        auto matrix=flatten(f);f.clear();auto expanded=mul(coefficients,matrix,pool);matrix.clear();coefficients.clear();Mat next(n*suffix*alphabet,previous);
        for(size_t r=0;r<expanded.nrow;++r)for(auto [c,x]:expanded[r])for(auto [l,y]:letter[c%letters])next[(r/suffix)*(suffix*alphabet)+size_t(l)*suffix+r%suffix].push_back(index(c/letters),x*y);
        expanded.clear();for(auto& row:next.rows){normalize(row);row.compress();}suffix*=alphabet;coefficients=std::move(next);
        std::cout<<"SMALL_PROJECTION weight="<<w<<" suffix_words="<<suffix<<" nnz="<<coefficients.nnz()<<std::endl;
    }
    Mat result(n,suffix);for(size_t r=0;r<coefficients.nrow;++r)if(coefficients[r].nnz()){require(coefficients[r].nnz()==1&&coefficients[r](0)==0,"invalid projected terminal");result[r/suffix].push_back(index(r%suffix),coefficients[r][0]);}return result;
}
COO project_words(const Tensor& input,const Mat& projection){
    COO source(input),out(std::vector<size_t>(input.rank(),projection.ncol));std::map<size_t,Q> sums;
    // The discovery projection maps each original letter to at most one
    // letter; no word expansion is needed to project lower-loop references.
    for(const auto& row:projection.rows)require(row.nnz()<=1,"word projection must be monomial");
    for(size_t p=0;p<source.nnz();++p){size_t key=0;Q x=source.val(p);for(auto l:source.index_vector(p)){if(!projection[l].nnz()){x=0;break;}key=key*projection.ncol+projection[l](0);x*=projection[l][0];}if(x!=0)sums[key]+=x;}
    for(auto& [key,x]:sums)if(x!=0){std::vector<I> word(input.rank());size_t k=key;for(size_t i=word.size();i-->0;){word[i]=index(k%projection.ncol);k/=projection.ncol;}out.push_back(word,x);}return out;
}
Mat fit_coefficients(const Mat& basis,const COO& previous,const Mat& first,const Mat& projection,const std::vector<I>& div,size_t weight,rref_option_t opt){
    Mat equations(0,basis.nrow+1);size_t alphabet=projection.ncol,words=power(alphabet,weight-1);auto bt=transpose(basis);
    for(size_t d=0;d<div.size();++d){COO derivative(std::vector<size_t>{alphabet});std::vector<I> letter_index(derivative.rank());require(letter_index.size()==1,"derivative must be a one-letter tensor");
        for(size_t i=0;i<first.nrow;++i)if(auto x=first[i].find(div[d]))for(auto [j,y]:projection[i]){letter_index[0]=j;derivative.push_back(letter_index,(*x)*y);}
        auto rhs=tensor_shuffle_product_parallel(previous,derivative,QQ,nullptr);std::vector<Q> values(words);
        for(size_t p=0;p<rhs.nnz();++p){size_t key=0;for(auto j:rhs.index_vector(p))key=key*alphabet+j;values[key]+=rhs.val(p);}
        for(size_t w=0;w<words;++w){Vec row;for(auto [j,x]:bt[w*alphabet+d])row.push_back(j,x);if(values[w]!=0)row.push_back(index(basis.nrow),values[w]);if(row.nnz())equations.rows.push_back(std::move(row));}
    }
    equations.nrow=equations.rows.size();auto pivots=sparse_mat_rref_reconstruct(equations,opt);Mat solution(1,basis.nrow);size_t rank=0;
    for(const auto& batch:pivots)for(auto [r,c]:batch){require(size_t(c)<basis.nrow,"collinear discovery equations are inconsistent");++rank;if(auto x=equations[r].find(index(basis.nrow)))solution[0].push_back(c,*x);}
    if(rank!=basis.nrow)return Mat(0,basis.nrow);normalize(solution[0]);return solution;
}
Tensor boundary_from_lower(const fs::path& data,const fs::path& lower,size_t loop,thread_pool* pool){
    std::map<size_t,COO> e,r;e.emplace(1,COO(read(data/"E1.wxf")));
    for(size_t l=2;l<loop;++l){e.emplace(l,COO(read(lower/(std::to_string(l)+"loop")/("E"+std::to_string(l)+".wxf"))));r.emplace(l,COO(read(lower/(std::to_string(l)+"loop")/("R"+std::to_string(l)+".wxf"))));}
    packed_shuffle::Accumulator sum(2*loop,e.at(1).dim(0));
    for(size_t k=1;k<loop;++k){sum.add(k==1?e.at(1):r.at(k),e.at(loop-k),Q(long(k),long(loop)));
        std::cout<<"PACKED_BOUNDARY term="<<k<<" stored_words="<<sum.stored_words()<<std::endl;}
    e.clear();r.clear();auto rhs=sum.finish();std::cout<<"PACKED_BOUNDARY_RESULT nnz="<<rhs.nnz()<<std::endl;return Tensor(std::move(rhs),pool);
}
int main(int argc,char** argv){try{
    require(argc==8,"early_mhv_bootstrap boundary|coefficients|verify DATA FEC_CHAIN EARLY_DIR LOWER_OUTPUT OUTPUT LOOP");
    std::string mode=argv[1];fs::path data=fs::absolute(argv[2]),chain=fs::absolute(argv[3]),early=fs::absolute(argv[4]),lower=fs::absolute(argv[5]),out=fs::absolute(argv[6]);
    size_t loop=std::stoul(argv[7]),weight=2*loop,fw=weight-2;rref_option_t opt;opt->pool.reset(8);fs::create_directories(out);
    require(loop>=2&&loop<=5,"this validation driver supports loops two through five");
    if(mode=="boundary"){auto boundary=boundary_from_lower(data,lower,loop,&opt->pool);write(out/("boundary_"+std::to_string(loop)+"L.wxf"),boundary);return 0;}
    auto s=read(early/("EARLY_"+std::to_string(fw)+"p2.wxf")),last2=read(early/"LEC_2.wxf"),last1=read(data/"LEC_1.wxf");
    auto col=read_matrix(data/"colmat42.wxf"),divproj=read_matrix(data/"colprojdiv.wxf");std::vector<I> div;
    for(size_t i=0;i<divproj.nrow;++i)if(divproj[i].nnz())div.push_back(index(i));require(!div.empty(),"empty divergent-letter space");
    auto first=flatten(read(data/"E1.wxf"));auto previous_path=loop==2?data/"E1.wxf":lower/(std::to_string(loop-1)+"loop")/("E"+std::to_string(loop-1)+".wxf");
    if(mode=="coefficients"){
        Mat solution;for(size_t extra=0;extra<2&&solution.nrow==0;++extra){Mat projection(col.ncol,div.size()+extra);for(size_t i=0;i<div.size();++i)projection[div[i]].push_back(index(i),Q(1));
            if(extra)for(size_t i=0;i<projection.nrow;++i)if(!projection[i].nnz())projection[i].push_back(index(div.size()),Q(long(i+1)));
            auto basis=project_basis(s,last2,last1,mul(col,projection),data,chain,fw,&opt->pool);auto previous=project_words(read(previous_path),projection);
            solution=fit_coefficients(basis,previous,first,projection,div,weight,opt);std::cout<<"DISCOVERY alphabet="<<projection.ncol<<" unique="<<(solution.nrow==1)<<std::endl;
        }
        require(solution.nrow==1,"small projections do not determine a unique coefficient vector");write(out/"candidate_coefficients.wxf",solution);
        auto physical=mul(solution,flatten(s),&opt->pool);write_tensor_view(out/"candidate_recursive.wxf",MatrixTensorView(physical,{1,s.dim(1),s.dim(2)}));
        std::cout<<"CANDIDATE saved; full collinear derivative verification is still required"<<std::endl;return 0;
    }
    require(mode=="verify","unknown amplitude mode");
    auto physical=read(out/"candidate_recursive.wxf");
    require(equal(mul(read_matrix(out/"candidate_coefficients.wxf"),flatten(s),&opt->pool),flatten(physical)),"candidate coefficient/basis mismatch");
    auto projection_root=out/"projection";fs::create_directories(projection_root);
    for(size_t w=2;w<fw;++w){auto target=projection_root/("FEC_"+std::to_string(w)+".wxf");if(!fs::exists(target))fs::create_symlink(chain/target.filename(),target);}
    compute_rhs_threads()=8;
    run_bootstrap_cmd(shell_quote(bootstrap_exe_path())+" --project --symmetry collinear --target FEC_"+std::to_string(fw-1),data,projection_root);
    auto prefix=read_matrix(projection_root/"collinear"/("first_w"+std::to_string(fw-1)+".wxf"));
    auto right=mul(flatten(last2),kron(col,mul(flatten(last1),col)),&opt->pool);
    auto f=read(chain/("FEC_"+std::to_string(fw)+".wxf"));auto projected=restricted_projection::sew_one_compact_right(physical,f,prefix,col,right,opt);f.clear();
    auto projected_matrix=flatten(projected);projected.clear();auto start=tensor(projected_matrix,{1,prefix.ncol,col.ncol,col.ncol,col.ncol},&opt->pool);projected_matrix.clear();prefix.clear();
    std::vector<fs::path> bases;for(size_t w=fw-1;w>=2;--w)bases.push_back(projection_root/"collinear"/("first_w"+std::to_string(w)+"_basis.wxf"));
    auto expanded=expand_tensor<Q,I>(std::move(start),bases,QQ,&opt->pool);COO expression(std::move(expanded));expression.reshape(std::vector<size_t>(weight,col.ncol));
    auto actual_path=out/("E"+std::to_string(loop)+"_candidate.wxf");
    {Tensor actual(std::move(expression),&opt->pool);write(actual_path,actual);std::cout<<"FULL_COLLINEAR_EXPRESSION nnz="<<actual.nnz()<<std::endl;}
    Tensor boundary;
    auto boundary_name="boundary_"+std::to_string(loop)+"L.wxf";
    auto reference_boundary=lower/(std::to_string(loop)+"loop")/boundary_name;
    if(fs::exists(reference_boundary)){boundary=read(reference_boundary);std::cout<<"Using independent reference boundary "<<reference_boundary<<std::endl;}
    else if(fs::exists(out/boundary_name)){boundary=read(out/boundary_name);std::cout<<"Using separately computed boundary "<<out/boundary_name<<std::endl;}
    else boundary=boundary_from_lower(data,lower,loop,&opt->pool);
    auto actual=read(actual_path);require(actual.dims()==boundary.dims(),"collinear boundary dimension mismatch");size_t checked=0;
    auto divergent=[&](size_t first,const Tensor& t,size_t p){if(std::find(div.begin(),div.end(),index(first))!=div.end())return true;
        for(auto l:t.index_vector(p))if(std::find(div.begin(),div.end(),l)!=div.end())return true;return false;};
    for(size_t first=0;first<actual.dim(0);++first){size_t i=actual.rowptr()[first],j=boundary.rowptr()[first],ie=actual.rowptr()[first+1],je=boundary.rowptr()[first+1];
        while(true){while(i<ie&&!divergent(first,actual,i))++i;while(j<je&&!divergent(first,boundary,j))++j;if(i==ie&&j==je)break;
            require(i<ie&&j<je,"collinear divergent support mismatch");
            require(actual.index_vector(i)==boundary.index_vector(j)&&actual.val(i)==boundary.val(j),"collinear divergent coefficient mismatch");++i;++j;++checked;
        }
    }
    std::cout<<"PASS complete collinear divergent part; checked_words="<<checked<<std::endl;
    write(out/("hepMHV_"+std::to_string(loop)+"L_recursive.wxf"),physical);
    std::cout<<"CERTIFIED unique "<<loop<<"-loop MHV symbol with full collinear divergent-part verification"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
