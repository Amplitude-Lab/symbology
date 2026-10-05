#include "symrep_bootstrap.hpp"
#include <random>
using namespace symrep;
int main(){try{
    std::mt19937 rng(48271);rref_option_t opt;opt->pool.reset(2);size_t cases=0;
    bool integral=false;auto random_matrix=[&](size_t nr,size_t nc){Mat a(nr,nc);for(size_t i=0;i<nr;++i)for(size_t j=0;j<nc;++j)if(rng()%3==0){Q x=Q(long(rng()%9)-4)/Q(integral?1:1+rng()%5);if(x!=0)a[i].push_back(index(j),x);}return a;};
    for(size_t trial=0;trial<40;++trial)for(bool backward:{false,true}){
        integral=trial%2==0;
        size_t letters=2+rng()%3,previous=1+rng()%5,older=1+rng()%3;
        auto d=condition_tensor(random_matrix(1+rng()%5,letters*letters),letters);
        auto f=tensor(random_matrix(previous,older*letters),backward?std::vector<size_t>{previous,letters,older}:std::vector<size_t>{previous,older,letters});
        auto reference=assemble(f,d,backward),k=tensor_kernel::extension(f,d,backward,opt);
        require(!mul(reference,transpose(k)).nnz(),"factored kernel has nonzero exact residual");
        require(k.nrow==reference.ncol-row_basis(reference).nrow,"factored kernel dimension mismatch");
        require(row_basis(k).nrow==k.nrow,"factored kernel is dependent");++cases;
    }
    for(size_t trial=0;trial<20;++trial){
        integral=trial%2==0;size_t letters=2+rng()%3,previous=1+rng()%4,older=1+rng()%3,last=1+rng()%4,terminal=1+rng()%4;
        auto dm=random_matrix(1+rng()%4,letters*letters),fm=random_matrix(previous,older*letters),lm=random_matrix(last,letters*terminal);
        auto d=condition_tensor(dm,letters),f=tensor(fm,{previous,older,letters}),l=tensor(lm,{last,letters,terminal});
        // Direct four-index contraction, independent of the local-row builder.
        Mat a(older*dm.nrow*terminal,previous*last);
        for(size_t b=0;b<previous;++b)for(auto [ol,x]:fm[b])for(size_t j=0;j<last;++j)for(auto [mt,y]:lm[j])
            for(size_t q=0;q<dm.nrow;++q)for(auto [ll,z]:dm[q])if(size_t(ll)==(ol%letters)*letters+mt/terminal)
                a[(ol/letters*dm.nrow+q)*terminal+mt%terminal].push_back(index(b*last+j),x*y*z);
        for(auto& row:a.rows)normalize(row);
        auto dc=d,fc=f,lc=l;
        auto raw=flatten(tensor_kernel::sew_one(std::move(dc),std::move(fc),std::move(lc),opt,tensor_kernel::Strategy::staged,{},tensor_kernel::BoundaryReduction::raw));
        require(!mul(a,transpose(raw)).nnz()&&raw.nrow==a.ncol-row_basis(a).nrow&&row_basis(raw).nrow==raw.nrow,"unreduced local interface changed the exact kernel");
        auto k=flatten(tensor_kernel::sew_one(std::move(d),std::move(f),std::move(l),opt));
        require(!mul(a,transpose(k)).nnz(),"multi-terminal sewing residual");
        require(k.nrow==a.ncol-row_basis(a).nrow&&row_basis(k).nrow==k.nrow,"multi-terminal sewing rank");++cases;
    }
    // A residual bound larger than one prime must use the complete CRT
    // modulus, never mistake one-prime congruence for an exact certificate.
    Mat large(1,4);large[0].push_back(0,Q(1));large[0].push_back(1,Q(int_t(2).pow(80ul)+int_t(17)));large[0].push_back(2,Q(3));
    auto bigd=condition_tensor(large,2),bigf=tensor(identity(2),{2,1,2});
    auto biga=assemble(bigf,bigd,false),bigk=tensor_kernel::extension(bigf,bigd,false,opt);
    require(!mul(biga,transpose(bigk)).nnz()&&bigk.nrow==3&&row_basis(bigk).nrow==3,"multi-prime factored certificate");++cases;
    {
        integral=false;auto fm=random_matrix(6,300*4),cm=random_matrix(16,4*4);auto f=tensor(fm,{6,300,4});
        tensor_kernel::ConstraintRows rows(f,cm,4);structured_kernel::Relations relations(24);std::vector<I> columns(24);std::iota(columns.begin(),columns.end(),I(0));
        field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));sparse_mat<ulong,I> serial(0,24),parallel(0,24);
        rows.modular(relations,columns,fp,[&](auto&& row){serial.rows.push_back(std::move(row));});
        rows.modular(relations,columns,fp,[&](auto&& row){parallel.rows.push_back(std::move(row));},&opt->pool);
        require(serial.rows.size()==parallel.rows.size(),"parallel assembly changed row count");
        for(size_t r=0;r<serial.rows.size();++r)require(serial.rows[r]==parallel.rows[r],"parallel assembly changed an equation or its order");++cases;
    }
    std::cout<<"PASS "<<cases<<" rational factored extension/sewing operators against independent exact ranks and residuals"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
