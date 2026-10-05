#include "early_sew.hpp"
using namespace symrep;
int main(){try{
    rref_option_t opt;opt->pool.reset(3);
    Mat p(3,3),s(3,3),change=identity(3);p[0].push_back(1,Q(1));p[1].push_back(2,Q(1));p[2].push_back(0,Q(1));
    s[0].push_back(1,Q(1));s[1].push_back(0,Q(1));s[2].push_back(2,Q(1));change[0].push_back(1,Q(1,3));
    auto ci=chart_inverse(change,opt);std::vector<Mat> a{mul(mul(change,p),ci),mul(mul(change,s),ci)},b{p,s};
    for(size_t trial=0;trial<16;++trial){
        Mat fm(3,4);for(size_t i=0;i<3;++i)for(size_t j=0;j<4;++j){long x=long((i*7+j*11+trial*5)%9)-4;if(x)fm[i].push_back(index(j),Q(x,long(j%3+1)));}
        auto f=tensor(fm,{3,2,2});Mat local(trial?2:0,6);for(size_t i=0;i<local.nrow;++i)for(size_t j=0;j<6;++j){long x=long((i*5+j*3+trial)%7)-3;if(x)local[i].push_back(index(j),Q(x,long(i+1)));}
        tensor_kernel::ConstraintRows eq(f,local,3);early_sew::ProductInvarianceRows inv(a,b);early_sew::Rows source{eq,inv};
        Mat explicit_a(0,9);structured_kernel::Relations rel(9);eq(rel,SIZE_MAX,[&](Vec&& v){explicit_a.rows.push_back(std::move(v));});
        for(size_t g=0;g<a.size();++g){auto h=transpose(kron(a[g],b[g]));for(size_t c=0;c<9;++c){h[c].push_back(index(c),Q(-1));normalize(h[c]);explicit_a.rows.push_back(std::move(h[c]));}}
        explicit_a.nrow=explicit_a.rows.size();auto expected=kernel(explicit_a,opt);
        auto actual=structured_kernel::solve_stream(9,source,opt,nullptr,[&](const Mat& k,ulong prime,const int_t& modulus){return inv.verify(k,modulus)&&eq.verify(k,opt,true,prime,&modulus);});
        require(actual.nrow==expected.nrow&&mul(actual,transpose(explicit_a)).nnz()==0&&row_basis(actual).nrow==actual.nrow,"early invariance kernel differs from explicit product");
        early_sew::RefinedRows staged{source,opt,trial%2?1ul:2ul,trial%2?64ul:2ul,true,trial%3==0,trial%4?0ul:8ul};
        auto staged_result=structured_kernel::solve_stream(9,staged,opt,nullptr,[&](const Mat& k,ulong prime,const int_t& modulus){return inv.verify(k,modulus)&&eq.verify(k,opt,true,prime,&modulus);});
        require(staged_result.nrow==expected.nrow&&mul(staged_result,transpose(explicit_a)).nnz()==0&&row_basis(staged_result).nrow==staged_result.nrow,"staged invariant kernel differs from explicit product");
        if(!trial)require(actual.nrow==2,"lost paired nontrivial S3 representations");
        // Bounded discovery is only a subset; its rows must be identical to
        // complete rows, never truncated equations.
        field_t fp(FIELD_Fp,n_nextprime(1ULL<<60,0));std::vector<I> columns(9);std::iota(columns.begin(),columns.end(),0);
        sparse_mat<ulong,I> candidate(3,9);for(size_t i=0;i<3;++i)for(size_t j=0;j<9;++j)candidate[i].push_back(index(j),ulong((i*11+j*7+trial)%19+1));
        std::vector<sparse_vec<ulong,I>> residuals;eq.modular_residual(candidate,fp,[&](auto&& row){residuals.push_back(std::move(row));},&opt->pool);
        for(size_t g=0;g<a.size();++g)inv.modular_residual(candidate,fp,g,[&](auto&& row){residuals.push_back(std::move(row));},&opt->pool);
        size_t ri=0;for(const auto& row:explicit_a.rows){auto mr=row%fp.mod;sparse_vec<ulong,I> want;for(size_t i=0;i<candidate.nrow;++i){ulong x=0;for(auto [j,y]:mr)x=nmod_add(x,nmod_mul(y,*candidate[i].find(j),fp.mod),fp.mod);if(x)want.push_back(index(i),x);}
            if(want.nnz()){require(ri<residuals.size()&&want.nnz()==residuals[ri].nnz(),"factored residual support mismatch");for(size_t p=0;p<want.nnz();++p)require(want(p)==residuals[ri](p)&&want[p]==residuals[ri][p],"factored residual differs from explicit operator");++ri;}}
        require(ri==residuals.size(),"extra factored residual rows");
        std::set<std::string> full;source.modular(rel,columns,fp,[&](auto&& v){std::ostringstream o;for(auto [j,x]:v)o<<j<<':'<<x<<',';full.insert(o.str());},&opt->pool);
        source.max_terms=4;source.modular(rel,columns,fp,[&](auto&& v){std::ostringstream o;for(auto [j,x]:v)o<<j<<':'<<x<<',';require(full.count(o.str()),"selected equation differs from complete source");},&opt->pool);
        std::vector<unsigned char> seen_equations(eq.nrow()),seen_invariants(inv.nrow());source.equation_seen=&seen_equations;source.invariant_seen=&seen_invariants;std::set<std::string> collected;
        auto collect=[&](auto&& v){std::ostringstream o;for(auto [j,x]:v)o<<j<<':'<<x<<',';collected.insert(o.str());};
        source.max_terms=2;source.modular(rel,columns,fp,collect,&opt->pool);source.max_terms=SIZE_MAX;source.modular(rel,columns,fp,collect,&opt->pool);
        require(collected==full,"equation discovery mask omitted a complete operator row");size_t repeated=0;source.modular(rel,columns,fp,[&](auto&&){++repeated;},&opt->pool);require(!repeated,"discovery mask repeated consumed equations");
    }
    // A loose global operator bound must not force another rank solve when
    // the small reconstructed kernel is already exact. Also skip a prime
    // dividing an action denominator, and reject a spurious congruent lift.
    for(size_t trial=0;trial<2;++trial){ulong prime=n_nextprime(1ULL<<60,0),next=n_nextprime(prime,0);
        Mat g=identity(3);g[1].push_back(2,trial?Q(1)/Q(next):Q(int_t(2).pow(100ul)));g[2][0]=Q(-1);
        Mat f(3,1);f[1].push_back(0,Q(1));auto ft=tensor(f,{3,1,1});Mat local(1,1);local[0].push_back(0,Q(1));
        tensor_kernel::ConstraintRows eq(ft,local,1);early_sew::ProductInvarianceRows inv({g},{identity(1)});early_sew::Rows rows{eq,inv};size_t calls=0;
        auto k=structured_kernel::solve_stream(3,rows,opt,nullptr,[&](const Mat& candidate,ulong p,const int_t& modulus){++calls;return early_sew::verify_with_extra_primes(candidate,rows,opt,p,modulus);});
        require(calls==1&&k.nrow==1&&k[0].nnz()==1&&k[0](0)==0,"extra residue certificate repeated elimination or changed the kernel");
        Mat bad(1,3);bad[0].push_back(0,Q(1));bad[0].push_back(2,Q(prime));require(!early_sew::verify_with_extra_primes(bad,rows,opt,prime,int_t(prime)),"extra residue certificate accepted a spurious congruent lift");
    }
    std::cout<<"PASS 16 exact product-invariance intersections, nontrivial irrep pairing, and bounded equation subsets"<<std::endl;
    std::cout<<"PASS extra-prime residual certificates, denominator-prime rejection and false-lift rejection"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
