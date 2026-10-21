#include "irr_physics_pipeline.hpp"
#include "irr_bulk.hpp"
using irr::require;
using MM=sparse_mat<ulong,int32_t>;

int_t factor_norm(const irr::matrices& blocks,int_t& den) {
    den=1;for(const auto& b:blocks)for(const auto& row:b.rows)for(auto [j,x]:row)den=LCM(den,x.den());
    int_t sum=0;for(const auto& b:blocks)for(const auto& row:b.rows)for(auto [j,x]:row){int_t v=den;v/=x.den();v*=x.num().abs();sum+=v;}return sum;
}
int main(int argc,char**argv)try {
    require(argc==7,"joint_certificate H K QQ_ROOT FEC LEC DATA");
    Flint::set_memory_functions();size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]);
    fs::path root=argv[3],fec=argv[4],lec=argv[5],data=argv[6];
    rref_option_t opt;irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;auto* pool=&opt->pool;
    irr_stage::verify_complete(root);irr_stage::verify_basis(root,irr_stage::basis_identity(fec,lec,h,k),h,k);
    auto object=irr_io::read_object(root,h,k,fec,lec,pool);irr::channels_t cg(pool);
    irr_bulk::check_lec_convention(lec,cg,pool);std::ofstream discard("/dev/null");
    auto mf=irr_bulk::read_copies(fec,h,true),ml=irr_bulk::read_copies(lec,k,false);
    auto of=irr_bulk::read_copies(fec,h-1,true),ol=irr_bulk::read_copies(lec,k-1,false);
    auto first=irr_bulk::read_chain(fec,h,true,mf,discard,pool),last=irr_bulk::read_chain(lec,k,false,ml,discard,pool);
    auto W=irr::read_tensor(data/"dlogmat_full_ca_pd3.wxf",pool);
    auto t=irr::make_seam_template(cg,W,nmhv_sew_read_irrep_axis(data/"dlogmat_full_ca_pd3_blocks.tsv"),opt);
    auto layout=irr::seam_layout(t,of,ol,mf,ml);int_t df,dl,dr=1,dc=1,maximum=0;
    auto sf=factor_norm(first,df),sl=factor_norm(last,dl);int_t sr=0;
    for(const auto& b:t.routes)for(const auto& l:b)for(const auto& r:l)for(const auto& v:r)dr=LCM(dr,v.value.den());
    for(const auto& b:t.routes)for(const auto& l:b)for(const auto& r:l)for(const auto& v:r){int_t x=dr;x/=v.value.den();x*=v.value.num().abs();sr+=x;}
    for(const auto& b:object.S)for(const auto& row:b.rows)for(auto [j,x]:row)dc=LCM(dc,x.den());
    for(const auto& b:object.S)for(const auto& row:b.rows)for(auto [j,x]:row){int_t v=dc;v/=x.den();v*=x.num().abs();maximum=std::max(maximum,v);}
    int_t bound=sf*sl*sr*maximum,modulus=1;size_t primes=0;
    std::vector<fs::path> files;for(const auto& file:fs::directory_iterator(root/"cache"))if(file.path().extension()==".bin")files.push_back(file.path());std::sort(files.begin(),files.end());
    size_t native_primes=0,extra_primes=0;ulong maximum_prime=0;
    auto verify_prime=[&](ulong p,bool native){field_t F(FIELD_Fp,p);
        require(df%F.mod!=0&&dl%F.mod!=0&&dr%F.mod!=0&&dc%F.mod!=0,"Bad certificate prime");
        MM cached;if(native){cached=read_cache(root,p,h,k);require(cached.nrow==object.directions+1&&cached.ncol==layout.columns,"QQ/cache dimensions");}
        size_t affine_rows=object.directions+1;
        size_t offset=0;std::array<MM,6> C;
        for(size_t b=0;b<6;++b){C[b]=nmhv_sew_modular_reduce(object.S[b],F,pool);
            for(size_t r=0;r<affine_rows;++r){std::vector<ulong> expected(mf[b]*ml[b]),actual(expected.size());
                for(auto [j,x]:C[b][r])expected[j]=x;
                if(native)for(auto [j,x]:cached[r])if(size_t(j)>=offset&&size_t(j)<offset+expected.size())actual[j-offset]=x;
                if(native)require(actual==expected,"QQ/cache congruence mismatch");}offset+=mf[b]*ml[b];}
        std::vector<std::vector<ulong>> residual(affine_rows,std::vector<ulong>(layout.rows));
        for(size_t b=0;b<6;++b){
            auto A=nmhv_sew_modular_reduce(first[b],F,pool);A.transpose_replace();
            auto B=nmhv_sew_modular_reduce(last[b],F,pool);
            auto left=irr::decode_channels(cg,b,of),right=irr::decode_channels(cg,b,ol);
            for(size_t r=0;r<affine_rows;++r){
                MM coefficients(mf[b],ml[b]);for(auto [j,x]:C[b][r])coefficients[j/ml[b]].push_back(j%ml[b],x);
                auto middle=sparse_mat_mul(A,coefficients,F,pool);
                auto image=sparse_mat_mul(middle,B,F,pool);
                for(size_t i=0;i<image.nrow;++i)for(auto [j,x]:image[i]){
                    auto l=left[i],z=right[j];for(const auto& route:t.routes[b][l.channel][z.channel]){
                        size_t row=layout.rowbase[route.row]+l.copy*ol[z.block]+z.copy;
                        residual[r][row]=nmod_add(residual[r][row],nmod_mul(x,route.value%F.mod,F.mod),F.mod);
                    }
                }
            }
        }
        for(const auto& row:residual)for(auto x:row)require(!x,"Original unprereduced sewing equation failed");
        modulus*=p;++primes;std::cout<<"ORIGINAL_SEAM_PRIME prime="<<p<<" rows="<<layout.rows<<" affine_rows="<<affine_rows<<" native_congruence="<<native<<" passed=1"<<std::endl;
    };
    for(const auto& file:files){auto name=file.filename().string();require(name.starts_with("p"),"Cache filename");
        ulong p=std::stoull(name.substr(1,name.size()-5));verify_prime(p,true);++native_primes;maximum_prime=std::max(maximum_prime,p);}
    require(native_primes,"No native reconstruction primes");
    // Further residual-only primes strengthen the integer height certificate.
    // They are not native symbol solves and are never written as solve caches.
    while(bound>=modulus){maximum_prime=n_nextprime(maximum_prime,0);verify_prime(maximum_prime,false);++extra_primes;}
    require(primes&&bound<modulus,"Exact original residual certificate");
    std::cout<<"JOINT_CERTIFICATE {\"exact_QQ_original_seam\":true,\"directions\":"<<object.directions
        <<",\"original_rows\":"<<layout.rows<<",\"columns\":"<<layout.columns
        <<",\"residual_primes\":"<<primes<<",\"native_congruence_primes\":"<<native_primes<<",\"extra_residual_primes\":"<<extra_primes<<",\"bound_bits\":"<<bound.bits()
        <<",\"modulus_bits\":"<<modulus.bits()<<"}"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
