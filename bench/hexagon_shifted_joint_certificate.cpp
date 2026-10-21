#include "irr_physics_pipeline.hpp"
#include "irr_bulk.hpp"
#include "irr_shift_cache.hpp"
using irr::require;
using MM=sparse_mat<ulong,int32_t>;

// This helper accepts a direct shift-cache child of a native joint solve.
// Its exact QQ conclusion comes from the independent original target seam;
// the parent hash receipts establish provenance, not native target solves.
int_t factor_norm(const irr::matrices& blocks,int_t& den) {
    den=1;for(const auto& b:blocks)for(const auto& row:b.rows)for(auto [j,x]:row)den=LCM(den,x.den());
    int_t sum=0;for(const auto& b:blocks)for(const auto& row:b.rows)for(auto [j,x]:row){int_t v=den;v/=x.den();v*=x.num().abs();sum+=v;}return sum;
}

struct shift_provenance {
    fs::path parent;size_t h,k,rows;uint64_t fingerprint;
    std::string source_algorithm;irr_stage::fields reconstruction;
};

shift_provenance verify_shift_provenance(const fs::path& root,const fs::path& fec,const fs::path& lec,
    size_t h,size_t k,size_t rows) {
    auto stage=irr_stage::read(root/"stage.tsv");
    require(stage.at("schema")=="irr-shift-cache-v1"&&stage.at("operator")=="shift-cache","Certificate requires a shifted cache stage");
    require(irr_io::number(stage.at("h"))==h&&irr_io::number(stage.at("k"))==k,"Shift target split mismatch");
    require(stage.at("kind")=="unique"||stage.at("kind")=="affine","Certificate requires affine origin-last caches");
    require(irr_io::number(stage.at("rows"))==rows,"Shift target affine row count mismatch");
    shift_provenance out;out.parent=fs::absolute(stage.at("parent")).lexically_normal();out.rows=rows;
    irr_stage::no_links(out.parent);require(out.parent!=fs::absolute(root).lexically_normal(),"Shift cannot be its own parent");
    auto parent=irr_stage::read(out.parent/"stage.tsv");
    require(parent.at("operator")=="bulk"&&parent.at("schema")=="irr-physics-v2-candidate","Expected a direct native joint-solve parent");
    out.source_algorithm=parent.at("algorithm");
    require(out.source_algorithm.starts_with("joint-discovery-physical-affine-"),"Parent is not a joint physical solve");
    require(parent.at("object_kind")=="affine"&&parent.at("helicity")==stage.at("helicity"),"Source affine/helicity convention mismatch");
    require(parent.at("history")==stage.at("history"),"Shift physical history mismatch");
    out.h=irr_io::number(parent.at("h"));out.k=irr_io::number(parent.at("k"));
    require(out.h>h&&out.k>0&&out.h+out.k==h+k,"Invalid source-to-target left shift");
    auto source_basis=irr_stage::basis_identity(fec,lec,out.h,out.k),target_basis=irr_stage::basis_identity(fec,lec,h,k);
    irr_stage::verify_basis(out.parent,source_basis,out.h,out.k);
    std::string identity="irr-shift-cache-v1\n"+hex_sha256_file(out.parent/"input_identity.txt")+"\n"+
        stage.at("helicity")+" "+stage.at("kind")+" "+std::to_string(out.h)+" "+std::to_string(out.k)+" "+
        std::to_string(h)+" "+std::to_string(k)+" "+std::to_string(rows)+"\n";
    for(const auto& fields:{source_basis,target_basis})for(auto [key,value]:fields)identity+=key+"\t"+value+"\n";
    require(irr_stage::text(root/"input_identity.txt")==identity,"Shift map/source identity mismatch");
    out.fingerprint=irr_shift_cache::fingerprint(out.parent);
    out.reconstruction=irr_stage::read(root/"reconstruction_primes.tsv");
    auto summary=irr_stage::read(root/"summary.tsv");size_t crt=0,validation=0;
    for(const auto& [name,role]:out.reconstruction){require(name.starts_with("p"),"Invalid reconstruction prime name");
        size_t consumed=0;ulong p=std::stoull(name.substr(1),&consumed);
        require(consumed==name.size()-1&&p>2&&n_is_prime(p),"Invalid reconstruction prime");
        if(role=="CRT")++crt;else if(role=="validation")++validation;else require(false,"Unknown reconstruction prime role");}
    require(crt>0&&validation==2&&crt==irr_io::number(summary.at("crt_primes"))&&validation==irr_io::number(summary.at("validation_primes")),"Reconstruction role receipt mismatch");
    std::cout<<"SHIFT_PROVENANCE parent="<<out.parent<<" source_h="<<out.h<<" source_k="<<out.k
        <<" target_h="<<h<<" target_k="<<k<<" parent_algorithm="<<out.source_algorithm
        <<" reconstruction_receipt_sha256="<<hex_sha256_file(root/"reconstruction_primes.tsv")
        <<" native_target_solves=0 verified=1"<<std::endl;
    return out;
}

int main(int argc,char**argv)try {
    require(argc==7,"shifted_joint_certificate H K QQ_ROOT FEC LEC DATA");
    Flint::set_memory_functions();size_t h=std::stoul(argv[1]),k=std::stoul(argv[2]);
    fs::path root=argv[3],fec=argv[4],lec=argv[5],data=argv[6];
    rref_option_t opt;irr_run::threads(opt);opt->method=1;opt->is_back_sub=true;opt->verbose=false;auto* pool=&opt->pool;
    irr_stage::no_links(root);irr_stage::verify_complete(root);irr_stage::verify_basis(root,irr_stage::basis_identity(fec,lec,h,k),h,k);
    auto object=irr_io::read_object(root,h,k,fec,lec,pool);irr::channels_t cg(pool);
    require(object.kind=="affine"||object.kind=="unique","Expected an affine QQ object");
    auto provenance=verify_shift_provenance(root,fec,lec,h,k,object.directions+1);
    irr_bulk::check_lec_convention(lec,cg,pool);std::ofstream discard("/dev/null");
    auto mf=irr_bulk::read_copies(fec,h,true),ml=irr_bulk::read_copies(lec,k,false);
    auto source_mf=irr_bulk::read_copies(fec,provenance.h,true),source_ml=irr_bulk::read_copies(lec,provenance.k,false);
    size_t source_columns=irr_shift_cache::width(source_mf,source_ml);
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
    std::map<ulong,fs::path> files;for(const auto& file:fs::directory_iterator(root/"cache"))if(file.path().extension()==".bin"){
        auto name=file.path().filename().string();require(name.starts_with("p"),"Cache filename");size_t consumed=0;
        auto number=name.substr(1,name.size()-5);ulong p=std::stoull(number,&consumed);
        require(consumed==number.size()&&p>2&&n_is_prime(p),"Invalid cache prime filename");require(files.emplace(p,file.path()).second,"Duplicate derived prime");}
    for(const auto& [name,role]:provenance.reconstruction)require(files.count(std::stoull(name.substr(1))),"Reconstruction receipt references absent derived cache");
    size_t derived_primes=0,extra_primes=0;ulong maximum_prime=0;
    auto verify_prime=[&](ulong p,bool derived){field_t F(FIELD_Fp,p);
        require(df%F.mod!=0&&dl%F.mod!=0&&dr%F.mod!=0&&dc%F.mod!=0,"Bad certificate prime");
        MM cached;if(derived){
            auto source=provenance.parent/"cache"/("p"+std::to_string(p)+".bin");
            auto source_header=irr_shift_cache::header(source,provenance.h,provenance.k,provenance.fingerprint);
            require(source_header.prime==p&&source_header.nrow==provenance.rows&&source_header.ncol==source_columns,"Source cache dimensions/prime mismatch");
            require(irr_stage::text(root/("source_p"+std::to_string(p)+".sha256"))==hex_sha256_file(source)+"\n","Shifted cache source receipt mismatch");
            cached=read_cache(root,p,h,k);require(cached.nrow==object.directions+1&&cached.ncol==layout.columns,"QQ/derived cache dimensions");}
        size_t affine_rows=object.directions+1,offset=0;std::array<MM,6> C;
        for(size_t b=0;b<6;++b){C[b]=nmhv_sew_modular_reduce(object.S[b],F,pool);
            for(size_t r=0;r<affine_rows;++r){std::vector<ulong> expected(mf[b]*ml[b]),actual(expected.size());
                for(auto [j,x]:C[b][r])expected[j]=x;
                if(derived)for(auto [j,x]:cached[r])if(size_t(j)>=offset&&size_t(j)<offset+expected.size())actual[j-offset]=x;
                if(derived)require(actual==expected,"QQ/derived cache congruence mismatch");}offset+=mf[b]*ml[b];}
        std::vector<std::vector<ulong>> residual(affine_rows,std::vector<ulong>(layout.rows));
        for(size_t b=0;b<6;++b){
            auto A=nmhv_sew_modular_reduce(first[b],F,pool);A.transpose_replace();auto B=nmhv_sew_modular_reduce(last[b],F,pool);
            auto left=irr::decode_channels(cg,b,of),right=irr::decode_channels(cg,b,ol);
            for(size_t r=0;r<affine_rows;++r){
                MM coefficients(mf[b],ml[b]);for(auto [j,x]:C[b][r])coefficients[j/ml[b]].push_back(j%ml[b],x);
                auto middle=sparse_mat_mul(A,coefficients,F,pool);auto image=sparse_mat_mul(middle,B,F,pool);
                for(size_t i=0;i<image.nrow;++i)for(auto [j,x]:image[i]){
                    auto l=left[i],z=right[j];for(const auto& route:t.routes[b][l.channel][z.channel]){
                        size_t row=layout.rowbase[route.row]+l.copy*ol[z.block]+z.copy;
                        residual[r][row]=nmod_add(residual[r][row],nmod_mul(x,route.value%F.mod,F.mod),F.mod);}}}}
        for(const auto& row:residual)for(auto x:row)require(!x,"Original unprereduced target sewing equation failed");
        modulus*=p;++primes;std::cout<<"ORIGINAL_TARGET_SEAM_PRIME prime="<<p<<" rows="<<layout.rows<<" affine_rows="<<affine_rows
            <<" derived_cache_congruence="<<derived<<" native_target_solve=0 passed=1"<<std::endl;
    };
    for(const auto& [p,file]:files){verify_prime(p,true);++derived_primes;maximum_prime=std::max(maximum_prime,p);}
    require(derived_primes,"No derived reconstruction caches");
    while(bound>=modulus){maximum_prime=n_nextprime(maximum_prime,0);verify_prime(maximum_prime,false);++extra_primes;}
    require(primes&&bound<modulus,"Exact original target residual certificate");
    std::cout<<"SHIFTED_JOINT_CERTIFICATE {\"exact_QQ_original_target_seam\":true,\"shift_source_hash_receipts_verified\":true,\"directions\":"<<object.directions
        <<",\"original_rows\":"<<layout.rows<<",\"columns\":"<<layout.columns<<",\"source_h\":"<<provenance.h<<",\"source_k\":"<<provenance.k
        <<",\"target_h\":"<<h<<",\"target_k\":"<<k<<",\"residual_primes\":"<<primes<<",\"native_target_solve_primes\":0"
        <<",\"derived_cache_congruence_primes\":"<<derived_primes<<",\"extra_residual_primes\":"<<extra_primes<<",\"bound_bits\":"<<bound.bits()
        <<",\"modulus_bits\":"<<modulus.bits()<<"}"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}
