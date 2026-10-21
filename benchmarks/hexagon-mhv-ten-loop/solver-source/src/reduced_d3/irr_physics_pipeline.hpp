#pragma once
#include "irr_projection.hpp"
#include "irr_stage.hpp"
#include "irr_shift_cache.hpp"
#include "irr_profile.hpp"
#include <functional>
#include <memory>

namespace fs=std::filesystem;
using irr::require;
using irr_stage::fields;
using irr_physics::equations;

// Keep at most n independent equations while slot projections arrive in batches.
// This is only the small affine solve, never a presentation RREF of sewing.
struct affine_rows {
    field_t F;size_t n;std::vector<std::optional<std::vector<ulong>>> pivots;thread_pool* pool;
    affine_rows(size_t count,const field_t& field,thread_pool* p=nullptr):F(field),n(count),pivots(n),pool(p){}
    void add(const equations& rows,bool origin_column=true) {
        std::vector<std::vector<ulong>> pending;pending.reserve(rows.size());
        for(const auto& eq:rows){auto r=eq.coefficients;require(r.size()==n+size_t(origin_column),"Projection width mismatch");
            ulong rhs=eq.rhs;if(origin_column){rhs=nmod_sub(rhs,r.back(),F.mod);r.pop_back();}r.push_back(rhs);pending.push_back(std::move(r));}
        // Existing pivots are read-only across the whole incoming batch.
        // Admit new pivots afterwards in the original sample order.
        if(pool&&n>=128)irr_parallel::rows(pending.size(),pool,[&](size_t i){auto& r=pending[i];
            for(size_t j=0;j<n;++j)if(pivots[j]&&r[j]){auto x=r[j];for(size_t k=j;k<=n;++k)r[k]=nmod_sub(r[k],nmod_mul(x,(*pivots[j])[k],F.mod),F.mod);}});
        for(size_t i=0;i<rows.size();++i){auto& r=pending[i];const auto& eq=rows[i];
            for(size_t j=0;j<n;++j)if(pivots[j]&&r[j]){auto x=r[j];for(size_t k=j;k<=n;++k)r[k]=nmod_sub(r[k],nmod_mul(x,(*pivots[j])[k],F.mod),F.mod);}
            size_t j=0;while(j<n&&!r[j])++j;
            if(j==n){require(!r[n],"Inconsistent "+eq.condition+" at sample "+std::to_string(eq.sample));continue;}
            auto inv=n_invmod(r[j],F.mod.n);for(size_t k=j;k<=n;++k)r[k]=nmod_mul(r[k],inv,F.mod);
            pivots[j]=std::move(r);
        }
    }
    auto solve() {
        // Back-substitute only once. Each worker owns an earlier pivot row;
        // the normalized source row is immutable until all workers finish.
        for(size_t j=n;j-->0;)if(pivots[j]){
            auto reduce=[&](size_t i){auto& old=pivots[i];if(old&&(*old)[j]){auto x=(*old)[j];
                for(size_t k=j;k<=n;++k)(*old)[k]=nmod_sub((*old)[k],nmod_mul(x,(*pivots[j])[k],F.mod),F.mod);}};
            if(pool&&j>=128&&n-j>=128)irr_parallel::rows(j,pool,reduce);else for(size_t i=0;i<j;++i)reduce(i);
        }
        nmhv_physical_affine_solution_t<int32_t> result;result.origin.resize(n);std::vector<size_t> free;
        for(size_t j=0;j<n;++j)if(pivots[j]){result.pivot_columns.push_back(j);result.origin[j]=pivots[j]->back();}else free.push_back(j);
        result.rank=result.pivot_columns.size();result.directions=sparse_mat<ulong,int32_t>(free.size(),n);
        for(size_t r=0;r<free.size();++r)for(size_t j=0;j<n;++j){ulong x=j==free[r]?1:(pivots[j]?nmod_neg((*pivots[j])[free[r]],F.mod):0);
            if(x)result.directions[r].push_back(j,x);}
        return result;
    }
};

inline bool affine_satisfies(const equations& rows,const nmhv_physical_affine_solution_t<int32_t>& x,const field_t& F){
    for(const auto& row:rows){
        require(row.coefficients.size()==x.origin.size()+1,"Affine validation width mismatch");
        ulong value=row.coefficients.back();
        for(size_t j=0;j<x.origin.size();++j)value=nmod_add(value,nmod_mul(row.coefficients[j],x.origin[j],F.mod),F.mod);
        if(value!=row.rhs)return false;
        for(const auto& direction:x.directions.rows){value=0;for(auto [j,v]:direction)value=nmod_add(value,nmod_mul(row.coefficients[j],v,F.mod),F.mod);if(value)return false;}
    }
    return true;
}

// Project once, then admit nested endpoint truncations without repeating jets.
// Continue after full rank: unused higher powers must still be consistent.
inline void ope_rank_trace(affine_rows& reduced,equations rows,std::ostream& report,ulong prime) {
    for(const auto& row:rows)require(row.ope_coordinate[0]>=0,"Missing OPE power label");
    std::stable_sort(rows.begin(),rows.end(),[](const auto& a,const auto& b){return a.ope_coordinate[0]<b.ope_coordinate[0];});
    report<<"prime\tpower\tnew_equations\ttotal_equations\tnew_rank\trank\tremaining\tfirst_unique\n";
    size_t rank=std::count_if(reduced.pivots.begin(),reduced.pivots.end(),[](const auto& p){return bool(p);});
    bool unique=rank==reduced.n;
    report<<prime<<"\t-1\t0\t0\t0\t"<<rank<<'\t'<<reduced.n-rank<<'\t'<<unique<<'\n';
    for(size_t i=0;i<rows.size();){size_t end=i+1;int power=rows[i].ope_coordinate[0];
        while(end<rows.size()&&rows[end].ope_coordinate[0]==power)++end;
        equations batch;batch.reserve(end-i);for(size_t j=i;j<end;++j)batch.push_back(std::move(rows[j]));
        reduced.add(batch);
        size_t next=std::count_if(reduced.pivots.begin(),reduced.pivots.end(),[](const auto& p){return bool(p);});
        report<<prime<<'\t'<<power<<'\t'<<end-i<<'\t'<<end<<'\t'<<next-rank<<'\t'<<next<<'\t'<<reduced.n-next<<'\t'<<(!unique&&next==reduced.n)<<std::endl;
        std::cout<<"OPE_RANK prime="<<prime<<" power="<<power<<" rank="<<next<<" remaining="<<reduced.n-next<<" new_rank="<<next-rank<<std::endl;
        unique=next==reduced.n;rank=next;i=end;
    }
}

struct settings {
    fs::path fec,lec,data,output,input,plan,shifted_cache;
    std::string helicity,conditions;size_t h=0,k=0,samples=0,batch=32,checks=32,primes=0;
    uint64_t seed=0x193a4527;bool all=false;
    std::string sampling="adaptive",projection="scalar";
    bool shoup=true,resume_legacy=false,ope_rank_only=false;
    size_t qq_cache_mib=0;
};
struct operation {std::string label,name;size_t samples,checks;uint64_t seed;std::string coverage,projection;fs::path data;bool data_override=false;};
operation parse_operation(const settings& s,std::string label,std::string name,const fs::path& config={}) {
    require(name=="collinear"||name=="mrk"||name=="ope"||name=="ope_small_s","Unknown condition: "+name);
    require(name!="ope_small_s"||s.helicity=="MHV","ope_small_s supports MHV only");
    operation op{label,name,s.samples,s.checks,s.seed,name=="ope_small_s"?"MHV-small-S-T2-all":"",name=="collinear"?s.projection:"scalar",s.data};
    if(!config.empty()&&config!="-")for(auto [key,value]:irr_stage::read(config)){
        if(key=="samples")op.samples=irr_io::number(value);else if(key=="validation_samples")op.checks=irr_io::number(value);
        else if(key=="seed")op.seed=irr_io::number(value);
        else if(key=="projection")op.projection=value;
        else if(key=="data"){op.data=fs::path(value);if(op.data.is_relative())op.data=config.parent_path()/op.data;op.data=fs::absolute(op.data).lexically_normal();op.data_override=true;}
        else if(key=="coverage"){
            if(name=="ope_small_s"){require(value=="MHV-small-S-T2-all"||value=="MHV-small-S-T1T2-all","Unsupported small-S OPE coverage");op.coverage=value;}
            else require(value==(s.helicity=="MHV"?"T1T2-top3":"NMHV-T1T2-all"),"Unsupported OPE coverage for helicity");
        }else throw std::runtime_error("Unknown condition configuration: "+key);
    }
    require(op.checks&&(op.samples||s.sampling=="adaptive"),"Sampling and validation counts must be positive");
    require(op.projection=="scalar"||op.projection=="tail-exact","Choose scalar or tail-exact projection");
    require(op.projection!="tail-exact"||s.k<=3,"Tail-exact requires k<=3");
    require(op.projection!="tail-exact"||name=="collinear","Tail-exact is only implemented for collinear");
    return op;
}
bool legacy_recipe(const fields& r){return r.at("schema")=="irr-physics-v1";}
std::vector<operation> operations(const settings& s) {
    std::vector<operation> out;auto append=[&](std::string label,std::string name,fs::path config){
        auto op=parse_operation(s,label,name,config);
        for(auto& old:out)require(old.label!=label,"Duplicate condition label");out.push_back(op);
    };
    if(!s.plan.empty()){
        std::ifstream in(s.plan);require(bool(in),"Cannot read condition plan");std::string line;
        while(std::getline(in,line))if(!line.empty()&&line[0]!='#'){std::istringstream row(line);std::string label,name,config;row>>label>>name>>config;require(bool(row),"Plan requires label operator config");append(label,name,config=="-"?fs::path{}:s.plan.parent_path()/config);}
    }else {std::istringstream in(s.conditions);std::string name;size_t i=0;while(std::getline(in,name,'+'))append(std::to_string(++i),name,{});}
    require(!out.empty(),"CONDITIONS or CONDITION_PLAN is required");return out;
}
uint64_t fingerprint(const fs::path& root){return std::stoull(hex_sha256_file(root/"input_identity.txt").substr(0,16),nullptr,16);}
fs::path cache_file(const fs::path& root,ulong p){return root/"cache"/("p"+std::to_string(p)+".bin");}
sparse_mat<ulong,int32_t> read_cache(const fs::path& root,ulong prime,size_t h,size_t k) {
    irr_run::phase timer("cache_read");
    auto path=cache_file(root,prime);require(irr_stage::text(path.string()+".sha256")==hex_sha256_file(path)+"\n","Damaged modular cache: "+path.string());
    std::ifstream in(path,std::ios::binary);nmhv_sew_cache_header_t header{};in.read(reinterpret_cast<char*>(&header),sizeof(header));
    require(bool(in)&&nmhv_sew_cache_header_matches(header,(h+k)/2,h,k,header.nrow,header.ncol,prime,fingerprint(root)),"Cache identity mismatch");
    require(header.nrow&&header.nrow<=INT32_MAX&&header.ncol<=INT32_MAX,"Cache shape exceeds index range");
    require(fs::file_size(path)==sizeof(header)+header.nrow*header.ncol*sizeof(ulong),"Cache length mismatch");
    sparse_mat<ulong,int32_t> out(header.nrow,header.ncol);std::vector<ulong> row(header.ncol);
    for(size_t i=0;i<out.nrow;++i){in.read(reinterpret_cast<char*>(row.data()),row.size()*sizeof(ulong));for(size_t j=0;j<out.ncol;++j){require(row[j]<prime,"Invalid residue");if(row[j])out[i].push_back(j,row[j]);}}
    require(bool(in),"Truncated cache");return out;
}
void save_cache(const fs::path& root,const sparse_mat<ulong,int32_t>& mat,ulong prime,size_t h,size_t k) {
    irr_run::phase timer("cache_write");
    fs::create_directories(root/"cache");auto p=cache_file(root,prime);
    nmhv_sew_write_modular_basis(p,mat,prime,(h+k)/2,h,k,fingerprint(root));
    irr_stage::fixed_text(p.string()+".sha256",hex_sha256_file(p)+"\n");
}

// A source may be an exact object or a stage recipe with modular checkpoints.
// Recipes refer to their actual parent; no implicit collinear dependency exists.
struct pipeline {
    std::function<sparse_mat<ulong,int32_t>(const fs::path&,ulong)> parent_evaluator;
    std::function<sparse_mat<ulong,int32_t>(const fs::path&,ulong)> special_evaluator;
    settings s;irr_stage::lock fec_guard,lec_guard;irr_physics::context c;fields basis;std::map<std::string,std::unique_ptr<irr_stage::lock>> locks;
    std::map<std::string,std::string> input_hashes,source_digests;
    std::map<std::string,fields> bound;
    std::set<std::string> verified_qq;
    std::set<std::string> checked_legacy;
    bool shifted_verified=false;
    fs::path shifted_original,shifted_parent;
    size_t shifted_h=0,shifted_k=0;
    pipeline(settings value):s(std::move(value)),fec_guard(s.fec/".active.lock",true),lec_guard(s.lec/".active.lock",true),c(s.fec,s.lec,s.output/"diagnostics"/std::to_string(std::chrono::high_resolution_clock::now().time_since_epoch().count()),s.data,s.h,s.k,s.batch,s.helicity=="MHV",n_nextprime(1ULL<<60,0)),basis(irr_stage::basis_identity(s.fec,s.lec,s.h,s.k)){c.shoup=s.shoup;}
    void lock_stage(const fs::path& root){auto key=fs::absolute(root).lexically_normal().generic_string();if(!locks.count(key)){irr_stage::no_links(root);locks[key]=std::make_unique<irr_stage::lock>(root/".active.lock");}}
    void protect_source(const fs::path& root){auto key=fs::absolute(root).lexically_normal().generic_string();if(!locks.count(key)){irr_stage::no_links(root);locks[key]=std::make_unique<irr_stage::lock>(root/".active.lock",true);}}
    std::string input_hash(const fs::path& path){auto key=fs::weakly_canonical(path).generic_string();
        if(auto it=input_hashes.find(key);it!=input_hashes.end())return it->second;
        return input_hashes[key]=hex_sha256_file(path);}
    std::vector<fs::path> inputs(const std::string& op,size_t h=0,size_t k=0,fs::path data={}) {
        if(!h){h=s.h;k=s.k;}
        if(data.empty())data=s.data;
        std::vector<fs::path> paths{s.lec/"LEC_0_components.wxf",s.lec/"LEC_0_blocks.tsv",s.lec/"LEC_1_physical.wxf",s.lec/"common/CG.tsv"};
        for(bool first:{true,false})for(size_t w=1;w<=(first?h:k);++w){auto root=first?s.fec:s.lec;auto dir=root/joint_eet_weight_name(w);
            paths.push_back(dir/(first?"multiplicities.tsv":"blocks.tsv"));for(size_t b=0;b<6;++b)paths.push_back(irr::tensor_path(dir/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf")));
        }
        if(op=="mrk")for(size_t l=1;l<=(h+k)/2;++l)for(size_t k=0;k<l;++k)paths.push_back(data/(s.helicity=="MHV"?"mrk_mhv":"mrk")/("bfkl_f_L"+std::to_string(l)+"_k"+std::to_string(k)+".tsv"));
        if(op=="ope_small_s"){
            for(auto f:std::vector<std::string>{"projection.tsv","generator_manifest.json","kernels_F2.tsv","kernels_F3.tsv","kernels_F4.tsv","calE_L"+std::to_string((h+k)/2)+".tsv"})paths.push_back(data/"ope_mhv_small_s"/f);
        }
        if(op=="ope"){
            if(s.helicity=="MHV")for(auto f:std::vector<std::string>{"kernels_F2.tsv","kernels_F3.tsv","calE_ope_L"+std::to_string((h+k)/2)+".tsv"})paths.push_back(data/"ope"/f);
            else {paths.push_back(data/"ope_nmhv/projection.tsv");paths.push_back(data/"ope_nmhv"/("E1111_L"+std::to_string((h+k)/2)+".tsv"));
                for(int f=2;f<=6;++f)for(auto prefix:{"kernels_F","prefactors_F"})paths.push_back(data/"ope_nmhv"/(prefix+std::to_string(f)+".tsv"));}
        }
        if(op=="ope_small_s"||(op=="ope"&&s.helicity=="NMHV")){
            auto coverage=data/(op=="ope_small_s"?"ope_mhv_small_s":"ope_nmhv")/"coverage_identity.tsv";
            if(fs::exists(coverage))paths.push_back(coverage);
        }
        return paths;
    }
    void preflight(const operation& op) {
        for(const auto& path:inputs(op.name,0,0,op.data))require(fs::is_regular_file(path),"Missing stage input before computation: "+path.string());
        if(op.name!="ope_small_s"&&!(op.name=="ope"&&s.helicity=="NMHV"))return;
        auto root=op.data/(op.name=="ope_small_s"?"ope_mhv_small_s":"ope_nmhv");
        fields config;std::ifstream meta(root/"projection.tsv");std::string line;
        while(std::getline(meta,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();line.erase(std::remove(line.begin(),line.end(),'"'),line.end());
            auto row=nmhv_limit_split_tabs(line);require(row.size()==2,"Invalid OPE projection metadata");require(config.emplace(row[0],row[1]).second,"Duplicate projection setting");}
        require(config.at("max_t")=="2","Only complete low-twist projection is supported");
        size_t grade=irr_io::number(config.at("grade")),loop=(s.h+s.k)/2;
        auto coverage=root/"coverage_identity.tsv";
        // New libraries carry a compact, hash-bound projection receipt. Legacy
        // <=9-loop libraries retain their original cache identities.
        if(loop>=10||fs::exists(coverage)){
            auto receipt=irr_stage::read(coverage);
            require(receipt.at("status")=="complete"&&irr_io::number(receipt.at("max_loop"))>=loop,"Incomplete OPE coverage");
            require(receipt.at("generator_sha256")==input_hash(root/"generator_manifest.json"),"OPE generator manifest changed");
            for(const auto& [key,value]:receipt)if(key.starts_with("file:"))require(input_hash(root/key.substr(5))==value,"OPE coverage file changed: "+key);
        }
        auto path=root/((op.name=="ope_small_s"?"calE_L":"E1111_L")+std::to_string(loop)+".tsv");
        std::ifstream in(path);std::getline(in,line);size_t solve=0,check=0;std::set<std::array<int,6>> seen;
        while(std::getline(in,line))if(!line.empty()){
            auto row=nmhv_ope_row(line);require(row.size()==9,"Invalid OPE RHS row");std::array<int,6> key;
            for(size_t i=0;i<6;++i)key[i]=std::stoi(row[i]);
            require(key[0]==int(loop)&&key[1]>=1&&key[1]<=2&&std::abs(key[2])==key[1]&&key[3]>=key[1]&&size_t(key[1]+key[3])<=grade&&key[4]>=0&&key[5]>=0&&size_t(key[4]+key[5])<=2*loop,"OPE RHS exceeds kernel/log coverage");
            require(seen.insert(key).second,"Duplicate OPE coefficient");require(row[8]=="solve"||row[8]=="check","Invalid OPE role");
            require(int_t(row[7])!=0,"Zero OPE denominator");rat_t value(row[6]+"/"+row[7]);
            bool held=row[8]=="check"||(op.name=="ope_small_s"&&key[1]==1&&op.coverage!="MHV-small-S-T1T2-all");
            if(held)++check;else ++solve;
        }
        require(in.eof()&&solve&&check,"OPE requires nonempty solve and held-out rows");
        std::cout<<"PREFLIGHT operator="<<op.name<<" loop="<<loop<<" solve="<<solve<<" check="<<check<<" grade="<<grade<<std::endl;
    }
    std::string input_digest(const fs::path& root) {
        irr_run::phase timer("input_identity");
        if(!fs::exists(root/"stage.tsv"))protect_source(root);
        auto key=fs::weakly_canonical(root).generic_string();
        if(auto it=source_digests.find(key);it!=source_digests.end())return it->second;
        if(fs::exists(root/"stage.tsv"))return hex_sha256_file(root/"input_identity.txt");
        std::vector<fs::path> files;for(auto& e:fs::directory_iterator(root))if(e.is_regular_file()){
            auto name=e.path().filename().string();if(name=="columns.tsv"||name=="object.tsv"||name.ends_with(".wxf")||name.ends_with(".wxf.shape"))files.push_back(e.path());}
        std::sort(files.begin(),files.end());require(!files.empty(),"Empty sewing source");std::string value;for(auto p:files)value+=p.filename().string()+":"+input_hash(p)+";";return source_digests[key]=value;
    }
    void prepare_shifted(const fs::path& root) {
        if(s.shifted_cache.empty())return;
        irr_stage::no_links(s.shifted_cache);irr_stage::no_links(s.shifted_cache/"cache");
        auto shifted=irr_stage::read(s.shifted_cache/"stage.tsv"),r=irr_stage::read(root/"stage.tsv");
        require(shifted.at("schema")=="irr-shift-cache-v1"&&shifted.at("kind")=="unique"&&shifted.at("rows")=="1","Shifted OPE resume currently requires unique per-prime results");
        for(auto key:{"helicity","h","k","history"})require(shifted.at(key)==r.at(key),std::string("Shifted stage mismatch: ")+key);
        irr_stage::verify_basis(s.shifted_cache,basis,s.h,s.k);
        shifted_original=shifted.at("parent");if(shifted_original.is_relative())shifted_original=s.shifted_cache/shifted_original;
        auto old=irr_stage::read(shifted_original/"stage.tsv");
        require(old.at("operator")=="ope"&&r.at("operator")=="ope","Shifted resume supports OPE only");
        for(auto key:{"schema","history","operator","helicity","seed","samples","sampling","projection","validation_samples"})require(old.at(key)==r.at(key),std::string("Source OPE recipe differs: ")+key);
        shifted_h=irr_io::number(old.at("h"));shifted_k=irr_io::number(old.at("k"));
        require(shifted_h>s.h&&shifted_h+shifted_k==s.h+s.k,"Not a left shift of this weight");
        irr_stage::verify_basis(shifted_original,irr_stage::basis_identity(s.fec,s.lec,shifted_h,shifted_k),shifted_h,shifted_k);
        // Rebuild the old producer identity using today's recurrence and RHS
        // files, so a different OPE table cannot silently reuse old residues.
        std::string label="irr-physics-stream-v1";for(auto [key,value]:old)if(key!="parent"&&key!="label")label+=" "+key+"="+value;
        std::string identity="irr-bulk-v1 "+label+"\n";for(auto file:inputs("ope",shifted_h,shifted_k))identity+=hex_sha256_file(file)+"\n";
        require(irr_stage::text(shifted_original/"input_identity.txt")==identity,"Original OPE input identity/RHS changed");
        require(irr_stage::text(s.shifted_cache/"input_identity.txt").starts_with("irr-shift-cache-v1\n"+hex_sha256_file(shifted_original/"input_identity.txt")+"\n"),"Shifted cache has another source identity");
        shifted_parent=old.at("parent");if(shifted_parent.is_relative())shifted_parent=shifted_original/shifted_parent;
        require(input_digest(shifted_parent)==old.at("parent_identity"),"Original pre-OPE parent changed");
        if(fs::exists(shifted_parent/"qq_complete.tsv"))irr_stage::verify_complete(shifted_parent);
        irr_stage::fixed_text(root/"shifted_resume.tsv","shifted_identity\t"+hex_sha256_file(s.shifted_cache/"input_identity.txt")+"\noriginal_identity\t"+hex_sha256_file(shifted_original/"input_identity.txt")+"\n");
    }
    std::optional<sparse_mat<ulong,int32_t>> shifted_prime(ulong p) {
        if(s.shifted_cache.empty())return std::nullopt;auto path=cache_file(s.shifted_cache,p);if(!fs::exists(path))return std::nullopt;
        auto h=irr_shift_cache::header(path,s.h,s.k,fingerprint(s.shifted_cache));
        require(h.nrow==1&&h.ncol==irr_shift_cache::width(c.sew.first,c.sew.last),"Shifted cache shape/unique convention mismatch");
        auto source=cache_file(shifted_original,p);auto original=irr_shift_cache::header(source,shifted_h,shifted_k,fingerprint(shifted_original));require(original.nrow==1,"Original OPE cache is not unique");
        require(irr_stage::text(s.shifted_cache/("source_p"+std::to_string(p)+".sha256"))==hex_sha256_file(source)+"\n","Original prime differs from shifted provenance");
        return irr_shift_cache::read(path,h);
    }
    void verify_shifted_parent(const sparse_mat<ulong,int32_t>& current,ulong p) {
        auto obj=irr_io::read_object(shifted_parent,shifted_h,shifted_k,s.fec,s.lec,&c.opt->pool);field_t F(FIELD_Fp,p);
        // A homogeneous parent has an implicit zero origin in the physics API.
        sparse_mat<ulong,int32_t> matrix(obj.directions+1,irr_shift_cache::width(obj.first,obj.last));size_t offset=0;
        for(const auto& block:obj.S){auto reduced=irr_shift_cache::reduce(block,F);for(size_t i=0;i<reduced.nrow;++i)for(auto [j,x]:reduced[i])matrix[i].push_back(offset+j,x);offset+=block.ncol;}
        irr::seam_recoupling_t recouple(c.cg,c.opt);
        for(size_t h=shifted_h;h>s.h;--h)matrix=irr_shift_cache::shift_one(matrix,h,shifted_h+shifted_k-h,s.fec,s.lec,c.cg,recouple,F,c.opt);
        require(nmhv_sew_matrices_equal(matrix,current),"Target pre-OPE parent is not the same coordinate-preserving shift");
        std::cout<<"SHIFTED PARENT VERIFIED prime="<<p<<std::endl;
    }
    fields recipe(const fs::path& parent,const operation& op,const std::string& history) {
        fields result={{"schema","irr-physics-v2-candidate"},{"sampling",s.sampling},{"projection",op.projection},{"history",history},{"operator",op.name},{"label",op.label},{"helicity",s.helicity},{"h",std::to_string(s.h)},{"k",std::to_string(s.k)},
            {"samples",std::to_string(op.samples)},{"validation_samples",std::to_string(op.checks)},{"seed",std::to_string(op.seed)},{"parent",fs::relative(parent,s.output/history).generic_string()},{"parent_identity",input_digest(parent)}};
        if(op.data_override)result["data"]=fs::relative(op.data,s.output/history).generic_string();
        if(op.name=="ope_small_s")result["coverage"]=op.coverage;return result;
    }
    fields resume_recipe(const fs::path& root,fields wanted){
        if(!s.resume_legacy||!fs::exists(root/"stage.tsv"))return wanted;
        auto old=irr_stage::read(root/"stage.tsv");if(!legacy_recipe(old))return wanted;
        for(const auto& key:{"history","operator","helicity","h","k","parent_identity","seed","validation_samples"})
            require(old.at(key)==wanted.at(key),std::string("Legacy resume mismatch: ")+key);
        auto parent=fs::path(old.at("parent"));if(parent.is_relative())parent=root/parent;
        require(fs::weakly_canonical(parent)==fs::weakly_canonical(root/wanted.at("parent")),"Legacy parent path mismatch");
        // Keep the original identity and cache headers. Record the new execution
        // policy separately; it does not change the sewing coordinate system.
        bind(root,old);
        irr_stage::fixed_text(root/"resume_v2.tsv","sampling\t"+s.sampling+"\nprojection\t"+s.projection+"\nsamples\t"+std::to_string(s.samples)+"\nlegacy_identity_sha256\t"+hex_sha256_file(root/"input_identity.txt")+"\n");
        return old;
    }
    void bind(const fs::path& root,const fields& r) {
        auto key=fs::weakly_canonical(root).generic_string();
        if(auto it=bound.find(key);it!=bound.end()){require(it->second==r,"Stage configuration changed in this process");return;}
        irr_run::phase timer("bind_identity");
        auto files=inputs(r.at("operator"),0,0,r.count("data")?root/r.at("data"):s.data);std::string label="irr-physics-stream-v1";
        for(auto [key,value]:r)if(key!="parent"&&key!="label")label+=" "+key+"="+value;
        fs::create_directories(root);lock_stage(root);
        std::string identity="irr-bulk-v1 "+label+"\n";for(auto file:files)identity+=input_hash(file)+"\n";
        if(!fs::exists(root/"input_identity.txt"))for(const auto& entry:fs::directory_iterator(root))require(entry.path().filename()==".active.lock","Unmanaged output directory; retained");
        irr_stage::fixed_text(root/"input_identity.txt",identity);
        if(fs::exists(root/"basis_identity.tsv"))require(irr_stage::read(root/"basis_identity.tsv")==basis,"Sewing recurrence basis changed");else irr_stage::write(root/"basis_identity.tsv",basis);
        if(fs::exists(root/"stage.tsv")){
            auto old=irr_stage::read(root/"stage.tsv");for(auto [key,value]:r)if(key!="parent"&&key!="label")require(old.at(key)==value,"Stage configuration mismatch");
        }else irr_stage::write(root/"stage.tsv",r);
        bound[key]=r;
        if(!fs::exists(root/"manifest.json")){std::ofstream out(root/"manifest.json");out<<"{\"schema\":\"hexagon-irr-physics-v1\",\"complete\":false,\"status\":\"cache-only\",\"history\":"<<irr_io::quote(r.at("history"))<<"}\n";}
    }
    void unpack(const sparse_mat<ulong,int32_t>& packed) {
        size_t off=0;c.rows=packed.nrow;
        for(size_t b=0;b<6;++b){size_t cols=c.sew.first[b]*c.sew.last[b];c.S[b]=sparse_mat<ulong,int32_t>(c.rows,cols);
            for(size_t i=0;i<c.rows;++i)for(auto [j,x]:packed[i])if(j>=off&&size_t(j)<off+cols)c.S[b][i].push_back(j-off,x);off+=cols;}
        require(off==packed.ncol,"Packed sewing column mismatch");
    }
    sparse_mat<ulong,int32_t> exact_source(const fs::path& root,ulong p) {
        irr_run::phase timer("exact_source");protect_source(root);
        auto key=fs::weakly_canonical(root).generic_string();
        if(fs::exists(root/"qq_complete.tsv")&&!verified_qq.count(key)){irr_stage::verify_complete(root);verified_qq.insert(key);}
        irr_stage::verify_basis(root,basis,s.h,s.k);
        auto obj=irr_io::read_object(root,s.h,s.k,s.fec,s.lec,&c.opt->pool);field_t F(FIELD_Fp,p);size_t cols=0;for(auto& m:obj.S)cols+=m.ncol;
        sparse_mat<ulong,int32_t> out(obj.directions+1,cols);size_t off=0;
        for(auto& m:obj.S){auto mm=nmhv_sew_modular_reduce(m,F,&c.opt->pool);for(size_t i=0;i<mm.nrow;++i)for(auto [j,x]:mm[i])out[i].push_back(off+j,x);off+=m.ncol;}return out;
    }
    void batches(const fields& r,bool validation,const std::function<void(const equations&)>& consume) {
        c.seed=std::stoull(r.at("seed"));if(validation)c.seed^=0x6a09e667f3bcc909ULL;
        if(r.at("operator")=="ope_small_s"){consume(irr_physics::ope_small_s(c,validation,r.at("coverage")=="MHV-small-S-T1T2-all"));c.cache_summary();return;}
        if(r.at("operator")=="ope"){consume(irr_physics::ope(c,validation));c.cache_summary();return;}
        size_t count=irr_io::number(r.at(validation?"validation_samples":"samples"));
        for(size_t offset=0;offset<count;offset+=s.batch){c.samples=std::min(s.batch,count-offset);c.sample_offset=offset;
            std::cout<<"PROJECTION "<<r.at("operator")<<" prime="<<c.F.mod.n<<" validation="<<validation<<" samples="<<offset<<".."<<offset+c.samples<<"/"<<count<<std::endl;
            consume(!validation&&r.at("operator")=="collinear"&&r.count("projection")&&r.at("projection")=="tail-exact"?
                irr_physics::collinear_tail_exact(c):irr_physics::slots(c,r.at("operator")=="mrk"));c.cache_summary();}
    }
    auto train(const fields& r){
        affine_rows reduced(c.rows-1,c.F,&c.opt->pool);
        if(r.at("operator")=="ope"||r.at("operator")=="ope_small_s"){
            auto directory=c.out/"ope_rank"/r.at("history");fs::create_directories(directory);
            std::ofstream report(directory/("prime_"+std::to_string(c.F.mod.n)+".tsv"));require(bool(report),"Cannot write OPE rank report");
            std::cout<<"OPE_RANK_REPORT "<<directory<<" operator="<<r.at("operator")<<" initial_directions="<<reduced.n<<std::endl;
            batches(r,false,[&](const equations& eq){auto start=c.now();ope_rank_trace(reduced,eq,report,c.F.mod.n);c.record("affine_forward",eq.size(),c.rows,0,0,start);});
            require(bool(report),"Failed writing OPE rank report");auto start=c.now();auto child=reduced.solve();
            c.record("affine_backward",child.rank,c.rows,child.directions.nnz(),0,start);return child;
        }
        if(!r.count("sampling")||r.at("sampling")=="fixed"||r.at("operator")=="ope"||r.at("operator")=="ope_small_s"){
            batches(r,false,[&](const equations& eq){auto start=c.now();reduced.add(eq);c.record("affine_forward",eq.size(),c.rows,0,0,start);});
            auto start=c.now();auto child=reduced.solve();c.record("affine_backward",child.rank,c.rows,child.directions.nnz(),0,start);return child;
        }
        size_t maximum=irr_io::number(r.at("samples"));if(!maximum)maximum=std::max<size_t>(384,2*reduced.n+32);
        size_t offset=0,last_rank=0,stagnant=0,attempt=0,total_rows=0;
        for(;;){
            c.seed=std::stoull(r.at("seed"));c.samples=std::min(s.batch,maximum-offset);c.sample_offset=offset;
            auto eq=r.at("operator")=="collinear"&&r.count("projection")&&r.at("projection")=="tail-exact"?
                irr_physics::collinear_tail_exact(c):irr_physics::slots(c,r.at("operator")=="mrk");
            auto start=c.now();reduced.add(eq);c.record("affine_forward",eq.size(),c.rows,0,0,start);
            size_t rank=std::count_if(reduced.pivots.begin(),reduced.pivots.end(),[](const auto& p){return bool(p);});
            stagnant=rank==last_rank?stagnant+eq.size():0;last_rank=rank;offset+=c.samples;total_rows+=eq.size();
            std::cout<<"SAMPLING "<<r.at("operator")<<" prime="<<c.F.mod.n<<" points="<<offset<<" equations="<<total_rows
                <<" rank="<<rank<<" unknowns="<<reduced.n<<" limit="<<maximum<<std::endl;c.cache_summary();
            if(rank==reduced.n||stagnant>=32||offset==maximum){
                start=c.now();auto child=reduced.solve();c.record("affine_backward",child.rank,c.rows,child.directions.nnz(),0,start);
                // A fresh validation stream for every attempted stop. If it
                // exposes a missing constraint, add it and continue training.
                c.seed=std::stoull(r.at("seed"))^0xbb67ae8584caa73bULL^(++attempt*0x9e3779b97f4a7c15ULL);
                bool valid=true;size_t checks=irr_io::number(r.at("validation_samples"));
                for(size_t i=0;i<checks;i+=s.batch){c.sample_offset=i;c.samples=std::min(s.batch,checks-i);
                    auto test=irr_physics::slots(c,r.at("operator")=="mrk");
                    if(!affine_satisfies(test,child,c.F)){valid=false;reduced.add(test);}}
                std::cout<<"SAMPLING_CHECK "<<r.at("operator")<<" attempt="<<attempt<<" passed="<<valid<<std::endl;
                if(valid)return child;
                require(offset<maximum,"Adaptive sample limit exhausted; checkpoints retained, increase --samples");stagnant=0;
            }
        }
    }
    void observe_ope(const fs::path& root,ulong p) {
        auto r=irr_stage::read(root/"stage.tsv");
        require(r.at("operator")=="ope"||r.at("operator")=="ope_small_s","Rank-only requires OPE");
        fs::path parent=r.at("parent");if(parent.is_relative())parent=root/parent;
        require(fs::exists(parent/"qq_complete.tsv")||fs::exists(parent/"object.tsv"),"Rank-only needs an existing rational pre-OPE object");
        auto current=exact_source(parent,p);c.F=field_t(FIELD_Fp,p);unpack(current);current.clear();
        c.data=r.count("data")?fs::weakly_canonical(root/r.at("data")):fs::weakly_canonical(s.data);
        auto child=train(r);
        batches(r,true,[&](const equations& eq){require(affine_satisfies(eq,child,c.F),"OPE rank experiment failed held-out check");});
        std::cout<<"OPE_RANK_CHECK prime="<<p<<" remaining="<<child.directions.nrow<<" passed=1"<<std::endl;
    }
    sparse_mat<ulong,int32_t> evaluate(const fs::path& root,ulong p,size_t depth=0) {
        irr_run::phase timer("evaluate_"+root.filename().string()+"_"+std::to_string(p));
        require(depth<64,"Cyclic/overlong stage ancestry");
        if(!fs::exists(root/"stage.tsv"))return exact_source(root,p);
        auto special=irr_stage::read(root/"stage.tsv");
        if(special_evaluator&&(special.at("operator")=="bulk"||special.at("operator")=="shift"))return special_evaluator(root,p);
        auto r=irr_stage::read(root/"stage.tsv");require(r.at("helicity")==s.helicity&&irr_io::number(r.at("h"))==s.h&&irr_io::number(r.at("k"))==s.k,"Incompatible cached stage");
        fs::path parent=r.at("parent");if(parent.is_relative())parent=root/parent;
        // Completed QQ is self-contained even after upstream caches are removed.
        if(fs::exists(root/"qq_complete.tsv"))return exact_source(root,p);
        require(input_digest(parent)==r.at("parent_identity"),"Parent sewing changed");bind(root,r);
        bool cached=fs::exists(cache_file(root,p));
        auto imported=s.shifted_cache.empty()?std::optional<sparse_mat<ulong,int32_t>>{}:shifted_prime(p);
        if(imported&&shifted_verified&&!cached){save_cache(root,*imported,p,s.h,s.k);std::cout<<"SHIFTED CACHE REUSED prime="<<p<<" nnz="<<imported->nnz()<<std::endl;return std::move(*imported);}
        auto key=fs::weakly_canonical(root).generic_string();
        bool check=s.resume_legacy&&legacy_recipe(r)&&!checked_legacy.count(key);
        if(cached&&!check&&(!imported||shifted_verified))return read_cache(root,p,s.h,s.k);
        auto execution=r;
        if(s.resume_legacy&&legacy_recipe(r)){
            execution["sampling"]=s.sampling;execution["projection"]=s.projection;execution["samples"]=std::to_string(s.samples);
            if(s.sampling=="fixed"&&!s.samples)execution["samples"]=r.at("samples");
        }
        auto current=parent_evaluator?parent_evaluator(parent,p):evaluate(parent,p,depth+1);c.F=field_t(FIELD_Fp,p);unpack(current);
        if(imported&&!shifted_verified)verify_shifted_parent(current,p);
        auto data=r.count("data")?fs::weakly_canonical(root/r.at("data")):fs::weakly_canonical(s.data);
        if(c.data!=data){c.data=data;c.mrk_symbols.reset();}
        auto child=train(execution);
        std::ostringstream profile;profile<<current.nrow<<' '<<current.ncol<<' '<<child.directions.nrow<<'\n';for(auto j:child.pivot_columns)profile<<j<<' ';profile<<'\n';
        irr_stage::fixed_text(root/"pivot_profile.txt",profile.str());
        sparse_mat<ulong,int32_t> transform(child.directions.nrow+1,current.nrow);
        for(size_t i=0;i<child.directions.nrow;++i)transform[i]=child.directions[i];
        for(size_t j=0;j<child.origin.size();++j)if(child.origin[j])transform[child.directions.nrow].push_back(j,child.origin[j]);
        transform[child.directions.nrow].push_back(current.nrow-1,1);
        auto result=sparse_mat_mul(transform,current,c.F,&c.opt->pool);current.clear();
        nmhv_sew_modular_canonicalize(result,&c.opt->pool);unpack(result);
        // Fresh slot samples check both the origin and every remaining direction.
        batches(r,true,[&](const equations& eq){for(const auto& row:eq){require(row.coefficients.back()==row.rhs,"Affine representative failed independent condition check");for(size_t j=0;j+1<row.coefficients.size();++j)require(!row.coefficients[j],"Direction failed independent condition check; increase samples");}});
        if(imported){require(nmhv_sew_matrices_equal(result,*imported),"Shifted cache differs from native target OPE result; not mixed");shifted_verified=true;
            std::cout<<"SHIFTED OPE VERIFIED prime="<<p<<" rows="<<result.nrow<<" nnz="<<result.nnz()<<std::endl;}
        if(cached){auto previous=read_cache(root,p,s.h,s.k);
            require(nmhv_sew_matrices_equal(result,previous),"Legacy cache differs from newly computed sewing; retained, not mixed");
            checked_legacy.insert(key);
            std::cout<<"LEGACY VERIFIED "<<root<<" prime="<<p<<" rows="<<result.nrow<<" nnz="<<result.nnz()<<std::endl;
            return previous;
        }
        if(check)checked_legacy.insert(key);
        save_cache(root,result,p,s.h,s.k);
        std::ofstream status(root/"manifest.json");status<<"{\"schema\":\"hexagon-irr-physics-v1\",\"complete\":false,\"status\":\"cache-only\",\"history\":"<<irr_io::quote(r.at("history"))<<",\"h\":"<<s.h<<",\"k\":"<<s.k<<",\"directions\":"<<result.nrow-1<<",\"columns\":"<<result.ncol<<",\"origin_row_in_cache\":\"last\"}\n";
        std::cout<<"STAGE "<<r.at("history")<<" prime="<<p<<" directions="<<result.nrow-1<<" nnz="<<result.nnz()<<std::endl;return result;
    }
    void export_qq(const fs::path& root,sparse_mat<rat_t,int32_t> matrix,size_t primes,size_t bits) {
        irr_run::phase timer("export_qq");
        irr::sewing_object_t obj;obj.h=s.h;obj.k=s.k;obj.first=c.sew.first;obj.last=c.sew.last;obj.directions=matrix.nrow-1;obj.kind=obj.directions?"affine":"unique";
        auto descriptor=irr_stage::read(root/"stage.tsv");
        if(descriptor.count("object_kind")&&descriptor.at("object_kind")=="homogeneous"){
            require(matrix[matrix.nrow-1].nnz()==0,"Homogeneous stage has nonzero origin");
            obj.kind="homogeneous";matrix.rows.pop_back();--matrix.nrow;
        }
        auto stats=joint_eet_rational_stats(matrix);size_t nnz=matrix.nnz();
        irr_io::write_reduced_packed(root,obj,std::move(matrix),&c.opt->pool);
        auto memory=hex_platform::current_process_memory();
        std::ofstream summary(root/"summary.tsv");summary<<"directions\tnnz\theight_bits\tcrt_primes\tvalidation_primes\tmodulus_bits\tpeak_rss_bytes\tpeak_private_bytes\n"<<obj.directions<<'\t'<<nnz<<'\t'<<stats.height_bits<<'\t'<<primes<<"\t2\t"<<bits<<'\t'<<memory.peak_working_set_bytes<<'\t'<<memory.peak_private_commit_bytes<<'\n';summary.close();
        auto r=irr_stage::read(root/"stage.tsv");std::ofstream manifest(root/"manifest.json");manifest<<"{\"schema\":\"hexagon-irr-physics-v1\",\"complete\":true,\"status\":\"QQ-verified\",\"history\":"<<irr_io::quote(r.at("history"))<<",\"helicity\":"<<irr_io::quote(s.helicity)<<",\"h\":"<<s.h<<",\"k\":"<<s.k<<",\"kind\":"<<irr_io::quote(obj.kind)<<",\"directions\":"<<obj.directions<<",\"crt_primes\":"<<primes<<",\"validation_primes\":2,\"modulus_bits\":"<<bits<<",\"final_presentation_rref\":false,\"origin_row_in_cache\":\"last\"}\n";manifest.close();
        irr_stage::complete(root);std::cout<<"QQ COMPLETE "<<root<<" directions="<<obj.directions<<" height_bits="<<stats.height_bits<<" modulus_bits="<<bits<<std::endl;
    }
    void reconstruct(const fs::path& root) {
        if(fs::exists(root/"qq_complete.tsv")){irr_stage::verify_complete(root);std::cout<<"QQ already complete "<<root<<std::endl;return;}
        std::vector<fs::path> paths;std::vector<ulong> primes;ulong p=n_nextprime(1ULL<<60,0);size_t nr=0,nc=0;int_t modulus=1;
        for(size_t n=1;;++n){auto mat=evaluate(root,p);if(n==1){nr=mat.nrow;nc=mat.ncol;}require(mat.nrow==nr&&mat.ncol==nc,"Unlucky prime: affine shape changed");mat.clear();
            paths.push_back(cache_file(root,p));primes.push_back(p);modulus*=p;p=n_nextprime(p,0);
            if(s.primes&&n<s.primes)continue;
            std::vector<fs::path> validation;std::vector<ulong> vp;ulong q=p;
            for(size_t i=0;i<2;++i){auto test=evaluate(root,q);require(test.nrow==nr&&test.ncol==nc,"Unlucky validation prime");validation.push_back(cache_file(root,q));vp.push_back(q);q=n_nextprime(q,0);}
            sparse_mat<rat_t,int32_t> result;
            bool ok;{irr_run::phase timer("CRT_and_validation");ok=nmhv_sew_reconstruct_basis(result,paths,primes,validation,vp,(s.h+s.k)/2,s.h,s.k,nr,nc,fingerprint(root),&c.opt->pool);}
            std::cout<<"RECONSTRUCT primes="<<n<<" modulus_bits="<<modulus.bits()<<" success="<<ok<<std::endl;
            if(ok){export_qq(root,std::move(result),n,modulus.bits());return;}
            require(!s.primes,"Manual prime budget insufficient; all checkpoints retained");
        }
    }
};
