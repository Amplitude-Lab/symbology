#pragma once
#include "../reduced_d3/irr_projection.hpp"

// Direct low-grade MHV OPE rows in the original invariant copy coordinates.
// This is a prototype: validate its adjoint and small-control projections
// before using it in a solve. Existing grade-12 physical holdouts are retained.
namespace irr_physics {

inline void require_ope_grade_filtration(const nmhv_ope_germ_t& germ){
    for(const auto& terms:germ.kernels)for(const auto& term:terms)
        require(term.t==-1?term.z>=0:term.t+term.z>=0,
            "OPE kernel can feed discarded higher grade into a retained coefficient");
}

// Exact transpose of nmhv_ope_germ_t::integrate. Reverse its three phases:
// origin log integration, triangular log solve, then kernel convolution.
// No stride-by-stride integration matrix is materialized.
inline void integrate_ope_adjoint(const nmhv_ope_germ_t& germ,const ulong* output,
    ulong* input,size_t letter,size_t w){
    require(letter<9&&w>0&&w<=germ.weight,"OPE adjoint integration index");
    std::fill(input,input+germ.stride,0);std::vector<ulong> adjoint(output,output+germ.stride);
    for(size_t a=0;a<w;++a)for(size_t b=0;b<w-a;++b){
        ulong value=germ.add(germ.mul(adjoint[germ.index(0,0,a+1,b)],germ.at[letter]),
            germ.mul(adjoint[germ.index(0,0,a,b+1)],germ.az[letter]));
        input[germ.index(0,0,a,b)]=germ.add(input[germ.index(0,0,a,b)],germ.mul(value,germ.inverse[w]));
    }
    // Forward integration traverses descending a,b. Its transpose traverses
    // ascending a,b and scatters into the larger-log cotangent coordinates.
    for(int t=0;t<=germ.max_t;++t)for(int z=-t;z<=germ.grade-t;++z){
        if(t==0&&z==0)continue;
        for(size_t a=0;a<=w;++a)for(size_t b=0;b<=w-a;++b){
            size_t position=germ.index(t,z,a,b);ulong value=germ.mul(adjoint[position],germ.inverse[t?t:z]);
            adjoint[position]=value;
            if(t&&a+b<w){size_t next=germ.index(t,z,a+1,b);
                adjoint[next]=nmod_sub(adjoint[next],germ.mul(a+1,value),germ.F.mod);}
            if(!t&&a+b<w){size_t next=germ.index(t,z,a,b+1);
                adjoint[next]=nmod_sub(adjoint[next],germ.mul(b+1,value),germ.F.mod);}
        }
    }
    for(const auto& term:germ.kernels[letter])for(int t=0;t<=germ.max_t;++t)for(int z=-t;z<=germ.grade-t;++z){
        if((t==0&&z==0)||(t==0)!=(term.t==-1))continue;
        int ot=t-(term.t==-1?0:term.t),oz=z-term.z;if(!germ.valid(ot,oz))continue;
        for(size_t a=0;a<w;++a)for(size_t b=0;b<w-a;++b){
            ulong value=adjoint[germ.index(t,z,a,b)];if(value){size_t position=germ.index(ot,oz,a,b);
                input[position]=germ.add(input[position],germ.mul(value,term.value));}
        }
    }
}

inline std::shared_ptr<const std::vector<ulong>> direct_ope_head(context& c,const nmhv_ope_germ_t& germ,
    const std::string& head_identity){
    if(c.ope_head_cache_prime!=c.F.mod.n){c.clear_ope_head_cache();c.ope_head_cache_prime=c.F.mod.n;}
    auto key=head_identity+":"+std::to_string(c.h)+":"+std::to_string(germ.stride);
    if(auto it=c.ope_head_cache.find(key);it!=c.ope_head_cache.end())return it->second;
    std::vector<ulong> head(germ.stride);head[germ.index(0,0,0,0)]=1;
    for(size_t w=1;w<=c.h;++w){auto K=c.modular_recurrence(w,true);
        axis old(irr::chain_copies(c.fec,w-1,true)),next(irr::chain_copies(c.fec,w,true));
        head=step(c,irr::matrices{},old,next,head,germ.stride,[&](const ulong* a,ulong* b,size_t letter){
            germ.integrate(a,b,letter,w);},"ope_head_"+std::to_string(w),K.get());
    }
    auto result=std::make_shared<const std::vector<ulong>>(std::move(head));size_t bytes=result->size()*sizeof(ulong);
    if(!head_identity.empty()&&bytes<=c.ope_head_cache_limit-c.ope_head_cache_bytes){
        c.ope_head_cache.emplace(key,result);c.ope_head_cache_bytes+=bytes;
        std::cout<<"DIRECT_OPE_HEAD_CACHE stored=1 prime="<<c.F.mod.n<<" bytes="<<bytes<<" total="<<c.ope_head_cache_bytes<<std::endl;
    }
    return result;
}

struct direct_ope_target {int t,z,a,b;ulong rhs;};
inline std::vector<direct_ope_target> direct_ope_lowgrade_targets(context& c,bool include_twist1){
    require(c.mhv,"Direct low-grade OPE prototype is MHV only");auto root=c.data/"ope_mhv_small_s";
    c.input(root/"projection.tsv");std::map<std::string,std::string> config;
    std::ifstream settings(root/"projection.tsv");std::string line;
    while(std::getline(settings,line)){if(!line.empty()&&line.back()=='\r')line.pop_back();
        line.erase(std::remove(line.begin(),line.end(),'"'),line.end());auto fields=nmhv_limit_split_tabs(line);
        require(fields.size()==2,"Invalid direct OPE projection configuration");require(config.emplace(fields[0],fields[1]).second,"Duplicate direct OPE configuration");}
    require(config.at("schema")=="mhv-small-s-v1"&&config.at("endpoint")=="SmallS"&&std::stoi(config.at("max_t"))==2&&std::stoi(config.at("grade"))>=4,"Incompatible direct OPE library");
    auto path=root/("calE_L"+std::to_string((c.h+c.k)/2)+".tsv");c.input(path);
    std::ifstream in(path);require(bool(in),"Missing direct OPE RHS");std::getline(in,line);
    std::vector<direct_ope_target> targets;
    while(std::getline(in,line))if(!line.empty()){
        if(line.back()=='\r')line.pop_back();line.erase(std::remove(line.begin(),line.end(),'"'),line.end());auto fields=nmhv_limit_split_tabs(line);
        require(fields.size()==9&&(fields[8]=="solve"||fields[8]=="check"),"Invalid direct OPE RHS row");
        int t=std::stoi(fields[1]),z=std::stoi(fields[3]),a=std::stoi(fields[4]),b=std::stoi(fields[5]);
        require(std::stoul(fields[0])==(c.h+c.k)/2&&(t==1||t==2)&&std::stoi(fields[2])==t&&z>=t&&a>=0&&b>=0&&size_t(a+b)<=c.h+c.k,"Unsupported direct OPE RHS coefficient");
        if(fields[8]=="solve"&&((include_twist1&&t==1&&z==1)||(t==2&&z==2)))
            targets.push_back({t,z,a,b,rat_t(fields[6]+"/"+fields[7])%c.F.mod});
    }
    require(in.eof()&&!targets.empty(),"No existing low-grade OPE training coefficients");return targets;
}

// A range of existing T1*S1/T2*S2 solve coefficients can be selected to bound
// row storage. No existing check row is used to fit; RHS entries are never
// synthesized for missing coefficients. Origin remains a separate RHS.
inline equations direct_ope_small_s_lowgrade(context& c,bool include_twist1=true,
    size_t first_target=0,size_t target_count=SIZE_MAX){
    auto terminal=joint_eet_flatten_tensor(irr::read_tensor(c.lec/"LEC_0_components.wxf",&c.opt->pool),&c.opt->pool);
    require(terminal.nrow==1&&terminal.ncol==1&&terminal[0].nnz()==1&&terminal[0](0)==0&&terminal[0][0]==rat_t(1),
        "Direct MHV OPE requires the verified unit scalar terminal");
    auto all_targets=direct_ope_lowgrade_targets(c,include_twist1);
    require(first_target<all_targets.size()&&target_count>0,"Invalid direct OPE target range");
    size_t count=std::min(target_count,all_targets.size()-first_target);
    std::vector<direct_ope_target> targets(all_targets.begin()+first_target,all_targets.begin()+first_target+count);
    size_t columns=0;for(size_t b=0;b<6;++b)columns+=c.sew.first[b]*c.sew.last[b];
    c.rows=columns;for(auto& block:c.S)block.clear();equations out(targets.size());
    for(size_t i=0;i<targets.size();++i){auto x=targets[i];out[i].condition="MHV small-S OPE";out[i].sample=first_target+i;
        out[i].rhs=x.rhs;out[i].coefficients.resize(columns);out[i].ope_coordinate={x.z,x.t,x.t,x.a,x.b};}
    auto root=c.data/"ope_mhv_small_s";std::array<ulong,3> q;
    for(size_t j=0;j<3;++j)q[j]=nmod_add(j+2,n_invmod(j+2,c.F.mod.n),c.F.mod);
    axis fa(c.sew.first),la(c.sew.last);size_t direct_nnz=0;auto start=c.now();
    std::vector<std::shared_ptr<const std::array<mat,6>>> tail_chain;
    for(size_t k=1;k<=c.k;++k)tail_chain.push_back(c.modular_recurrence(k,false));
    for(size_t Findex=0;Findex<3;++Findex){
        size_t l=(Findex+1)%3,r=(Findex+2)%3;
        auto den=nmod_mul(nmod_sub(q[Findex],q[l],c.F.mod),nmod_sub(q[Findex],q[r],c.F.mod),c.F.mod);
        require(den,"Degenerate physical-helicity interpolation prime");
        std::array<ulong,3> coeff{0,nmod_neg(nmod_add(q[l],q[r],c.F.mod),c.F.mod),1};
        for(auto& x:coeff)x=nmod_mul(x,n_invmod(den,c.F.mod.n),c.F.mod);
        auto kernel=root/("kernels_F"+std::to_string(Findex+2)+".tsv");c.input(kernel);
        nmhv_ope_germ_t germ(c.h+c.k,c.F.mod.n,kernel,2,4);require_ope_grade_filtration(germ);
        std::string identity="MHV-SmallS:"+hex_sha256_file(kernel)+":2:4";
        auto head_storage=direct_ope_head(c,germ,identity);const auto& head=*head_storage;
        require(head.size()==fa.dim*germ.stride,"Direct OPE head dimensions");
        for(size_t target=0;target<targets.size();++target){auto x=targets[target];
            // MHV has a scalar LEC0 terminal. Reverse the ordered tail map:
            // its forward weights are h+1,...,h+k; the adjoint uses h+k,...,h+1.
            std::vector<ulong> cotangent(germ.stride);cotangent[germ.index(x.t,x.z,x.a,x.b)]=coeff[x.t];
            for(size_t k=1;k<=c.k;++k){axis old(irr::chain_copies(c.lec,k-1,false)),next(irr::chain_copies(c.lec,k,false));
                size_t w=c.h+c.k-k+1;
                cotangent=step(c,irr::matrices{},old,next,cotangent,germ.stride,[&](const ulong* a,ulong* b,size_t letter){
                    integrate_ope_adjoint(germ,a,b,letter,w);},"ope_adjoint_tail_"+std::to_string(k),tail_chain[k-1].get());
            }
            require(cotangent.size()==la.dim*germ.stride,"Direct OPE tail adjoint dimensions");
            // Compress only the small cotangent support, never the FEC head or
            // a 275k identity candidate. One endpoint row is formed at a time.
            std::vector<std::vector<std::pair<size_t,ulong>>> support(la.dim);
            for(size_t i=0;i<la.dim;++i)for(size_t z=0;z<germ.stride;++z)if(cotangent[i*germ.stride+z])
                support[i].push_back({z,cotangent[i*germ.stride+z]});
            size_t offset=0;
            for(size_t b=0;b<6;++b){size_t width=fa.m[b]*la.m[b];
                irr_parallel::rows(width,&c.opt->pool,[&](size_t pair){ulong value=0;
                    for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){
                        size_t hi=fa.off[b]+mu*fa.m[b]+pair/la.m[b],ti=la.off[b]+mu*la.m[b]+pair%la.m[b];ulong image=0;
                        for(auto [z,v]:support[ti])image=germ.add(image,germ.mul(head[hi*germ.stride+z],v));
                        value=germ.add(value,germ.mul(image,mu?3:1));}
                    auto& entry=out[target].coefficients[offset+pair];entry=germ.add(entry,value);
                });offset+=width;
            }
            require(offset==columns,"Direct OPE invariant layout mismatch");
        }
        std::cout<<"DIRECT_OPE_LOWGRADE F="<<Findex+2<<" grade=4 targets="<<targets.size()<<" total_available="<<all_targets.size()
            <<" first_target="<<first_target<<" original_columns="<<columns<<" ordered_adjoint=1 grade_filtration=1"<<std::endl;
    }
    for(const auto& row:out)for(auto x:row.coefficients)direct_nnz+=bool(x);
    c.record("ope_direct_lowgrade",out.size(),columns,direct_nnz,out.size()*columns*sizeof(ulong),start);return out;
}
} // namespace irr_physics
