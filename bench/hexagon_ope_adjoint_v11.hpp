#pragma once
#include "shared/direct_physical_rows.hpp"

namespace hexagon_ope_adjoint {
// Transpose the exact ordered germ integration operator. This is an adjoint
// of the linear map, without complex conjugation. It preserves the original
// integration constants and ordered log recurrence.
inline void integrate(const nmhv_ope_germ_t& germ,const ulong* cotangent,ulong* previous,
                      size_t letter,size_t weight){
    irr::require(weight&&weight<=germ.weight&&letter<9,"OPE adjoint integration range");
    std::fill(previous,previous+germ.stride,0);
    std::vector<ulong> work(cotangent,cotangent+germ.stride);
    // Reverse the pure-log contribution (the final forward operation).
    for(size_t a=0;a<weight;++a)for(size_t b=0;b<weight-a;++b){
        auto value=germ.add(germ.mul(work[germ.index(0,0,a+1,b)],germ.at[letter]),
                            germ.mul(work[germ.index(0,0,a,b+1)],germ.az[letter]));
        previous[germ.index(0,0,a,b)]=germ.add(previous[germ.index(0,0,a,b)],germ.mul(value,germ.inverse[weight]));
    }
    // Forward antiderivative back-substitution goes in descending log order;
    // its transpose goes in ascending order, accumulating into higher logs.
    for(int t=0;t<=germ.max_t;++t)for(int z=-t;z<=germ.grade-t;++z){
        if(t==0&&z==0)continue;
        for(size_t a=0;a<=weight;++a)for(size_t b=0;b<=weight-a;++b){
            auto index=germ.index(t,z,a,b);
            auto value=germ.mul(work[index],germ.inverse[t?t:z]);work[index]=value;
            if(a+b<weight){
                auto next=t?germ.index(t,z,a+1,b):germ.index(t,z,a,b+1);
                auto scale=t?a+1:b+1;
                work[next]=nmod_sub(work[next],germ.mul(scale,value),germ.F.mod);
            }
        }
    }
    // Transpose kernel convolution into the preceding weight's germ.
    for(auto term:germ.kernels[letter])for(int t=0;t<=germ.max_t;++t)
        for(int z=-t;z<=germ.grade-t;++z){
            if((t==0&&z==0)||(t==0)!=(term.t==-1))continue;
            int ot=t-(term.t==-1?0:term.t),oz=z-term.z;if(!germ.valid(ot,oz))continue;
            for(size_t a=0;a<weight;++a)for(size_t b=0;b<weight-a;++b){
                auto index=germ.index(ot,oz,a,b);
                previous[index]=germ.add(previous[index],germ.mul(work[germ.index(t,z,a,b)],term.value));
            }
        }
}
struct Target {int t,z,a,b;ulong rhs;};
inline std::vector<Target> selected_targets(irr_physics::context& c,int grade,size_t maximum){
    auto path=c.data/"ope_mhv_small_s"/("calE_L"+std::to_string((c.h+c.k)/2)+".tsv");c.input(path);
    std::ifstream input(path);std::string line;std::getline(input,line);std::vector<Target> targets;
    while(std::getline(input,line))if(!line.empty()){
        auto fields=nmhv_ope_row(line);int t=std::stoi(fields[1]),z=std::stoi(fields[3]),a=std::stoi(fields[4]),b=std::stoi(fields[5]);
        irr::require(std::stoul(fields[0])==(c.h+c.k)/2,"Direct OPE loop mismatch");
        if(fields[8]!="solve"||t+z>grade||a+b>=int(c.h+c.k))continue;
        irr::require((t==1||t==2)&&std::stoi(fields[2])==t&&z>=t&&a>=0&&b>=0,"Unsupported direct OPE coordinate");
        targets.push_back({t,z,a,b,rat_t(fields[6]+"/"+fields[7])%c.F.mod});
    }
    irr::require(input.eof()&&!targets.empty(),"No known solve-role OPE targets in the selected grade");
    std::stable_sort(targets.begin(),targets.end(),[](const auto& a,const auto& b){
        return std::tuple(a.t+a.z,a.t,a.z,-a.a-a.b,-a.a)<std::tuple(b.t+b.z,b.t,b.z,-b.a-b.b,-b.a);});
    if(maximum&&targets.size()>maximum)targets.resize(maximum);return targets;
}
inline std::vector<ulong> head(irr_physics::context& c,const nmhv_ope_germ_t& germ){
    std::vector<ulong> value(germ.stride);value[germ.index(0,0,0,0)]=1;
    for(size_t w=1;w<=c.h;++w){auto recurrence=c.modular_recurrence(w,true);
        irr_physics::axis old(irr::chain_copies(c.fec,w-1,true)),next(irr::chain_copies(c.fec,w,true));
        value=irr_physics::step(c,irr::matrices{},old,next,value,germ.stride,
            [&](const ulong* a,ulong* b,size_t letter){germ.integrate(a,b,letter,w);},
            "ope_dual_head_"+std::to_string(w),recurrence.get());
    }return value;
}
inline std::vector<ulong> tail(irr_physics::context& c,const nmhv_ope_germ_t& germ,
                             const std::vector<Target>& targets,size_t begin,size_t count){
    const size_t stride=count*germ.stride;std::vector<ulong> value(stride);
    auto saved_samples=c.samples;c.samples=1;auto terminal=irr_physics::terminal(c,false,1);c.samples=saved_samples;
    irr::require(c.mhv&&terminal.size()==1,"Direct OPE adjoint requires the MHV scalar terminal");
    for(size_t j=0;j<count;++j){auto target=targets[begin+j];
        value[j*germ.stride+germ.index(target.t,target.z,target.a,target.b)]=terminal[0];}
    for(size_t k=1;k<=c.k;++k){auto recurrence=c.modular_recurrence(k,false);
        irr_physics::axis old(irr::chain_copies(c.lec,k-1,false)),next(irr::chain_copies(c.lec,k,false));
        size_t weight=c.h+c.k-k+1;
        value=irr_physics::step(c,irr::matrices{},old,next,value,stride,
            [&](const ulong* a,ulong* b,size_t letter){for(size_t j=0;j<count;++j)
                integrate(germ,a+j*germ.stride,b+j*germ.stride,letter,weight);},
            "ope_dual_tail_"+std::to_string(k),recurrence.get());
    }return value;
}
// Original-coordinate OPE rows from one head germ and small batches of
// endpoint cotangents. Never enumerate a candidate or identity basis. Each
// F=2,3,4 contribution is interpolated exactly as in ope_small_s(). Select
// only known solve-role RHS rows; held-out rows and unknown T^2 F^0 stay out.
inline irr_physics::equations selected_rows(irr_physics::context& c,int grade=4,
                                            size_t maximum_targets=32,size_t target_batch=4){
    irr::require(c.mhv&&grade>=2&&target_batch,"Invalid direct OPE projection policy");
    irr_physics::direct_coordinate_space(c);auto targets=selected_targets(c,grade,maximum_targets);
    irr_physics::equations rows(targets.size());
    for(size_t i=0;i<rows.size();++i){auto target=targets[i];rows[i].condition="MHV small-S OPE adjoint";
        rows[i].sample=i;rows[i].rhs=target.rhs;rows[i].ope_coordinate={target.z,target.t,target.t,target.a,target.b};
        rows[i].coefficients.resize(c.rows);}
    std::array<ulong,3> q;for(size_t j=0;j<3;++j)q[j]=nmod_add(j+2,n_invmod(j+2,c.F.mod.n),c.F.mod);
    irr_physics::axis fa(c.sew.first),la(c.sew.last);
    for(size_t f=0;f<3;++f){size_t l=(f+1)%3,r=(f+2)%3;
        ulong den=nmod_mul(nmod_sub(q[f],q[l],c.F.mod),nmod_sub(q[f],q[r],c.F.mod),c.F.mod);
        std::array<ulong,3> weights{0,nmod_neg(nmod_add(q[l],q[r],c.F.mod),c.F.mod),1};
        for(auto& x:weights)x=nmod_mul(x,n_invmod(den,c.F.mod.n),c.F.mod);
        auto kernel=c.data/"ope_mhv_small_s"/("kernels_F"+std::to_string(f+2)+".tsv");c.input(kernel);
        nmhv_ope_germ_t germ(c.h+c.k,c.F.mod.n,kernel,2,grade);
        // A smaller grade is exact only if integration never lowers t+z.
        for(const auto& terms:germ.kernels)for(auto term:terms)
            irr::require(term.t==-1?term.z>=0:term.t+term.z>=0,"OPE kernel lowers the selected jet grade");
        auto first=head(c,germ);
        for(size_t begin=0;begin<targets.size();begin+=target_batch){size_t count=std::min(target_batch,targets.size()-begin);
            auto last=tail(c,germ,targets,begin,count);using Jet=std::vector<std::pair<size_t,ulong>>;
            std::vector<Jet> sparse(la.dim*count);
            for(size_t j=0;j<la.dim;++j)for(size_t sample=0;sample<count;++sample)
                for(size_t z=0;z<germ.stride;++z){auto x=last[(j*count+sample)*germ.stride+z];if(x)sparse[j*count+sample].push_back({z,x});}
            last.clear();last.shrink_to_fit();size_t offset=0;
            for(size_t b=0;b<6;++b){size_t width=c.sew.first[b]*c.sew.last[b];
                irr_parallel::rows(width,&c.opt->pool,[&](size_t coordinate){
                    for(size_t sample=0;sample<count;++sample){ulong value=0;
                        for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){
                            size_t fi=fa.off[b]+mu*fa.m[b]+coordinate/la.m[b],li=la.off[b]+mu*la.m[b]+coordinate%la.m[b];
                            ulong dot=0;for(auto [z,x]:sparse[li*count+sample])dot=germ.add(dot,germ.mul(first[fi*germ.stride+z],x));
                            value=germ.add(value,germ.mul(dot,mu?3:1));
                        }
                        auto scale=weights[targets[begin+sample].t];auto& entry=rows[begin+sample].coefficients[offset+coordinate];
                        entry=germ.add(entry,germ.mul(value,scale));
                    }
                });offset+=width;
            }
            std::cout<<"DIRECT_OPE_ADJOINT F="<<f+2<<" grade="<<grade<<" targets="<<begin<<".."<<begin+count
                     <<" columns="<<c.rows<<std::endl;
        }
    }return rows;
}
} // namespace hexagon_ope_adjoint
