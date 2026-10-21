#pragma once
#include "mhv_ope_germ.hpp"

inline std::vector<std::string> nmhv_ope_row(std::string line) {
    // getline strips LF, but Linux text streams retain the CR of Windows TSV.
    if(!line.empty()&&line.back()=='\r')line.pop_back();
    line.erase(std::remove(line.begin(),line.end(),'"'),line.end());
    auto fields=nmhv_limit_split_tabs(line);
    if(fields.size()!=9)throw std::runtime_error("Invalid NMHV OPE table: expected nine columns");
    if(fields[8]!="solve"&&fields[8]!="check")throw std::runtime_error("Invalid NMHV OPE role: expected solve or check");
    return fields;
}

// Ordered T,z jets. t+z>=0 is preserved by every kernel and rational prefactor.
struct nmhv_ope_germ_t {
    size_t weight,logs,stride;
    int max_t,grade;
    field_t F;
    struct term_t {int t,z;ulong value;};
    std::array<std::vector<term_t>,9> kernels;
    std::array<ulong,9> at{},az{};
    std::vector<ulong> inverse;
    nmhv_ope_germ_t(size_t w,ulong p,const std::filesystem::path& path,int mt,int gd,int cycle=0)
        :weight(w),logs(w+1),stride((mt+1)*(gd+1)*logs*logs),max_t(mt),grade(gd),F(FIELD_Fp,p),inverse(std::max<size_t>(w,gd)+2) {
        for(size_t i=1;i<inverse.size();++i)inverse[i]=n_invmod(i,p);
        std::ifstream in(path);int l,t,z;std::string value;
        while(in>>l>>t>>z>>value){
            if(!value.empty()&&value.front()=='"')value=value.substr(1,value.size()-2);
            if(l<0||l>=9||t< -1||t>mt)throw std::runtime_error("Invalid NMHV OPE kernel");
            ulong x=rat_t(value)%F.mod;kernels[l].push_back({t,z,x});
            if(t==0&&z==0)at[l]=x;if(t== -1&&z==0)az[l]=x;
        }
        if(!in.eof()||kernels[0].empty())throw std::runtime_error("Missing NMHV OPE kernel: "+path.string());
        // E(w,u,v) and E(v,w,u), with the same (not parity-flipped) square-root branch.
        auto old=kernels;auto old_at=at,old_az=az;
        for(int i=0;i<9;++i){int j=(i/3)*3+(i%3+cycle)%3;kernels[i]=old[j];at[i]=old_at[j];az[i]=old_az[j];}
    }
    ulong add(ulong a,ulong b)const{return nmod_add(a,b,F.mod);}
    ulong mul(ulong a,ulong b)const{return nmod_mul(a,b,F.mod);}
    bool valid(int t,int z)const{return t>=0&&t<=max_t&&z>=-t&&z+t<=grade;}
    size_t index(int t,int z,size_t a,size_t b)const{return ((t*(grade+1)+z+t)*logs+a)*logs+b;}
    void integrate(const ulong* previous,ulong* next,size_t letter,size_t w)const {
        std::fill(next,next+stride,0);
        for(auto term:kernels[letter])for(int t=0;t<=max_t;++t)for(int z=-t;z<=grade-t;++z){
            if((t==0&&z==0)||(t==0)!=(term.t== -1))continue;
            int ot=t-(term.t== -1?0:term.t),oz=z-term.z;if(!valid(ot,oz))continue;
            for(size_t a=0;a<w;++a)for(size_t b=0;b<w-a;++b){auto x=previous[index(ot,oz,a,b)];if(x){auto& y=next[index(t,z,a,b)];y=add(y,mul(x,term.value));}}
        }
        for(int t=0;t<=max_t;++t)for(int z=-t;z<=grade-t;++z){if(t==0&&z==0)continue;
            for(int a=int(w);a>=0;--a)for(int b=int(w)-a;b>=0;--b){auto& x=next[index(t,z,a,b)];
                if(t&&size_t(a+1)<logs)x=nmod_sub(x,mul(a+1,next[index(t,z,a+1,b)]),F.mod);
                if(!t&&size_t(b+1)<logs)x=nmod_sub(x,mul(b+1,next[index(t,z,a,b+1)]),F.mod);
                x=mul(x,inverse[t?t:z]);
            }
        }
        for(size_t a=0;a<w;++a)for(size_t b=0;b<w-a;++b){auto x=mul(previous[index(0,0,a,b)],inverse[w]);
            auto& xt=next[index(0,0,a+1,b)];xt=add(xt,mul(x,at[letter]));
            auto& xz=next[index(0,0,a,b+1)];xz=add(xz,mul(x,az[letter]));}
    }
};
