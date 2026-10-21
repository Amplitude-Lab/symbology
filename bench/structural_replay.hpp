#pragma once
#include "symrep.hpp"
#include <fstream>
#include <numeric>
#include <stdexcept>

namespace structural_replay {
using I=int32_t;
using M=sparse_mat<ulong,I>;
using Groups=std::vector<std::vector<pivot_t<I>>>;
struct Mismatch:std::runtime_error{using std::runtime_error::runtime_error;};
inline void guard(bool ok,const char* why){if(!ok)throw Mismatch(why);}
struct Batch {size_t inputs=0;uint64_t support=0;std::vector<size_t> selected;Groups groups;};
struct Trace {
    size_t columns=0,cursor=0;bool replay=false;
    std::vector<Batch> batches;
    void complete()const{guard(!replay||cursor==batches.size(),"Replay batch count changed");}
    void write(std::ostream& out,const std::string& identity)const{
        out<<"sparse-structural-replay-v1\n"<<identity<<'\n'<<columns<<' '<<batches.size()<<'\n';
        for(const auto& b:batches){out<<b.inputs<<' '<<b.support<<' '<<b.selected.size()<<' '<<b.groups.size()<<'\n';
            for(auto r:b.selected)out<<r<<' ';out<<'\n';
            for(const auto& g:b.groups){out<<g.size();for(auto[r,c]:g)out<<' '<<r<<' '<<c;out<<'\n';}}
        guard(bool(out),"Replay plan write failed");
    }
    static Trace read(std::istream& in,const std::string& identity){
        std::string schema,key;std::getline(in,schema);std::getline(in,key);
        guard(schema=="sparse-structural-replay-v1"&&key==identity,"Replay identity changed");
        Trace t;t.replay=true;size_t count;in>>t.columns>>count;
        guard(bool(in)&&t.columns<=INT32_MAX&&count<1000000,"Replay plan dimensions");t.batches.resize(count);
        for(auto& b:t.batches){size_t nr,ng;in>>b.inputs>>b.support>>nr>>ng;
            guard(bool(in)&&b.inputs<=INT32_MAX&&nr<=t.columns&&nr<=b.inputs&&ng<=nr+1,"Replay batch dimensions");
            b.selected.resize(nr);std::vector<unsigned char> rows(b.inputs),cols(t.columns),mapped(nr);
            for(auto& r:b.selected){in>>r;guard(bool(in)&&r<b.inputs&&!rows[r],"Replay selected row invalid");rows[r]=1;}
            b.groups.resize(ng);size_t total=0;
            for(auto& g:b.groups){size_t n;in>>n;guard(bool(in)&&n<=nr-total,"Replay group dimensions");total+=n;g.resize(n);
                for(auto& p:g){in>>p.r>>p.c;guard(bool(in)&&p.r>=0&&size_t(p.r)<nr&&p.c>=0&&size_t(p.c)<t.columns&&!mapped[p.r]&&!cols[p.c],"Replay pivot invalid");mapped[p.r]=1;cols[p.c]=1;}}
            guard(total==nr,"Replay missing pivot");}
        std::string extra;guard(!(in>>extra),"Replay plan trailing data");return t;
    }
};
inline uint64_t support_hash(const M& a){
    uint64_t h=1469598103934665603ULL;
    auto mix=[&](uint64_t x){h^=x;h*=1099511628211ULL;};
    mix(a.ncol);mix(a.nrow);for(const auto& row:a.rows){mix(row.nnz());for(auto c:row.index_span())mix(uint32_t(c));}return h;
}
// Replay only the learned pivot groups. Numerical coefficients are freshly
// supplied. The caller independently reduces ALL discarded input equations
// against the resulting basis; a rank gain or changed RHS triggers fallback.
inline Groups forward(M& a,const Batch& plan,const field_t& F,rref_option* opt){
    const auto& groups=plan.groups;
    guard(a.nrow==plan.selected.size(),"Replay selected shape changed");
    if(groups.empty()){guard(!a.nrow,"Replay empty groups");return groups;}
    std::vector<unsigned char> solved(a.ncol),done(groups[0].size());size_t remaining=done.size();
    while(remaining){size_t before=remaining;
        for(size_t i=0;i<done.size();++i)if(!done[i]){auto[r,c]=groups[0][i];size_t count=0;I column=-1;
            for(auto[j,x]:a[r])if(x&&!solved[j]){++count;column=j;}
            if(count==1){guard(column==c,"Replay singleton pivot changed");done[i]=1;solved[c]=1;--remaining;}}
        guard(remaining<before,"Replay singleton rank lost");}
    for(auto[r,c]:groups[0]){a[r].zero();a[r].push_back(c,1);}
    schur_helper_buffer<ulong,I> buffer(opt->pool.get_thread_count(),a.ncol);
    for(size_t g=0;g<groups.size();++g){const auto& group=groups[g];if(group.empty())continue;
        // SparseRREF groups form a triangular pivot block, possibly with
        // upper-group dependencies; schur_complete respects their order.
        for(auto[r,c]:group){auto x=a[r].find(c);guard(x&&*x,"Replay zero pivot");}
        opt->pool.detach_loop(size_t(0),group.size(),[&](size_t i){auto[r,c]=group[i];sparse_vec_rescale(a[r],n_invmod(*a[r].find(c),F.mod.n),F);});opt->pool.wait();
        std::vector<size_t> rest;for(size_t j=g+1;j<groups.size();++j)for(auto p:groups[j])rest.push_back(p.r);
        opt->pool.detach_blocks<size_t>(0,rest.size(),[&](size_t begin,size_t end){schur_helper<ulong,I> helper(buffer,SparseRREF::thread_id());for(size_t j=begin;j<end;++j)schur_complete(a,rest[j],group,F,helper);});opt->pool.wait();
    }
    return groups;
}
}
