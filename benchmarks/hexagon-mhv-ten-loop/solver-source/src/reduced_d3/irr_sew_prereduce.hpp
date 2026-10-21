#pragma once
#include "irr_sewing.hpp"
#include "nmhv_irrep_sew.hpp"

namespace irr {
void merge_mod(sparse_vec<ulong,int32_t>& row,const field_t& F) {
    row.sort_indices();
    size_t n=0;
    for(size_t p=0;p<row.nnz();) {
        auto col=row(p);ulong value=0;
        do {value=nmod_add(value,row[p++],F.mod);} while(p<row.nnz()&&row(p)==col);
        if(value){row(n)=col;row[n++]=value;}
    }
    row.resize(n);
}
template<class T> void stats(const std::string& label,const sparse_mat<T,int32_t>& mat) {
    std::cout<<"STATS "<<label<<" rows="<<mat.nrow<<" cols="<<mat.ncol<<" nnz="<<mat.nnz()<<std::endl;
}

// Contract one side of the existing singlet routes, retaining the open channel.
// For left reduction, group by the old-last irrep; for right, by old-first.
// All target irreps sharing that open irrep stay in the SAME partial matrix.
sparse_mat<ulong,int32_t> partial_seam(bool left,const irr::channels_t& cg,
    const irr::seam_template_t& t,const irr::seam_layout_t& layout,
    const irr::copies& of,const irr::copies& ol,const irr::matrices& first,
    const irr::matrices& last,const field_t& F,rref_option_t opt) {
    irr_parallel::phase timing(left?"left_assembly_and_prereduce":"right_assembly_and_prereduce");
    const auto& contracted=left?first:last;const auto& remaining=left?last:first;
    const auto& old=left?of:ol;const auto& open=left?ol:of;
    std::array<sparse_mat<ulong,int32_t>,6> chain,other;
    std::array<std::vector<irr::decoded_channel_t>,6> decoded;
    std::array<std::vector<std::vector<size_t>>,6> open_indices;
    for(size_t b=0;b<6;++b) {
        chain[b]=nmhv_sew_modular_reduce(contracted[b],F,&opt->pool);
        other[b]=nmhv_sew_modular_reduce(remaining[b],F,&opt->pool);other[b].transpose_replace();
        decoded[b]=irr::decode_channels(cg,b,old);
        auto channels=irr::decode_channels(cg,b,open);open_indices[b].resize(cg.cg[b].size());
        for(size_t j=0;j<channels.size();++j)open_indices[b][channels[j].channel].push_back(j);
    }
    sparse_mat<ulong,int32_t> result(0,layout.columns);
    struct column_t{size_t block,copy,channel;};
    for(size_t a=0;a<6;++a) {
        if(!open[a])continue;
        std::vector<size_t> base(t.rows.size(),SIZE_MAX);size_t rows=0;
        for(size_t r=0;r<t.rows.size();++r) {
            const auto& row=t.rows[r];
            if((left?row.right:row.left)==a){base[r]=rows;rows+=old[left?row.left:row.right];}
        }
        std::vector<column_t> columns;
        for(size_t b=0;b<6;++b)for(size_t i=0;i<chain[b].nrow;++i)
            for(size_t g:cg.channels[b][a])columns.push_back({b,i,g});
        if(!rows||columns.empty())continue;
        sparse_mat<ulong,int32_t> mat(columns.size(),rows);
        opt->pool.detach_loop(0,columns.size(),[&](size_t c){
            const auto [b,i,g]=columns[c];
            for(auto [j,x]:chain[b][i]) {
                const auto& d=decoded[b][j];
                const auto& routes=left?t.routes[b][d.channel][g]:t.routes[b][g][d.channel];
                for(const auto& route:routes)if(base[route.row]!=SIZE_MAX) {
                    ulong value=nmod_mul(x,route.value%F.mod,F.mod);
                    if(value)mat[c].push_back(int32_t(base[route.row]+d.copy),value);
                }
            }
            merge_mod(mat[c],F);
        });opt->pool.wait();mat.transpose_replace();
        auto label=std::string(left?"left_":"right_")+irr::specs()[a].name;
        stats(label+"_before",mat);mat.sort_rows_by_nnz();
        auto pivots=sparse_mat_rref(mat,F,opt);
        auto reduced=nmhv_sew_modular_pivot_rows(std::move(mat),pivots);mat.clear();
        stats(label+"_after",reduced);
        size_t offset=result.rows.size(),count=reduced.nrow*open[a];result.rows.resize(offset+count);result.nrow=result.rows.size();
        opt->pool.detach_loop(0,count,[&](size_t i){
            size_t row=i/open[a],copy=i%open[a];auto& out=result[offset+i];
            for(auto [c,x]:reduced[row]) {
                const auto [b,j,g]=columns[c];
                for(auto [v,y]:other[b][open_indices[b][g][copy]]) {
                    size_t col=layout.colbase[b]+(left?j*last[b].nrow+v:size_t(v)*last[b].nrow+j);
                    out.push_back(int32_t(col),nmod_mul(x,y,F.mod));
                }
            }
            merge_mod(out,F);
        });opt->pool.wait();
    }
    return result;
}


}

