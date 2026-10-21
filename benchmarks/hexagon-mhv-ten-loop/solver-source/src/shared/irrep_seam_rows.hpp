#pragma once
#include "reduced_d3/irr_sew_prereduce.hpp"
#include "staged_kernel.hpp"

// The right-side row basis is computed in this field only. Its contraction
// with the left recurrence is replayed in bounded tiles instead of stored.
class IrrepSeamRows {
    using MM=staged_kernel::MM;using MV=staged_kernel::MV;
    struct Column {size_t block,copy,channel;};
    struct Partial {MM local;std::vector<Column> columns;size_t copies,offset;};
    const field_t& F;
    bool contract_left=false;
    std::array<MM,6> first;
    std::array<std::vector<std::vector<size_t>>,6> open_indices;
    std::vector<Partial> partials;
    irr::copies colbase,last_copies;
    size_t rows=0,columns=0;
    const Partial& part(size_t id)const{
        for(const auto& p:partials)if(id>=p.offset&&id<p.offset+p.local.nrow*p.copies)return p;
        throw std::runtime_error("Irrep seam row outside layout");
    }
    size_t cost(size_t id)const{
        const auto& p=part(id);size_t row=(id-p.offset)/p.copies,copy=(id-p.offset)%p.copies,n=0;
        for(auto [c,x]:p.local[row]){const auto& v=p.columns[c];n+=first[v.block][open_indices[v.block][v.channel][copy]].nnz();}return n;
    }
    MV build(size_t id,size_t limit,const streaming_echelon::Substitution* substitution=nullptr,const std::vector<int32_t>* coordinate=nullptr)const{
        if(limit!=SIZE_MAX&&cost(id)>limit)return {};
        const auto& p=part(id);size_t row=(id-p.offset)/p.copies,copy=(id-p.offset)%p.copies;MV out;
        for(auto [c,x]:p.local[row]){const auto& v=p.columns[c];
            for(auto [j,y]:first[v.block][open_indices[v.block][v.channel][copy]]){
                size_t col=colbase[v.block]+(contract_left?v.copy*last_copies[v.block]+size_t(j):size_t(j)*last_copies[v.block]+v.copy);
                ulong value=nmod_mul(x,y,F.mod);int32_t mapped=coordinate?(*coordinate)[col]:int32_t(col);if(mapped<0)continue;
                if(substitution){int32_t previous=mapped;mapped=substitution->image[previous];if(mapped<0)continue;value=nmod_mul(value,substitution->scale[previous],F.mod);}
                if(value)out.push_back(mapped,value);
            }}
        staged_kernel::normalize_mod(out,F);return out;
    }
public:
    IrrepSeamRows(const irr::channels_t& cg,const irr::seam_template_t& t,const irr::seam_layout_t& layout,
        const irr::copies& of,const irr::copies& ol,const irr::matrices& left,const irr::matrices& right,
        const field_t& field,rref_option_t& opt,bool use_left=false):F(field),contract_left(use_left),colbase(layout.colbase),columns(layout.columns){
        const auto& contracted=contract_left?left:right;const auto& remaining=contract_left?right:left;
        const auto& old=contract_left?of:ol;const auto& open=contract_left?ol:of;
        std::array<MM,6> chain;std::array<std::vector<irr::decoded_channel_t>,6> decoded;
        for(size_t b=0;b<6;++b){chain[b]=nmhv_sew_modular_reduce(contracted[b],F,&opt->pool);
            first[b]=nmhv_sew_modular_reduce(remaining[b],F,&opt->pool);first[b].transpose_replace();
            decoded[b]=irr::decode_channels(cg,b,old);auto channels=irr::decode_channels(cg,b,open);
            open_indices[b].resize(cg.cg[b].size());for(size_t j=0;j<channels.size();++j)open_indices[b][channels[j].channel].push_back(j);
            last_copies[b]=right[b].nrow;}
        for(size_t a=0;a<6;++a){if(!open[a])continue;std::vector<size_t> base(t.rows.size(),SIZE_MAX);size_t nr=0;
            for(size_t r=0;r<t.rows.size();++r)if((contract_left?t.rows[r].right:t.rows[r].left)==a){base[r]=nr;nr+=old[contract_left?t.rows[r].left:t.rows[r].right];}
            Partial p;p.copies=open[a];p.offset=rows;
            for(size_t b=0;b<6;++b)for(size_t i=0;i<chain[b].nrow;++i)for(size_t g:cg.channels[b][a])p.columns.push_back({b,i,g});
            if(!nr||p.columns.empty())continue;MM image(p.columns.size(),nr);
            opt->pool.detach_loop(size_t(0),p.columns.size(),[&](size_t c){auto [b,i,g]=p.columns[c];
                for(auto [j,x]:chain[b][i]){const auto& d=decoded[b][j];for(const auto& route:(contract_left?t.routes[b][d.channel][g]:t.routes[b][g][d.channel]))if(base[route.row]!=SIZE_MAX){
                    ulong value=nmod_mul(x,route.value%F.mod,F.mod);if(value)image[c].push_back(int32_t(base[route.row]+d.copy),value);}}
                irr::merge_mod(image[c],F);});opt->pool.wait();image.transpose_replace();image.sort_rows_by_nnz();
            auto pivots=sparse_mat_rref(image,F,opt);p.local=nmhv_sew_modular_pivot_rows(std::move(image),pivots);
            rows+=p.local.nrow*p.copies;partials.push_back(std::move(p));
        }
        std::cout<<"FACTORED_SEAM rows="<<rows<<" columns="<<columns<<" contraction="<<(contract_left?"left":"right")<<" full_product_stored=0"<<std::endl;
    }
    size_t ncol()const{return columns;}size_t nrow()const{return rows;}
    template<class Emit> void modular(structured_kernel::Relations&,const std::vector<int32_t>& coordinate,
        const field_t& field,Emit&& emit,thread_pool* pool,size_t limit=SIZE_MAX,size_t=0,
        const streaming_echelon::Substitution* substitution=nullptr,std::vector<unsigned char>* seen=nullptr,unsigned shift=0)const{
        symrep::require(field.mod.n==F.mod.n&&coordinate.size()==columns,"Factored seam field/layout");
        for(size_t c=0;c<coordinate.size();++c)symrep::require(coordinate[c]>=-1&&coordinate[c]<int32_t(columns),"Factored seam coordinate map");
        for(size_t begin=0;begin<rows;){size_t end=begin,terms=0;
            while(end<rows&&end-begin<4096){size_t n=std::min(cost(end),limit==SIZE_MAX?SIZE_MAX:limit+1);
                if(end>begin&&n>2*1024*1024-std::min(terms,size_t(2*1024*1024)))break;terms+=n;++end;}
            std::vector<MV> tile(end-begin);auto fill=[&](size_t i){size_t id=begin+i;
                if((seen&&(*seen)[id])||!streaming_echelon::sampled_equation(id,shift,UINT64_C(0x9e3779b97f4a7c15)))return;
                if(limit!=SIZE_MAX&&cost(id)>limit)return;
                tile[i]=build(id,limit,substitution,&coordinate);if(seen)(*seen)[id]=1;};
            if(pool){pool->detach_loop(size_t(0),tile.size(),fill);pool->wait();}else for(size_t i=0;i<tile.size();++i)fill(i);
            for(auto& row:tile)if(row.nnz())emit(std::move(row));begin=end;
        }
    }
    template<class Emit> void residual_transpose(const MM& kt,size_t count,size_t block,Emit&& emit,thread_pool* pool)const{
        size_t begin=block*128,end=std::min(rows,begin+128);std::vector<MV> tile(end-begin);
        auto fill=[&](size_t i){auto equation=build(begin+i,SIZE_MAX);std::vector<ulong> values(count);
            for(auto [c,x]:equation)for(auto [j,y]:kt[c])values[j]=nmod_add(values[j],nmod_mul(x,y,F.mod),F.mod);
            for(size_t j=0;j<count;++j)if(values[j])tile[i].push_back(int32_t(j),values[j]);};
        if(pool){pool->detach_loop(size_t(0),tile.size(),fill);pool->wait();}else for(size_t i=0;i<tile.size();++i)fill(i);
        for(auto& row:tile)emit(std::move(row));
    }
    size_t residual_blocks()const{return (rows+127)/128;}
    template<class Emit> void residual(const MM& k,const field_t& field,size_t block,Emit&& emit,thread_pool*)const{
        symrep::require(field.mod.n==F.mod.n&&k.ncol==columns&&block<residual_blocks(),"Factored seam residual layout");auto kt=k.transpose();
        for(size_t id=block*128;id<std::min(rows,block*128+128);++id){auto equation=build(id,SIZE_MAX);std::vector<ulong> values(k.nrow);
            for(auto [c,x]:equation)for(auto [j,y]:kt[c])values[j]=nmod_add(values[j],nmod_mul(x,y,F.mod),F.mod);
            MV row;for(size_t j=0;j<values.size();++j)if(values[j])row.push_back(int32_t(j),values[j]);if(row.nnz())emit(std::move(row));}
    }
};
inline staged_kernel::MM factored_staged_kernel(IrrepSeamRows& rows,const field_t& F,rref_option_t& opt){
    auto source=staged_kernel::stack(rows);structured_kernel::Relations relations(rows.ncol());std::vector<int32_t> coordinate(rows.ncol());std::iota(coordinate.begin(),coordinate.end(),0);
    staged_kernel::Options options;staged_kernel::RefinedRows<decltype(source)> refined{source,opt,options};staged_kernel::MM reduced(0,rows.ncol());
    refined.modular(relations,coordinate,F,[&](staged_kernel::MV&& row){reduced.rows.push_back(std::move(row));},&opt->pool);reduced.nrow=reduced.rows.size();
    std::vector<std::vector<pivot_t<int32_t>>> pivots;if(!certified_kernel::private_unit_pivots(reduced,pivots))pivots=sparse_mat_rref(reduced,F,opt);
    return staged_kernel::modular_kernel(reduced,pivots,F);
}
