#pragma once
#include "reduced_d3/irr_physics_pipeline.hpp"
#include <numeric>

// CONDITIONAL empirical MHV tail restriction from Dixon--Liu,
// arXiv:2308.08199, equations (76)--(79), with the full D3 orbit.
// No all-orders proof or ten-loop amplitude result is claimed by this helper.
namespace mhv_double_bonus {
using Q=rat_t;using QM=sparse_mat<Q,int32_t>;using PM=sparse_mat<ulong,int32_t>;
using irr::require;
inline size_t power9(size_t w){size_t n=1;for(size_t j=0;j<w;++j)n*=9;return n;}
struct Axis{irr::copies m,off{};size_t dimension=0;explicit Axis(irr::copies counts):m(counts){
    for(size_t b=0;b<6;++b){off[b]=dimension;dimension+=m[b]*irr::specs()[b].irrep_dimension;}}};

// New letters: a=u/(vw), b=v/(wu), c=w/(uv), m_i=(1-u_i)/u_i, y_i.
// B[new,old] obeys dlog(new)=B*dlog(old), with old order
// (u,v,w,1-u,1-v,1-w,y_u,y_v,y_w). Thus coefficient functionals
// transform as R_new=B*R_old*B^T, NOT with the inverse alphabet matrix.
inline QM alphabet_matrix(){
    QM B(9,9);for(size_t i=0;i<3;++i){for(size_t j=0;j<3;++j)B[i].push_back(j,Q(i==j?1:-1));
        B[3+i].push_back(i,Q(-1));B[3+i].push_back(3+i,Q(1));B[6+i].push_back(6+i,Q(1));}
    return B;
}
inline QM published_bonus_relations(){
    QM R(4,81);auto add=[&](size_t r,size_t i,size_t j,int x){R[r].push_back(i*9+j,Q(x));};
    // Eq76: E^{yv,yu}-E^{yu,yv}-E^{u,u}-E^{1-u,u}+E^{v,v}+E^{1-v,v}=0.
    add(0,7,6,1);add(0,6,7,-1);add(0,0,0,-1);add(0,3,0,-1);add(0,1,1,1);add(0,4,1,1);
    // Eq77 and Eq78 are the two cyclic images of Eq76.
    add(1,8,7,1);add(1,7,8,-1);add(1,1,1,-1);add(1,4,1,-1);add(1,2,2,1);add(1,5,2,1);
    add(2,6,8,1);add(2,8,6,-1);add(2,2,2,-1);add(2,5,2,-1);add(2,0,0,1);add(2,3,0,1);
    // Eq79: E^{1-u,yu}+E^{yu,u}+E^{v,yw}-E^{yw,v}=0.
    add(3,3,6,1);add(3,6,0,1);add(3,1,8,1);add(3,8,1,-1);
    for(auto& row:R.rows)fec_irrep_merge_entries(row);return R;
}
inline QM published_qbar_calibration(){
    // Eq70--75 are expected to be automatic in generic LEC2. They calibrate
    // alphabet duality and chronological last-slot orientation independently
    // of imposing any bonus relation.
    QM R(6,81);auto add=[&](size_t r,size_t i,size_t j,int x){R[r].push_back(i*9+j,Q(x));};
    add(0,7,6,1);add(0,8,6,-1);add(0,7,8,1);add(0,7,7,-1);add(0,2,0,-1);
    add(1,8,7,1);add(1,6,7,-1);add(1,8,6,1);add(1,8,8,-1);add(1,0,1,-1);
    add(2,6,8,1);add(2,7,8,-1);add(2,6,7,1);add(2,6,6,-1);add(2,1,2,-1);
    add(3,3,6,1);add(3,7,1,-1);add(3,8,0,1);add(3,1,8,1);
    add(4,4,7,1);add(4,8,2,-1);add(4,6,1,1);add(4,2,6,1);
    add(5,5,8,1);add(5,6,0,-1);add(5,7,2,1);add(5,0,7,1);
    for(auto& row:R.rows)fec_irrep_merge_entries(row);return R;
}
inline QM full_d3_orbit(const QM& base){
    // Possible uniform y inversion under physical dihedral generators changes
    // these parity-homogeneous relations only by an overall sign: Eq76--78
    // have zero/two y letters, while every Eq79 term has exactly one. Thus
    // all six permutations generate the same relation space in either branch.
    QM R(base.nrow*6,81);std::array<size_t,3> permutation{0,1,2};size_t orbit=0;
    do{for(size_t r=0;r<base.nrow;++r)for(auto [pair,value]:base[r]){
        size_t i=pair/9,j=pair%9;i=(i/3)*3+permutation[i%3];j=(j/3)*3+permutation[j%3];
        R[orbit*base.nrow+r].push_back(i*9+j,value);
    }++orbit;}while(std::next_permutation(permutation.begin(),permutation.end()));
    require(orbit==6,"D3 orbit size");for(auto& row:R.rows)fec_irrep_merge_entries(row);return R;
}
inline QM coefficient_dual_transform(const QM& old,const QM& B){
    QM R(old.nrow,81);for(size_t r=0;r<old.nrow;++r)for(auto [pair,value]:old[r]){
        size_t i=pair/9,j=pair%9;for(size_t a=0;a<9;++a){auto x=B[a].find(i);if(!x)continue;
            for(size_t b=0;b<9;++b){auto y=B[b].find(j);if(y)R[r].push_back(a*9+b,value*(*x)*(*y));}}
    }for(auto& row:R.rows)fec_irrep_merge_entries(row);return R;
}

inline std::vector<QM> expand_lec_words(const fs::path& lec,size_t maximum,const irr::channels_t& cg,rref_option_t opt){
    require(maximum>=2&&maximum<=4,"Bonus audit expands only tiny LEC2--4 word spaces");
    std::vector<QM> words;words.emplace_back(1,1);words[0][0].push_back(0,Q(1));
    require(irr::chain_copies(lec,0,false)==irr::copies{1,0,0,0,0,0},"Bonus audit requires scalar MHV terminal");
    for(size_t w=1;w<=maximum;++w){auto old=irr::chain_copies(lec,w-1,false),current=irr::chain_copies(lec,w,false);
        Axis axis(current);auto K=irr::read_recurrence(lec,w,false,cg,&opt->pool);QM next(axis.dimension,power9(w));
        for(size_t b=0;b<6;++b){auto expanded=irr::expand(cg,b,old,K[b]);
            require(expanded.nrow==current[b]*irr::specs()[b].irrep_dimension&&expanded.ncol==words[w-1].nrow*9,"Expanded LEC/word layout");
            for(size_t row=0;row<expanded.nrow;++row)for(auto [column,value]:expanded[row])
                for(auto [word,previous]:words[w-1][column/9])
                    // LEC growth prepends the new letter to the existing suffix.
                    next[axis.off[b]+row].push_back((column%9)*power9(w-1)+word,value*previous);
        }
        for(auto& row:next.rows)fec_irrep_merge_entries(row);words.push_back(std::move(next));
    }return words;
}

inline std::array<QM,6> constraint_images(const QM& words,const irr::copies& copies,size_t w,const QM& R){
    require(w>=2&&words.ncol==power9(w)&&R.ncol==81,"Bonus constraint word dimensions");Axis axis(copies);
    std::array<std::vector<std::pair<size_t,Q>>,81> routes;
    for(size_t r=0;r<R.nrow;++r)for(auto [pair,value]:R[r])routes[pair].push_back({r,value});
    std::array<QM,6> images;size_t prefixes=power9(w-2);
    for(size_t b=0;b<6;++b){size_t components=irr::specs()[b].irrep_dimension;
        images[b]=QM(copies[b],components*prefixes*R.nrow);
        for(size_t mu=0;mu<components;++mu)for(size_t copy=0;copy<copies[b];++copy)
            for(auto [word,value]:words[axis.off[b]+mu*copies[b]+copy])for(const auto& [r,scale]:routes[word%81])
                images[b][copy].push_back((mu*prefixes+word/81)*R.nrow+r,value*scale);
        for(auto& row:images[b].rows)fec_irrep_merge_entries(row);
    }return images;
}
inline void calibrate_single_final_entry(const QM& words,const QM& B){
    require(words.ncol==9,"Single final-entry calibration width");
    for(const auto& row:words.rows)for(size_t i=0;i<3;++i){Q residual=0;
        for(auto [letter,value]:row){Q functional=0;if(auto x=B[letter].find(i))functional+=*x;
            if(auto x=B[letter].find(3+i))functional+=*x;residual+=value*functional;}
        require(residual==0,"Published single final entry fails coefficient-dual alphabet calibration");
    }
}
struct Restriction{size_t weight;irr::copies original{},reduced{};std::array<QM,6> images,lift,coordinates;};
inline Restriction restrict_tail(const QM& words,const irr::copies& copies,size_t weight,const QM& constraints,rref_option_t& opt){
    Restriction result;result.weight=weight;result.original=copies;result.images=constraint_images(words,copies,weight,constraints);
    const field_t QQ(FIELD_QQ);
    for(size_t b=0;b<6;++b){auto equations=result.images[b].transpose();
        result.lift[b]=nmhv_sew_kernel_rows(std::move(equations),QQ,opt);
        // The generic kernel helper clears both dimensions for a zero
        // nullspace. Preserve its known ambient copy dimension explicitly.
        if(!result.lift[b].nrow)result.lift[b]=QM(0,copies[b]);
        result.reduced[b]=result.lift[b].nrow;
        require(result.lift[b].ncol==copies[b],"Bonus copy-kernel dimensions");
        require(!sparse_mat_mul(result.lift[b],result.images[b],QQ,&opt->pool).nnz(),"Exact bonus copy-kernel residual");
        result.coordinates[b]=irr::right_inverse(result.lift[b],opt);
        require(nmhv_sew_matrices_equal(sparse_mat_mul(result.lift[b],result.coordinates[b],QQ,&opt->pool),joint_eet_identity<Q,int32_t>(result.reduced[b])),"Bonus copy coordinate chart");
        std::cout<<"BONUS_COPY_KERNEL weight="<<weight<<" block="<<irr::specs()[b].name<<" old="<<copies[b]<<" new="<<result.reduced[b]
            <<" irrep_components="<<irr::specs()[b].irrep_dimension<<" common_copy_transform=1 exact_residual=0"<<std::endl;
    }return result;
}
inline size_t dimension(const irr::copies& copies){return Axis(copies).dimension;}
inline void export_matrix(const fs::path& path,const QM& matrix,rref_option_t opt){
    sparse_tensor<Q,int32_t,SPARSE_CSR> tensor(matrix,&opt->pool);irr::write_tensor(path,tensor);
    auto readback=joint_eet_flatten_tensor(irr::read_tensor(path,&opt->pool),&opt->pool);
    require(nmhv_sew_matrices_equal(readback,matrix),"Bonus map WXF readback mismatch");
}
inline irr_stage::fields native_check(const Restriction& restriction,size_t h,const fs::path& fec,const fs::path& lec,
    const fs::path& native,ulong prime,rref_option_t opt){
    size_t k=restriction.weight;require(prime>2&&n_is_prime(prime),"Bonus native-control prime");field_t F(FIELD_Fp,prime);
    auto first=irr::chain_copies(fec,h,true),last=irr::chain_copies(lec,k,false);require(last==restriction.original,"Bonus native-control LEC convention mismatch");
    auto cache=read_cache(native,prime,h,k);size_t width=0;for(size_t b=0;b<6;++b)width+=first[b]*last[b];require(cache.ncol==width,"Bonus native control columns");
    size_t offset=0,visited=0,reduced_columns=0;
    for(size_t b=0;b<6;++b){PM C(first[b],last[b]);
        auto A=nmhv_sew_modular_reduce(restriction.images[b],F,&opt->pool);
        auto K=nmhv_sew_modular_reduce(restriction.lift[b],F,&opt->pool),Qmap=nmhv_sew_modular_reduce(restriction.coordinates[b],F,&opt->pool);
        require(nmhv_sew_matrices_equal(sparse_mat_mul(K,Qmap,F,&opt->pool),joint_eet_identity<ulong,int32_t>(K.nrow)),"Bonus coordinate chart has bad control prime");
        for(size_t solution=0;solution<cache.nrow;++solution){C=PM(first[b],last[b]);
            for(auto [column,value]:cache[solution])if(size_t(column)>=offset&&size_t(column)<offset+first[b]*last[b]){
                size_t local=size_t(column)-offset;C[local/last[b]].push_back(local%last[b],value);}
            auto residual=sparse_mat_mul(C,A,F,&opt->pool);require(!residual.nnz(),"Known native amplitude violates empirical double-bonus relations: block="+irr::specs()[b].name);
            auto reduced=sparse_mat_mul(C,Qmap,F,&opt->pool);
            require(nmhv_sew_matrices_equal(sparse_mat_mul(reduced,K,F,&opt->pool),C),"Bonus control coordinate round-trip mismatch");
            visited+=first[b]*A.ncol;
        }offset+=first[b]*last[b];reduced_columns+=first[b]*restriction.reduced[b];
    }
    std::cout<<"BONUS_NATIVE_CHECK h="<<h<<" k="<<k<<" prime="<<prime<<" rows="<<cache.nrow<<" columns="<<cache.ncol
        <<" reduced_columns="<<reduced_columns<<" original_constraint_scalar_rows="<<visited<<" residual=0 copy_roundtrip=PASS"<<std::endl;
    return {{"h",std::to_string(h)},{"k",std::to_string(k)},{"prime",std::to_string(prime)},{"native_root",fs::absolute(native).generic_string()},
        {"source_input_sha256",hex_sha256_file(native/"input_identity.txt")},{"source_cache_sha256",hex_sha256_file(cache_file(native,prime))},
        {"native_rows",std::to_string(cache.nrow)},{"original_columns",std::to_string(cache.ncol)},{"reduced_columns",std::to_string(reduced_columns)},
        {"all_affine_rows_checked","true"},{"double_bonus_residual","zero"},{"coordinate_roundtrip","passed"}};
}
} // namespace mhv_double_bonus
