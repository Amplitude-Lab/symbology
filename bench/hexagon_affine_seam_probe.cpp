#include "hexagon_affine_seam_discovery_v10.hpp"
using symrep::require;
struct ToyRows {
    staged_kernel::MM rows;
    size_t ncol()const{return rows.ncol;}
    auto masks(bool remember)const{return std::vector<std::vector<unsigned char>>{std::vector<unsigned char>(remember?rows.nrow:0)};}
    template<class Emit> void discover(structured_kernel::Relations&,const std::vector<int32_t>& columns,
        const field_t& F,const staged_kernel::Selection& select,std::vector<std::vector<unsigned char>>& seen,
        Emit&& emit,thread_pool*)const{
        for(size_t r=0;r<rows.nrow;++r){auto& mask=seen[0];if(!mask.empty()&&mask[r])continue;
            if(rows[r].nnz()>select.max_terms)continue;staged_kernel::MV row;
            for(auto [j,x]:rows[r]){int32_t c=columns[j];if(c<0)continue;ulong v=x;
                if(select.substitution){auto old=c;c=select.substitution->image[old];if(c<0)continue;v=nmod_mul(v,select.substitution->scale[old],F.mod);}
                row.push_back(c,v);}
            staged_kernel::normalize_mod(row,F);if(!mask.empty())mask[r]=1;if(row.nnz())emit(std::move(row));
        }
    }
};
ulong coefficient(const staged_kernel::MV& row,int32_t j){auto p=row.find(j);return p?*p:0;}
int main()try{
    Flint::set_memory_functions();field_t F(FIELD_Fp,101);rref_option_t opt;opt->pool.reset(2);opt->method=1;opt->is_back_sub=true;opt->verbose=false;
    for(bool free_direction:{false,true}){
        size_t n=free_direction?3:2;ToyRows source{staged_kernel::MM(1,n)};source.rows[0].push_back(0,1);source.rows[0].push_back(1,99);
        std::vector<int32_t> order=free_direction?std::vector<int32_t>{2,0,1}:std::vector<int32_t>{1,0};
        staged_kernel::Options settings;settings.target_nullity=free_direction?2:1;settings.max_residual_nullity=settings.target_nullity;
        auto physical=[&](streaming_echelon::Reducer& reducer,const std::vector<int32_t>& columns,size_t round){if(round)return;
            staged_kernel::MV row;row.push_back(columns[0],1);row.push_back(int32_t(n),95);row.sort_indices();reducer.consume(std::move(row));};
        auto result=staged_kernel::discover_affine_candidate(source,F,opt,settings,physical,order);
        require(result.nrow==size_t(free_direction?2:1)&&result.ncol==n,"Wrong affine dimension");
        auto& origin=result[result.nrow-1];require(coefficient(origin,0)==6&&coefficient(origin,1)==3,"Wrong affine origin");
        for(size_t r=0;r+1<result.nrow;++r)require(!coefficient(result[r],0)&&!coefficient(result[r],1)&&coefficient(result[r],2),"Wrong homogeneous direction");
        std::cout<<"AFFINE_SEAM_PROBE free_direction="<<free_direction<<" passed=1"<<std::endl;
    }
    {
        staged_kernel::MM kernel(2,3);kernel[0].push_back(0,12);kernel[0].push_back(1,6);kernel[0].push_back(2,2);
        kernel[1].push_back(0,19);kernel[1].push_back(1,9);kernel[1].push_back(2,3);
        auto result=staged_kernel::affine_slice_augmented_kernel(std::move(kernel),F);
        require(result.nrow==2&&coefficient(result[1],0)==6&&coefficient(result[1],1)==3&&coefficient(result[0],0)==1&&!coefficient(result[0],1),"Nonzero t direction slicing failed");
        std::cout<<"AFFINE_SEAM_PROBE nonzero_t_direction=1 passed=1"<<std::endl;
    }
    {
        ToyRows source{staged_kernel::MM(1,1)};source.rows[0].push_back(0,1);staged_kernel::Options settings;settings.target_nullity=1;settings.max_residual_nullity=1;
        auto physical=[&](streaming_echelon::Reducer& reducer,const std::vector<int32_t>& columns,size_t round){if(round)return;staged_kernel::MV row;row.push_back(columns[0],1);row.push_back(1,100);reducer.consume(std::move(row));};
        bool rejected=false;try{staged_kernel::discover_affine_candidate(source,F,opt,settings,physical);}catch(const std::exception&){rejected=true;}
        require(rejected,"Inconsistent affine equations were accepted");std::cout<<"AFFINE_SEAM_PROBE inconsistent_rhs_rejected=1 passed=1"<<std::endl;
    }
}catch(const std::exception& e){std::cerr<<"AFFINE_SEAM_PROBE_ERROR "<<e.what()<<std::endl;return 1;}
