#pragma once
#include "symrep_io.hpp"

namespace symrep {
// Concatenate saved row partitions while retaining at most one input tensor.
// Works with the existing bounded WXF writer, including exact byte readback.
// Metadata is captured from each certified matrix before its storage is freed.
class PartitionedTensorView {
public:
    struct Part {fs::path path;std::vector<size_t> row_nonzeros;};
private:
    std::vector<Part> parts;
    std::vector<size_t> dimensions,offsets,part_offsets;
    mutable size_t loaded=SIZE_MAX;
    mutable Tensor cache;
    size_t locate(size_t p)const{
        require(p<offsets.back(),"partitioned tensor index out of range");
        size_t s=std::upper_bound(part_offsets.begin(),part_offsets.end(),p)-part_offsets.begin()-1;
        if(loaded!=s){
            cache.clear();cache=read(parts[s].path);loaded=s;
            require(cache.rank()==rank()&&cache.dim(0)==parts[s].row_nonzeros.size(),"partitioned tensor row shape mismatch");
            for(size_t j=1;j<rank();++j)require(cache.dim(j)==dim(j),"partitioned tensor inner shape mismatch");
            require(cache.nnz()==part_offsets[s+1]-part_offsets[s],"partitioned tensor nonzero count mismatch");
            for(size_t r=0;r<cache.dim(0);++r)require(cache.rowptr()[r+1]-cache.rowptr()[r]==parts[s].row_nonzeros[r],"partitioned tensor row counts changed");
        }
        return p-part_offsets[s];
    }
public:
    PartitionedTensorView(std::vector<Part> files,std::vector<size_t> inner_dimensions):parts(std::move(files)),dimensions{0},offsets{0},part_offsets{0}{
        require(!inner_dimensions.empty(),"partitioned tensor requires inner dimensions");
        dimensions.insert(dimensions.end(),inner_dimensions.begin(),inner_dimensions.end());
        for(const auto& part:parts){dimensions[0]+=part.row_nonzeros.size();for(auto n:part.row_nonzeros)offsets.push_back(offsets.back()+n);part_offsets.push_back(offsets.back());}
    }
    const auto& dims()const{return dimensions;}size_t dim(size_t i)const{return dimensions.at(i);}
    size_t rank()const{return dimensions.size();}size_t nnz()const{return offsets.back();}
    const size_t* rowptr()const{return offsets.data();}
    const Q& val(size_t p)const{auto i=locate(p);return cache.val(i);}
    const I* index(size_t p)const{auto i=locate(p);return cache.index(i);}
};
}
