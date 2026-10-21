#pragma once
#include <cstdint>
#include <stdexcept>

namespace workload {
// Exact delayed reduction. If unsigned addition wraps, the mathematical sum
// is represented by wrapped + 2^128. Reduce that exceptional case immediately;
// otherwise retain the full sum until its coefficient is actually needed.
// This works for any nonzero 64-bit modulus. Dispatch can separately restrict
// it to fields where infrequent overflow makes it profitable.
struct WideModularAccumulator {
    using Wide=__uint128_t;
    uint64_t modulus,wrap_residue;
    explicit WideModularAccumulator(uint64_t p):modulus(p){
        if(!p)throw std::runtime_error("Zero accumulator modulus");
        Wide half=(Wide(1)<<127)%p;wrap_residue=uint64_t((half*2)%p);
    }
    uint64_t reduce(Wide value)const{return uint64_t(value%modulus);}
    bool add(Wide& value,Wide product)const{
        Wide previous=value;value+=product;
        if(value>=previous)return false;
        uint64_t r=reduce(value);
        value=r>=modulus-wrap_residue?r-(modulus-wrap_residue):r+wrap_residue;
        return true;
    }
};
}
