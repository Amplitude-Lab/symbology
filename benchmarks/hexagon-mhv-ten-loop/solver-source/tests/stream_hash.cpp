#include "hex_hash.hpp"
#include <iostream>
int main(int argc,char** argv) {
    if(argc!=2)return 2;
    std::ifstream in(argv[1],std::ios::binary);
    std::cout<<hex_sha256_stream(in)<<'\n';
    return 0;
}
