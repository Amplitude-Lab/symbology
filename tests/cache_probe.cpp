#include "native_cache.hpp"
#include <iostream>
int main(int argc, char** argv) {
    for(int i=1;i<argc;++i) std::cout<<native_cache::digest(argv[i])<<'\n';
}
