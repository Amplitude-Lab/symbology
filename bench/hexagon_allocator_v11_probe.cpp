#include "hexagon_mapped_buffers_v11.hpp"
#include <iostream>
#include <sys/stat.h>

void require(bool value,const char* message){if(!value)throw std::runtime_error(message);}
size_t allocated_bytes(){
    struct stat info{};const char* path=std::getenv("HEXAGON_MAPPED_FILE");
    require(path&&::stat(path,&info)==0,"Cannot inspect probe scratch");return size_t(info.st_blocks)*512;
}
void check(const unsigned char* bytes,size_t count,unsigned char value){
    for(size_t i=0;i<count;++i)require(bytes[i]==value,"Live coefficient bytes changed");
}
int main()try{
    auto& arena=hexagon_mapped_buffers::pool();require(arena.enabled(),"Probe requires mapped storage");
    constexpr size_t full=4*1024*1024+600,prefix=8192+17,edge=6000;
    auto* left=static_cast<unsigned char*>(arena.allocate(edge));
    auto* middle=static_cast<unsigned char*>(arena.allocate(full));
    auto* right=static_cast<unsigned char*>(arena.allocate(edge));
    std::memset(left,0xa5,edge);std::memset(middle,0x3c,full);std::memset(right,0x5a,edge);
    arena.maybe_evict(true);
    size_t before=allocated_bytes();auto initial=arena.statistics();
    auto* shrunk=static_cast<unsigned char*>(arena.reallocate(middle,prefix));
    require(shrunk==middle,"Shrinking block moved unexpectedly");
    check(left,edge,0xa5);check(middle,prefix,0x3c);check(right,edge,0x5a);
    auto reduced=arena.statistics();
    require(initial.live_capacity_bytes-reduced.live_capacity_bytes>3*1024*1024,"Shrinking did not return unused capacity");
    if(reduced.reclaim_enabled){require(reduced.reclaim_calls>0,"Freed range was not punched");require(allocated_bytes()<before,"Freed file extents did not shrink");}
    auto* reused=static_cast<unsigned char*>(arena.allocate(2*1024*1024));
    require(reinterpret_cast<uintptr_t>(reused)>reinterpret_cast<uintptr_t>(middle)&&
            reinterpret_cast<uintptr_t>(reused)<reinterpret_cast<uintptr_t>(right),"Returned tail was not reused");
    std::memset(reused,0x7e,2*1024*1024);
    arena.maybe_evict(true);check(left,edge,0xa5);check(middle,prefix,0x3c);check(right,edge,0x5a);check(reused,2*1024*1024,0x7e);
    arena.release(reused);arena.release(middle);check(left,edge,0xa5);check(right,edge,0x5a);
    // The freed middle and its returned tail must coalesce for a request
    // larger than either of the fragments remaining after tail reuse.
    auto* joined=static_cast<unsigned char*>(arena.allocate(full));
    require(joined==middle,"Adjacent free ranges did not coalesce");
    std::memset(joined,0x19,full);check(left,edge,0xa5);check(right,edge,0x5a);
    arena.release(joined);arena.release(right);arena.release(left);
    auto end=arena.statistics();require(end.live_capacity_bytes==0&&end.live_requested_bytes==0&&end.used_address_bytes==0,"Empty arena retained live allocation state");
    std::cout<<"ALLOCATOR_V11_PROBE PASS prefix_preserved=1 neighbours_preserved=1 tail_reused=1 free_coalesced=1 arena_suffix_trimmed=1 reclaim_calls="
             <<end.reclaim_calls<<" reclaim_requested_bytes="<<end.reclaim_requested_bytes<<" file_allocated_bytes="<<allocated_bytes()<<std::endl;
}catch(const std::exception& error){std::cerr<<error.what()<<std::endl;return 1;}
