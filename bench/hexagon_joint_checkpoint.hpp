#pragma once
#include "irr_physics_pipeline.hpp"
#include <array>
#include <cstring>

// Private intermediate checkpoints use sparse rows. Final native solve caches
// retain the package's original dense format and origin-last convention.
namespace hexagon_joint_checkpoint {
using MM = sparse_mat<ulong, int32_t>;
struct Header {
    char magic[8];
    uint64_t version, loop, h, k, rows, columns, prime, identity, nnz;
};
struct Entry { uint64_t column, value; };
static_assert(sizeof(Header) == 80 && sizeof(Entry) == 16 && sizeof(ulong) == 8);
inline fs::path file(const fs::path& root, ulong p) {
    return root / "cache" / ("p" + std::to_string(p) + ".csr");
}
inline Header header(const fs::path& root, const MM& a, ulong p, size_t h, size_t k) {
    Header out{}; std::memcpy(out.magic, "HXJCSR1", 7);
    out.version=1; out.loop=(h+k)/2; out.h=h; out.k=k;
    out.rows=a.nrow; out.columns=a.ncol; out.prime=p;
    out.identity=fingerprint(root); out.nnz=a.nnz(); return out;
}
inline Header read_header(std::ifstream& in, const fs::path& path,
                          const fs::path& root, ulong p, size_t h, size_t k) {
    Header out{}; in.read(reinterpret_cast<char*>(&out), sizeof(out));
    irr::require(bool(in) && !std::memcmp(out.magic,"HXJCSR1\0",8)
        && out.version==1 && out.loop==(h+k)/2 && out.h==h && out.k==k
        && out.prime==p && out.identity==fingerprint(root), "Sparse checkpoint identity mismatch");
    irr::require(out.rows && out.rows<=INT32_MAX && out.columns<=INT32_MAX
        && out.nnz<=out.rows*out.columns, "Sparse checkpoint dimensions");
    irr::require(fs::file_size(path)==sizeof(Header)+out.rows*sizeof(uint64_t)
        +out.nnz*sizeof(Entry), "Sparse checkpoint length mismatch");
    return out;
}
inline void save(const fs::path& root, const MM& a, ulong p, size_t h, size_t k) {
    irr_run::phase timer("sparse_checkpoint_write_verify");
    fs::create_directories(root/"cache"); auto path=file(root,p);
    irr::require(!fs::exists(path), "Do not overwrite a completed sparse checkpoint");
    auto temporary=fs::path(path.string()+".tmp"); auto expected=header(root,a,p,h,k);
    std::array<Entry,2048> buffer;
    {
        std::ofstream out(temporary,std::ios::binary);
        out.write(reinterpret_cast<const char*>(&expected),sizeof(expected));
        for(const auto& row:a.rows) {
            uint64_t count=row.nnz(); out.write(reinterpret_cast<const char*>(&count),sizeof(count));
            size_t used=0; int32_t previous=-1;
            for(auto [column,value]:row) {
                irr::require(column>previous && size_t(column)<a.ncol && value && value<p,
                    "Noncanonical sparse checkpoint row"); previous=column;
                buffer[used++]={uint64_t(column),uint64_t(value)};
                if(used==buffer.size()) { out.write(reinterpret_cast<const char*>(buffer.data()),used*sizeof(Entry)); used=0; }
            }
            if(used)out.write(reinterpret_cast<const char*>(buffer.data()),used*sizeof(Entry));
        }
        out.flush(); irr::require(bool(out),"Cannot write sparse checkpoint");
    }
    // Exact streaming readback against the live matrix avoids a second large
    // matrix allocation while checking every index, residue and affine row.
    {
        std::ifstream in(temporary,std::ios::binary); auto actual=read_header(in,temporary,root,p,h,k);
        irr::require(!std::memcmp(&actual,&expected,sizeof(actual)),"Sparse checkpoint header readback");
        for(const auto& row:a.rows) {
            uint64_t count=0; in.read(reinterpret_cast<char*>(&count),sizeof(count));
            irr::require(count==row.nnz(),"Sparse checkpoint row readback"); auto it=row.begin();
            for(uint64_t done=0;done<count;) {
                size_t n=std::min<uint64_t>(buffer.size(),count-done);
                in.read(reinterpret_cast<char*>(buffer.data()),n*sizeof(Entry));
                for(size_t i=0;i<n;++i,++it) { auto [column,value]=*it;
                    irr::require(buffer[i].column==uint64_t(column) && buffer[i].value==uint64_t(value),
                        "Sparse checkpoint coefficient readback"); }
                done+=n;
            }
        }
        irr::require(bool(in),"Sparse checkpoint readback truncated");
    }
    auto digest=hex_sha256_file(temporary);
    fs::rename(temporary,path); irr_stage::fixed_text(path.string()+".sha256",digest+"\n");
    std::cout<<"SPARSE_CHECKPOINT rows="<<a.nrow<<" columns="<<a.ncol<<" nnz="<<a.nnz()
        <<" bytes="<<fs::file_size(path)<<" exact_readback=1 origin_row=last"<<std::endl;
}
inline MM load(const fs::path& root, ulong p, size_t h, size_t k) {
    irr_run::phase timer("sparse_checkpoint_read"); auto path=file(root,p);
    irr::require(irr_stage::text(path.string()+".sha256")==hex_sha256_file(path)+"\n",
        "Damaged sparse checkpoint");
    std::ifstream in(path,std::ios::binary); auto info=read_header(in,path,root,p,h,k);
    MM a(info.rows,info.columns); std::array<Entry,2048> buffer; uint64_t total=0;
    for(size_t r=0;r<a.nrow;++r) {
        uint64_t count=0; in.read(reinterpret_cast<char*>(&count),sizeof(count));
        irr::require(bool(in) && count<=a.ncol && count<=info.nnz-total,"Sparse checkpoint row length");
        a[r].reserve(count); int64_t previous=-1;
        for(uint64_t done=0;done<count;) {
            size_t n=std::min<uint64_t>(buffer.size(),count-done);
            in.read(reinterpret_cast<char*>(buffer.data()),n*sizeof(Entry));
            irr::require(bool(in),"Truncated sparse checkpoint row");
            for(size_t i=0;i<n;++i) { auto entry=buffer[i];
                irr::require(entry.column<info.columns && int64_t(entry.column)>previous
                    && entry.value && entry.value<p,"Invalid sparse checkpoint coefficient");
                a[r].push_back(int32_t(entry.column),ulong(entry.value)); previous=entry.column;
            }
            done+=n;
        }
        total+=count;
    }
    irr::require(bool(in) && total==info.nnz,"Sparse checkpoint nonzero count"); return a;
}
}
