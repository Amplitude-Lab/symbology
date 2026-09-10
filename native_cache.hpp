#pragma once
// Content provenance for automatically generated intermediates. A legacy file
// without a matching receipt is recomputed; file existence is never proof.
#include <array>
#include <bit>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace native_cache {
class sha256 {
    std::array<uint32_t, 8> h{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
    std::array<unsigned char, 64> buf{};
    uint64_t length = 0;
    size_t used = 0;
    void block() {
        static constexpr uint32_t k[64] = {
            0x428a2f98,0x71374491,0xb5c0fbcf,0xe9b5dba5,0x3956c25b,0x59f111f1,0x923f82a4,0xab1c5ed5,
            0xd807aa98,0x12835b01,0x243185be,0x550c7dc3,0x72be5d74,0x80deb1fe,0x9bdc06a7,0xc19bf174,
            0xe49b69c1,0xefbe4786,0x0fc19dc6,0x240ca1cc,0x2de92c6f,0x4a7484aa,0x5cb0a9dc,0x76f988da,
            0x983e5152,0xa831c66d,0xb00327c8,0xbf597fc7,0xc6e00bf3,0xd5a79147,0x06ca6351,0x14292967,
            0x27b70a85,0x2e1b2138,0x4d2c6dfc,0x53380d13,0x650a7354,0x766a0abb,0x81c2c92e,0x92722c85,
            0xa2bfe8a1,0xa81a664b,0xc24b8b70,0xc76c51a3,0xd192e819,0xd6990624,0xf40e3585,0x106aa070,
            0x19a4c116,0x1e376c08,0x2748774c,0x34b0bcb5,0x391c0cb3,0x4ed8aa4a,0x5b9cca4f,0x682e6ff3,
            0x748f82ee,0x78a5636f,0x84c87814,0x8cc70208,0x90befffa,0xa4506ceb,0xbef9a3f7,0xc67178f2};
        uint32_t w[64];
        for (size_t i=0; i<16; ++i) w[i]=uint32_t(buf[4*i])<<24|uint32_t(buf[4*i+1])<<16|uint32_t(buf[4*i+2])<<8|buf[4*i+3];
        for (size_t i=16; i<64; ++i) {
            auto x=w[i-15], y=w[i-2];
            w[i]=w[i-16]+(std::rotr(x,7)^std::rotr(x,18)^(x>>3))+w[i-7]+(std::rotr(y,17)^std::rotr(y,19)^(y>>10));
        }
        auto a=h[0],b=h[1],c=h[2],d=h[3],e=h[4],f=h[5],g=h[6],v=h[7];
        for (size_t i=0; i<64; ++i) {
            auto t=v+(std::rotr(e,6)^std::rotr(e,11)^std::rotr(e,25))+((e&f)^(~e&g))+k[i]+w[i];
            auto u=(std::rotr(a,2)^std::rotr(a,13)^std::rotr(a,22))+((a&b)^(a&c)^(b&c));
            v=g;g=f;f=e;e=d+t;d=c;c=b;b=a;a=t+u;
        }
        h[0]+=a;h[1]+=b;h[2]+=c;h[3]+=d;h[4]+=e;h[5]+=f;h[6]+=g;h[7]+=v;
    }
public:
    void update(const char* bytes, size_t n) {
        length += n;
        for (size_t i=0;i<n;++i) {
            buf[used++]=static_cast<unsigned char>(bytes[i]);
            if (used==64) { block();used=0; }
        }
    }
    std::string finish() {
        uint64_t bits=length*8;
        buf[used++]=0x80;
        if (used>56) { while(used<64) buf[used++]=0;block();used=0; }
        while(used<56) buf[used++]=0;
        for (int i=7;i>=0;--i) buf[used++]=static_cast<unsigned char>(bits>>(8*i));
        block();
        std::ostringstream out;out<<std::hex<<std::setfill('0');
        for(auto x:h) out<<std::setw(8)<<x;
        return out.str();
    }
};
inline std::string digest(const std::filesystem::path& p) {
    std::ifstream f(p, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot fingerprint input: " + p.string());
    sha256 hash; std::vector<char> block(1<<20);
    while(f) { f.read(block.data(), block.size());hash.update(block.data(), f.gcount()); }
    if (!f.eof()) throw std::runtime_error("Cannot finish fingerprint: " + p.string());
    return hash.finish();
}
// Resolve invocation through PATH as well as ./relative and absolute paths.
inline std::filesystem::path executable_path(const std::string& arg) {
    auto p = std::filesystem::path(arg);
    if (p.has_parent_path()) return std::filesystem::absolute(p);
    if (const char* env = std::getenv("PATH")) {
        std::istringstream parts(env); std::string dir;
        while (std::getline(parts, dir,
#ifdef _WIN32
            ';'
#else
            ':'
#endif
        )) {
            auto candidate = (dir.empty() ? std::filesystem::current_path() : std::filesystem::path(dir)) / p;
            if (std::filesystem::is_regular_file(candidate)) return std::filesystem::absolute(candidate);
        }
    }
    return std::filesystem::absolute(p);
}
// Set once by each executable; embedded library callers get a build identity.
inline std::string executable_identity = std::string(__DATE__) + " " + __TIME__;
inline void set_executable(const std::filesystem::path& p) { executable_identity = digest(p); }

class receipt {
    std::filesystem::path stamp;
    std::vector<std::filesystem::path> inputs, outputs;
    std::string recipe, before;
    bool allow_absent;
    std::string input_state() const {
        std::ostringstream out;
        out << "symbology-native-cache-v1\n" << std::quoted(executable_identity) << '\n' << std::quoted(recipe) << '\n';
        for (const auto& p:inputs) out << digest(p) << '\n';
        return out.str();
    }
    std::string output_state() const {
        std::ostringstream out;
        for (const auto& p:outputs) out << (std::filesystem::exists(p) ? digest(p) : "absent") << '\n';
        return out.str();
    }
public:
    receipt(std::filesystem::path marker, std::string tag,
            std::vector<std::filesystem::path> in, std::vector<std::filesystem::path> out, bool absent = false)
        :stamp(std::move(marker)), inputs(std::move(in)), outputs(std::move(out)), recipe(std::move(tag)), before(input_state()), allow_absent(absent) {}
    bool valid() const {
        for (const auto& p:outputs) if ((!allow_absent && !std::filesystem::exists(p)) || (std::filesystem::exists(p) && (!std::filesystem::is_regular_file(p) || !std::filesystem::file_size(p)))) return false;
        std::ifstream f(stamp, std::ios::binary);
        if (!f) return false;
        return std::string(std::istreambuf_iterator<char>(f), {}) == before + output_state();
    }
    void save() const {
        for (const auto& p:outputs) if ((!allow_absent && !std::filesystem::exists(p)) || (std::filesystem::exists(p) && (!std::filesystem::is_regular_file(p) || !std::filesystem::file_size(p)))) throw std::runtime_error("Missing native output: " + p.string());
        if(input_state()!=before) throw std::runtime_error("Inputs changed during native calculation; cache receipt not saved");
        auto tmp=stamp; tmp += ".tmp";
        { std::ofstream f(tmp, std::ios::binary); f<<before<<output_state();f.flush();
          if(!f) throw std::runtime_error("Cannot write cache receipt: " + tmp.string()); }
        std::filesystem::rename(tmp, stamp);
    }
};
}
