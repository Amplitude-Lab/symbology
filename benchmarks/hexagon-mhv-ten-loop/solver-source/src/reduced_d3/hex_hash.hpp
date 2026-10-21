#pragma once

#include <array>
#include <algorithm>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

inline uint32_t hex_sha256_rotate_right(uint32_t value, unsigned int shift) {
	return (value >> shift) | (value << (32 - shift));
}

struct hex_sha256_state {
	std::array<uint32_t,8> hash{0x6a09e667,0xbb67ae85,0x3c6ef372,0xa54ff53a,0x510e527f,0x9b05688c,0x1f83d9ab,0x5be0cd19};
	std::array<unsigned char,64> pending{};
	uint64_t length=0;size_t used=0;
	void block(const unsigned char* message) {
	static constexpr std::array<uint32_t, 64> constants{
		0x428a2f98, 0x71374491, 0xb5c0fbcf, 0xe9b5dba5,
		0x3956c25b, 0x59f111f1, 0x923f82a4, 0xab1c5ed5,
		0xd807aa98, 0x12835b01, 0x243185be, 0x550c7dc3,
		0x72be5d74, 0x80deb1fe, 0x9bdc06a7, 0xc19bf174,
		0xe49b69c1, 0xefbe4786, 0x0fc19dc6, 0x240ca1cc,
		0x2de92c6f, 0x4a7484aa, 0x5cb0a9dc, 0x76f988da,
		0x983e5152, 0xa831c66d, 0xb00327c8, 0xbf597fc7,
		0xc6e00bf3, 0xd5a79147, 0x06ca6351, 0x14292967,
		0x27b70a85, 0x2e1b2138, 0x4d2c6dfc, 0x53380d13,
		0x650a7354, 0x766a0abb, 0x81c2c92e, 0x92722c85,
		0xa2bfe8a1, 0xa81a664b, 0xc24b8b70, 0xc76c51a3,
		0xd192e819, 0xd6990624, 0xf40e3585, 0x106aa070,
		0x19a4c116, 0x1e376c08, 0x2748774c, 0x34b0bcb5,
		0x391c0cb3, 0x4ed8aa4a, 0x5b9cca4f, 0x682e6ff3,
		0x748f82ee, 0x78a5636f, 0x84c87814, 0x8cc70208,
		0x90befffa, 0xa4506ceb, 0xbef9a3f7, 0xc67178f2
	};
		std::array<uint32_t, 64> words{};
		for (size_t i = 0; i < 16; ++i) {
			const size_t offset = 4 * i;
			words[i] = static_cast<uint32_t>(message[offset]) << 24
				| static_cast<uint32_t>(message[offset + 1]) << 16
				| static_cast<uint32_t>(message[offset + 2]) << 8
				| static_cast<uint32_t>(message[offset + 3]);
		}
		for (size_t i = 16; i < words.size(); ++i) {
			const uint32_t s0 = hex_sha256_rotate_right(words[i - 15], 7)
				^ hex_sha256_rotate_right(words[i - 15], 18)
				^ (words[i - 15] >> 3);
			const uint32_t s1 = hex_sha256_rotate_right(words[i - 2], 17)
				^ hex_sha256_rotate_right(words[i - 2], 19)
				^ (words[i - 2] >> 10);
			words[i] = words[i - 16] + s0 + words[i - 7] + s1;
		}
		uint32_t a = hash[0], b = hash[1], c = hash[2], d = hash[3];
		uint32_t e = hash[4], f = hash[5], g = hash[6], h = hash[7];
		for (size_t i = 0; i < words.size(); ++i) {
			const uint32_t sum1 = hex_sha256_rotate_right(e, 6)
				^ hex_sha256_rotate_right(e, 11) ^ hex_sha256_rotate_right(e, 25);
			const uint32_t choice = (e & f) ^ (~e & g);
			const uint32_t temporary1 = h + sum1 + choice + constants[i] + words[i];
			const uint32_t sum0 = hex_sha256_rotate_right(a, 2)
				^ hex_sha256_rotate_right(a, 13) ^ hex_sha256_rotate_right(a, 22);
			const uint32_t majority = (a & b) ^ (a & c) ^ (b & c);
			const uint32_t temporary2 = sum0 + majority;
			h = g; g = f; f = e; e = d + temporary1;
			d = c; c = b; b = a; a = temporary1 + temporary2;
		}
		hash[0] += a; hash[1] += b; hash[2] += c; hash[3] += d;
		hash[4] += e; hash[5] += f; hash[6] += g; hash[7] += h;
	}
	void update(const unsigned char* data,size_t size) {
		length+=size;
		if(used){auto n=std::min(size,64-used);std::copy_n(data,n,pending.data()+used);used+=n;data+=n;size-=n;if(used==64){block(pending.data());used=0;}}
		while(size>=64){block(data);data+=64;size-=64;}
		if(size){std::copy_n(data,size,pending.data());used=size;}
	}
	std::string finish() {
		const uint64_t bits=length*8;pending[used++]=0x80;
		if(used>56){std::fill(pending.begin()+used,pending.end(),0);block(pending.data());used=0;}
		std::fill(pending.begin()+used,pending.begin()+56,0);
		for(size_t i=0;i<8;++i)pending[56+i]=static_cast<unsigned char>(bits>>(56-8*i));
		block(pending.data());
	std::ostringstream result;
	result << std::hex << std::setfill('0');
	for (const auto word : hash)
		result << std::setw(8) << word;
	return result.str();
	}
};

inline std::string hex_sha256(const std::vector<unsigned char>& input) {
	hex_sha256_state state;state.update(input.data(),input.size());return state.finish();
}

inline std::string hex_sha256_stream(std::istream& in) {
	hex_sha256_state state;std::array<unsigned char,65536> buffer;
	for(;;){in.read(reinterpret_cast<char*>(buffer.data()),buffer.size());state.update(buffer.data(),static_cast<size_t>(in.gcount()));
		if(in.bad()||(!in.eof()&&in.fail()))throw std::runtime_error("SHA256 read failed");if(in.eof())break;}
	return state.finish();
}

inline std::vector<unsigned char> hex_read_binary_file(
	const std::filesystem::path& path) {
	std::ifstream stream(path, std::ios::binary);
	if (!stream)
		throw std::runtime_error("Cannot open file for hashing: " + path.string());
	stream.seekg(0, std::ios::end);
	const auto length = stream.tellg();
	if (length < 0)
		throw std::runtime_error("Cannot determine file size: " + path.string());
	stream.seekg(0, std::ios::beg);
	std::vector<unsigned char> bytes(static_cast<size_t>(length));
	if (!bytes.empty())
		stream.read(reinterpret_cast<char*>(bytes.data()), length);
	if (!stream)
		throw std::runtime_error("Cannot read file for hashing: " + path.string());
	return bytes;
}

inline std::string hex_sha256_file(const std::filesystem::path& path) {
	std::ifstream in(path,std::ios::binary);if(!in)throw std::runtime_error("Cannot open file for hashing: "+path.string());
	return hex_sha256_stream(in);
}
