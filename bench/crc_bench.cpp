// micro-bench: old istreambuf_iterator CRC vs new chunked-buffer CRC
#include <array>
#include <chrono>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

uint32_t crc_old(const std::filesystem::path& filename) {
	std::ifstream file(filename, std::ios::binary);
	constexpr auto crc_table = [] {
		std::array<uint32_t, 256> table{};
		for (uint32_t i = 0; i < 256; ++i) {
			uint32_t c = i;
			for (size_t j = 0; j < 8; ++j) c = (c & 1) ? (0xEDB88320 ^ (c >> 1)) : (c >> 1);
			table[i] = c;
		}
		return table;
	}();
	return std::accumulate(std::istreambuf_iterator<char>(file), std::istreambuf_iterator<char>(),
		uint32_t(0xFFFFFFFF), [&](uint32_t crc, char c) {
			return crc_table[(crc ^ c) & 0xFF] ^ (crc >> 8);
		}) ^ 0xFFFFFFFF;
}

uint32_t crc_new(const std::filesystem::path& filename) {
	std::ifstream file(filename, std::ios::binary);
	std::vector<char> buf(1 << 20);
	static const auto table = [] {
		std::array<uint32_t, 256> t{};
		for (uint32_t i = 0; i < 256; ++i) {
			uint32_t c = i;
			for (size_t j = 0; j < 8; ++j) c = (c & 1) ? (0xEDB88320 ^ (c >> 1)) : (c >> 1);
			t[i] = c;
		}
		return t;
	}();
	uint32_t crc = 0xFFFFFFFF;
	while (file) {
		file.read(buf.data(), (std::streamsize)buf.size());
		std::streamsize got = file.gcount();
		if (got <= 0) break;
		for (std::streamsize i = 0; i < got; i++)
			crc = table[(crc ^ (uint8_t)buf[i]) & 0xFF] ^ (crc >> 8);
	}
	return crc ^ 0xFFFFFFFF;
}

int main(int argc, char** argv) {
	auto path = argv[1];
	auto sz = std::filesystem::file_size(path);
	auto t0 = std::chrono::steady_clock::now();
	uint32_t a = crc_old(path);
	auto t1 = std::chrono::steady_clock::now();
	uint32_t b = crc_new(path);
	auto t2 = std::chrono::steady_clock::now();
	auto ms = [](auto d) { return std::chrono::duration<double, std::milli>(d).count(); };
	printf("size=%.1f MB  old=%.0f ms (%.0f MB/s)  new=%.1f ms (%.0f MB/s)  equal=%s\n",
		sz / 1e6, ms(t1 - t0), sz / 1e6 / (ms(t1 - t0) / 1e3), ms(t2 - t1), sz / 1e6 / (ms(t2 - t1) / 1e3),
		a == b ? "yes" : "NO");
	return a == b ? 0 : 1;
}
