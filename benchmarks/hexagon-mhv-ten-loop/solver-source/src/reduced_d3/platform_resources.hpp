#pragma once

#include <cstdint>
#include <fstream>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#include <psapi.h>
#else
#include <sys/resource.h>
#include <unistd.h>
#endif

namespace hex_platform {

struct process_memory_t {
	std::uint64_t working_set_bytes = 0;
	std::uint64_t peak_working_set_bytes = 0;
	std::uint64_t private_commit_bytes = 0;
	std::uint64_t peak_private_commit_bytes = 0;
	bool valid = false;
};

inline process_memory_t current_process_memory() {
	process_memory_t result;
#ifdef _WIN32
	PROCESS_MEMORY_COUNTERS_EX counters{};
	counters.cb = sizeof(counters);
	if (GetProcessMemoryInfo(
		GetCurrentProcess(),
		reinterpret_cast<PROCESS_MEMORY_COUNTERS*>(&counters),
		sizeof(counters))) {
		result.working_set_bytes = counters.WorkingSetSize;
		result.peak_working_set_bytes = counters.PeakWorkingSetSize;
		result.private_commit_bytes = counters.PrivateUsage;
		result.peak_private_commit_bytes = counters.PeakPagefileUsage;
		result.valid = true;
	}
#else
	long resident_pages = 0;
	long total_pages = 0;
	std::ifstream status("/proc/self/statm");
	status >> total_pages >> resident_pages;
	const long page_size = sysconf(_SC_PAGESIZE);
	if (status && resident_pages >= 0 && page_size > 0) {
		result.working_set_bytes = static_cast<std::uint64_t>(resident_pages)
			* static_cast<std::uint64_t>(page_size);
		result.valid = true;
	}
	rusage usage{};
	if (getrusage(RUSAGE_SELF, &usage) == 0) {
		result.peak_working_set_bytes = static_cast<std::uint64_t>(usage.ru_maxrss) * 1024;
		result.valid = true;
	}
#endif
	return result;
}

}
