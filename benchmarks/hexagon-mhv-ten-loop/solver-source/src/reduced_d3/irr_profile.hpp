#pragma once
#include "platform_resources.hpp"
#include <chrono>
#include <iostream>

namespace irr_run {
// Process peak is cumulative; phase RSS is a boundary observation, not a peak.
struct phase {
    std::string name;
    std::chrono::steady_clock::time_point start=std::chrono::steady_clock::now();
    explicit phase(std::string value):name(std::move(value)){}
    ~phase(){auto m=hex_platform::current_process_memory();
        std::cout<<"PROFILE phase="<<name<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()
            <<" rss_bytes="<<m.working_set_bytes<<" process_peak_bytes="<<m.peak_working_set_bytes<<std::endl;}
};
}
