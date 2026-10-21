#include "nmhv_ope_germ.hpp"
#include <iostream>
int main(int argc,char** argv)try{
    if(argc!=2)throw std::runtime_error("test_nmhv_ope_rows DATA");
    // istringstream and binary input preserve CR on Windows too.
    const std::string row="9\t1\t-1\t11\t0\t0\t1\t2\t\"check\"";
    std::istringstream in(row+"\r\n");std::string line;std::getline(in,line);
    if(nmhv_ope_row(line)!=nmhv_ope_row(row))throw std::runtime_error("CRLF differs from LF");
    for(const std::string bad:{"typo","","check extra"}){
        bool rejected=false;try{nmhv_ope_row("9\t1\t-1\t11\t0\t0\t1\t2\t"+bad);}catch(const std::exception&){rejected=true;}
        if(!rejected)throw std::runtime_error("Unknown role was accepted");
    }
    size_t solve=0,check=0;
    std::ifstream data(std::filesystem::path(argv[1])/"ope_nmhv/E1111_L9.tsv",std::ios::binary);
    if(!std::getline(data,line))throw std::runtime_error("Missing L9 RHS");
    while(std::getline(data,line))if(!line.empty()){
        auto f=nmhv_ope_row(line);if(f.size()!=9)throw std::runtime_error("Wrong columns");
        if(f[8]=="solve")++solve;else if(f[8]=="check")++check;else throw std::runtime_error("Unrecognized role");
    }
    if(solve!=3420||check!=1140)throw std::runtime_error("Wrong L9 solve/check counts");
    std::cout<<"NMHV OPE ROWS PASS solve="<<solve<<" check="<<check<<'\n';
}catch(const std::exception& e){std::cerr<<e.what()<<'\n';return 1;}
