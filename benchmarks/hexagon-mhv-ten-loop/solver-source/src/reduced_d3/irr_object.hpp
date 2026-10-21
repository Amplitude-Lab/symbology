#pragma once
#include "irr_conversion.hpp"
namespace irr_io {
namespace fs=std::filesystem;
using irr::require;
std::string quote(const std::string& s) {
	std::string out="\"";for(unsigned char c:s){if(c=='\\'||c=='"'){out+='\\';out+=c;}else if(c=='\n')out+="\\n";else if(c=='\r')out+="\\r";else if(c=='\t')out+="\\t";else{require(c>=32,"Control character in text");out+=c;}}return out+'"';
}
size_t number(const std::string& s) {
	require(!s.empty()&&std::all_of(s.begin(),s.end(),[](char c){return c>='0'&&c<='9';}),"Expected a nonnegative integer");return std::stoull(s);
}

std::string stem(size_t h,size_t k,size_t b){return "irrSEW_"+std::to_string(h)+"p"+std::to_string(k)+"_"+irr::specs()[b].name;}


irr::sewing_object_t read_object(const fs::path& root,size_t h,size_t k,const fs::path& fec,const fs::path& lec,thread_pool* pool) {
	irr::sewing_object_t out;out.h=h;out.k=k;out.first=irr::chain_copies(fec,h,true);out.last=irr::chain_copies(lec,k,false);
	if(fs::exists(root/"object.tsv")) {
		std::ifstream meta(root/"object.tsv");std::string key,value;std::map<std::string,std::string> fields;
		while(meta>>key>>value)fields[key]=value;
		require(fields.count("kind")&&fields.count("h")&&fields.count("k")&&fields.count("directions"),"Incomplete sewing object metadata");
		out.kind=fields.at("kind");require(out.kind=="homogeneous"||out.kind=="affine"||out.kind=="unique","Unknown sewing object kind");
		require(number(fields.at("h"))==h&&number(fields.at("k"))==k,"Sewing metadata weight mismatch");out.directions=number(fields.at("directions"));
	}
	std::ifstream columns(root/"columns.tsv");require(bool(columns),"Missing sewing columns.tsv");std::string line;std::getline(columns,line);size_t offset=0;
	for(size_t b=0;b<6;++b) {
		std::string block,file;size_t first,last,base;columns>>block>>first>>last>>base>>file;
		require(bool(columns)&&block==irr::specs()[b].name&&first==out.first[b]&&last==out.last[b]&&base==offset,"Sewing block metadata mismatch");offset+=first*last;
		require(file==stem(h,k,b)+".wxf"||(out.kind=="unique"&&file==stem(h,k,b)+"_origin.wxf"),"Unsupported reduced sewing filename");
		auto t=irr::read_tensor(root/file,pool);
		if(out.kind=="unique") {
			require(out.directions==0&&t.rank()==2&&t.dim(0)==first&&t.dim(1)==last,"Invalid unique sewing block");
			sparse_tensor<rat_t,int32_t,SPARSE_COO> c(std::move(t));c.reshape({1,first*last});out.S[b]=c.to_sparse_mat(pool);
		}else {
			require(t.rank()==3&&t.dim(1)==first&&t.dim(2)==last,"Invalid homogeneous sewing block");
			if(b==0&&!fs::exists(root/"object.tsv"))out.directions=t.dim(0);
			require(t.dim(0)==out.directions,"Solution axis differs between blocks");out.S[b]=joint_eet_flatten_tensor(t,pool);
			if(out.kind=="affine") {
				auto origin=irr::read_tensor(root/(stem(h,k,b)+"_origin.wxf"),pool);
				require(origin.rank()==2&&origin.dim(0)==first&&origin.dim(1)==last,"Invalid affine origin block");
				sparse_tensor<rat_t,int32_t,SPARSE_COO> c(std::move(origin));c.reshape({1,first*last});out.S[b].append(c.to_sparse_mat(pool));
			}
		}
	}return out;
}


void write_tensor(const fs::path& path,sparse_tensor<rat_t,int32_t,SPARSE_CSR> t,std::ostream* stats=nullptr) {
	auto actual=irr::write_tensor(path,t,stats!=nullptr);
	if(stats) {
		require(actual==path,"Release tensor has a zero-dimensional shape: "+path.string());
		auto reload=irr::read_tensor(actual,nullptr);
		require(reload.dims()==t.dims()&&reload.nnz()==t.nnz()&&reload.rowptr()==t.rowptr(),"Release WXF roundtrip shape mismatch: "+path.string());
		for(size_t p=0;p<t.nnz();++p) {
			for(size_t i=1;i<t.rank();++i)require(reload.index(p)[i-1]==t.index(p)[i-1],"Release WXF roundtrip index mismatch: "+path.string());
			require(reload.val(p)==t.val(p),"Release WXF roundtrip value mismatch: "+path.string());
		}
		std::cout<<"WXF_ROUNDTRIP_PASS "<<path.string()<<std::endl;
		size_t height=0;for(size_t p=0;p<t.nnz();++p)height=std::max<size_t>(height,t.val(p).height_bits());
		*stats<<actual.filename().string()<<'\t';for(size_t i=0;i<t.rank();++i)*stats<<(i?",":"")<<t.dim(i);
		*stats<<'\t'<<t.nnz()<<'\t'<<height<<'\t'<<fs::file_size(actual)<<'\n';
	}
}
sparse_tensor<rat_t,int32_t,SPARSE_CSR> shaped(sparse_mat<rat_t,int32_t> matrix,const std::vector<size_t>& dims,thread_pool* pool) {
	sparse_tensor<rat_t,int32_t,SPARSE_COO> t(std::move(matrix),pool);t.reshape(dims);t.sort_indices(pool);return sparse_tensor<rat_t,int32_t,SPARSE_CSR>(std::move(t),pool);
}
void write_reduced(const fs::path& root,const irr::sewing_object_t& obj,thread_pool* pool) {
	std::ofstream columns(root/"columns.tsv"),meta(root/"object.tsv");columns<<"block\tfirst_copies\tlast_copies\tcolumn_offset\tfilename\n";
	meta<<"kind\t"<<obj.kind<<"\nh\t"<<obj.h<<"\nk\t"<<obj.k<<"\ndirections\t"<<obj.directions<<'\n';size_t off=0;
	for(size_t b=0;b<6;++b) {
		size_t n=obj.kind=="unique"?1:obj.directions;sparse_mat<rat_t,int32_t> S(n,obj.S[b].ncol);
		for(size_t i=0;i<n;++i)S[i]=obj.S[b][i];std::vector<size_t> dims{obj.first[b],obj.last[b]};if(obj.kind!="unique")dims.insert(dims.begin(),n);
		std::string name=stem(obj.h,obj.k,b)+(obj.kind=="unique"?"_origin.wxf":".wxf");
		write_tensor(root/name,shaped(std::move(S),dims,pool));
		if(obj.kind=="affine") {
			sparse_mat<rat_t,int32_t> O(1,obj.S[b].ncol);O[0]=obj.S[b][obj.directions];
			write_tensor(root/(stem(obj.h,obj.k,b)+"_origin.wxf"),shaped(std::move(O),{obj.first[b],obj.last[b]},pool));
		}
		columns<<irr::specs()[b].name<<'\t'<<obj.first[b]<<'\t'<<obj.last[b]<<'\t'<<off<<'\t'<<name<<'\n';off+=obj.first[b]*obj.last[b];
	}
}

// Consume packed QQ coefficients without materializing all six blocks plus a reload.
// Each tensor is checked while its exact values are still available in memory.
void write_reduced_packed(const fs::path& root,const irr::sewing_object_t& obj,
    sparse_mat<rat_t,int32_t>&& packed,thread_pool* pool) {
    require(packed.nrow==obj.directions+size_t(obj.kind!="homogeneous"),"Packed solution count mismatch");
    std::ofstream columns(root/"columns.tsv"),meta(root/"object.tsv");
    columns<<"block\tfirst_copies\tlast_copies\tcolumn_offset\tfilename\n";
    meta<<"kind\t"<<obj.kind<<"\nh\t"<<obj.h<<"\nk\t"<<obj.k<<"\ndirections\t"<<obj.directions<<'\n';
    size_t off=0;
    for(size_t b=0;b<6;++b){size_t width=obj.first[b]*obj.last[b];
        auto output=[&](const std::string& name,size_t begin,size_t count,const std::vector<size_t>& dims){
            sparse_mat<rat_t,int32_t> block(count,width);
            irr_parallel::rows(count,pool,[&](size_t i){auto& row=packed[begin+i];size_t lo=0,hi=row.nnz();
                while(lo<hi){size_t mid=(lo+hi)/2;if(size_t(row(mid))<off)lo=mid+1;else hi=mid;}
                for(size_t j=lo;j<row.nnz()&&size_t(row(j))<off+width;++j)block[i].push_back(row(j)-off,std::move(row[j]));
            });
            auto tensor=shaped(std::move(block),dims,pool);irr::write_tensor(root/name,tensor);
            auto reload=irr::read_tensor(root/name,pool);
            require(tensor.dims()==reload.dims()&&tensor.rowptr()==reload.rowptr()&&tensor.nnz()==reload.nnz(),"QQ export shape mismatch");
            for(size_t p=0;p<tensor.nnz();++p){require(tensor.val(p)==reload.val(p),"QQ export value mismatch");
                for(size_t j=1;j<tensor.rank();++j)require(tensor.index(p)[j-1]==reload.index(p)[j-1],"QQ export index mismatch");}
        };
        auto name=stem(obj.h,obj.k,b)+(obj.kind=="unique"?"_origin.wxf":".wxf");
        if(obj.kind=="unique")output(name,0,1,{obj.first[b],obj.last[b]});
        else output(name,0,obj.directions,{obj.directions,obj.first[b],obj.last[b]});
        if(obj.kind=="affine")output(stem(obj.h,obj.k,b)+"_origin.wxf",obj.directions,1,{obj.first[b],obj.last[b]});
        columns<<irr::specs()[b].name<<'\t'<<obj.first[b]<<'\t'<<obj.last[b]<<'\t'<<off<<'\t'<<name<<'\n';off+=width;
    }
    require(off==packed.ncol,"Packed column count mismatch");packed.clear();
    columns.close();meta.close();require(bool(columns)&&bool(meta),"Cannot finish QQ object metadata");
}


}
