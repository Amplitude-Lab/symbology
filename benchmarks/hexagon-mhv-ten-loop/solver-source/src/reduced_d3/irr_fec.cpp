#include "irr_recurrence.hpp"
#include "joint_eet_io.hpp"
#include "irr_runtime.hpp"
#include "k_crt_cache.hpp"
#include "irr_stage.hpp"
#include "irr_block_cache.hpp"

namespace fs=std::filesystem;
const field_t QQ(FIELD_QQ);

int main(int argc,char** argv){try {
	if(argc<3)throw std::runtime_error("irr_fec MAX_WEIGHT OUTPUT [--condition WXF] [--condition-blocks TSV] [--clean-block-cache]");
	size_t maximum=std::stoul(argv[1]);if(!maximum)throw std::runtime_error("MAX_WEIGHT must be positive");
	fs::path root=argv[2],condition="data/dlogmat_ca_1term_pd3.wxf",blocks="data/dlogmat_ca_1term_pd3_blocks.tsv";
	bool clean_blocks=false;
	for(int i=3;i<argc;++i){std::string key=argv[i];if(key=="--clean-block-cache"){clean_blocks=true;continue;}if(i+1==argc)throw std::runtime_error("Missing option value");
		if(key=="--condition")condition=argv[++i];else if(key=="--condition-blocks")blocks=argv[++i];else throw std::runtime_error("Unknown option "+key);}
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=true;irr_run::threads(opt);
	if(clean_blocks)irr_stage::no_links(root);
	irr_run::bind(root,"FEC primitive method1",{condition,blocks});
	if(clean_blocks)irr_stage::no_links(root/".active.lock");
	irr_stage::lock active(root/".active.lock");
	auto* pool=&opt->pool;irr::channels_t cg(pool);
	std::ofstream terms(root/"CG.tsv");terms<<"target\tinput\tchannel\ttarget_component\tinput_component\tletter\tcoefficient\n";
	for(size_t b=0;b<6;++b)for(size_t a=0;a<6;++a)for(size_t g:cg.channels[b][a])
		for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)for(const auto& term:cg.cg[b].vectors[g].component[s])
			terms<<irr::specs()[b].name<<'\t'<<irr::specs()[a].name<<'\t'<<g<<'\t'<<s<<'\t'<<term.coordinate/9-cg.axis[a].offset<<'\t'<<term.coordinate%9<<'\t'<<term.coefficient<<'\n';
	terms.close();
	auto W=irr::read_tensor(condition,pool);auto cb=nmhv_sew_read_irrep_axis(blocks);
	std::array<irr::template_t,6> templates;for(size_t b=0;b<6;++b)templates[b]=irr::make_template(cg,b,W,cb,opt);
	irr::copies older{},previous{1,0,0,0,0,0};irr::matrices K;
	std::ofstream summary(root/"summary.tsv");summary<<"weight\tdimension\tirr_nnz\theight_bits\n";
	for(size_t w=1;w<=maximum;++w) {
		auto dir=root/joint_eet_weight_name(w);fs::create_directories(dir);irr::matrices next;bool legacy_verified=false;
		for(size_t b=0;b<6;++b) {
			auto file=dir/("irrFEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf");
			if(clean_blocks){irr_stage::no_links(file);irr_stage::no_links(fs::path(file.string()+".ready"));irr_stage::no_links(dir/"cache"/irr::specs()[b].name);}
			if(irr_run::ready(file)) {
				next[b]=joint_eet_flatten_tensor(irr::read_tensor(file,pool),pool);
				if(clean_blocks)irr_block_cache::clean(dir,irr::specs()[b].name,file,active,legacy_verified);
				std::cout<<"REUSE weight="<<w<<" block="<<irr::specs()[b].name<<std::endl;continue;
			}
			size_t nr=0,nc=0,nnz=0;
			if(w==1) {
				auto seed=irr::skeleton(previous,pool);auto c=make_hhex_candidates(fec_symmetry_t::parity_dihedral,seed,irr::specs()[b].name);
				next[b]=sparse_mat<rat_t,int32_t>(b==0||b==2?1:0,c.size());if(next[b].nrow)next[b][0].push_back(0,rat_t(1));
			}else {
				auto matrix=irr::assemble(cg,templates[b],b,cb,older,previous,K,opt);
				nr=matrix.nrow;nc=matrix.ncol;nnz=matrix.nnz();
				std::cout<<"IRR_FEC weight="<<w<<" block="<<irr::specs()[b].name<<" rows="<<nr<<" cols="<<nc<<" nnz="<<nnz<<std::endl;
				if(!nr||!nc)next[b]=irr::solve(matrix,opt);
				else {
					matrix.sort_rows_by_nnz();auto pivots=k_crt::rref(matrix,opt,dir/"cache"/irr::specs()[b].name);
					size_t rank=0;for(const auto& group:pivots)rank+=group.size();
					if(rank==nc)next[b]=sparse_mat<rat_t,int32_t>(0,nc);
					else {next[b]=sparse_mat_rref_kernel(matrix,pivots,QQ,opt);next[b].transpose_replace();}
					matrix.clear();pool->detach_loop(0,next[b].nrow,[&](size_t r){next[b][r].canonicalize();fec_irrep_cancel_divisor(next[b][r]);});pool->wait();next[b].sort_rows_by_nnz();
				}
			}
			sparse_tensor<rat_t,int32_t,SPARSE_CSR> tensor(next[b],pool);irr::write_tensor(file,tensor,false,clean_blocks);
			auto stats=joint_eet_rational_stats(next[b]);std::ofstream report(dir/(irr::specs()[b].name+".tsv"));
			report<<"rows\tcolumns\tmatrix_nnz\tcopies\tirr_nnz\theight_bits\tpeak_bytes\n"<<nr<<'\t'<<nc<<'\t'<<nnz<<'\t'<<next[b].nrow<<'\t'<<next[b].nnz()<<'\t'<<stats.height_bits<<'\t'<<hex_platform::current_process_memory().peak_working_set_bytes<<'\n';
			report.close();if(!report)throw std::runtime_error("Cannot write block report");
			irr_run::mark(file);
			if(clean_blocks)irr_block_cache::clean(dir,irr::specs()[b].name,file,active,legacy_verified);
		}
		irr::copies current{};size_t dim=0,nnz=0,height=0;std::ofstream mult(dir/"multiplicities.tsv"),hashes(dir/"checksums.sha256"),cols(dir/"columns.tsv");
		cols<<"target\tinput\told_copies\tcolumn_offset\tchannels_per_copy\n";
		for(size_t b=0;b<6;++b) {
			current[b]=next[b].nrow;dim+=current[b]*irr::specs()[b].irrep_dimension;nnz+=next[b].nnz();height=std::max(height,joint_eet_rational_stats(next[b]).height_bits);
			mult<<irr::specs()[b].name<<'\t'<<current[b]<<'\n';size_t offset=0;
			for(size_t a=0;a<6;++a){cols<<irr::specs()[b].name<<'\t'<<irr::specs()[a].name<<'\t'<<previous[a]<<'\t'<<offset<<'\t'<<cg.channels[b][a].size()<<'\n';offset+=previous[a]*cg.channels[b][a].size();}
			if(next[b].ncol!=offset)throw std::runtime_error("FEC channel dimensions mismatch");
			auto file=irr::tensor_path(dir/("irrFEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf"));hashes<<hex_sha256_file(file)<<"  "<<file.filename().string()<<'\n';
		}
		mult.close();cols.close();hashes<<hex_sha256_file(dir/"multiplicities.tsv")<<"  multiplicities.tsv\n"<<hex_sha256_file(dir/"columns.tsv")<<"  columns.tsv\n";hashes.close();
		std::ofstream manifest(dir/"manifest.json");manifest<<"{\"schema\":\"hexagon-irr-fec-v1\",\"complete\":true,\"weight\":"<<w<<",\"dimension\":"<<dim<<",\"irr_nnz\":"<<nnz<<",\"height_bits\":"<<height<<",\"method\":1,\"checksums_sha256\":\""<<hex_sha256_file(dir/"checksums.sha256")<<"\"}\n";
		summary<<w<<'\t'<<dim<<'\t'<<nnz<<'\t'<<height<<std::endl;
		manifest.close();irr_stage::complete_chain(dir);
		older=previous;previous=current;K=std::move(next);std::cout<<"COMPLETE FEC weight="<<w<<" dimension="<<dim<<" nnz="<<nnz<<" height_bits="<<height<<std::endl;
	}
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}}
