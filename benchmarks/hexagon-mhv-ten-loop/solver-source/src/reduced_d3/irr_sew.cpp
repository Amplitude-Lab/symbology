#include "irr_bulk.hpp"
namespace fs=std::filesystem;
using namespace irr_bulk;

int main(int argc,char** argv){try{
	require(argc>=5,"irr_sew MHV|NMHV FIRST_WEIGHT LAST_WEIGHT OUTPUT [--fec-root DIR] [--lec-root DIR] [--condition WXF] [--condition-blocks TSV] [--prereduce left|right|none] [--verify] [--expect DIM]");
	std::string helicity=argv[1];require(helicity=="MHV"||helicity=="NMHV","Unknown helicity");
	size_t h=std::stoul(argv[2]),k=std::stoul(argv[3]);require(h>0&&k>0&&(h+k)%2==0,"Require h,k>0 and even total weight; k=0 is not implemented here");
	fs::path root=argv[4],fec="output/FEC",lec=fs::path("output/LEC")/helicity;
	fs::path condition="data/dlogmat_full_ca_pd3.wxf",blocks="data/dlogmat_full_ca_pd3_blocks.tsv";
	std::string prereduce=k==1?"right":"left";
	bool verify=false;size_t expected=SIZE_MAX;std::string solver="baseline";
	for(int i=5;i<argc;++i) {
		std::string key=argv[i];if(key=="--verify"){verify=true;continue;}require(i+1<argc,"Missing option value");std::string value=argv[++i];
		if(key=="--solver")solver=value;else if(key=="--prereduce")prereduce=value;
		else if(key=="--fec-root")fec=value;else if(key=="--lec-root")lec=value;else if(key=="--condition")condition=value;
		else if(key=="--condition-blocks")blocks=value;else if(key=="--expect")expected=std::stoul(value);else throw std::runtime_error("Unknown option "+key);
	}
	require(solver=="baseline"||solver=="staged-factored","--solver must be baseline or staged-factored");
	require(solver!="staged-factored"||prereduce=="right","factored staged sewing requires --prereduce right");
	require(prereduce=="left"||prereduce=="right"||prereduce=="none","--prereduce must be left, right or none");
	std::vector<fs::path> identity{condition,blocks,lec/"common"/"CG.tsv",lec/"LEC_0_components.wxf",lec/"LEC_1_physical.wxf"};
	identity.push_back(fec/joint_eet_weight_name(h)/"multiplicities.tsv");
	if(h>1)identity.push_back(fec/joint_eet_weight_name(h-1)/"multiplicities.tsv");
	identity.push_back(lec/joint_eet_weight_name(k)/"blocks.tsv");
	identity.push_back(k>1?lec/joint_eet_weight_name(k-1)/"blocks.tsv":lec/"LEC_0_blocks.tsv");
	for(size_t b=0;b<6;++b) {
		identity.push_back(irr::tensor_path(fec/joint_eet_weight_name(h)/("irrFEC_"+std::to_string(h)+"_"+irr::specs()[b].name+".wxf")));
		identity.push_back(irr::tensor_path(lec/joint_eet_weight_name(k)/("irrLEC_"+std::to_string(k)+"_"+irr::specs()[b].name+".wxf")));
	}
	irr_run::bind(root,helicity+" SEW "+std::to_string(h)+"p"+std::to_string(k)+(prereduce=="none"?"":" prereduce="+prereduce)+" solver="+solver,identity);
	irr_stage::lock active(root/".active.lock");
	if(irr_run::ready(root/"manifest.json")) {irr_run::verify_checksums(root);std::cout<<"COMPLETE SEW already exists: "<<root<<std::endl;return 0;}
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=true;irr_run::threads(opt);auto* pool=&opt->pool;
	auto start=std::chrono::steady_clock::now();irr::channels_t cg(pool);check_lec_convention(lec,cg,pool);
	auto terminal=nmhv_sew_read_tensor<rat_t,int32_t>(lec/"LEC_0_components.wxf",QQ,pool);
	require(terminal.rank()==2&&terminal.dim(0)==size_t(helicity=="MHV"?1:5),"LEC helicity/terminal mismatch");
	auto mf=read_copies(fec,h,true),ml=read_copies(lec,k,false),of=read_copies(fec,h-1,true),ol=read_copies(lec,k-1,false);
	fs::create_directories(root);std::ofstream inputs(root/"inputs.sha256");
	auto first=read_chain(fec,h,true,mf,inputs,pool),last=read_chain(lec,k,false,ml,inputs,pool);
	for(const auto& path:{condition,blocks,lec/"common"/"CG.tsv",lec/"LEC_1_physical.wxf",lec/"seed_scale.tsv",lec/"LEC_0_components.wxf",lec/joint_eet_weight_name(k)/"manifest.json"})
		inputs<<hex_sha256_file(path)<<"  "<<fs::absolute(path).generic_string()<<'\n';
	for(const auto& path:{fec/joint_eet_weight_name(h)/"multiplicities.tsv",lec/joint_eet_weight_name(k)/"blocks.tsv",
		k>1?lec/joint_eet_weight_name(k-1)/"blocks.tsv":lec/"LEC_0_blocks.tsv"})
		inputs<<hex_sha256_file(path)<<"  "<<fs::absolute(path).generic_string()<<'\n';
	if(h>1){auto path=fec/joint_eet_weight_name(h-1)/"multiplicities.tsv";inputs<<hex_sha256_file(path)<<"  "<<fs::absolute(path).generic_string()<<'\n';}
	inputs.close();
	auto W=nmhv_sew_read_tensor<rat_t,int32_t>(condition,QQ,pool);auto cb=nmhv_sew_read_irrep_axis(blocks);
	auto t=irr::make_seam_template(cg,W,cb,opt);auto layout=irr::seam_layout(t,of,ol,mf,ml);
	// Keep the QQ matrix only for the original path or explicit expanded verification.
	sparse_mat<rat_t,int32_t> matrix,reference;
	size_t input_nnz=0,input_rows=0,input_height=0;
	if(prereduce=="none"||verify) {
		matrix=irr::assemble_seam(cg,t,layout,of,ol,first,last,opt);
		if(verify)verify_assembly(cg,t,layout,of,ol,first,last,W,matrix,pool);
		pool->detach_loop(0,matrix.nrow,[&](size_t r){fec_irrep_cancel_divisor(matrix[r]);});pool->wait();
		if(prereduce=="none") {
			input_nnz=matrix.nnz();input_rows=matrix.nrow;
			input_height=joint_eet_rational_stats(matrix).height_bits;
		}
		if(verify)reference=matrix;
		if(prereduce!="none")matrix.clear();
	}
	std::cout<<"IRR_SEW h="<<h<<" k="<<k<<" prereduce="<<prereduce<<" columns="<<layout.columns<<std::endl;
	auto assemble=[&](const field_t& F,rref_option_t& options) {
		if(prereduce=="none")return nmhv_sew_modular_reduce(matrix,F,&options->pool);
		auto reduced=irr::partial_seam(prereduce=="left",cg,t,layout,of,ol,first,last,F,options);
		input_nnz=std::max(input_nnz,reduced.nnz());input_rows=std::max(input_rows,reduced.nrow);
		std::cout<<"IRR_SEW_MOD prime="<<F.mod.n<<" rows="<<reduced.nrow<<" nnz="<<reduced.nnz()<<std::endl;
		return reduced;
	};
	size_t primes=0;
    std::function<staged_kernel::MM(const field_t&,rref_option_t&)> direct;
    if(solver=="staged-factored")direct=[&](const field_t& F,rref_option_t& options){
        IrrepSeamRows rows(cg,t,layout,of,ol,first,last,F,options);
        input_rows=std::max(input_rows,rows.nrow());
        return factored_staged_kernel(rows,F,options);
    };
    auto reconstructed=solve_seam(assemble,layout.columns,h,k,root,opt,primes,direct);
	auto kernel=std::move(reconstructed.basis);matrix.clear();
	pool->detach_loop(0,kernel.nrow,[&](size_t r){kernel[r].canonicalize();fec_irrep_cancel_divisor(kernel[r]);});pool->wait();kernel.sort_rows_by_nnz();
	require(expected==SIZE_MAX||kernel.nrow==expected,"Unexpected bulk solution dimension: "+std::to_string(kernel.nrow));
	if(verify) {
		auto exact=sparse_mat_mul(reference,kernel.transpose(),QQ,pool);for(auto& r:exact.rows)r.canonicalize();
		require(exact.nnz()==0,"Exact rational kernel residual is nonzero");exact.clear();
		for(ulong prime:{1000003UL,1000033UL}) {
			field_t F(FIELD_Fp,prime);auto A=nmhv_sew_modular_reduce(reference,F,pool),K=nmhv_sew_modular_reduce(kernel,F,pool);
			auto residual=sparse_mat_mul(A,K.transpose(),F,pool);for(auto& r:residual.rows)r.canonicalize();
			require(residual.nnz()==0,"Sewing kernel residual is nonzero");
			auto rank_opt=opt;rank_opt->verbose=false;auto p=sparse_mat_rref(K,F,rank_opt);size_t rank=0;for(const auto& group:p)rank+=group.size();
			require(rank==kernel.nrow,"Sewing solution rows are dependent");
			auto constraint_pivots=sparse_mat_rref(A,F,rank_opt);
			require(A.ncol-nmhv_sew_modular_rank(constraint_pivots)==kernel.nrow,"Kernel is not the full rational solution space");
		}
		std::cout<<"VERIFY exact_QQ_kernel_residual=PASS independent_solution_rows=PASS full_kernel_dimension=PASS primes=1000003,1000033"<<std::endl;
	}
	std::ofstream columns(root/"columns.tsv"),hashes(root/"checksums.sha256"),routes(root/"routes.tsv"),rows(root/"rows.tsv");
	columns<<"block\tfirst_copies\tlast_copies\tcolumn_offset\tfilename\n";
	rows<<"template_row\told_first_block\told_last_block\tcondition_block\tcondition_copy\trow_offset\tfirst_copies\tlast_copies\n";
	for(size_t r=0;r<t.rows.size();++r)rows<<r<<'\t'<<irr::specs()[t.rows[r].left].name<<'\t'<<irr::specs()[t.rows[r].right].name<<'\t'<<cb[t.rows[r].condition].spec.name<<'\t'<<t.rows[r].copy<<'\t'<<layout.rowbase[r]<<'\t'<<of[t.rows[r].left]<<'\t'<<ol[t.rows[r].right]<<'\n';
	routes<<"block\tfirst_channel\tlast_channel\ttemplate_row\tcoefficient\n";
	for(size_t b=0;b<6;++b) {
		for(size_t l=0;l<t.routes[b].size();++l)for(size_t r=0;r<t.routes[b][l].size();++r)
			for(const auto& route:t.routes[b][l][r])routes<<irr::specs()[b].name<<'\t'<<l<<'\t'<<r<<'\t'<<route.row<<'\t'<<route.value<<'\n';
		sparse_mat<rat_t,int32_t> block(kernel.nrow,mf[b]*ml[b]);
		irr_parallel::rows(kernel.nrow,pool,[&](size_t i){
			size_t lo=0,hi=kernel[i].nnz();
			while(lo<hi){size_t mid=(lo+hi)/2;if(size_t(kernel[i](mid))<layout.colbase[b])lo=mid+1;else hi=mid;}
			for(size_t p=lo;p<kernel[i].nnz()&&size_t(kernel[i](p))<layout.colbase[b]+block.ncol;++p)
				block[i].push_back(int32_t(kernel[i](p)-layout.colbase[b]),kernel[i][p]);
		});
		sparse_tensor<rat_t,int32_t,SPARSE_COO> coo(std::move(block),pool);coo.reshape({kernel.nrow,mf[b],ml[b]});
		sparse_tensor<rat_t,int32_t,SPARSE_CSR> tensor(std::move(coo),pool);
		std::string name="irrSEW_"+std::to_string(h)+"p"+std::to_string(k)+"_"+irr::specs()[b].name+".wxf";
		auto actual=irr::tensor_path(root/name);
		if(!irr_run::ready(root/name)){actual=irr::write_tensor(root/name,tensor);irr_run::mark(root/name);}
		auto reload=irr::read_tensor(root/name,pool);
		require(nmhv_sew_matrices_equal(block,joint_eet_flatten_tensor(reload,pool)),"Sewing WXF roundtrip mismatch");
		columns<<irr::specs()[b].name<<'\t'<<mf[b]<<'\t'<<ml[b]<<'\t'<<layout.colbase[b]<<'\t'<<name<<'\n';hashes<<hex_sha256_file(actual)<<"  "<<actual.filename().string()<<'\n';
	}
	columns.close();rows.close();routes.close();
	if(!fs::exists(root/"CG.tsv"))fs::copy_file(lec/"common"/"CG.tsv",root/"CG.tsv");
	std::ofstream chart(root/"row_chart.tsv");chart<<"template_row\told_first_component\tcondition_coordinate\told_last_component\tcoefficient\n";
	size_t prototype_dimension=fec_irrep_axis_dimension(cg.axis);
	for(size_t r=0;r<t.rows.size();++r)for(const auto& [coord,value]:t.rows[r].embedding)
		chart<<r<<'\t'<<coord/prototype_dimension/W.dim(2)<<'\t'<<coord/prototype_dimension%W.dim(2)<<'\t'<<coord%prototype_dimension<<'\t'<<value<<'\n';chart.close();
	for(const std::string file:{"inputs.sha256","rows.tsv","columns.tsv","routes.tsv","CG.tsv","row_chart.tsv"})hashes<<hex_sha256_file(root/file)<<"  "<<file<<'\n';hashes.close();
	auto stats=joint_eet_rational_stats(kernel);auto memory=hex_platform::current_process_memory();
	std::ofstream summary(root/"summary.tsv");summary<<"helicity\th\tk\trows\tcolumns\tconstraint_nnz\tconstraint_height_bits\tsolutions\tsewing_nnz\theight_bits\tmodulus_bits\tprimes\tpeak_bytes\tverified\n"
		<<helicity<<'\t'<<h<<'\t'<<k<<'\t'<<input_rows<<'\t'<<layout.columns<<'\t'<<input_nnz<<'\t'<<(prereduce=="none"?std::to_string(input_height):"NA")<<'\t'<<kernel.nrow<<'\t'<<kernel.nnz()<<'\t'<<stats.height_bits<<'\t'<<reconstructed.modulus_bits<<'\t'<<primes<<'\t'<<memory.peak_working_set_bytes<<'\t'<<verify<<'\n';
	std::ofstream manifest(root/"manifest.json");manifest<<"{\"schema\":\"hexagon-irr-sew-v1\",\"complete\":true,\"helicity\":\""<<helicity<<"\",\"h\":"<<h<<",\"k\":"<<k
		<<",\"solver\":\""<<solver<<"\",\"full_product_stored\":"<<(solver=="baseline"?"true":"false")<<",\"prereduce\":\""<<prereduce<<"\",\"constraint_arithmetic\":\""<<(prereduce=="none"?"QQ":"Fp")<<"\",\"solutions\":"<<kernel.nrow<<",\"method\":1,\"reconstruction_primes\":"<<primes<<",\"validation_primes\":2,\"modulus_bits\":"<<reconstructed.modulus_bits<<",\"seed_normalization\":\"physical-copy-primitive-v1\",\"layout\":\"solution,first_copy,last_copy\",\"singlet_E_metric\":[1,3],\"shared_solution_axis\":true,\"verified_expanded\":"<<(verify?"true":"false")
		<<",\"checksums_sha256\":\""<<hex_sha256_file(root/"checksums.sha256")<<"\",\"stage\":\"bulk\"}\n";
	manifest.close();require(bool(manifest),"Cannot write sewing manifest");irr_run::mark(root/"manifest.json");
	irr_stage::write(root/"basis_identity.tsv",irr_stage::basis_identity(fec,lec,h,k));irr_stage::complete_chain(root,false);
	std::cout<<"COMPLETE solutions="<<kernel.nrow<<" nnz="<<kernel.nnz()<<" height_bits="<<stats.height_bits<<" peak_bytes="<<memory.peak_working_set_bytes
		<<" seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}}
