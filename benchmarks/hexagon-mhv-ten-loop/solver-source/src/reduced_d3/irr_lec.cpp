#include "irr_recurrence.hpp"
#include "irr_stage.hpp"
#include "irr_runtime.hpp"
#include "mhv_sewing.hpp"
#include <chrono>

namespace fs=std::filesystem;
const field_t QQ(FIELD_QQ);
const char* seed_convention="physical-copy-primitive-v1";

void require(bool ok,const std::string& message) {
	if(!ok)throw std::runtime_error(message);
}
sparse_mat<rat_t,int32_t> read_rows(const fs::path& path,thread_pool* pool) {
	auto tensor=irr::read_tensor(path,pool);
	return joint_eet_flatten_tensor(tensor,pool);
}
void write_matrix(const fs::path& path,const sparse_mat<rat_t,int32_t>& m,thread_pool* pool) {
	if(irr_run::ready(path)) {require(nmhv_sew_matrices_equal(m,read_rows(path,pool)),"Existing matrix changed: "+path.string());return;}
	sparse_tensor<rat_t,int32_t,SPARSE_CSR> t(m,pool);
	irr::write_tensor(path,t);
	irr_run::mark(path);
}
std::vector<fec_irrep_axis_block_t> axis(const irr::copies& m,thread_pool* pool) {
	return fec_irrep_axis_from_blocks(irr::skeleton(m,pool));
}
void write_blocks(const fs::path& path,const std::vector<fec_irrep_axis_block_t>& blocks) {
	if(!fs::exists(path)){joint_eet_write_blocks(path,blocks);return;}
	auto old=nmhv_sew_read_irrep_axis(path);require(old.size()==blocks.size(),"Existing blocks changed");
	for(size_t b=0;b<blocks.size();++b)require(old[b].spec.name==blocks[b].spec.name&&old[b].offset==blocks[b].offset
		&&old[b].dimension==blocks[b].dimension&&old[b].multiplicity==blocks[b].multiplicity,"Existing block metadata changed");
}
irr::copies multiplicities(const std::vector<fec_irrep_axis_block_t>& blocks) {
	irr::copies m{};
	for(const auto& b:blocks)for(size_t i=0;i<6;++i)if(b.spec.name==irr::specs()[i].name)m[i]=b.multiplicity;
	return m;
}
sparse_mat<rat_t,int32_t> full_action(const irr::copies& m,const std::string& name) {
	size_t n=0;for(size_t b=0;b<6;++b)n+=m[b]*irr::specs()[b].irrep_dimension;
	sparse_mat<rat_t,int32_t> out(n,n);size_t off=0;
	for(size_t b=0;b<6;++b) {
		auto D=irr::action(b,name);
		for(size_t s=0;s<D.nrow;++s)for(size_t a=0;a<m[b];++a)
			for(auto [t,c]:D[s])out[off+s*m[b]+a].push_back(int32_t(off+t*m[b]+a),c);
		off+=m[b]*D.nrow;
	}
	return out;
}

// Convert only the small weight-one physical seed to copy/channel coordinates.
irr::matrices seed(bool mhv,const fs::path& data,const irr::channels_t& cg,
	irr::copies& zero,sparse_mat<rat_t,int32_t>& components,
	sparse_tensor<rat_t,int32_t,SPARSE_CSR>& physical_seed,std::vector<rat_t>& seed_scale,rref_option_t& opt) {
	auto cyc=nmhv_sew_read_matrix<rat_t,int32_t>(data/"cycrepmat.wxf",QQ);
	auto flip=nmhv_sew_read_matrix<rat_t,int32_t>(data/"fliprepmat.wxf",QQ);
	joint_eet_lec_t<int32_t> s;
	if(mhv)s=mhv_lec_seed(cyc,flip,opt);
	else s=make_joint_eet_lec_seed(nmhv_sew_read_matrix<rat_t,int32_t>(data/"nmhv_qbar_basis.wxf",QQ),
		nmhv_sew_read_matrix<rat_t,int32_t>(data/"nmhv_cyclic_R.wxf",QQ),cyc,
		nmhv_sew_read_matrix<rat_t,int32_t>(data/"nmhv_flip_R.wxf",QQ),flip,QQ,opt);
	components=joint_eet_identity<rat_t,int32_t>(mhv?1:5);zero={1,0,0,0,0,0};
	if(!mhv) {
		tail_state_t<rat_t,int32_t> state;state.weight=0;state.r_count=1;state.letter_count=5;
		state.basis=components;state.cyclic=s.cyclic_component.transpose();state.flip=s.flip_component.transpose();
		auto dec=nmhv_sew_irrep_tail_decomposition(state,QQ,opt);components=std::move(dec.basis);zero=multiplicities(dec.blocks);
	}
	auto rows=joint_eet_flatten_tensor(s.LEC,&opt->pool);
	seed_scale.assign(rows.nrow,rat_t(1));
	// Normalize in the PHYSICAL terminal basis, jointly across each irrep copy.
	// Independent rescaling of E components would change the fixed CG/actions.
	for(const auto& block:s.blocks)for(size_t a=0;a<block.multiplicity;++a) {
		int_t den=1,gcd=0;
		for(size_t c=0;c<block.spec.irrep_dimension;++c)
			for(auto [j,x]:rows[block.offset+c*block.multiplicity+a])den=LCM(den,x.den());
		for(size_t c=0;c<block.spec.irrep_dimension;++c)
			for(auto [j,x]:rows[block.offset+c*block.multiplicity+a]) {
				int_t value=den;value/=x.den();value*=x.num();gcd=GCD(gcd,value);
			}
		require(gcd!=0,"Zero physical seed copy");
		for(size_t c=0;c<block.spec.irrep_dimension;++c) {
			size_t r=block.offset+c*block.multiplicity+a;
			seed_scale[r]=rat_t(den,gcd);sparse_vec_rescale(rows[r],seed_scale[r],QQ);
		}
	}
	sparse_tensor<rat_t,int32_t,SPARSE_COO> physical(rows,&opt->pool);
	physical.reshape({rows.nrow,9,components.nrow});
	physical_seed=sparse_tensor<rat_t,int32_t,SPARSE_CSR>(std::move(physical),&opt->pool);
	if(!mhv)rows=nmhv_sew_transform_product_basis(rows,joint_eet_identity<rat_t,int32_t>(9),sparse_mat_inverse(components,QQ,opt),QQ,&opt->pool);
	sparse_tensor<rat_t,int32_t,SPARSE_COO> coo(rows,&opt->pool);
	coo.reshape({rows.nrow,9,components.nrow});coo.transpose_replace({0,2,1});
	coo.flatten({{0},{1,2}});rows=coo.to_sparse_mat(&opt->pool);
	auto zero_blocks=irr::skeleton(zero,&opt->pool);irr::matrices K;
	for(size_t b=0;b<6;++b) {
		auto candidate=make_hhex_candidates(fec_symmetry_t::parity_dihedral,zero_blocks,irr::specs()[b].name);
		sparse_mat<rat_t,int32_t> C(candidate.size(),rows.ncol);
		for(size_t c=0;c<candidate.size();++c)for(const auto& term:candidate.vectors[c].component[0])C[c].push_back(term.coordinate,term.coefficient);
		for(auto& row:C.rows)fec_irrep_merge_entries(row);
		auto m=multiplicities(s.blocks)[b];sparse_mat<rat_t,int32_t> target(m,rows.ncol);
		for(const auto& block:s.blocks)if(block.spec.name==irr::specs()[b].name)
			for(size_t r=0;r<m;++r)target[r]=rows[block.offset+r];
		if(!m){K[b]=sparse_mat<rat_t,int32_t>(0,C.nrow);continue;}
		auto pivots=dlogmat_irrep_pivot_columns(C,opt);
		require(pivots.size()==C.nrow,"Dependent seed CG channels");
		auto inverse=sparse_mat_inverse(dlogmat_irrep_select_columns(C,pivots),QQ,opt);
		K[b]=sparse_mat_mul(dlogmat_irrep_select_columns(target,pivots),inverse,QQ,&opt->pool);
		auto expanded=irr::expand(cg,b,zero,K[b]);
		for(const auto& block:s.blocks)if(block.spec.name==irr::specs()[b].name)
			for(size_t r=0;r<expanded.nrow;++r) {
				sparse_mat<rat_t,int32_t> a(1,rows.ncol),v(1,rows.ncol);a[0]=expanded[r];v[0]=rows[block.offset+r];
				require(nmhv_sew_matrices_equal(a,v),"Seed CG roundtrip mismatch");
			}
	}
	return K;
}

// Independent old-coordinate assembler used ONLY in --verify mode.
sparse_mat<rat_t,int32_t> reference_matrix(const fs::path& ref,size_t w,size_t target,
	const irr::copies& old,const irr::copies& previous,const sparse_tensor<rat_t,int32_t,SPARSE_CSR>& W,
	const std::vector<fec_irrep_axis_block_t>& condition,const std::vector<rat_t>& old_scale,
	const std::vector<rat_t>& previous_scale,rref_option_t opt) {
	auto prev=read_rows(ref/joint_eet_weight_name(w-1)/"LEC_internal.wxf",&opt->pool);
	for(size_t r=0;r<prev.nrow;++r)for(size_t p=0;p<prev[r].nnz();++p)
		prev[r][p]*=previous_scale.at(r)/old_scale.at(prev[r](p)/9);
	auto blocks=irr::skeleton(previous,&opt->pool);auto pa=axis(previous,&opt->pool);
	auto candidates=make_hhex_candidates(fec_symmetry_t::parity_dihedral,blocks,irr::specs()[target].name);
	auto codomain=make_hhex_codomain_basis(axis(old,&opt->pool),condition,irr::specs()[target],QQ,opt);
	std::vector<int32_t> selected(codomain.raw_dimension,-1);
	for(size_t r=0;r<codomain.multiplicity;++r)selected[codomain.pivot_coordinates[r]]=int32_t(r);
	sparse_mat<rat_t,int32_t> matrix(codomain.multiplicity,candidates.size());
	for(size_t a=0;a<prev.nrow;++a)for(auto [col,x]:prev[a]) {
		size_t right=col/9,letter=col%9;
		for(size_t p=W.rowptr()[letter];p<W.rowptr()[letter+1];++p) {
			auto idx=W.index(p);int32_t r=selected[right*W.dim(2)+idx[1]];if(r<0)continue;
			for(const auto& c:candidates.first_component_lookup[a*9+idx[0]])matrix[r].push_back(c.column,x*W.val(p)*c.coefficient);
		}
	}
	for(auto& r:matrix.rows)fec_irrep_merge_entries(r);
	return sparse_mat_mul(codomain.transpose_inverse_chart,matrix,QQ,&opt->pool);
}

std::vector<rat_t> verify(const fs::path& ref,size_t w,const irr::channels_t& cg,const irr::copies& old,
	const irr::copies& current,const irr::matrices& K,const sparse_mat<rat_t,int32_t>& components,
	const std::vector<rat_t>& old_scale,rref_option_t opt) {
	auto dir=ref/joint_eet_weight_name(w);
	auto expected=read_rows(dir/"LEC_internal.wxf",&opt->pool);
	sparse_mat<rat_t,int32_t> all(0,expected.ncol);
	for(size_t b=0;b<6;++b)all.append(irr::expand(cg,b,old,K[b]));
	require(all.nrow==expected.nrow&&all.ncol==expected.ncol,"Reference dimensions mismatch");
	std::vector<rat_t> scale(all.nrow,rat_t(1));
	// Infer the output-copy rescaling, then check EVERY coefficient exactly.
	for(size_t r=0;r<expected.nrow;++r) {
		for(size_t p=0;p<expected[r].nnz();++p)expected[r][p]/=old_scale.at(expected[r](p)/9);
		require(expected[r].nnz()&&all[r].nnz(),"Zero basis row in reference check");
		scale[r]=all[r][0]/expected[r][0];sparse_vec_rescale(expected[r],scale[r],QQ);
	}
	require(nmhv_sew_matrices_equal(all,expected),"Expanded LEC mismatch at weight "+std::to_string(w));
	auto expected_blocks=nmhv_sew_read_irrep_axis(dir/("LEC_"+std::to_string(w)+"_blocks.tsv"));
	require(current==multiplicities(expected_blocks),"Multiplicity metadata mismatch");
	for(const auto& b:expected_blocks)for(size_t a=0;a<b.multiplicity;++a)
		for(size_t c=1;c<b.spec.irrep_dimension;++c)
			require(scale[b.offset+a]==scale[b.offset+c*b.multiplicity+a],"Reference scaling changes irrep convention");
	sparse_tensor<rat_t,int32_t,SPARSE_COO> tensor(all,&opt->pool);
	tensor.reshape({all.nrow,all.ncol/9,9});tensor.transpose_replace({0,2,1});
	tensor.flatten({{0},{1,2}});auto published=tensor.to_sparse_mat(&opt->pool);
	if(w==1)published=nmhv_sew_transform_product_basis(published,joint_eet_identity<rat_t,int32_t>(9),components,QQ,&opt->pool);
	auto reference_published=read_rows(dir/("LEC_"+std::to_string(w)+".wxf"),&opt->pool);
	for(size_t r=0;r<reference_published.nrow;++r)for(size_t p=0;p<reference_published[r].nnz();++p)
		reference_published[r][p]*=scale[r]/old_scale.at(reference_published[r](p)%old_scale.size());
	require(nmhv_sew_matrices_equal(published,reference_published),"Published LEC mismatch after reference basis conversion");
	for(const std::string name:{"C","Cyclic","Flip","Parity"}) {
		auto D=full_action(current,name);
		if(name=="Cyclic"||name=="Flip")require(nmhv_sew_matrices_equal(D.transpose(),
			nmhv_sew_read_matrix<rat_t,int32_t>(dir/(name+"_LEC_"+std::to_string(w)+".wxf"),QQ)),"Old action mismatch");
		sparse_mat<rat_t,int32_t> transformed(0,all.ncol);
		for(size_t b=0;b<6;++b)transformed.append(irr::apply_action(b,name,irr::expand(cg,b,old,K[b])));
		require(nmhv_sew_matrices_equal(transformed,sparse_mat_mul(D,all,QQ,&opt->pool)),"Compact action mismatch");
		sparse_mat<rat_t,int32_t> letter(9,9);
		for(size_t i=0;i<9;++i) {
			size_t j=i;int sign=1;
			if(name=="C"||name=="Cyclic")j=3*(i/3)+(i+1)%3;
			if(name=="Flip")j=3*(i/3)+2-i%3;
			if((name=="Parity"||name=="Cyclic")&&i>=6)sign=-1;
			letter[i].push_back(int32_t(j),rat_t(sign));
		}
		auto image=nmhv_sew_transform_product_basis(all,full_action(old,name),letter,QQ,&opt->pool);
		require(nmhv_sew_matrices_equal(image,transformed),"Recurrence action intertwiner mismatch: "+name);
	}
	std::cout<<"VERIFY weight="<<w<<" reference_diagonal_change=PASS published_basis_change=PASS actions=PASS intertwiner=PASS expanded_nnz="<<all.nnz()<<std::endl;
	return scale;
}

void write_common(const fs::path& root,const irr::channels_t& cg,const std::array<irr::template_t,6>& templates,
	const fs::path& condition,const fs::path& condition_blocks,thread_pool* pool) {
	fs::create_directories(root);
	if(!fs::exists(root/"condition_blocks.tsv"))fs::copy_file(condition_blocks,root/"condition_blocks.tsv");
	std::ofstream channels(root/"channels.tsv"),terms(root/"CG.tsv"),rows(root/"condition_rows.tsv"),cols(root/"condition_columns.tsv");
	channels<<"target\tinput\tchannel\tletter_triplet\tletter_irrep\n";
	terms<<"target\tinput\tchannel\ttarget_component\tinput_component\tletter\tcoefficient\n";
	rows<<"target\trow\told_irrep\tcondition_block\tcondition_copy\traw_to_primitive\n";
	cols<<"target\tcolumn\tintermediate_irrep\tprevious_channel\tnew_local_channel\n";
	for(size_t b=0;b<6;++b) {
		for(size_t input=0;input<6;++input)for(size_t g:cg.channels[b][input]) {
			const auto& v=cg.cg[b].vectors[g];size_t triplet=(v.component[0][0].coordinate%9)/3;
			std::map<size_t,rat_t> sums;for(const auto& term:v.component[0])sums[term.coordinate/9]+=term.coefficient;
			bool A1=false;for(auto& [i,x]:sums)if(x!=0)A1=true;
			channels<<irr::specs()[b].name<<'\t'<<irr::specs()[input].name<<'\t'<<g<<'\t'<<triplet<<'\t'<<(A1?"A1":"E")<<'\n';
			for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)for(const auto& term:v.component[s])
				terms<<irr::specs()[b].name<<'\t'<<irr::specs()[input].name<<'\t'<<g<<'\t'<<s<<'\t'
					<<term.coordinate/9-cg.axis[input].offset<<'\t'<<term.coordinate%9<<'\t'<<term.coefficient<<'\n';
		}
		const auto& t=templates[b];size_t n=0;
		for(size_t input=0;input<6;++input)n+=cg.cg[input].size()*cg.channels[b][input].size();
		sparse_mat<rat_t,int32_t> matrix(t.rows.size(),n);size_t col=0;
		for(size_t input=0;input<6;++input)for(size_t g=0;g<cg.cg[input].size();++g) {
			for(size_t h=0;h<cg.channels[b][input].size();++h)
				cols<<irr::specs()[b].name<<'\t'<<col+h<<'\t'<<irr::specs()[input].name<<'\t'<<g<<'\t'<<h<<'\n';
			for(const auto& r:t.routes[input][g])matrix[r.row].push_back(int32_t(col+r.channel),r.coefficient*t.row_scale[r.row]);
			col+=cg.channels[b][input].size();
		}
		for(size_t r=0;r<t.rows.size();++r) {
			fec_irrep_merge_entries(matrix[r]);const auto& row=t.rows[r];
			rows<<irr::specs()[b].name<<'\t'<<r<<'\t'<<irr::specs()[row.left].name<<'\t'<<row.condition<<'\t'<<row.copy<<'\t'<<t.row_scale[r]<<'\n';
		}
		write_matrix(root/("dlogmatIRR_backward_"+irr::specs()[b].name+".wxf"),matrix,pool);
		for(const std::string name:{"C","Cyclic","Flip","Parity"})write_matrix(root/(name+"_"+irr::specs()[b].name+".wxf"),irr::action(b,name),pool);
	}
	std::ofstream meta(root/"convention.json");
	meta<<"{\n  \"schema\": \"hexagon-irr-lec-v2\",\n  \"seed_normalization\": \""<<seed_convention<<"\",\n  \"index_base\": 0,\n"
		<<"  \"block_order\": [\"pe_A1\",\"pe_A2\",\"pe_E\",\"po_A1\",\"po_A2\",\"po_E\"],\n"
		<<"  \"column_order\": [\"old_irrep\",\"old_copy\",\"channel\"],\n"
		<<"  \"cg_axes\": [\"target_component\",\"old_component\",\"physical_letter\"],\n"
		<<"  \"physical_order\": \"new_letter before old_tail; terminal component last\",\n"
		<<"  \"action\": \"row basis D tensor I_m; C order 3; Cyclic=Parity*C\",\n"
		<<"  \"condition_axes\": \"source W[new_letter,old_letter,q]; assembly swaps first two axes\",\n"
		<<"  \"condition_sha256\": \""<<hex_sha256_file(condition)<<"\",\n"
		<<"  \"condition_blocks_sha256\": \""<<hex_sha256_file(condition_blocks)<<"\"\n}\n";
}

int main(int argc,char** argv){try{
	require(argc>=4,"irr_lec MHV|NMHV MAX_WEIGHT OUTPUT [--verify REF] [--data DIR] [--condition WXF] [--condition-blocks TSV]");
	bool mhv=std::string(argv[1])=="MHV";require(mhv||std::string(argv[1])=="NMHV","Unknown helicity");
	size_t maximum=std::stoul(argv[2]);require(maximum>0,"MAX_WEIGHT must be positive");fs::path root=argv[3],ref,data="data";
	fs::path condition="data/dlogmat_full_ca_pd3.wxf",condition_blocks="data/dlogmat_full_ca_pd3_blocks.tsv";
	for(int i=4;i<argc;++i){std::string key=argv[i];require(i+1<argc,"Missing option value");fs::path value=argv[++i];
		if(key=="--verify")ref=value;else if(key=="--data")data=value;else if(key=="--condition")condition=value;
		else if(key=="--condition-blocks")condition_blocks=value;else throw std::runtime_error("Unknown option "+key);}
	std::vector<fs::path> identity{condition,condition_blocks,data/"cycrepmat.wxf",data/"fliprepmat.wxf"};
	if(!mhv)for(const auto& file:{"nmhv_qbar_basis.wxf","nmhv_cyclic_R.wxf","nmhv_flip_R.wxf"})identity.push_back(data/file);
	irr_run::bind(root,std::string("LEC physical-copy-primitive-v1 ")+(mhv?"MHV":"NMHV"),identity);
	irr_stage::lock active(root/".active.lock");
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=true;irr_run::threads(opt);
	auto start=std::chrono::steady_clock::now();auto*pool=&opt->pool;irr::channels_t cg(pool);
	auto W=sparse_tensor_read_wxf<rat_t,int32_t>(condition,QQ,pool);
	sparse_tensor<rat_t,int32_t,SPARSE_COO> dc(std::move(W));dc.transpose_replace({1,0,2});W=sparse_tensor<rat_t,int32_t,SPARSE_CSR>(std::move(dc),pool);
	auto cb=nmhv_sew_read_irrep_axis(condition_blocks);std::array<irr::template_t,6> templates;
	for(size_t b=0;b<6;++b)templates[b]=irr::make_template(cg,b,W,cb,opt);
	irr::copies older{},previous{};sparse_mat<rat_t,int32_t> components;
	sparse_tensor<rat_t,int32_t,SPARSE_CSR> physical_seed;std::vector<rat_t> seed_scale;
	auto K=seed(mhv,data,cg,previous,components,physical_seed,seed_scale,opt);
	write_common(root/"common",cg,templates,condition,condition_blocks,pool);
	if(!irr_run::ready(root/"LEC_1_physical.wxf")){irr::write_tensor(root/"LEC_1_physical.wxf",physical_seed);irr_run::mark(root/"LEC_1_physical.wxf");}
	std::ofstream seed_meta(root/"seed_scale.tsv");seed_meta<<"basis_index\tnew_over_original_seed\n";
	for(size_t r=0;r<seed_scale.size();++r)seed_meta<<r<<'\t'<<seed_scale[r]<<'\n';seed_meta.close();
	write_matrix(root/"LEC_0_components.wxf",components,pool);
	write_blocks(root/"LEC_0_blocks.tsv",axis(previous,pool));
	std::ofstream inputs(root/"inputs.sha256");
	for(const std::string name:{"cycrepmat.wxf","fliprepmat.wxf","nmhv_qbar_basis.wxf","nmhv_cyclic_R.wxf","nmhv_flip_R.wxf"}) {
		if(mhv&&name.substr(0,5)=="nmhv_")continue;
		inputs<<hex_sha256_file(data/name)<<"  "<<fs::absolute(data/name).generic_string()<<'\n';
	}
	inputs<<hex_sha256_file(condition)<<"  "<<fs::absolute(condition).generic_string()<<'\n';
	inputs<<hex_sha256_file(condition_blocks)<<"  "<<fs::absolute(condition_blocks).generic_string()<<'\n';inputs.close();
	std::ofstream common_hashes(root/"common"/"checksums.sha256");
	for(const auto& entry:fs::directory_iterator(root/"common"))if(entry.path().filename()!="checksums.sha256")
		common_hashes<<hex_sha256_file(entry.path())<<"  "<<entry.path().filename().string()<<'\n';
	common_hashes.close();
	std::vector<rat_t> older_reference_scale(components.nrow,rat_t(1)),previous_reference_scale=older_reference_scale;
	std::ofstream summary(root/"summary.tsv");summary<<"weight\tdimension\tirr_nnz\theight_bits\tfractions\tmax_rows\tmax_columns\tmax_matrix_nnz\tpeak_bytes\n";
	for(size_t w=1;w<=maximum;++w) {
		auto dir=root/joint_eet_weight_name(w);fs::create_directories(dir);
		auto block_path=[&](size_t b){return dir/("irrLEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf");};
		size_t maxrows=0,maxcols=0,maxnnz=0;
		if(w>1) {
			irr::matrices next;
			for(size_t b=0;b<6;++b) {
				if(irr_run::ready(block_path(b))) {next[b]=read_rows(block_path(b),pool);std::cout<<"REUSE LEC weight="<<w<<" block="<<irr::specs()[b].name<<std::endl;continue;}
				auto matrix=irr::assemble(cg,templates[b],b,cb,older,previous,K,opt);
				maxrows=std::max(maxrows,matrix.nrow);maxcols=std::max(maxcols,matrix.ncol);maxnnz=std::max(maxnnz,matrix.nnz());
				if(!ref.empty()) {
					auto raw=irr::assemble(cg,templates[b],b,cb,older,previous,K,opt,false);
					auto old=reference_matrix(ref,w,b,older,previous,W,cb,older_reference_scale,previous_reference_scale,opt);
					require(nmhv_sew_matrices_equal(raw,old),"Raw assembly differs from reference");
					auto normalized=matrix;
					for(auto* m:{&raw,&normalized})for(auto& row:m->rows) {
						fec_irrep_cancel_divisor(row);if(row.nnz()&&row[0]<0)sparse_vec_rescale(row,rat_t(-1),QQ);
					}
					require(nmhv_sew_matrices_equal(raw,normalized),"Primitive template changes equations");
				}
				std::cout<<"IRR_LEC weight="<<w<<" block="<<irr::specs()[b].name<<" rows="<<matrix.nrow<<" cols="<<matrix.ncol<<" nnz="<<matrix.nnz()<<std::endl;
				next[b]=irr::solve(matrix,opt,dir/"cache"/irr::specs()[b].name);
				// Keep each completed block even if a subsequent solve/check fails.
				write_matrix(block_path(b),next[b],pool);
			}
			K=std::move(next);
		}
		irr::copies current{};size_t nnz=0,dim=0,height=0,fractions=0;
		for(size_t b=0;b<6;++b){current[b]=K[b].nrow;nnz+=K[b].nnz();dim+=current[b]*irr::specs()[b].irrep_dimension;
			auto stats=joint_eet_rational_stats(K[b]);height=std::max(height,stats.height_bits);fractions+=stats.fractions;}
		if(!ref.empty()) {
			auto scale=verify(ref,w,cg,previous,current,K,components,previous_reference_scale,opt);
			if(w==1)require(scale==seed_scale,"Reference seed normalization mismatch");
			std::ofstream basis_change(dir/"reference_scale.tsv");basis_change<<"basis_index\tnew_over_reference\n";
			for(size_t r=0;r<scale.size();++r)basis_change<<r<<'\t'<<scale[r]<<'\n';
			older_reference_scale=previous_reference_scale;previous_reference_scale=std::move(scale);
		}
		std::ofstream columns(dir/"columns.tsv"),hashes(dir/"checksums.sha256");columns<<"target\tinput\told_copies\tcolumn_offset\tchannels_per_copy\n";
		for(size_t b=0;b<6;++b) {
			auto file=block_path(b);if(w==1)write_matrix(file,K[b],pool);
		hashes<<hex_sha256_file(irr::tensor_path(file))<<"  "<<irr::tensor_path(file).filename().string()<<'\n';
			auto reload=read_rows(file,pool);require(nmhv_sew_matrices_equal(K[b],reload),"WXF roundtrip mismatch");K[b]=std::move(reload);
			size_t off=0;for(size_t input=0;input<6;++input){columns<<irr::specs()[b].name<<'\t'<<irr::specs()[input].name<<'\t'<<previous[input]<<'\t'<<off<<'\t'<<cg.channels[b][input].size()<<'\n';off+=previous[input]*cg.channels[b][input].size();}
		}
		write_blocks(dir/"blocks.tsv",axis(current,pool));
		columns.close();hashes<<hex_sha256_file(dir/"columns.tsv")<<"  columns.tsv\n"<<hex_sha256_file(dir/"blocks.tsv")<<"  blocks.tsv\n";
		if(!ref.empty())hashes<<hex_sha256_file(dir/"reference_scale.tsv")<<"  reference_scale.tsv\n";hashes.close();
		std::ofstream completed(dir/"manifest.json");
		completed<<"{\"schema\":\"hexagon-irr-lec-v2\",\"seed_normalization\":\""<<seed_convention<<"\",\"helicity\":\""<<(mhv?"MHV":"NMHV")<<"\",\"weight\":"<<w<<",\"complete\":true,\"dimension\":"<<dim<<",\"irr_nnz\":"<<nnz<<",\"height_bits\":"<<height<<",\"fractions\":"<<fractions
			<<",\"method\":1,\"verified_reference\":"<<(!ref.empty()?"true":"false")
			<<",\"inputs_sha256\":\""<<hex_sha256_file(root/"inputs.sha256")
			<<"\",\"common_sha256\":\""<<hex_sha256_file(root/"common"/"checksums.sha256")
			<<"\",\"seed_components_sha256\":\""<<hex_sha256_file(root/"LEC_0_components.wxf")
			<<"\",\"physical_seed_sha256\":\""<<hex_sha256_file(root/"LEC_1_physical.wxf")
			<<"\",\"seed_scale_sha256\":\""<<hex_sha256_file(root/"seed_scale.tsv")
			<<"\",\"checksums_sha256\":\""<<hex_sha256_file(dir/"checksums.sha256")
			<<"\",\"previous_manifest_sha256\":\""<<(w>1?hex_sha256_file(root/joint_eet_weight_name(w-1)/"manifest.json"):"")<<"\"}\n";
		summary<<w<<'\t'<<dim<<'\t'<<nnz<<'\t'<<height<<'\t'<<fractions<<'\t'<<maxrows<<'\t'<<maxcols<<'\t'<<maxnnz<<'\t'<<hex_platform::current_process_memory().peak_working_set_bytes<<std::endl;
		older=previous;previous=current;
		completed.close();irr_stage::complete_chain(dir);
		std::cout<<"COMPLETE weight="<<w<<" dimension="<<dim<<" irr_nnz="<<nnz<<std::endl;
	}
	std::cout<<"DONE seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
