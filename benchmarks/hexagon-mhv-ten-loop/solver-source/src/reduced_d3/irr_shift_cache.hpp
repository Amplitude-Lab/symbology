#pragma once
#include "irr_object.hpp"
#include "irr_runtime.hpp"
#include "irr_stage.hpp"
#include "irr_recurrence_store.hpp"

// Cache conversion never evaluates physical conditions or reconstructs the
// source sewing. All large contractions and target charts are over F_p.
namespace irr_shift_cache {
namespace fs=std::filesystem;
using irr::require;
struct bad_prime:std::runtime_error {using std::runtime_error::runtime_error;};
inline ulong residue(const rat_t& x,const field_t& F) {
	if(x.den()%F.mod==0)throw bad_prime("Recurrence/CG denominator vanishes modulo prime");
	return x%F.mod;
}
inline sparse_mat<ulong,int32_t> reduce(const sparse_mat<rat_t,int32_t>& A,const field_t& F) {
	sparse_mat<ulong,int32_t> out(A.nrow,A.ncol);
	for(size_t i=0;i<A.nrow;++i)for(auto [j,x]:A[i]){auto y=residue(x,F);if(y)out[i].push_back(j,y);}return out;
}
inline auto recurrence(const fs::path& root,size_t w,bool first,const irr::channels_t& cg,const field_t& F,thread_pool* pool) {
	std::array<sparse_mat<ulong,int32_t>,6> out;
	auto old=irr::chain_copies(root,w-1,first),m=irr::chain_copies(root,w,first);
	for(size_t b=0;b<6;++b){auto file=root/joint_eet_weight_name(w)/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf");
		auto tensor=irr::read_tensor(file,pool);require(tensor.rank()==2,"Recurrence must be a matrix");
		auto matrix=joint_eet_flatten_tensor(tensor,pool);size_t columns=0;
		for(size_t a=0;a<6;++a)columns+=old[a]*cg.channels[b][a].size();
		require(matrix.nrow==m[b]&&matrix.ncol==columns,"Recurrence dimensions differ");out[b]=reduce(matrix,F);
	}return out;
}
inline sparse_mat<ulong,int32_t> right_inverse(const sparse_mat<ulong,int32_t>& K,const field_t& F,rref_option_t opt) {
	sparse_mat<ulong,int32_t> inverse(K.ncol,K.nrow);if(!K.nrow)return inverse;
	std::vector<size_t> count(K.ncol),pivot(K.nrow,SIZE_MAX);
	for(const auto& row:K.rows)for(auto [j,x]:row)++count[j];
	bool private_chart=true;
	for(size_t i=0;i<K.nrow;++i){for(auto [j,x]:K[i])if(count[j]==1){pivot[i]=j;inverse[j].push_back(i,n_invmod(x,F.mod.n));break;}
		if(pivot[i]==SIZE_MAX){private_chart=false;break;}}
	if(!private_chart){auto reduced=K;auto groups=sparse_mat_rref(reduced,F,opt);pivot.clear();
		for(const auto& group:groups)for(const auto& p:group)pivot.push_back(p.c);
		if(pivot.size()!=K.nrow)throw bad_prime("Target recurrence loses row rank modulo prime");
		std::sort(pivot.begin(),pivot.end());std::vector<int32_t> selected(K.ncol,-1);
		for(size_t i=0;i<pivot.size();++i)selected[pivot[i]]=i;
		sparse_mat<ulong,int32_t> square(K.nrow,K.nrow);
		for(size_t i=0;i<K.nrow;++i)for(auto [j,x]:K[i])if(selected[j]>=0)square[i].push_back(selected[j],x);
		auto chart=sparse_mat_inverse(square,F,opt);inverse=sparse_mat<ulong,int32_t>(K.ncol,K.nrow);
		for(size_t i=0;i<pivot.size();++i)inverse[pivot[i]]=std::move(chart[i]);
	}
	require(nmhv_sew_matrices_equal(sparse_mat_mul(K,inverse,F,&opt->pool),joint_eet_identity<ulong,int32_t>(K.nrow)),"Invalid modular recurrence chart");return inverse;
}
inline void check_moving_rank(const sparse_mat<ulong,int32_t>& K,const field_t& F,rref_option_t opt) {
	std::vector<size_t> count(K.ncol);for(const auto& row:K.rows)for(auto [j,x]:row)++count[j];
	bool private_chart=true;for(const auto& row:K.rows){bool found=false;for(auto [j,x]:row)if(count[j]==1){found=true;break;}if(!found){private_chart=false;break;}}
	if(private_chart)return;
	auto reduced=K;auto groups=sparse_mat_rref(reduced,F,opt);size_t rank=0;for(const auto& g:groups)rank+=g.size();
	if(rank!=K.nrow)throw bad_prime("Moving recurrence loses row rank modulo prime");
}
inline size_t width(const irr::copies& first,const irr::copies& last) {
	size_t n=0;for(size_t b=0;b<6;++b){require(!last[b]||first[b]<=size_t(INT32_MAX)/last[b],"Block exceeds int32 range");
		n+=first[b]*last[b];require(n<=INT32_MAX,"Packed sewing exceeds int32 range");}return n;
}
inline sparse_mat<ulong,int32_t> shift_one(const sparse_mat<ulong,int32_t>& input,size_t h,size_t k,
	const fs::path& fec,const fs::path& lec,const irr::channels_t& cg,const irr::seam_recoupling_t& recouple,const field_t& F,rref_option_t opt,irr::recurrence_store* store=nullptr) {
	require(h>1&&k>0,"Left shift must retain positive first weight");auto* pool=&opt->pool;
	auto first=irr::chain_copies(fec,h,true),old=irr::chain_copies(fec,h-1,true),last=irr::chain_copies(lec,k,false),next=irr::chain_copies(lec,k+1,false);
	require(input.ncol==width(first,last),"Source cache column layout mismatch");
	irr::conversion_progress("phase=cache_shift_read h="+std::to_string(h)+" k="+std::to_string(k));
	auto moving_ptr=store?store->modular(fec,h,true,cg,F,pool):std::make_shared<const irr::recurrence_store::residues>(recurrence(fec,h,true,cg,F,pool));
	auto target_ptr=store?store->modular(lec,k+1,false,cg,F,pool):std::make_shared<const irr::recurrence_store::residues>(recurrence(lec,k+1,false,cg,F,pool));
	const auto& moving=*moving_ptr;const auto& target=*target_ptr;
	for(const auto& block:moving)check_moving_rank(block,F,opt);
	std::array<sparse_mat<ulong,int32_t>,6> inverse;
	std::array<std::vector<std::tuple<size_t,size_t,size_t>>,6> decode;
	std::array<std::unique_ptr<irr_parallel::grouped_entries>,6> index;
	std::array<irr::copies,6> base{};irr::copies oldbase{},offset{},outoffset{},columns{};size_t groups=0,n=0,m=0;
	for(size_t a=0;a<6;++a){oldbase[a]=groups;groups+=old[a];offset[a]=n;n+=first[a]*last[a];outoffset[a]=m;m+=old[a]*next[a];
		inverse[a]=right_inverse(target[a],F,opt);
		for(size_t b=0;b<6;++b){base[a][b]=columns[a];columns[a]+=last[b]*cg.channels[a][b].size();}
		require(columns[a]==target[a].ncol,"Target channel ordering mismatch");}
	std::array<std::vector<std::vector<ulong>>,6> values;
	for(size_t b=0;b<6;++b){values[b].resize(recouple.routes[b].size());
		for(size_t g=0;g<values[b].size();++g)for(const auto& route:recouple.routes[b][g])values[b][g].push_back(residue(route.value,F));
		std::vector<size_t> group;
		for(size_t a=0;a<6;++a)for(size_t c=0;c<old[a];++c)for(size_t g:cg.channels[b][a]){decode[b].push_back({a,c,g});group.push_back(oldbase[a]+c);}
		index[b]=std::make_unique<irr_parallel::grouped_entries>(moving[b],group,groups);}
	sparse_mat<ulong,int32_t> out(input.nrow,width(old,next));
	for(size_t sol=0;sol<input.nrow;++sol){
		std::array<std::vector<std::vector<std::pair<size_t,ulong>>>,6> source;
		for(size_t b=0;b<6;++b)source[b].resize(first[b]);
		for(auto [col,x]:input[sol]){size_t b=0;while(b<6&&size_t(col)>=offset[b]+first[b]*last[b])++b;
			require(b<6,"Invalid packed cache column");size_t ij=col-offset[b];source[b][ij/last[b]].push_back({ij%last[b],x});}
		for(size_t a=0;a<6;++a){sparse_mat<ulong,int32_t> ambient(old[a],columns[a]);
			irr_parallel::rows(old[a],pool,[&](size_t copy){std::map<int32_t,ulong> sums;
				for(size_t b=0;b<6;++b){const auto& idx=*index[b];size_t group=oldbase[a]+copy;
					for(size_t p=idx.offsets[group];p<idx.offsets[group+1];++p){auto [r,pos]=idx.entries[p];auto [block,c,g]=decode[b][moving[b][r](pos)];
						for(auto [stationary,x]:source[b][r]){auto xy=nmod_mul(x,moving[b][r][pos],F.mod);
							for(size_t v=0;v<recouple.routes[b][g].size();++v){const auto& route=recouple.routes[b][g][v];require(route.block==a,"Recoupling irrep mismatch");
								auto col=int32_t(base[a][b]+stationary*cg.channels[a][b].size()+route.channel);
								sums[col]=nmod_add(sums[col],nmod_mul(xy,values[b][g][v],F.mod),F.mod);}}}}
				for(auto [j,x]:sums)if(x)ambient[copy].push_back(j,x);
			});
			auto coefficients=sparse_mat_mul(ambient,inverse[a],F,pool);
			require(nmhv_sew_matrices_equal(sparse_mat_mul(coefficients,target[a],F,pool),ambient),"Modular shifted tensor is outside target recurrence span");
			for(size_t i=0;i<coefficients.nrow;++i)for(auto [j,x]:coefficients[i])out[sol].push_back(outoffset[a]+i*next[a]+j,x);
		}
		irr::conversion_progress("phase=cache_shift_solution row="+std::to_string(sol+1)+"/"+std::to_string(input.nrow),true);
	}return out;
}

inline uint64_t fingerprint(const fs::path& root){return std::stoull(hex_sha256_file(root/"input_identity.txt").substr(0,16),nullptr,16);}
inline nmhv_sew_cache_header_t header(const fs::path& path,size_t h,size_t k,uint64_t id) {
	irr_stage::no_links(path);irr_stage::no_links(path.string()+".sha256");
	require(irr_stage::text(path.string()+".sha256")==hex_sha256_file(path)+"\n","Cache checksum mismatch: "+path.string());
	std::ifstream in(path,std::ios::binary);nmhv_sew_cache_header_t v{};in.read(reinterpret_cast<char*>(&v),sizeof(v));
	require(bool(in)&&nmhv_sew_cache_header_matches(v,(h+k)/2,h,k,v.nrow,v.ncol,v.prime,id),"Cache header/identity mismatch");
	require(v.nrow<=INT32_MAX&&v.ncol<=INT32_MAX&&v.prime>2&&n_is_prime(v.prime),"Invalid cache shape/prime");
	require(fs::file_size(path)==sizeof(v)+v.nrow*v.ncol*sizeof(ulong),"Cache size mismatch");return v;
}
inline sparse_mat<ulong,int32_t> read(const fs::path& path,const nmhv_sew_cache_header_t& v) {
	std::ifstream in(path,std::ios::binary);in.seekg(sizeof(v));sparse_mat<ulong,int32_t> out(v.nrow,v.ncol);std::vector<ulong> row(v.ncol);
	for(size_t i=0;i<v.nrow;++i){in.read(reinterpret_cast<char*>(row.data()),row.size()*sizeof(ulong));require(bool(in),"Truncated cache");
		for(size_t j=0;j<v.ncol;++j){require(row[j]<v.prime,"Invalid modular residue");if(row[j])out[i].push_back(j,row[j]);}}return out;
}
inline int run(int argc,char** argv) {
	require(argc>=7,"irr_convert shift-cache MHV|NMHV H K INPUT OUTPUT --head H [--fec-root DIR] [--lec-root DIR] [--kind homogeneous|affine|unique] [--reconstruct]");
	std::string helicity=argv[2],kind;size_t h=irr_io::number(argv[3]),k=irr_io::number(argv[4]),head=h;bool reconstruct=false;
	require((helicity=="MHV"||helicity=="NMHV")&&h&&k&&(h+k)%2==0,"Invalid helicity/split");
	fs::path input=fs::absolute(argv[5]).lexically_normal(),output=fs::absolute(argv[6]).lexically_normal(),fec="output/FEC",lec=fs::path("output/LEC")/helicity;
	for(int i=7;i<argc;++i){std::string key=argv[i];if(key=="--reconstruct"){reconstruct=true;continue;}require(i+1<argc,"Missing option value");std::string value=argv[++i];
		if(key=="--head")head=irr_io::number(value);else if(key=="--fec-root")fec=value;else if(key=="--lec-root")lec=value;else if(key=="--kind")kind=value;else throw std::runtime_error("Unknown option "+key);}
	require(head>0&&head<h,"shift-cache currently supports left shifts only");size_t tail=h+k-head;
	irr_stage::no_links(input);irr_stage::no_links(output);require(input!=output,"Source and target must differ");
	for(auto p=output;!p.empty();p=p.parent_path()){require(p!=input,"Target cannot be inside source");if(p==p.parent_path())break;}
	require(fs::is_directory(input/"cache"),"Missing source cache directory");irr_stage::no_links(input/"cache");
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=false;irr_run::threads(opt);irr::channels_t cg(&opt->pool);
	auto terminal=nmhv_sew_read_tensor<rat_t,int32_t>(lec/"LEC_0_components.wxf",field_t(FIELD_QQ),&opt->pool);
	require(terminal.rank()==2&&terminal.dim(0)==size_t(helicity=="MHV"?1:5),"LEC helicity mismatch");
	auto source_basis=irr_stage::basis_identity(fec,lec,h,k),basis=irr_stage::basis_identity(fec,lec,head,tail);
	irr_stage::verify_basis(input,source_basis,h,k);
	std::string declared;size_t declared_rows=SIZE_MAX;
	if(fs::exists(input/"object.tsv")){auto meta=irr_stage::read(input/"object.tsv");declared=meta.at("kind");declared_rows=irr_io::number(meta.at("directions"))+(declared!="homogeneous");
		require(irr_io::number(meta.at("h"))==h&&irr_io::number(meta.at("k"))==k,"Source object split mismatch");}
	else if(fs::exists(input/"stage.tsv")){auto stage=irr_stage::read(input/"stage.tsv");require(stage.at("helicity")==helicity&&irr_io::number(stage.at("h"))==h&&irr_io::number(stage.at("k"))==k,"Source stage mismatch");
		if(stage.at("operator")=="collinear"||stage.at("operator")=="mrk"||stage.at("operator")=="ope")declared="affine";
		else if(stage.at("operator")=="shift-cache"){declared=stage.at("kind");declared_rows=irr_io::number(stage.at("rows"));}}
	if(kind.empty())kind=declared;else require(declared.empty()||kind==declared,"Explicit kind conflicts with source metadata");
	require(kind=="homogeneous"||kind=="affine"||kind=="unique","Cannot infer origin convention; specify --kind for a homogeneous bulk cache");
	std::map<ulong,fs::path> sources;size_t rows=SIZE_MAX,cols=width(irr::chain_copies(fec,h,true),irr::chain_copies(lec,k,false));auto sid=fingerprint(input);
	for(const auto& e:fs::directory_iterator(input/"cache"))if(e.path().extension()==".bin"){
		if(!fs::exists(e.path().string()+".sha256")){std::cout<<"SKIP unfinished cache "<<e.path()<<std::endl;continue;}
		auto v=header(e.path(),h,k,sid);require(v.ncol==cols,"Source cache columns differ from recurrence layout");if(rows==SIZE_MAX)rows=v.nrow;
		require(v.nrow==rows,"Source cache solution counts differ");require(sources.emplace(v.prime,e.path()).second,"Duplicate source prime");}
	require(!sources.empty(),"No completed source primes");require(declared_rows==SIZE_MAX||declared_rows==rows,"Source object/cache row mismatch");
	if(kind!="homogeneous"){require(rows>0,"Affine cache must include origin");if(kind=="unique")require(rows==1,"Unique cache must have one row");kind=rows==1?"unique":"affine";}
	// Bind the coordinate map, not the current number of primes: later arrivals
	// can be resumed without changing the output's reconstruction identity.
	std::string identity="irr-shift-cache-v1\n"+hex_sha256_file(input/"input_identity.txt")+"\n"+helicity+" "+kind+" "+std::to_string(h)+" "+std::to_string(k)+" "+std::to_string(head)+" "+std::to_string(tail)+" "+std::to_string(rows)+"\n";
	for(const auto& fields:{source_basis,basis})for(auto [key,value]:fields)identity+=key+"\t"+value+"\n";
	if(fs::exists(output))require(fs::exists(output/"input_identity.txt")&&irr_stage::text(output/"input_identity.txt")==identity,"Unmanaged output or different shift identity");
	fs::create_directories(output/"cache");irr_stage::no_links(output/"cache");
	for(const auto& dir:{output,output/"cache"})for(const auto& e:fs::directory_iterator(dir))irr_stage::no_links(e.path());
	irr_stage::lock lock(output/".active.lock");irr_stage::fixed_text(output/"input_identity.txt",identity);
	irr_stage::write(output/"basis_identity.tsv",basis);auto id=fingerprint(output);
	irr_stage::write(output/"stage.tsv",{{"schema","irr-shift-cache-v1"},{"operator","shift-cache"},{"history",irr_stage::history(input)},{"helicity",helicity},{"h",std::to_string(head)},{"k",std::to_string(tail)},{"kind",kind},{"rows",std::to_string(rows)},{"parent",input.generic_string()}});
	irr::seam_recoupling_t recouple(cg,opt);std::vector<fs::path> paths;std::vector<ulong> primes;
	for(const auto& [p,source]:sources){auto file=output/"cache"/("p"+std::to_string(p)+".bin");auto provenance=output/("source_p"+std::to_string(p)+".sha256");
		if(fs::exists(file)){auto v=header(file,head,tail,id);require(v.nrow==rows&&v.ncol==width(irr::chain_copies(fec,head,true),irr::chain_copies(lec,tail,false)),"Shifted cache shape mismatch");
			require(irr_stage::text(provenance)==hex_sha256_file(source)+"\n","Source prime changed or missing conversion receipt");std::cout<<"RESUME prime="<<p<<std::endl;
		}else try{auto mat=read(source,header(source,h,k,sid));field_t F(FIELD_Fp,p);
			for(size_t current=h;current>head;--current)mat=shift_one(mat,current,h+k-current,fec,lec,cg,recouple,F,opt);
			irr_stage::fixed_text(provenance,hex_sha256_file(source)+"\n");nmhv_sew_write_modular_basis(file,mat,p,(h+k)/2,head,tail,id);
			std::cout<<"SHIFT_CACHE prime="<<p<<" rows="<<mat.nrow<<" columns="<<mat.ncol<<" nnz="<<mat.nnz()<<" verified=PASS"<<std::endl;
		}catch(const bad_prime& e){irr_stage::write(output/("bad_prime_"+std::to_string(p)+".tsv"),{{"reason",e.what()},{"source_sha256",hex_sha256_file(source)}});std::cout<<"BAD_PRIME "<<p<<" "<<e.what()<<std::endl;continue;}
		paths.push_back(file);primes.push_back(p);
	}
	if(fs::exists(output/"qq_complete.tsv")){irr_stage::verify_complete(output);irr_stage::complete(output);return 0;}
	if(fs::exists(output/"qq_cache_safe.tsv")){irr_stage::verify_chain(output);irr_stage::complete_chain(output,false);return 0;}
	bool success=false;size_t used=0;int_t modulus=1;sparse_mat<rat_t,int32_t> result;
	if(reconstruct)for(size_t n=1;n+2<=paths.size();++n){modulus*=primes[n-1];std::cout<<"RECONSTRUCT crt_primes="<<n<<" validation_primes=2"<<std::endl;
		if(nmhv_sew_reconstruct_basis(result,{paths.begin(),paths.begin()+n},{primes.begin(),primes.begin()+n},{paths.begin()+n,paths.begin()+n+2},{primes.begin()+n,primes.begin()+n+2},(h+k)/2,head,tail,rows,width(irr::chain_copies(fec,head,true),irr::chain_copies(lec,tail,false)),id,&opt->pool)){success=true;used=n;break;}}
	std::ofstream manifest(output/"manifest.json");manifest<<"{\"schema\":\"irr-shift-cache-v1\",\"complete\":"<<(success?"true":"false")<<",\"h\":"<<head<<",\"k\":"<<tail<<",\"kind\":"<<irr_io::quote(kind)<<",\"crt_primes\":"<<used<<",\"validation_primes\":"<<(success?2:0)<<",\"preserves_solution_coordinates\":true,\"final_presentation_rref\":false}\n";manifest.close();
	if(success){irr::sewing_object_t obj;obj.h=head;obj.k=tail;obj.kind=kind;obj.directions=rows-(kind!="homogeneous");obj.first=irr::chain_copies(fec,head,true);obj.last=irr::chain_copies(lec,tail,false);size_t off=0;
		for(size_t b=0;b<6;++b){obj.S[b]=sparse_mat<rat_t,int32_t>(rows,obj.first[b]*obj.last[b]);for(size_t i=0;i<rows;++i)for(auto [j,x]:result[i])if(size_t(j)>=off&&size_t(j)<off+obj.S[b].ncol)obj.S[b][i].push_back(j-off,x);off+=obj.S[b].ncol;}
		irr_io::write_reduced(output,obj,&opt->pool);
		auto check=irr_io::read_object(output,head,tail,fec,lec,&opt->pool);for(size_t b=0;b<6;++b)require(nmhv_sew_matrices_equal(check.S[b],obj.S[b]),"Export readback mismatch");
		irr_stage::fields receipt;for(size_t i=0;i<used+2;++i)receipt["p"+std::to_string(primes[i])]=i<used?"CRT":"validation";irr_stage::write(output/"reconstruction_primes.tsv",receipt);
		irr_stage::write(output/"summary.tsv",{{"crt_primes",std::to_string(used)},{"validation_primes","2"},{"modulus_bits",std::to_string(modulus.sizeinbase(2))},{"directions",std::to_string(obj.directions)}});
		if(kind!="homogeneous")irr_stage::complete(output);else {std::ofstream hashes(output/"checksums.sha256");for(size_t b=0;b<6;++b){auto file=irr::tensor_path(output/(irr_io::stem(head,tail,b)+".wxf"));hashes<<hex_sha256_file(file)<<"  "<<file.filename().string()<<'\n';}hashes.close();irr_stage::complete_chain(output,false);}
	}
	std::cout<<(success?"QQ_VERIFIED":"CACHE_ONLY")<<" completed_primes="<<paths.size()<<" (source caches retained; no new primes evaluated)"<<std::endl;return 0;
}
}
