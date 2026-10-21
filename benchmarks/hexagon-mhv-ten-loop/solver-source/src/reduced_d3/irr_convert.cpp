#include "irr_object.hpp"
#include "irr_runtime.hpp"
#include "irr_stage.hpp"
#include "irr_shift_cache.hpp"
#include <iomanip>

namespace fs=std::filesystem;
const field_t QQ(FIELD_QQ);
using irr::require; using namespace irr_io;

struct names_t {
	std::string first="hexagon_coproduct_entries_first_{w}.wxf";
	std::string last="hexagon_{helicity}_coproduct_entries_last_{w}.wxf";
	std::string sew="hexagon_{helicity}_sewing_{h}p{k}.wxf";
	std::string origin="hexagon_{helicity}_sewing_{h}p{k}_origin.wxf";
};
std::string filename(std::string pattern,const std::string& helicity,size_t w,size_t h,size_t k) {
	std::map<std::string,std::string> values={{"w",std::to_string(w)},{"h",std::to_string(h)},{"k",std::to_string(k)},
		{"weight",std::to_string(h+k)},{"loop",std::to_string((h+k)/2)},{"helicity",helicity=="MHV"?"mhv":"nmhv"}};
	std::string out;
	for(size_t i=0;i<pattern.size();) {
		if(pattern[i]!='{'){require(pattern[i]!='}',"Unmatched filename brace");out+=pattern[i++];continue;}
		auto end=pattern.find('}',i+1);require(end!=std::string::npos,"Unmatched filename brace");auto key=pattern.substr(i+1,end-i-1);
		require(values.count(key)!=0,"Unknown filename placeholder: "+key);out+=values.at(key);i=end+1;
	}
	require(!out.empty()&&out!="."&&out!=".."&&out.find_first_of("/\\:<>\"|?*\r\n\t")==std::string::npos&&out.back()!='.'&&out.back()!=' ',"Release filenames must be plain filenames, not paths");
	require(fs::path(out).extension()==".wxf","Tensor filename must end in .wxf");return out;
}
sparse_mat<rat_t,int32_t> read_coefficients(const fs::path& path,size_t n,thread_pool* pool) {
	auto t=nmhv_sew_read_tensor<rat_t,int32_t>(path,QQ,pool);require(t.rank()==2&&t.dim(1)==n,"Coefficient matrix must have one column per input direction");return joint_eet_flatten_tensor(t,pool);
}
void select_object(irr::sewing_object_t& out,const fs::path& directions,const fs::path& origin,bool unique,thread_pool* pool) {
	if(directions.empty()&&origin.empty()&&!unique)return;
	require(out.kind!="unique","Cannot select directions from a unique object");
	require(origin.empty()||out.kind=="homogeneous","New origin coefficients require a homogeneous source");
	require(!unique||(out.kind=="homogeneous"&&origin.empty()),"Unique and affine selections are mutually exclusive");
	size_t old=out.directions;auto D=directions.empty()?joint_eet_identity<rat_t,int32_t>(old):
		(directions=="none"?sparse_mat<rat_t,int32_t>(0,old):read_coefficients(directions,old,pool));
	sparse_mat<rat_t,int32_t> O;if(!origin.empty()){O=read_coefficients(origin,old,pool);require(O.nrow==1,"Origin coefficients must be a single row");}
	if(unique)require(D.nrow==1,"Unique selection must have exactly one row");
	for(size_t b=0;b<6;++b) {
		sparse_mat<rat_t,int32_t> source(old,out.S[b].ncol);for(size_t r=0;r<old;++r)source[r]=out.S[b][r];
		auto selected=sparse_mat_mul(D,source,QQ,pool);
		if(out.kind=="affine"){sparse_mat<rat_t,int32_t> keep(1,source.ncol);keep[0]=out.S[b][old];selected.append(std::move(keep));}
		else if(!origin.empty())selected.append(sparse_mat_mul(O,source,QQ,pool));
		out.S[b]=std::move(selected);
	}
	out.directions=unique?0:D.nrow;if(unique)out.kind="unique";else if(!origin.empty())out.kind="affine";
}

void export_recurrence(const fs::path& root,const fs::path& chain,size_t w,bool first,const irr::channels_t& cg,
	const std::string& name,std::ostream& table,std::ostream& stats,rref_option_t opt) {
	auto* pool=&opt->pool;auto old=irr::chain_copies(chain,w-1,first),m=irr::chain_copies(chain,w,first);size_t nold=0,offset=0;
	for(size_t b=0;b<6;++b)nold+=old[b]*irr::specs()[b].irrep_dimension;
	auto K=irr::read_recurrence(chain,w,first,cg,pool);sparse_mat<rat_t,int32_t> full(0,nold*9);
	for(size_t b=0;b<6;++b) {
		full.append(irr::expand(cg,b,old,K[b]));const auto& s=irr::specs()[b];
		table<<w<<'\t'<<s.name<<'\t'<<s.parity<<'\t'<<(s.irrep==fec_irrep_t::A1?"A1":s.irrep==fec_irrep_t::A2?"A2":"E")<<'\t'<<s.irrep_dimension<<'\t'<<m[b]<<'\t'<<offset<<'\t'<<m[b]*s.irrep_dimension<<'\n';offset+=m[b]*s.irrep_dimension;
	}
	sparse_tensor<rat_t,int32_t,SPARSE_COO> tensor(std::move(full),pool);tensor.reshape({offset,nold,9});
	// The no-pool transpose already sorts in place; a second parallel sort
	// would allocate another full index/value buffer without changing the order.
	tensor.transpose_replace(first?std::vector<size_t>{1,2,0}:std::vector<size_t>{0,2,1});
	if(first&&w==1)tensor.reshape({9,offset});
	if(!first&&w==1) {
		tensor.flatten({{0},{1,2}});auto rows=tensor.to_sparse_mat(pool);
		auto ct=nmhv_sew_read_tensor<rat_t,int32_t>(chain/"LEC_0_components.wxf",QQ,pool);auto components=joint_eet_flatten_tensor(ct,pool);
		rows=nmhv_sew_transform_product_basis(rows,joint_eet_identity<rat_t,int32_t>(9),components,QQ,pool);
		auto physical=nmhv_sew_read_tensor<rat_t,int32_t>(chain/"LEC_1_physical.wxf",QQ,pool);
		require(nmhv_sew_matrices_equal(rows,joint_eet_flatten_tensor(physical,pool)),"Expanded LEC1 differs from the normalized physical seed");
		write_tensor(root/name,shaped(std::move(rows),{offset,9,components.ncol},pool),&stats);
	}else write_tensor(root/name,sparse_tensor<rat_t,int32_t,SPARSE_CSR>(std::move(tensor),pool),&stats);
}

void export_sewing(const fs::path& root,const std::string& name,const irr::sewing_object_t& obj,
	size_t begin,size_t count,bool solution_axis,std::ostream& stats,thread_pool* pool) {
	irr::copies fo{},lo{};size_t fd=0,ld=0;
	for(size_t b=0;b<6;++b){fo[b]=fd;lo[b]=ld;fd+=obj.first[b]*irr::specs()[b].irrep_dimension;ld+=obj.last[b]*irr::specs()[b].irrep_dimension;}
	require(fd*count<=size_t(INT32_MAX)&&ld<=size_t(INT32_MAX),"Release sewing exceeds int32 range");
	sparse_mat<rat_t,int32_t> matrix(count*fd,ld);
	for(size_t b=0;b<6;++b)for(size_t row=0;row<count;++row)for(auto [ij,x]:obj.S[b][begin+row])
		for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)
			matrix[row*fd+fo[b]+s*obj.first[b]+ij/obj.last[b]].push_back(int32_t(lo[b]+s*obj.last[b]+ij%obj.last[b]),x*(s?rat_t(3):rat_t(1)));
	std::vector<size_t> dims{fd,ld};if(solution_axis)dims.insert(dims.begin(),count);
	write_tensor(root/name,shaped(std::move(matrix),dims,pool),&stats);
}

void write_wl(const fs::path& root,const std::string& helicity,const irr::sewing_object_t& obj,
	const std::vector<std::string>& first,const std::vector<std::string>& last,const std::string& sew,const std::string& origin) {
	std::ofstream wl(root/"hexagon_symbols.wl");std::string lc=helicity+"LCE",sc=helicity+"SEW",oc=helicity+"Origin";
	wl<<"(* Self-contained WXF access. No SparseRREF package is required. *)\nBeginPackage[\"HexagonSymbolRelease`\"];\n"
		<<"FCE::usage=\"FCE[w] gives {old,letter,new}; FCE[1] omits the trivial old index.\";\n"
		<<lc<<"::usage=\""<<lc<<"[w] gives {new,letter,old}; the physical terminal index is retained.\";\n"
		<<sc<<"::usage=\""<<sc<<"[L] gives "<<(obj.kind=="unique"?"{first,last}":"{solution,first,last}")<<".\";\n"
		<<oc<<"::usage=\""<<oc<<"[L] gives the affine representative {first,last}, when present.\";\n";
	std::vector<std::pair<std::string,size_t>> components=helicity=="MHV"?std::vector<std::pair<std::string,size_t>>{{"MHV",1}}:std::vector<std::pair<std::string,size_t>>{{"NMHVE",1},{"NMHVEt",4}};
	for(const auto& [name,component]:components) {
		wl<<name<<"coeff::usage=\""<<name<<"coeff"<<(obj.kind=="unique"?"[letters...]":"[solution, {letters...}]")<<" returns a word coefficient.\";\n";
		if(obj.kind=="affine")wl<<name<<"OriginCoeff::usage=\""<<name<<"OriginCoeff[letters...] returns a coefficient of the affine representative, not of a homogeneous direction.\";\n";
	}
	wl<<"Begin[\"`Private`\"];\nroot=DirectoryName[$InputFileName];\nread[file_]:=Import[FileNameJoin[{root,file}],\"WXF\"];\n"
		<<"bad::arg=\"Invalid weight, loop, solution, component or symbol word.\";\n";
	for(size_t w=1;w<=obj.h;++w)wl<<"FCE["<<w<<"]:=FCE["<<w<<"]=read["<<quote(first[w-1])<<"];\n";
	for(size_t w=1;w<=obj.k;++w)wl<<lc<<"["<<w<<"]:="<<lc<<"["<<w<<"]=read["<<quote(last[w-1])<<"];\n";
	size_t fd=0,ld=0;for(size_t b=0;b<6;++b){fd+=obj.first[b]*irr::specs()[b].irrep_dimension;ld+=obj.last[b]*irr::specs()[b].irrep_dimension;}
	wl<<sc<<"["<<(obj.h+obj.k)/2<<"] := "<<sc<<"["<<(obj.h+obj.k)/2<<"] = ";
	if(obj.kind!="unique"&&!obj.directions)wl<<"{}; (* No homogeneous directions. *)\n";else wl<<"read["<<quote(sew)<<"];\n";
	if(obj.kind=="affine")wl<<oc<<"["<<(obj.h+obj.k)/2<<"]:="<<oc<<"["<<(obj.h+obj.k)/2<<"]=read["<<quote(origin)<<"];\n";
	wl<<"FCE[___]:=(Message[bad::arg];$Failed);\n"<<lc<<"[___]:=(Message[bad::arg];$Failed);\n"<<sc<<"[___]:=(Message[bad::arg];$Failed);\n"<<oc<<"[___]:=(Message[bad::arg];$Failed);\n"
		<<"coeff[word_List,s_Integer,r_Integer,origin_]:=Module[{f,t,m},If[Length[word]=!="<<obj.h+obj.k<<"||!VectorQ[word,IntegerQ[#]&&1<=#<=9&]||!MemberQ["<<(helicity=="MHV"?"{1}":"{1,4}")<<",r],Message[bad::arg];Return[$Failed]];\n"
		<<"If[!origin&&!(1<=s<="<<(obj.kind=="unique"?1:obj.directions)<<"),Message[bad::arg];Return[$Failed]];\n"
		<<"m=If[origin,"<<oc<<"["<<(obj.h+obj.k)/2<<"],"<<sc<<"["<<(obj.h+obj.k)/2<<"]"<<(obj.kind=="unique"?"":"[[s]]")<<"];If[m===$Failed,Return[$Failed]];\n"
		<<"f=FCE[1][[word[[1]],All]];Do[f=f.FCE[i][[All,word[[i]],All]],{i,2,"<<obj.h<<"}];\n"
		<<"t="<<lc<<"[1][[All,Last[word],r]];Do[t="<<lc<<"[i][[All,word[[Length[word]-i+1]],All]].t,{i,2,"<<obj.k<<"}];f.m.t];\n";
	for(const auto& [name,component]:components) {
		if(obj.kind=="unique")wl<<name<<"coeff[word__Integer]:=coeff[{word},1,"<<component<<",False];\n";
		else wl<<name<<"coeff[s_Integer,word_List]:=coeff[word,s,"<<component<<",False];\n";
		wl<<name<<"coeff[___]:=(Message[bad::arg];$Failed);\n";
		if(obj.kind=="affine")wl<<name<<"OriginCoeff[word__Integer]:=coeff[{word},1,"<<component<<",True];\n"<<name<<"OriginCoeff[___]:=(Message[bad::arg];$Failed);\n";
	}wl<<"End[];\nEndPackage[];\n";
}

void release(const fs::path& root,const std::string& helicity,const irr::sewing_object_t& obj,const names_t& names,
	const fs::path& fec,const fs::path& lec,const irr::channels_t& cg,rref_option_t opt) {
	std::vector<std::string> first,last;std::set<std::string> files;
	auto claim=[&](const std::string& name){std::string key=name;std::transform(key.begin(),key.end(),key.begin(),[](unsigned char c){return char(std::tolower(c));});require(files.insert(key).second,"Colliding release filenames: "+name);};
	for(size_t w=1;w<=obj.h;++w){first.push_back(filename(names.first,helicity,w,obj.h,obj.k));claim(first.back());}
	for(size_t w=1;w<=obj.k;++w){last.push_back(filename(names.last,helicity,w,obj.h,obj.k));claim(last.back());}
	auto sew=filename(names.sew,helicity,obj.h+obj.k,obj.h,obj.k),origin=filename(names.origin,helicity,obj.h+obj.k,obj.h,obj.k);
	if(obj.kind=="unique"||obj.directions)claim(sew);if(obj.kind=="affine")claim(origin);
	fs::create_directories(root);std::ofstream ft(root/"FCE_blocks.tsv"),lt(root/(helicity+"_LCE_blocks.tsv")),stats(root/"TENSORS.tsv");
	for(auto* out:{&ft,&lt})*out<<"weight\tblock\tparity\tirrep\tirrep_dimension\tmultiplicity\toffset\tdimension\n";
	stats<<"file\tdimensions\tnnz\theight_bits\tbytes\n";
	for(size_t w=1;w<=obj.h;++w)export_recurrence(root,fec,w,true,cg,first[w-1],ft,stats,opt);
	for(size_t w=1;w<=obj.k;++w)export_recurrence(root,lec,w,false,cg,last[w-1],lt,stats,opt);
	if(obj.kind=="unique"||obj.directions)export_sewing(root,sew,obj,0,obj.kind=="unique"?1:obj.directions,obj.kind!="unique",stats,&opt->pool);
	if(obj.kind=="affine")export_sewing(root,origin,obj,obj.directions,1,false,stats,&opt->pool);
	write_wl(root,helicity,obj,first,last,sew,origin);
	std::ofstream readme(root/"README.md");
	readme<<"# Hexagon "<<helicity<<" Symbol Data at Weight "<<obj.h+obj.k<<"\n\n"
		<<"This self-contained directory stores a "<<obj.kind<<" object at the split "<<obj.h<<" + "<<obj.k<<". It does not assert that bulk constraints alone determine a physical amplitude.\n\n"
		<<"The alphabet is (a,b,c,m_u,m_v,m_w,y_u,y_v,y_w), with a=u/(v w), b=v/(w u), c=w/(u v), m_u=(1-u)/u and cyclic variants; the y letters use the usual hexagon square root. In particular a,b,c are not their square roots.\n\n"
		<<"FCE tensors have indices {old,letter,new}; FCE1 omits the trivial old index. LCE tensors have indices {new,letter,old}. LCE1 retains its "<<(helicity=="MHV"?"dimension-one terminal":"five-component physical terminal, with components 1 and 4 selecting E and Et")<<". The physical LEC1 is integer and primitive per full irrep copy, not separately per E component.\n\n"
		<<"Homogeneous directions have {solution,first_basis,last_basis}, including when there is just one direction. An affine representative or an explicitly selected unique object has {first_basis,last_basis}. An affine object means origin + arbitrary linear combination of directions, not an extra solution vector. Empty affine directions are omitted.\n\n"
		<<"There are "<<obj.directions<<" homogeneous directions"<<(obj.kind=="unique"?" (this is a unique object)":"")<<". All blocks share the same solution ordering; conversion preserves exact coefficients without renormalization.\n\n"
		<<"Load `Get[\"hexagon_symbols.wl\"]`. Access `FCE[w]`, `"<<helicity<<"LCE[w]`, `"<<helicity<<"SEW["<<(obj.h+obj.k)/2<<"]`"<<(obj.kind=="affine"?" and `"+helicity+"Origin["+std::to_string((obj.h+obj.k)/2)+"]`":"")<<". The script uses only built-in SparseArray/WXF support and resolves all data relative to its own directory.\n\n"
		<<"Word-coefficient access uses "<<(helicity=="MHV"?"`MHVcoeff`":"`NMHVEcoeff` and `NMHVEtcoeff`")<<" with "<<(obj.kind=="unique"?"integer letter arguments":"arguments `solution, {letters...}`")<<". Affine origin coefficients use the corresponding `OriginCoeff[letters...]` functions. Prefer these contractions over expanding entire high-weight symbols.\n\n"
		<<"Contract FCE1 ... FCE_h, the selected sewing matrix, LCE_k ... LCE1, then select the physical terminal. Metadata list only files present in this directory; `inputs.sha256` is an audit trail, not a runtime dependency.\n";
}

int main(int argc,char** argv){try{
	if(argc>1&&std::string(argv[1])=="shift-cache")return irr_shift_cache::run(argc,argv);
	require(argc>=7,"irr_convert shift|export MHV|NMHV H K INPUT OUTPUT [--head H] [--directions WXF] [--origin WXF] [--unique] [--fec-root DIR] [--lec-root DIR] [--first-name PATTERN] [--last-name PATTERN] [--sew-name PATTERN] [--origin-name PATTERN]");
	std::string mode=argv[1],helicity=argv[2];require(mode=="shift"||mode=="export","Unknown operation");require(helicity=="MHV"||helicity=="NMHV","Unknown helicity");
	size_t h=number(argv[3]),k=number(argv[4]),target=h;require(h&&k&&(h+k)%2==0,"Require h,k>=1 and even total weight");
	fs::path input=argv[5],output=argv[6],fec="output/FEC",lec=fs::path("output/LEC")/helicity,directions,origin;
	names_t names;bool unique=false;
	for(int i=7;i<argc;++i){std::string key=argv[i];if(key=="--unique"){unique=true;continue;}require(i+1<argc,"Missing option value");std::string value=argv[++i];
		if(key=="--head")target=number(value);else if(key=="--directions")directions=value;else if(key=="--origin")origin=value;else if(key=="--fec-root")fec=value;else if(key=="--lec-root")lec=value;
		else if(key=="--first-name")names.first=value;else if(key=="--last-name")names.last=value;else if(key=="--sew-name")names.sew=value;else if(key=="--origin-name")names.origin=value;else throw std::runtime_error("Unknown option "+key);}
	require(target&&target<h+k,"Target split must retain positive first/last weights");require(!fs::exists(output),"Output exists; use a new independent directory");
	Flint::set_memory_functions();rref_option_t opt;opt->method=1;opt->is_back_sub=true;opt->verbose=false;irr_run::threads(opt);irr::channels_t cg(&opt->pool);
	auto terminal=nmhv_sew_read_tensor<rat_t,int32_t>(lec/"LEC_0_components.wxf",QQ,&opt->pool);
	require(terminal.rank()==2&&terminal.dim(0)==size_t(helicity=="MHV"?1:5),"LEC helicity mismatch");
	irr_stage::verify_basis(input,irr_stage::basis_identity(fec,lec,h,k),h,k);
	irr::conversion_progress("phase=read_object begin");
	auto obj=read_object(input,h,k,fec,lec,&opt->pool);select_object(obj,directions,origin,unique,&opt->pool);
	irr::conversion_progress("phase=read_object end directions="+std::to_string(obj.directions));
	if(target!=h) {
		irr::seam_recoupling_t recouple(cg,opt);
		while(obj.h!=target){obj=irr::shift_one(obj,obj.h>target,fec,lec,cg,recouple,opt);std::cout<<"SHIFT exact_factorization=PASS h="<<obj.h<<" k="<<obj.k<<" kind="<<obj.kind<<" directions="<<obj.directions<<std::endl;}
	}
	irr::conversion_progress("phase=write_output begin");
	if(mode=="export")release(output,helicity,obj,names,fec,lec,cg,opt);
	else{fs::create_directories(output);write_reduced(output,obj,&opt->pool);}
	irr::conversion_progress("phase=write_output end");
	if(!irr_stage::history(input).empty())irr_stage::write(output/"condition_history.tsv",{{"history",irr_stage::history(input)}});
	if(mode=="shift")irr_stage::write(output/"basis_identity.tsv",irr_stage::basis_identity(fec,lec,obj.h,obj.k));
	std::ofstream inputs(output/"inputs.sha256");std::set<fs::path> paths;
	for(const auto& file:fs::directory_iterator(input))if(file.is_regular_file())paths.insert(file.path());
	for(size_t w=1;w<=std::max(h,obj.h);++w)for(size_t b=0;b<6;++b)paths.insert(fec/joint_eet_weight_name(w)/("irrFEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf"));
	for(size_t w=1;w<=std::max(k,obj.k);++w){paths.insert(lec/joint_eet_weight_name(w)/"manifest.json");for(size_t b=0;b<6;++b)paths.insert(lec/joint_eet_weight_name(w)/("irrLEC_"+std::to_string(w)+"_"+irr::specs()[b].name+".wxf"));}
	paths.insert(lec/"LEC_0_components.wxf");paths.insert(lec/"LEC_1_physical.wxf");paths.insert(lec/"common"/"CG.tsv");
	if(!directions.empty()&&directions!="none")paths.insert(directions);if(!origin.empty())paths.insert(origin);
	for(const auto& path:paths){auto actual=irr::tensor_path(path);inputs<<hex_sha256_file(actual)<<"  "<<fs::absolute(actual).generic_string()<<'\n';}inputs.close();
	std::ofstream manifest(output/"manifest.json");manifest<<"{\"schema\":\"hexagon-irr-conversion-v1\",\"complete\":true,\"operation\":"<<quote(mode)<<",\"helicity\":"<<quote(helicity)<<",\"kind\":"<<quote(obj.kind)<<",\"h\":"<<obj.h<<",\"k\":"<<obj.k<<",\"directions\":"<<obj.directions
		<<",\"normalization\":\"physical-copy-primitive-v1\",\"preserves_solution_coordinates\":true,\"files\":[";
	bool comma=false;for(const auto& file:fs::directory_iterator(output))if(file.is_regular_file()&&file.path().filename()!="manifest.json") {if(comma)manifest<<',';manifest<<quote(file.path().filename().string());comma=true;}manifest<<"]}\n";manifest.close();
	std::ofstream hashes(output/"checksums.sha256");for(const auto& file:fs::directory_iterator(output))if(file.is_regular_file()&&file.path().filename()!="checksums.sha256")hashes<<hex_sha256_file(file.path())<<"  "<<file.path().filename().string()<<'\n';
	std::cout<<"COMPLETE "<<mode<<" h="<<obj.h<<" k="<<obj.k<<" kind="<<obj.kind<<" directions="<<obj.directions<<std::endl;
}catch(const std::exception& e){std::cerr<<"Error: "<<e.what()<<std::endl;return 1;}}
