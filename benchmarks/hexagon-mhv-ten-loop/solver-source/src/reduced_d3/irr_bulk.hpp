#pragma once
#include "irr_sew_prereduce.hpp"
#include "shared/irrep_seam_rows.hpp"
#include "irr_runtime.hpp"
#include "irr_stage.hpp"
#include "joint_eet_io.hpp"
#include "nmhv_irrep_sew.hpp"
#include "reduced_d3/platform_resources.hpp"
#include "reduced_d3/hex_hash.hpp"
#include <chrono>
#include "irr_profile.hpp"
#include "irr_recurrence_store.hpp"

namespace irr_bulk {
namespace fs=std::filesystem;
const field_t QQ(FIELD_QQ);
using irr::require;
irr::copies read_copies(const fs::path& root,size_t w,bool first) {
	irr::copies m{};
	if(w==0&&first){m[0]=1;return m;}
	if(first) {
		std::ifstream in(root/joint_eet_weight_name(w)/"multiplicities.tsv");
		require(bool(in),"Missing irrFEC multiplicities");
		for(size_t b=0;b<6;++b){std::string name;in>>name>>m[b];require(bool(in)&&name==irr::specs()[b].name,"irrFEC block order mismatch");}
	}else {
		auto blocks=nmhv_sew_read_irrep_axis(w?root/joint_eet_weight_name(w)/"blocks.tsv":root/"LEC_0_blocks.tsv");
		require(blocks.size()==6,"irrLEC must have six block descriptors");
		for(size_t b=0;b<6;++b){require(blocks[b].spec.name==irr::specs()[b].name,"irrLEC block order mismatch");m[b]=blocks[b].multiplicity;}
	}
	return m;
}
irr::matrices read_chain(const fs::path& root,size_t w,bool first,const irr::copies& m,
	std::ofstream& inputs,thread_pool* pool) {
	irr::matrices K;
	for(size_t b=0;b<6;++b) {
		auto path=root/joint_eet_weight_name(w)/((first?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf");
		auto t=irr::read_tensor(path,pool);require(t.rank()==2,"irr recurrence must be rank two");
		K[b]=joint_eet_flatten_tensor(std::move(t),pool);require(K[b].nrow==m[b],"Recurrence multiplicity mismatch");
		auto actual=irr::tensor_path(path);inputs<<hex_sha256_file(actual)<<"  "<<fs::absolute(actual).generic_string()<<'\n';
	}
	return K;
}

void check_lec_convention(const fs::path& root,const irr::channels_t& cg,thread_pool* pool) {
	// Compare the explicit fixed CG table, not merely a convention name.
	std::ifstream in(root/"common"/"CG.tsv");require(bool(in),"Missing irrLEC CG convention");
	std::string line;std::getline(in,line);
	for(size_t b=0;b<6;++b)for(size_t input=0;input<6;++input)for(size_t g:cg.channels[b][input])
		for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)for(const auto& term:cg.cg[b].vectors[g].component[s]) {
			std::string target,source,value;size_t channel,component,old,letter;
			in>>target>>source>>channel>>component>>old>>letter>>value;
			std::ostringstream expected;expected<<term.coefficient;
			require(bool(in)&&target==irr::specs()[b].name&&source==irr::specs()[input].name&&channel==g&&component==s
				&&old==term.coordinate/9-cg.axis[input].offset&&letter==term.coordinate%9&&value==expected.str(),"irrLEC CG convention mismatch");
		}
	std::string extra;require(!(in>>extra),"Unexpected CG records");
	auto seed=nmhv_sew_read_tensor<rat_t,int32_t>(root/"LEC_1_physical.wxf",QQ,pool);
	auto m=read_copies(root,1,false);size_t off=0,total=0;
	for(size_t b=0;b<6;++b)total+=m[b]*irr::specs()[b].irrep_dimension;
	require(seed.rank()==3&&seed.dim(0)==total&&seed.dim(1)==9,"Physical seed layout mismatch");
	for(size_t b=0;b<6;++b) {
		for(size_t a=0;a<m[b];++a) {
			int_t gcd=0;
			for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)
				for(size_t p=seed.rowptr()[off+s*m[b]+a];p<seed.rowptr()[off+s*m[b]+a+1];++p) {
					require(seed.val(p).den()==1,"Physical LEC1 contains fractions");gcd=GCD(gcd,seed.val(p).num());
				}
			require(gcd==1,"Physical LEC1 is not copy-primitive");
		}off+=m[b]*irr::specs()[b].irrep_dimension;
	}
	require(seed.rank()==3&&seed.dim(0)==off&&seed.dim(1)==9,"Physical seed layout mismatch");
}

// Independent expanded contraction, used ONLY by --verify. One candidate at a
// time bounds memory; exact reconstruction checks the ENTIRE equation image.
void verify_assembly(const irr::channels_t& cg,const irr::seam_template_t& t,const irr::seam_layout_t& layout,
	const irr::copies& oldfirst,const irr::copies& oldlast,const irr::matrices& first,const irr::matrices& last,
	const sparse_tensor<rat_t,int32_t,SPARSE_CSR>& W,sparse_mat<rat_t,int32_t>& matrix,thread_pool* pool) {
	auto fa=fec_irrep_axis_from_blocks(irr::skeleton(oldfirst,pool));auto la=fec_irrep_axis_from_blocks(irr::skeleton(oldlast,pool));
	size_t ld=fec_irrep_axis_dimension(la),rawdim=fec_irrep_axis_dimension(fa)*W.dim(2)*ld,dim=fec_irrep_axis_dimension(cg.axis);
	require(rawdim<=size_t(INT32_MAX),"Expanded reference exceeds int32");
	auto image=matrix.transpose();
	std::vector<std::tuple<size_t,size_t,size_t>> rows;
	for(size_t r=0;r<t.rows.size();++r)for(size_t a=0;a<oldfirst[t.rows[r].left];++a)
		for(size_t b=0;b<oldlast[t.rows[r].right];++b)rows.push_back({r,a,b});
	for(size_t b=0;b<6;++b) {
		auto FEC=irr::expand(cg,b,oldfirst,first[b]),LEC=irr::expand(cg,b,oldlast,last[b]);
		for(size_t a=0;a<first[b].nrow;++a)for(size_t c=0;c<last[b].nrow;++c) {
			sparse_mat<rat_t,int32_t> raw(1,rawdim),restored(1,rawdim);
			for(size_t s=0;s<irr::specs()[b].irrep_dimension;++s)
				for(auto [i,x]:FEC[s*first[b].nrow+a])for(auto [j,y]:LEC[s*last[b].nrow+c])
					for(size_t p=W.rowptr()[i%9];p<W.rowptr()[i%9+1];++p)if(W.index(p)[0]==j%9)
						raw[0].push_back(int32_t(((i/9)*W.dim(2)+W.index(p)[1])*ld+j/9),x*y*W.val(p)*(s==1?rat_t(3):rat_t(1)));
			for(auto [r,x]:image[layout.colbase[b]+a*last[b].nrow+c]) {
				auto [row,i,j]=rows[r];const auto& tr=t.rows[row];
				for(const auto& [coord,value]:tr.embedding) {
					size_t left=coord/dim/W.dim(2),q=coord/dim%W.dim(2),right=coord%dim;
					left=fa[tr.left].offset+(left-cg.axis[tr.left].offset)*oldfirst[tr.left]+i;
					right=la[tr.right].offset+(right-cg.axis[tr.right].offset)*oldlast[tr.right]+j;
					restored[0].push_back(int32_t((left*W.dim(2)+q)*ld+right),x*value);
				}
			}
			fec_irrep_merge_entries(raw[0]);fec_irrep_merge_entries(restored[0]);
			require(nmhv_sew_matrices_equal(raw,restored),"Expanded sewing image differs from direct multiplicity assembly");
		}
	}
	std::cout<<"VERIFY exact_expanded_image=PASS injective_CG_chart=PASS"<<std::endl;
}

// Reconstruct the FINAL finite-field kernel, not an intermediate rational RREF.
// Existing cache/reconstruction routines retain every completed prime. This
// driver increases the reconstruction budget automatically and reuses caches.
template<class Assemble>
nmhv_sew_reconstruction_result_t<int32_t> solve_seam(Assemble assemble, size_t columns,
	size_t h,size_t k,const fs::path& root,rref_option_t& opt,size_t& primes,std::function<staged_kernel::MM(const field_t&,rref_option_t&)> direct={}) {
	uint64_t fingerprint=std::stoull(hex_sha256_file(root/"input_identity.txt").substr(0,16),nullptr,16);
	for(primes=1;;++primes) {
		try {
			opt->verbose=true;
			return nmhv_sew_reconstruct_with_evaluator<int32_t>((h+k)/2,h,k,primes,root/"cache",fingerprint,opt,
				[&](ulong prime,nmhv_sew_modular_cache_t<int32_t>&,rref_option_t& options) {
					field_t F(FIELD_Fp,prime);
					if(direct){nmhv_sew_prime_result_t<int32_t> result;result.basis=direct(F,options);result.dimensions={columns,result.basis.nrow};result.seam_nnz={0};return result;}
					auto reduced=assemble(F,options);
					size_t nnz=reduced.nnz();reduced.sort_rows_by_nnz();
					nmhv_sew_pivot_stage_t<int32_t> stage;nmhv_sew_prime_result_t<int32_t> result;
					result.basis=nmhv_sew_modular_kernel_rows(std::move(reduced),stage,F,options);
					result.dimensions={columns,result.basis.nrow};result.seam_nnz={nnz};return result;
				},true);
		}catch(const std::runtime_error& e) {
			if(std::string(e.what())!="Rational reconstruction failed validation; increase --modular-primes.")throw;
			std::cout<<"Final kernel reconstruction needs another prime; completed caches retained"<<std::endl;
		}
	}
}

// A single-prime source for the finite-field pipeline. No QQ sewing is needed.
struct source {
    size_t h,k;fs::path root;std::string mode;rref_option_t& opt;irr::channels_t cg;
    irr::copies mf,ml,of,ol;std::shared_ptr<const irr::matrices> first,last;
    sparse_tensor<rat_t,int32_t,SPARSE_CSR> W;
    irr::seam_template_t t;irr::seam_layout_t layout;
    source(size_t H,size_t K,fs::path output,const fs::path& fec,const fs::path& lec,
        const fs::path& condition,const fs::path& blocks,std::string direction,rref_option_t& options,irr::recurrence_store* store=nullptr)
        :h(H),k(K),root(output),mode(std::move(direction)),opt(options),cg(&opt->pool){
        irr_run::phase timer("bulk_prepare");opt->method=1;opt->is_back_sub=true;
        check_lec_convention(lec,cg,&opt->pool);
        mf=read_copies(fec,h,true);ml=read_copies(lec,k,false);of=read_copies(fec,h-1,true);ol=read_copies(lec,k-1,false);
        std::ofstream inputs(root/"inputs.sha256");
        first=store?store->rational(fec,h,true,cg,&opt->pool):std::make_shared<irr::matrices>(read_chain(fec,h,true,mf,inputs,&opt->pool));
        last=store?store->rational(lec,k,false,cg,&opt->pool):std::make_shared<irr::matrices>(read_chain(lec,k,false,ml,inputs,&opt->pool));
        if(store)for(bool isfirst:{true,false})for(size_t b=0;b<6;++b){auto w=isfirst?h:k;
            auto file=irr::tensor_path((isfirst?fec:lec)/joint_eet_weight_name(w)/((isfirst?"irrFEC_":"irrLEC_")+std::to_string(w)+"_"+irr::specs()[b].name+".wxf"));
            inputs<<store->file_hash(file)<<"  "<<fs::absolute(file).generic_string()<<'\n';}
        W=irr::read_tensor(condition,&opt->pool);auto cb=nmhv_sew_read_irrep_axis(blocks);
        t=irr::make_seam_template(cg,W,cb,opt);layout=irr::seam_layout(t,of,ol,mf,ml);
    }
    sparse_mat<ulong,int32_t> evaluate(ulong p){
        irr_run::phase timer("bulk_prime_"+std::to_string(p));field_t F(FIELD_Fp,p);
        for(const auto* chain:{first.get(),last.get()})for(const auto& block:*chain)for(const auto& row:block.rows)
            for(auto [j,x]:row)require(x.den()%F.mod!=0,"Bad bulk prime: recurrence denominator vanishes");
        sparse_mat<ulong,int32_t> A;
        {irr_run::phase assembly("bulk_assemble_prereduce");
            if(mode=="none"){auto exact=irr::assemble_seam(cg,t,layout,of,ol,*first,*last,opt);A=nmhv_sew_modular_reduce(exact,F,&opt->pool);}
            else A=irr::partial_seam(mode=="left",cg,t,layout,of,ol,*first,*last,F,opt);}
        std::cout<<"BULK_MATRIX prime="<<p<<" rows="<<A.nrow<<" columns="<<A.ncol<<" nnz="<<A.nnz()<<std::endl;
        A.sort_rows_by_nnz();nmhv_sew_pivot_stage_t<int32_t> stage;
        sparse_mat<ulong,int32_t> K;
        {irr_run::phase kernel("bulk_final_kernel");K=nmhv_sew_modular_kernel_rows(std::move(A),stage,F,opt);}
        std::vector<size_t> cols;for(const auto& group:stage.pivots)for(auto pivot:group)cols.push_back(pivot.c);std::sort(cols.begin(),cols.end());
        std::ostringstream chart;chart<<layout.columns<<' '<<K.nrow<<'\n';for(auto col:cols)chart<<col<<' ';chart<<'\n';
        irr_stage::fixed_text(root/"pivot_profile.txt",chart.str());
        // The last row is the explicit zero origin, shared by all six blocks.
        sparse_mat<ulong,int32_t> out(K.nrow+1,K.ncol);for(size_t i=0;i<K.nrow;++i)out[i]=std::move(K[i]);return out;
    }
};

}
