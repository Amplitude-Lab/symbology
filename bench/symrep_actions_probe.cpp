// Exact recursive transformations versus direct irreducible block actions.
// Legacy implementation is retained here only for controlled measurements.
#include "symrep_bootstrap.hpp"
using namespace symrep;
inline std::vector<Mat> legacy_carrier_actions(const Mat& carrier,const std::vector<Mat>& older,
    const std::vector<Mat>& letters,bool backward,rref_option_t opt){
    Chart chart(carrier,opt);std::vector<I> selected(carrier.ncol,-1);
    for(size_t i=0;i<chart.pivots.size();++i)selected[chart.pivots[i]]=index(i);
    std::vector<Mat> result;
    for(size_t h=0;h<letters.size();++h){const auto& left=backward?letters[h]:older[h];const auto& right=backward?older[h]:letters[h];
        require(left.nrow*right.nrow==carrier.ncol,"carrier action dimension mismatch");
        Mat restriction(carrier.ncol,carrier.nrow);
        for(size_t i=0;i<left.nrow;++i)for(size_t j=0;j<right.nrow;++j)for(auto [k,x]:left[i])for(auto [l,y]:right[j]){
            I c=selected[k*right.nrow+l];if(c>=0)restriction[i*right.nrow+j].push_back(c,x*y);}
        for(auto& row:restriction.rows)normalize(row);
        result.push_back(mul(mul(carrier,restriction,&opt->pool),chart.inverse,&opt->pool));
    }return result;
}
int main(int argc,char** argv){try{
    require(argc==8,"actions_probe PREPARED CHAIN MAX legacy|carrier|adapted OUTPUT THREADS verify|measure");
    rref_option_t opt;opt->pool.reset(std::stoul(argv[6]));opt->verbose=false;
    auto g=load_group(argv[1]);auto ms=load_models(argv[1],g,opt);auto p=load_prepared(argv[1],ms,false);
    fs::path chain=argv[2],out=argv[5];size_t maximum=std::stoul(argv[3]);std::string mode=argv[4];
    require(mode=="legacy"||mode=="carrier"||mode=="adapted","unknown action mode");fresh(out);
    bool factorized=fs::exists(chain/"w1_carrier.wxf"),verify=std::string(argv[7])=="verify";
    auto terminal=read_matrix(fs::path(argv[1])/"terminal_basis.wxf");auto inverse=sparse_mat_inverse(terminal,QQ,opt);
    std::vector<Mat> previous;for(auto h:g.generator_ids)previous.push_back(mul(mul(inverse,action(p.terminal,ms,h)),terminal));
    std::ofstream stats(out/"timings.tsv");stats<<"weight\tdimension\tconstruction_s\tapply_s\timplicit_apply_s\taction_nnz\n";
    for(size_t w=1;w<=maximum;++w){Mat basis;Layout layout;
        if(mode=="adapted")layout=load_layout(chain/(weight_name(w)+"_copies.tsv"),ms);
        else basis=flatten(read(chain/(weight_name(w)+(factorized?"_carrier.wxf":".wxf"))));
        auto start=std::chrono::steady_clock::now();std::vector<Mat> next;
        if(mode=="adapted")next=generators(g,ms,layout);
        else if(mode=="legacy")next=legacy_carrier_actions(basis,previous,g.generators,p.backward,opt);
        else next=carrier_actions(basis,previous,g.generators,p.backward,opt);
        double elapsed=seconds(start);size_t nnz=0;for(const auto& a:next)nnz+=a.nnz();
        if(verify&&mode!="adapted")check_carrier_action(basis,next,previous,g,p.backward,opt);
        // Apply all supplied generators to a deterministic dense coefficient row.
        // The same coordinate entries represent different physical vectors in
        // different bases; this measures operation cost, not vector equivalence.
        Mat vector(1,next[0].nrow);for(size_t i=0;i<vector.ncol;++i)vector[0].push_back(index(i),Q(long(i%13)-6));
        normalize(vector[0]);size_t sink=0;std::vector<double> samples,implicit;
        for(size_t repeat=0;repeat<12;++repeat){start=std::chrono::steady_clock::now();
            for(const auto& a:next)sink+=mul(vector,a,&opt->pool).nnz();if(repeat)samples.push_back(seconds(start));
            if(mode=="adapted"){start=std::chrono::steady_clock::now();for(auto h:g.generator_ids)sink+=apply_irrep_action(vector[0],layout,ms,h).nnz();if(repeat)implicit.push_back(seconds(start));}}
        std::sort(samples.begin(),samples.end());std::sort(implicit.begin(),implicit.end());double apply=samples[samples.size()/2],direct=implicit.empty()?0:implicit[implicit.size()/2];
        if(mode=="adapted")for(size_t j=0;j<next.size();++j)require(equal(rows({apply_irrep_action(vector[0],layout,ms,g.generator_ids[j])},vector.ncol),mul(vector,next[j],&opt->pool)),"implicit action differs from explicit action");
        std::cout<<"ACTIONS weight="<<w<<" dim="<<next[0].nrow<<" seconds="<<elapsed<<" nnz="<<nnz<<" apply_s="<<apply<<" implicit_s="<<direct<<" sink="<<sink<<std::endl;
        stats<<std::setprecision(9)<<w<<'\t'<<next[0].nrow<<'\t'<<elapsed<<'\t'<<apply<<'\t'<<direct<<'\t'<<nnz<<'\n';
        for(size_t j=0;j<next.size();++j)write(out/(weight_name(w)+"_g"+std::to_string(j)+".wxf"),next[j]);
        previous=std::move(next);
    }stats.close();seal(out);
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
