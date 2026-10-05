// Exact NMHV fixtures and independent comparisons. No physical amplitude is
// claimed by a last-entry basis comparison alone.
#include "symrep_bootstrap.hpp"
using namespace symrep;
bool private_rank(const Mat& a){
    std::vector<size_t> counts(a.ncol);for(const auto& row:a.rows)for(auto c:row.index_span())++counts[c];
    for(const auto& row:a.rows){bool found=false;for(auto [c,x]:row)if(counts[c]==1&&x!=0){found=true;break;}if(!found)return false;}return true;
}
int main(int argc,char** argv){try{
    require(argc>=2,"nmhv_kernel_check fixture|compare ...");
    rref_option_t opt;opt->pool.reset(8);
    std::string mode=argv[1];
    if(mode=="fixture") {
        require(argc==3,"fixture OUTPUT");fs::path out=argv[2];fresh(out);
        // Eq. (7),(8), arXiv:2511.09669v2. Zero-based cyclic labels.
        std::vector<std::pair<int,int>> pairs;
        for(int d=1;d<=3;++d)for(int c=0;c<7;++c)pairs.push_back(std::minmax(c,(c+d)%7));
        auto id=[&](int a,int b){std::pair<int,int> p=std::minmax((a+7)%7,(b+7)%7);auto it=std::find(pairs.begin(),pairs.end(),std::pair<int,int>(p));require(it!=pairs.end(),"pair not found");return size_t(it-pairs.begin());};
        Mat terminal(21,15);
        for(int c=0;c<7;++c)terminal[id(c,c+2)].push_back(0,Q(1));
        for(int c=0;c<7;++c){
            for(auto [a,b,s]:std::vector<std::tuple<int,int,int>>{{0,1,1},{0,5,-1},{1,3,-1},{3,5,-1}})terminal[id(a+c,b+c)].push_back(1+c,Q(s));
            for(auto [a,b,s]:std::vector<std::tuple<int,int,int>>{{0,3,1},{0,5,-1},{3,5,-1}})terminal[id(a+c,b+c)].push_back(8+c,Q(s));
        }
        for(auto& row:terminal.rows)normalize(row);
        require(row_basis(terminal).nrow==15,"R-invariant quotient rank");
        const std::array<std::array<int,7>,3> hns{{{{15,21,26,32,34,53,57}},{{21,23,31,33,41,43,62}},{{11,14,21,24,31,34,46}}}};
        Mat seed(147,42*15);size_t r=0;
        for(size_t d=0;d<3;++d)for(int c=0;c<7;++c)for(int a:hns[d]){
            int letter=7*(a/10-1)+(a%10-1+c)%7;
            for(auto [t,x]:terminal[d*7+c])seed[r].push_back(index(letter*15+t),x);++r;
        }
        seed=row_basis(seed);write(out/"LEC_1.wxf",tensor(seed,{seed.nrow,42,15}));write(out/"R_quotient.wxf",terminal);
        for(bool flip:{false,true}){
            Mat transformed(21,15);for(size_t i=0;i<21;++i){auto [a,b]=pairs[i];transformed[i]=terminal[id(flip?6-a:a+1,flip?6-b:b+1)];}
            auto action=transpose(Chart(transpose(terminal),opt).coordinates(transpose(transformed)));
            write(out/(flip?"flip_terminal.wxf":"cyclic_terminal.wxf"),action);
        }
        seal(out);std::cout<<"HEPTAGON_NMHV_SEED dimension="<<seed.nrow<<" terminal=15\n";
    } else if(mode=="sew-compare") {
        require(argc==4,"sew-compare ACTUAL REFERENCE");auto a=flatten(read(argv[2])),b=flatten(read(argv[3]));
        require(a.nrow==b.nrow,"sewing nullity mismatch");auto change=Chart(b,opt).coordinates(a);
        require(kernel(change,opt).nrow==0,"singular sewing basis change");std::cout<<"EXACT_SEW_PASS dimension="<<a.nrow<<std::endl;
    } else if(mode=="hept-symmetry") {
        require(argc==6,"hept-symmetry FIXTURE FEC_CHAIN FEC_WEIGHT SEW");fs::path fixture=argv[2],chain=argv[3];size_t weight=std::stoul(argv[4]);
        auto t=read_matrix(fixture/"R_quotient.wxf"),last=flatten(read(fixture/"LEC_1.wxf")),sew=flatten(read(argv[5]));
        std::vector<std::pair<int,int>> pairs;for(int d=1;d<=3;++d)for(int c=0;c<7;++c)pairs.push_back(std::minmax(c,(c+d)%7));
        auto id=[&](int a,int b){std::pair<int,int> p=std::minmax((a+14)%7,(b+14)%7);return size_t(std::find(pairs.begin(),pairs.end(),p)-pairs.begin());};
        Mat equations(0,sew.nrow);
        for(const auto& name:{"cyclic","flip"}) {
            auto g=read_matrix(std::string("data/")+(std::string(name)=="cyclic"?"cycrepmat.wxf":"fliprepmat.wxf"));
            Mat terminal;bool closed=false;int found=-1;
            for(int shift=0;shift<7&&!closed;++shift){Mat image(21,15);
                for(size_t i=0;i<21;++i){auto [a,b]=pairs[i];image[i]=t[id(std::string(name)=="cyclic"?a+1:shift-a,std::string(name)=="cyclic"?b+1:shift-b)];}
                terminal=transpose(Chart(transpose(t),opt).coordinates(transpose(image)));
                try{Chart(last,opt).coordinates(mul(last,kron(g,terminal)));closed=true;found=shift;}catch(const std::runtime_error&){}
            }
            require(closed,"NMHV seed not closed under "+std::string(name));std::cout<<"NMHV terminal "<<name<<" shift="<<found<<std::endl;
            auto previous=identity(1);for(size_t w=1;w<=weight;++w){auto f=flatten(read(chain/("w"+std::to_string(w)+".wxf")));previous=Chart(f,opt).coordinates(mul(f,kron(previous,g)));}
            auto dl=Chart(last,opt).coordinates(mul(last,kron(g,terminal)));
            auto ds=Chart(sew,opt).coordinates(mul(sew,kron(previous,dl)));
            for(size_t r=0;r<ds.nrow;++r){ds[r].push_back(index(r),Q(-1));normalize(ds[r]);}
            auto eq=transpose(ds);for(auto& row:eq.rows)equations.rows.push_back(std::move(row));equations.nrow=equations.rows.size();
        }
        auto invariant=kernel(std::move(equations),opt);std::cout<<"HEPTAGON_NMHV_INVARIANTS weight="<<weight+1<<" dimension="<<invariant.nrow<<std::endl;
        require((weight==1&&invariant.nrow==5)||(weight==3&&invariant.nrow==11)||(weight==5&&invariant.nrow==24),"published NMHV dimension mismatch");
    } else if(mode=="compare") {
        require(argc==6,"compare backward|forward ACTUAL REFERENCE MAX");bool backward=std::string(argv[2])=="backward";
        fs::path actual=argv[3],reference=argv[4];size_t maximum=std::stoul(argv[5]);
        auto seed=read(actual/"w1.wxf");auto previous=identity(seed.dim(backward?2:1)),letter=identity(seed.dim(backward?1:2));
        for(size_t w=1;w<=maximum;++w){auto stem="w"+std::to_string(w)+".wxf";
            auto a=flatten(read(actual/stem)),b=flatten(read(reference/stem));require(a.nrow==b.nrow,"nullity mismatch");
            auto physical=mul(a,backward?kron(letter,previous):kron(previous,letter),&opt->pool);
            auto change=Chart(b,opt).coordinates(physical);
            // A has full row rank if every row has a private nonzero column.
            // The previous change is invertible by induction. Together with
            // exact inclusion and equal dimensions this proves invertibility
            // of this change, without row reducing a large dense change matrix.
            if(!private_rank(a))require(kernel(change,opt).nrow==0,"singular recursive basis change");
            previous=std::move(change);std::cout<<"EXACT_CHAIN_PASS weight="<<w<<" dimension="<<b.nrow<<std::endl;
        }
    } else throw std::runtime_error("unknown mode");
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
