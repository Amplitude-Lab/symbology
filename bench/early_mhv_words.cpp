// Independent published word coefficients, used only AFTER construction.
// He et al., arXiv:2511.09669v2, equation (12) and its following paragraph.
#include "symrep_io.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==6,"early_mhv_words DATA FEC_CHAIN LEC2 SYMBOL LOOP");fs::path data=argv[1],chain=argv[2];size_t loop=std::stoul(argv[5]);
    require(loop>=2&&loop<=5,"published validation covers loops two through five");size_t weight=2*loop,fw=weight-2;
    std::vector<std::vector<Q>> previous(2,std::vector<Q>{Q(1)});
    for(size_t w=1;w<=fw;++w){auto f=read(w==1?data/"FEC_1.wxf":chain/("FEC_"+std::to_string(w)+".wxf"));std::vector<std::vector<Q>> next(2,std::vector<Q>(f.dim(0)));
        for(size_t b=0;b<f.dim(0);++b)for(size_t p=f.rowptr()[b];p<f.rowptr()[b+1];++p){auto i=f.index(p);
            if(i[1]==0&&previous[0][i[0]]!=0)next[0][b]+=f.val(p)*previous[0][i[0]];
            if(i[1]==(w==1?0:7)&&previous[1][i[0]]!=0)next[1][b]+=f.val(p)*previous[1][i[0]];}
        previous=std::move(next);
    }
    auto last=flatten(read(data/"LEC_1.wxf"));auto last2=read(argv[3]);auto symbol=read(argv[4]);require(symbol.dim(0)==1&&symbol.dim(1)==previous[0].size()&&symbol.dim(2)==last2.dim(0),"published check basis mismatch");
    long expected=(loop%2?-1:1)*long(1UL<<(2*(loop-1)));for(size_t j=1;j<=2*loop-3;j+=2)expected*=j;
    for(size_t test=0;test<4;++test){I penultimate=test<3?0:7,final=test==0?7:test==1?8:test==2?11:7;std::vector<Q> right(last2.dim(0));
        for(size_t b=0;b<last2.dim(0);++b)for(size_t p=last2.rowptr()[b];p<last2.rowptr()[b+1];++p){auto i=last2.index(p);if(i[0]==penultimate)if(auto x=last[i[1]].find(final))right[b]+=last2.val(p)*(*x);}
        Q value=0;for(size_t p=symbol.rowptr()[0];p<symbol.rowptr()[1];++p){auto i=symbol.index(p);value+=symbol.val(p)*previous[test<3?0:1][i[0]]*right[i[1]];}
        Q want(test==2?-expected:test==3?expected/long(1UL<<(loop-1)):expected);require(value==want,"published heptagon word coefficient mismatch");
        std::cout<<"loop="<<loop<<" published_case="<<test<<" coefficient="<<value<<std::endl;
    }
    std::cout<<"PASS four independent published original-alphabet word coefficients"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
