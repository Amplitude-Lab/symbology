// Published coefficients are independent checks, never solver inputs.
// He et al., arXiv:2511.09669v2, Eq. (12) and the following paragraph.
#include "recursive_word_chart.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==3,"heptagon_amplitude_check DATA OUTPUT");fs::path data=argv[1],out=argv[2];
    for(size_t loop=2;loop<=4;++loop){size_t weight=2*loop;std::vector<fs::path> chain{data/"FEC_1.wxf"};
        for(size_t w=2;w<weight;++w)chain.push_back(out/("FEC_"+std::to_string(w)+".wxf"));
        recursive_word_chart::Space space(chain,data/"LEC_1.wxf",out/(std::to_string(loop)+"loop")/("hepMHV_"+std::to_string(loop)+"L_recursive.wxf"));
        auto id=identity(42);long odd=1;for(size_t i=1;i<=2*loop-3;i+=2)odd*=i;
        long expected=(loop%2?-1:1)*long(1UL<<(2*(loop-1)))*odd;
        for(I last:{I(7),I(8),I(11)}){std::vector<I> word(weight,0);word.back()=last;
            auto value=space.evaluate(word,id);require(value.size()==1&&value[0]==Q(last==11?-expected:expected),"published repeated-first coefficient mismatch");
            std::cout<<"loop="<<loop<<" repeated_first_last="<<last+1<<" coefficient="<<value[0]<<std::endl;}
        std::vector<I> word(weight,7);word[0]=0;auto value=space.evaluate(word,id);
        require(value[0]==Q(expected/long(1UL<<(loop-1))),"published repeated-last coefficient mismatch");
        std::cout<<"loop="<<loop<<" repeated_last coefficient="<<value[0]<<std::endl;
    }
    std::cout<<"PASS 12 independent published word coefficients at two, three and four loops"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
