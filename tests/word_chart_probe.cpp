#include "recursive_word_chart.hpp"
#include "restricted_projection.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==2,"word_chart_probe OUTPUT");fs::path out=argv[1];fresh(out);
    auto f1=tensor(identity(3),{3,1,3}),f2=tensor(identity(9),{9,3,3}),last=tensor(identity(3),{3,3,1}),sew=tensor(identity(27),{27,9,3});
    write(out/"f1.wxf",f1);write(out/"f2.wxf",f2);write(out/"last.wxf",last);write(out/"sew.wxf",sew);
    recursive_word_chart::Space space({out/"f1.wxf",out/"f2.wxf"},out/"last.wxf",out/"sew.wxf");
    std::vector<Mat> generators;Mat cycle(3,3);for(size_t i=0;i<3;++i)cycle[i].push_back(index((i+1)%3),Q(1));generators.push_back(cycle);
    Mat dense(3,3);for(size_t i=0;i<3;++i)for(size_t j=0;j<3;++j)dense[i].push_back(index(j),Q(long(1+i+3*j))/Q(long(1+i+j)));generators.push_back(dense);
    generators.push_back(identity(3));
    for(const auto& g:generators)require(equal(space.action(g),kron(kron(g,g),g)),"private-word action differs from full tensor action");
    Mat sm(2,6),fm(3,6),prefix(2,2),letter(3,2),lastmap(2,2);
    auto fill=[](Mat& a){for(size_t i=0;i<a.nrow;++i)for(size_t j=0;j<a.ncol;++j){Q x=Q(long((i*3+j*2)%7)-3)/Q(1+i+j);if(x!=0)a[i].push_back(index(j),x);}};
    fill(sm);fill(fm);fill(prefix);fill(letter);fill(lastmap);
    auto st=tensor(sm,{2,3,2}),ft=tensor(fm,{3,2,3});
    auto restricted=restricted_projection::sew_one(st,ft,prefix,letter,lastmap,nullptr);Mat expected(2,8);
    for(size_t n=0;n<2;++n)for(auto [fl,s]:sm[n])for(auto [ol,f]:fm[fl/2])
        for(auto [o,p]:prefix[ol/3])for(auto [l,c]:letter[ol%3])for(auto [r,d]:lastmap[fl%2])expected[n].push_back(index((o*2+l)*2+r),s*f*p*c*d);
    for(auto& row:expected.rows)normalize(row);
    require(equal(flatten(restricted),expected),"restricted projection differs from explicit five-factor sum");
    rref_option_t projection_opt;projection_opt->pool.reset(2);
    require(equal(flatten(restricted_projection::sew_one_compact_right(st,ft,prefix,letter,lastmap,projection_opt)),expected),"full-rank right frame changed the projection");
    for(size_t rank=0;rank<3;++rank){Mat wide(2,5);
        if(rank)for(size_t j=0;j<5;++j){Q x(long(j)-2,long(j+1));if(x!=0)wide[0].push_back(index(j),x);Q y=rank==1?Q(2,3)*x:Q(long(3*j)-1,long(j+2));if(y!=0)wide[1].push_back(index(j),y);}
        Mat direct(2,20);for(size_t n=0;n<2;++n)for(auto [fl,s]:sm[n])for(auto [ol,f]:fm[fl/2])
            for(auto [o,p]:prefix[ol/3])for(auto [l,c]:letter[ol%3])for(auto [r,d]:wide[fl%2])direct[n].push_back(index((o*2+l)*5+r),s*f*p*c*d);
        for(auto& row:direct.rows)normalize(row);auto got=restricted_projection::sew_one_compact_right(st,ft,prefix,letter,wide,projection_opt);
        require(equal(flatten(got),direct),"compact right projection differs from independent five-factor sum");
    }
    for(size_t nc:{size_t(3),size_t(200),size_t(40000),size_t(1ULL<<32)}){
        // Large tensor dimensions exercise every packed-index width without
        // allocating their ambient dense space.
        size_t actual=std::min(nc,size_t(INT32_MAX));Mat m(2,actual);m[0].push_back(0,Q(-7)/Q(13));m[1].push_back(index(actual-1),Q(int_t(2).pow(90ul)+int_t(3)));
        auto t=tensor(m,{2,actual});MatrixTensorView view(m,{2,actual});auto path=out/("encoding-"+std::to_string(nc)+".wxf");write_tensor_view(path,view);
        require(file_to_ustr(path)==encode_tensor_view(t),"streamed WXF differs from SparseRREF encoder");
        require(equal(flatten(read(path)),m),"streamed WXF readback differs");
    }
    Mat empty(2,3);auto ep=out/"encoding-empty.wxf";write_tensor_view(ep,MatrixTensorView(empty,{2,3}));require(equal(flatten(read(ep)),empty),"empty streamed WXF differs");
    Mat many(5,20000);for(size_t r=0;r<many.nrow;++r)for(size_t c=0;c<many.ncol;++c)many[r].push_back(index(c),Q(1+r+c)/Q(7));
    MatrixTensorView many_view(many,{5,20000});auto mp=out/"encoding-chunks.wxf";write_tensor_view(mp,many_view);
    require(file_to_ustr(mp)==encode_tensor_view(tensor(many,{5,20000}))&&equal(flatten(read(mp)),many),"multi-chunk WXF differs from independent encoder/parser");
    {std::fstream damage(mp,std::ios::in|std::ios::out|std::ios::binary);damage.seekg(-1,std::ios::end);char c;damage.read(&c,1);c^=1;damage.seekp(-1,std::ios::end);damage.write(&c,1);}
    bool corrupt=false;try{verify_streamed_encoding(mp,many_view);}catch(const std::exception&){corrupt=true;}require(corrupt,"streamed readback accepted a corrupted byte");
    std::cout<<"PASS restricted projection against explicit rational five-factor sum"<<std::endl;
    std::cout<<"PASS recursive word charts against explicit tensor actions (cycle, dense rational transform, identity)"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
