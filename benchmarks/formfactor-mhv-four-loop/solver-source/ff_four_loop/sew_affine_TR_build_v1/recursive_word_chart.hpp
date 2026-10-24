#pragma once
// Exact transformations of a small COMPLETE sewn symbol space without forming
// transformation matrices for every intermediate basis. Each private kernel
// coordinate recursively selects a word on which just one basis element is
// nonzero. These words form an exact dual coordinate chart.
#include "symrep_bootstrap.hpp"
namespace recursive_word_chart {
using namespace symrep;
struct Words {std::vector<std::vector<I>> words;std::vector<Q> values;};
inline std::pair<std::vector<I>,std::vector<Q>> private_columns(const Tensor& t){
    require(t.rank()==3,"word charts require rank-three recursive tensors");size_t nc=t.dim(1)*t.dim(2);
    std::vector<I> owners(nc,-1),columns(t.dim(0),-1);std::vector<Q> values(t.dim(0));
    for(size_t b=0;b<t.dim(0);++b)for(size_t p=t.rowptr()[b];p<t.rowptr()[b+1];++p){auto i=t.index(p);size_t c=i[0]*t.dim(2)+i[1];owners[c]=owners[c]==-1?index(b):I(-2);}
    for(size_t b=0;b<t.dim(0);++b)for(size_t p=t.rowptr()[b];p<t.rowptr()[b+1];++p){auto i=t.index(p);size_t c=i[0]*t.dim(2)+i[1];
        if(owners[c]==index(b)&&columns[b]<0){columns[b]=index(c);values[b]=t.val(p);}}
    require(std::all_of(columns.begin(),columns.end(),[](I c){return c>=0;}),"recursive basis lacks a private-word chart");return {columns,values};
}
class Space {
    std::vector<Tensor> forward;
    Tensor last,sewing;
    Words chart;
public:
    Space(const std::vector<fs::path>& chain,const fs::path& last_path,const fs::path& sew_path):last(read(last_path)),sewing(read(sew_path)){
        require(last.dim(2)==1,"private-word sewing requires a one-step last boundary");Words prefix{{{}},{Q(1)}};
        for(const auto& path:chain){auto t=read(path);require(t.dim(1)==prefix.words.size(),"recursive chain dimension mismatch");auto [columns,values]=private_columns(t);Words next;
            for(size_t b=0;b<t.dim(0);++b){size_t old=columns[b]/t.dim(2),letter=columns[b]%t.dim(2);auto word=prefix.words[old];word.push_back(index(letter));next.words.push_back(std::move(word));next.values.push_back(values[b]*prefix.values[old]);}
            prefix=std::move(next);forward.push_back(std::move(t));}
        require(sewing.dim(1)==prefix.words.size()&&sewing.dim(2)==last.dim(0),"sewing dimension mismatch");
        auto [lc,lv]=private_columns(last);auto [sc,sv]=private_columns(sewing);
        for(size_t b=0;b<sewing.dim(0);++b){size_t f=sc[b]/last.dim(0),l=sc[b]%last.dim(0);auto word=prefix.words[f];word.push_back(lc[l]);chart.words.push_back(std::move(word));chart.values.push_back(sv[b]*prefix.values[f]*lv[l]);}
    }
    std::vector<Q> evaluate(const std::vector<I>& word,const Mat& generator)const{
        require(word.size()==forward.size()+1,"word weight mismatch");std::vector<Q> previous(1,Q(1));
        for(size_t w=0;w<forward.size();++w){const auto& t=forward[w];std::vector<Q> letters(generator.nrow),next(t.dim(0));
            for(size_t l=0;l<letters.size();++l)if(auto x=generator[l].find(word[w]))letters[l]=*x;
            for(size_t b=0;b<t.dim(0);++b)for(size_t p=t.rowptr()[b];p<t.rowptr()[b+1];++p){auto i=t.index(p);if(letters[i[1]]!=0&&previous[i[0]]!=0)next[b]+=t.val(p)*letters[i[1]]*previous[i[0]];}
            previous=std::move(next);}
        std::vector<Q> right(last.dim(0)),out(sewing.dim(0));
        for(size_t b=0;b<last.dim(0);++b)for(size_t p=last.rowptr()[b];p<last.rowptr()[b+1];++p)if(auto x=generator[last.index(p)[0]].find(word.back()))right[b]+=last.val(p)*(*x);
        for(size_t b=0;b<sewing.dim(0);++b)for(size_t p=sewing.rowptr()[b];p<sewing.rowptr()[b+1];++p){auto i=sewing.index(p);if(previous[i[0]]!=0&&right[i[1]]!=0)out[b]+=sewing.val(p)*previous[i[0]]*right[i[1]];}
        return out;
    }
    Mat action(const Mat& generator)const{Mat result(sewing.dim(0),sewing.dim(0));
        for(size_t j=0;j<chart.words.size();++j){auto values=evaluate(chart.words[j],generator);for(size_t i=0;i<values.size();++i)if(values[i]!=0)result[i].push_back(index(j),values[i]/chart.values[j]);}return result;}
    const Tensor& seed()const{return forward.front();}
    const Tensor& boundary()const{return last;}
    const Tensor& basis()const{return sewing;}
};
inline void save(const fs::path& path,const Mat& m){auto t=tensor(m,{m.nrow,m.ncol});auto bytes=sparse_tensor_write_wxf(t);std::ofstream f(path,std::ios::binary);f.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());require(bool(f),"cannot save word-chart result");}
inline void restrict_invariants(const fs::path& sewn,const std::vector<fs::path>& chain,const fs::path& last,
    const fs::path& conditions,const std::vector<std::pair<std::string,fs::path>>& generators,rref_option_t opt){
    if(generators.empty())return;auto start=std::chrono::steady_clock::now();Space space(chain,last,sewn);size_t n=space.basis().dim(0),letters=space.seed().dim(2);
    require(equal(space.action(identity(letters)),identity(n)),"recursive private-word chart failed identity check");
    auto rules=condition_rows(read(conditions)),first=flatten(space.seed()),boundary=flatten(space.boundary());Chart rule_chart(rules,opt),first_chart(first,opt),last_chart(boundary,opt);
    Mat equations(generators.size()*n,n);size_t count=0;std::vector<Mat> actions;
    for(const auto& [name,path]:generators){auto g=read_matrix(path);require(g.nrow==letters&&g.ncol==letters,"symmetry alphabet mismatch");
        // Local conditions and both endpoints imply invariance of every
        // complete recursively constructed symbol space, by induction.
        rule_chart.coordinates(mul(rules,transpose(kron(g,g)),&opt->pool));
        first_chart.coordinates(mul(flatten(space.seed()),g,&opt->pool));last_chart.coordinates(mul(flatten(space.boundary()),g,&opt->pool));
        auto action=space.action(g);save(sewn.string()+"."+name+".raw-action.wxf",action);
        for(size_t i=0;i<n;++i){for(auto [j,x]:action[i])equations[count*n+j].push_back(index(i),x);equations[count*n+i].push_back(index(i),Q(-1));}
        actions.push_back(std::move(action));++count;}
    for(auto& row:equations.rows)normalize(row);auto invariant=kernel(std::move(equations),opt);
    for(const auto& action:actions)require(equal(mul(invariant,action,&opt->pool),invariant),"invariant subspace failed exact generator check");
    std::cout<<"MHV symmetry dimension "<<n<<" -> "<<invariant.nrow<<" word_chart_seconds="<<std::chrono::duration<double>(std::chrono::steady_clock::now()-start).count()<<std::endl;
    if(invariant.nrow==n){save(sewn.string()+".invariant-map.wxf",invariant);return;}
    auto raw_bytes=sparse_tensor_write_wxf(space.basis());{std::ofstream f(sewn.string()+".unrestricted.wxf",std::ios::binary);f.write(reinterpret_cast<const char*>(raw_bytes.data()),raw_bytes.size());require(bool(f),"cannot save unrestricted sewing basis");}
    auto reduced=mul(invariant,flatten(space.basis()),&opt->pool);for(size_t r=0;r<reduced.nrow;++r){Q before=reduced[r][0];vec_cancel_divisor(reduced[r]);sparse_vec_rescale(invariant[r],reduced[r][0]/before,QQ);}
    save(sewn.string()+".invariant-map.wxf",invariant);
    auto result=tensor(reduced,{reduced.nrow,space.basis().dim(1),space.basis().dim(2)},&opt->pool);
    auto bytes=sparse_tensor_write_wxf(result);std::ofstream f(sewn,std::ios::binary);f.write(reinterpret_cast<const char*>(bytes.data()),bytes.size());require(bool(f),"cannot save invariant sewing basis");
}
}
