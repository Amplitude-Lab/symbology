#pragma once
#include "staged_kernel.hpp"

// Invariance on a tensor product, from arbitrary paired generator matrices.
namespace staged_kernel {
class ProductInvarianceRows {
    size_t left_dimension,right_dimension;
    std::vector<Mat> left_columns,right_columns;
    int_t height=0;
public:
    ProductInvarianceRows(std::vector<Mat> left,std::vector<Mat> right)
        :left_dimension(left.at(0).nrow),right_dimension(right.at(0).nrow){
        require(left.size()==right.size(),"invariant product generator count mismatch");
        for(size_t g=0;g<left.size();++g){
            require(left[g].nrow==left_dimension&&left[g].ncol==left_dimension&&right[g].nrow==right_dimension&&right[g].ncol==right_dimension,"invariant product generator dimensions mismatch");
            auto a=transpose(left[g]),b=transpose(right[g]);left[g].clear();right[g].clear();
            int_t da=1,db=1,ha=0,hb=0;
            // Each output coordinate is an independent equation and can
            // clear its own denominators. A global LCM across all coordinates
            // creates an unnecessarily enormous certificate bound.
            auto bound=[](const Mat& m,int_t& md,int_t& mh){for(const auto& row:m.rows){int_t d=1,h=0;for(auto [j,x]:row)d=LCM(d,x.den());
                    for(auto [j,x]:row){int_t v=d;v/=x.den();v*=x.num().abs();h+=v;}if(d>md)md=d;if(h>mh)mh=h;}};
            bound(a,da,ha);bound(b,db,hb);
            auto h=ha*hb+da*db;if(h>height)height=h;
            left_columns.push_back(std::move(a));right_columns.push_back(std::move(b));
        }
    }
    size_t ncol()const{return left_dimension*right_dimension;}
    size_t nrow()const{return generator_count()*ncol();}
    size_t residual_blocks()const{return generator_count();}
    template<class Emit> void residual(const sparse_mat<ulong,I>& k,const field_t& fp,size_t block,Emit&& emit,thread_pool* pool)const{
        modular_residual(k,fp,block,std::forward<Emit>(emit),pool);
    }
    size_t generator_count()const{return left_columns.size();}
    template<class Emit> void modular_residual(const sparse_mat<ulong,I>& k,const field_t& fp,size_t generator,Emit&& emit,thread_pool* pool)const{
        using MV=sparse_vec<ulong,I>;using MM=sparse_mat<ulong,I>;
        require(k.ncol==ncol()&&generator<left_columns.size(),"invariant residual dimensions mismatch");
        MM by_left(left_dimension,k.nrow*right_dimension);for(size_t s=0;s<k.nrow;++s)for(auto [c,x]:k[s])by_left[c/right_dimension].push_back(index(s*right_dimension+c%right_dimension),x);
        MM a(left_dimension,left_dimension),b(right_dimension,right_dimension);
        auto convert=[&](const Mat& q,MM& m){for(size_t i=0;i<q.nrow;++i){for(auto [j,x]:q[i])if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("invariant residual denominator");m[i]=q[i]%fp.mod;}};
        convert(left_columns[generator],a);convert(right_columns[generator],b);auto transformed=sparse_mat_mul(a,by_left,fp,pool);a.clear();
        for(size_t first=0;first<left_dimension;first+=128){size_t end=std::min(left_dimension,first+128);std::vector<std::vector<MV>> result(end-first);
            auto evaluate=[&](size_t i){size_t l=first+i;std::vector<ulong> v(k.nrow*right_dimension),original(v.size());
                for(auto [j,x]:transformed[l])v[j]=x;for(auto [j,x]:by_left[l])original[j]=x;
                for(size_t r=0;r<right_dimension;++r){MV row;for(size_t s=0;s<k.nrow;++s){ulong value=nmod_neg(original[s*right_dimension+r],fp.mod);
                        for(auto [j,x]:b[r])value=nmod_add(value,nmod_mul(x,v[s*right_dimension+j],fp.mod),fp.mod);
                        if(value)row.push_back(index(s),value);}
                    if(row.nnz())result[i].push_back(std::move(row));}
            };
            if(pool){pool->detach_loop(size_t(0),end-first,evaluate);pool->wait();}else for(size_t i=0;i<end-first;++i)evaluate(i);
            for(auto& block:result)for(auto& row:block)emit(std::move(row));
        }
    }
    template<class Emit> void operator()(structured_kernel::Relations& relations,size_t limit,Emit&& emit)const{
        for(size_t g=0;g<left_columns.size();++g)for(size_t l=0;l<left_dimension;++l)for(size_t r=0;r<right_dimension;++r){
            size_t diagonal=l*right_dimension+r;
            if(limit!=SIZE_MAX){size_t count=relations.is_zero(diagonal)?0:1;bool skip=false;
                for(auto i:left_columns[g][l].index_span()){
                    for(auto j:right_columns[g][r].index_span())if(!relations.is_zero(size_t(i)*right_dimension+j)&&++count>limit){skip=true;break;}
                    if(skip)break;
                }if(skip)continue;
            }
            Vec row;auto add=[&](size_t c,const Q& v){auto [root,scale]=relations.image(c);if(scale!=0)row.push_back(index(root),v*scale);};
            add(diagonal,Q(-1));for(auto [i,x]:left_columns[g][l])for(auto [j,y]:right_columns[g][r])add(size_t(i)*right_dimension+j,x*y);
            normalize(row);if(row.nnz())emit(std::move(row));
        }
    }
    template<class Emit> void modular(structured_kernel::Relations& relations,const std::vector<I>& columns,const field_t& fp,Emit&& emit,thread_pool* pool=nullptr,size_t max_terms=SIZE_MAX,size_t min_terms=0,const streaming_echelon::Substitution* substitution=nullptr,std::vector<unsigned char>* seen=nullptr,unsigned sample_shift=0)const{
        require(!seen||seen->size()==nrow(),"invariance discovery mask size mismatch");
        require(sample_shift<64,"invalid invariant sampling shift");
        using MV=sparse_vec<ulong,I>;std::vector<I> mapped(ncol(),-1);std::vector<ulong> scales(ncol());
        for(size_t c=0;c<ncol();++c){auto [r,x]=relations.image(c);if(x==0)continue;if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("invariant relation denominator");mapped[c]=columns[r];scales[c]=x%fp.mod;
            if(substitution&&mapped[c]>=0){I j=mapped[c];mapped[c]=substitution->image[j];scales[c]=nmod_mul(scales[c],substitution->scale[j],fp.mod);}}
        auto convert=[&](const Mat& q){sparse_mat<ulong,I> m(q.nrow,q.ncol);for(size_t i=0;i<q.nrow;++i){for(auto [j,x]:q[i])if(x.den()%fp.mod==0)throw certified_kernel::BadPrime("invariant action denominator");m[i]=q[i]%fp.mod;}return m;};
        for(size_t g=0;g<left_columns.size();++g){auto a=convert(left_columns[g]),b=convert(right_columns[g]);
            auto build=[&](size_t c){size_t id=g*ncol()+c;if((seen&&(*seen)[id])||!streaming_echelon::sampled_equation(id,sample_shift,UINT64_C(0xd1b54a32d192ed03)))return MV{};size_t l=c/right_dimension,r=c%right_dimension;MV row;
                if(max_terms!=SIZE_MAX){size_t support=mapped[c]>=0?1:0;
                    for(auto i:a[l].index_span())for(auto j:b[r].index_span())
                        if(mapped[size_t(i)*right_dimension+j]>=0&&++support>max_terms)return MV{};
                }
                if(mapped[c]>=0)row.push_back(mapped[c],nmod_neg(scales[c],fp.mod));
                for(auto [i,x]:a[l])for(auto [j,y]:b[r]){size_t k=size_t(i)*right_dimension+j;if(mapped[k]>=0){auto v=nmod_mul(nmod_mul(x,y,fp.mod),scales[k],fp.mod);if(v)row.push_back(mapped[k],v);if(row.nnz()>max_terms)return MV{};}}
                if(row.nnz()<=min_terms){if(seen&&!row.nnz())(*seen)[id]=1;return MV{};}row.sort_indices();size_t out=0;
                for(size_t first=0;first<row.nnz();){size_t end=first+1;I col=row(first);ulong value=row[first];while(end<row.nnz()&&row(end)==col)value=nmod_add(value,row[end++],fp.mod);
                    if(value){row(out)=col;row[out]=value;++out;}first=end;}row.resize(out);row.compress();if(seen)(*seen)[id]=1;return row;};
            for(size_t first=0;first<ncol();){size_t end=first,terms=0;
                while(end<ncol()&&end-first<1024){size_t cost=a[end/right_dimension].nnz()*b[end%right_dimension].nnz()+1;
                    if(max_terms!=SIZE_MAX)cost=std::min(cost,max_terms+1);
                    if(end>first&&cost>2*1024*1024-std::min(terms,size_t(2*1024*1024)))break;terms+=cost;++end;}
                std::vector<MV> tile(end-first);if(pool&&pool->get_thread_count()>1){pool->detach_loop(size_t(0),tile.size(),[&](size_t i){tile[i]=build(first+i);});pool->wait();}
                else for(size_t i=0;i<tile.size();++i)tile[i]=build(first+i);
                for(auto& row:tile)if(row.nnz())emit(std::move(row));first=end;
            }
        }
    }
    int_t integral_residual_bound(const Mat& candidate)const{
        int_t hk=0;for(const auto& row:candidate.rows)for(auto [j,x]:row){require(x.is_integer(),"invariant certificate requires primitive integral rows");auto h=x.num().abs();if(h>hk)hk=h;}
        return height*hk;
    }
    bool verify(const Mat& candidate,const int_t& modulus)const{
        auto bound=integral_residual_bound(candidate);std::cout<<"invariant_certificate_bound_bits="<<bound.bits()<<std::endl;
        return bound<modulus;
    }
};
}
