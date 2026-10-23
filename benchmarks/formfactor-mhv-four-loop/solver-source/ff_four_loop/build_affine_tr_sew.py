from pathlib import Path
import shutil
import os
_here=Path(__file__).resolve().parent
_root=Path(os.environ.get('FF_BUILD_ROOT',_here/'_build'))
D=_root/'sew_affine_TR_build';D.mkdir(parents=True,exist_ok=True);S=_root/'sew_masked_build'
for p in S.glob('*.hpp'):shutil.copy2(p,D/p.name)
p=D/'tensor_kernel.hpp';s=p.read_text()
def rep(a,b):
 global s
 assert s.count(a)==1,(a[:60],s.count(a));s=s.replace(a,b)
rep('    Mat local;','    Mat local;\n    Mat affine;\n    int_t affine_num=0,affine_den=1;')
rep('ConstraintRows(const Tensor& expansion,Mat conditions,size_t next_dimension,bool backward=false,bool order_short_rows=false)','ConstraintRows(const Tensor& expansion,Mat conditions,size_t next_dimension,Mat rhs,bool backward=false,bool order_short_rows=false)')
rep('local(std::move(conditions)),slices','local(std::move(conditions)),affine(std::move(rhs)),slices')
rep('        for(auto& row:local.rows)vec_cancel_divisor(row);','''        require(affine.nrow==older&&affine.ncol==local.nrow,"affine rhs dimensions");
        std::vector<Q> equation_scale(local.nrow);for(size_t q=0;q<local.nrow;++q){Q first=local[q][0];vec_cancel_divisor(local[q]);equation_scale[q]=local[q][0]/first;}''')
rep('    size_t ncol()const{return previous*next;}','    size_t ncol()const{return previous*next+1;}')
needle='''            slices[old*letters+i[backward?0:1]].push_back(index(b),std::move(value));
        }
    }'''
rep(needle,'''            slices[old*letters+i[backward?0:1]].push_back(index(b),std::move(value));
        }
        for(size_t o=0;o<older;++o)for(size_t p=0;p<affine[o].nnz();++p){auto q=affine[o](p);affine[o][p]*=Q(denominators[o])*equation_scale[q];auto x=affine[o][p];if(x.num().abs()>affine_num)affine_num=x.num().abs();if(x.den()>affine_den)affine_den=x.den();if(!ff_row_grades.empty())require(ff_old_grades[o]==ff_row_grades[q],"nonzero affine RHS in masked grade");}
    }''')
rep('return hf*hc*hk;','return (hf*hc*affine_den+affine_num)*hk;')
rep('for(const auto& row:k.rows)for(auto c:row.index_span())++sizes[c/next];','for(const auto& row:k.rows)for(auto c:row.index_span())if(size_t(c)<previous*next)++sizes[c/next];')
rep('for(size_t s=0;s<k.nrow;++s)for(auto [c,x]:k[s])by_basis[c/next].push_back(index(s*next+c%next),x);','std::vector<ulong> affine_value(k.nrow);for(size_t s=0;s<k.nrow;++s)for(auto [c,x]:k[s]){if(size_t(c)==previous*next)affine_value[s]=x;else by_basis[c/next].push_back(index(s*next+c%next),x);}')
rep('''                        if(value)row.push_back(index(s),value);}
                    if(row.nnz())result[o].push_back(std::move(row));}''','''                        size_t q=&equation-cm.data();if(auto h=affine[old+o].find(index(q))){if(h->den()%fp.mod==0)throw certified_kernel::BadPrime("affine denominator");value=nmod_add(value,nmod_mul((*h)%fp.mod,affine_value[s],fp.mod),fp.mod);}
                        if(value)row.push_back(index(s),value);}
                    if(row.nnz())result[o].push_back(std::move(row));}''')
# Add the affine coefficient to exact short-row scans; count it in limits.
rep('''            if(limit!=SIZE_MAX){size_t support=0;bool skip=false;''','''            auto hv=affine[old].find(index(&equation-local.rows.data()));
            if(limit!=SIZE_MAX){size_t support=hv&&!relations.is_zero(previous*next)?1:0;bool skip=false;''')
rep('''            Vec row;
            for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;''','''            Vec row;if(hv&&!relations.is_zero(previous*next)){auto[root,scale]=relations.image(previous*next);row.push_back(index(root),(*hv)*scale);}
            for(auto [lj,a]:equation){size_t l=lj/next,j=lj%next;''')
rep('''            if(max_terms!=SIZE_MAX){size_t support=0;''','''            auto hv=affine[old].find(index(id%cm.size()));
            if(max_terms!=SIZE_MAX){size_t support=hv&&mapped[previous*next]>=0?1:0;''')
rep('''            if(row.nnz()<=min_terms){if(seen&&!row.nnz())(*seen)[id]=1;return ModVec{};}''','''            if(hv&&mapped[previous*next]>=0){if(hv->den()%fp.mod==0)throw certified_kernel::BadPrime("affine source denominator");auto v=nmod_mul((*hv)%fp.mod,scales[previous*next],fp.mod);if(v)row.push_back(mapped[previous*next],v);}
            if(row.nnz()>max_terms)return ModVec{};
            if(row.nnz()<=min_terms){if(seen&&!row.nnz())(*seen)[id]=1;return ModVec{};}''')
begin=s.index('class ConstraintRows {');end=s.index('\n};',begin)+4
aff=s[begin:end].replace('ConstraintRows','AffineConstraintRows')
# These inherited convenience methods have homogeneous-only semantics. The
# production staged solver uses residual() and integral_residual_bound() above.
# Disable them explicitly rather than allow accidental affine certification.
def disable_method(source, signature):
 start=source.index(signature); brace=source.index('{',start); depth=1; j=brace+1
 while depth:
  depth += (source[j]=='{')-(source[j]=='}');j+=1
 return source[:brace]+"{throw std::logic_error(\"AffineConstraintRows requires staged residual and height certification\");}"+source[j:]
aff=disable_method(aff,'    void profile(')
aff=disable_method(aff,'    bool verify(')

original=(S/'tensor_kernel.hpp').read_text();pos=original.index('\n};',original.index('class ConstraintRows {'))+4
p.write_text(original[:pos]+'\n'+aff+original[pos:])
s=(S/'sew_masked.cpp').read_text();s=s.replace('size_t nf=f.dim(0),nl=175;','size_t nf=f.dim(0),nl=92;')
s=s.replace('out/("suffix_galois"+std::to_string(k)+"_w2.wxf")','out/"actions"/("suffix_TR_"+std::array<std::string,5>{"galois1a","galois1b","galois2a","galois2b","galois3"}[k-1]+"_w2.wxf")')
s=s.replace('"suffix_"+name+"_w2.wxf"','"suffix_TR_"+name+"_w2.wxf"')
s=s.replace('auto local=read_matrix(out/"boundary_local_rref.wxf");','auto local=read_matrix(out/"boundary_local_TR_rref.wxf");auto rhs=read_matrix(out/("affine_TR_rhs_w"+std::to_string(w+2)+"_residual.wxf"));')
s=s.replace('tensor_kernel::ConstraintRows seam(f,std::move(local),nl);f.clear();auto full=staged_kernel::stack(invariants,seam);','tensor_kernel::AffineConstraintRows seam(f,std::move(local),nl,std::move(rhs));f.clear();ZeroColumn inv(invariants);auto full=staged_kernel::stack(inv,seam);')
wrapper='''
struct ZeroColumn {
 const staged_kernel::ProductInvarianceRows& source;
 explicit ZeroColumn(const staged_kernel::ProductInvarianceRows&s):source(s){}
 size_t ncol()const{return source.ncol()+1;}size_t nrow()const{return source.nrow();}
 template<class Emit>void operator()(structured_kernel::Relations&r,size_t l,Emit&&e)const{source(r,l,std::forward<Emit>(e));}
 template<class Emit>void modular(structured_kernel::Relations&r,const std::vector<I>&c,const field_t&f,Emit&&e,thread_pool*p=nullptr,size_t mx=SIZE_MAX,size_t mn=0,const streaming_echelon::Substitution*s=nullptr,std::vector<unsigned char>*seen=nullptr,unsigned shift=0)const{source.modular(r,c,f,std::forward<Emit>(e),p,mx,mn,s,seen,shift);}
 size_t residual_blocks()const{return source.residual_blocks();}
 template<class Emit>void residual(const staged_kernel::MM&k,const field_t&f,size_t b,Emit&&e,thread_pool*p)const{auto q=k;--q.ncol;for(auto&r:q.rows)if(r.nnz()&&size_t(r(r.nnz()-1))==q.ncol)r.resize(r.nnz()-1);source.residual(q,f,b,std::forward<Emit>(e),p);}
 int_t integral_residual_bound(const Mat&k)const{return source.integral_residual_bound(k);}
};
'''
s=s.replace('int main(int argc',wrapper+'int main(int argc')
s=s.replace('if(mode=="full"||mode=="masked")k=staged_kernel::solve(full,opt,{},&stats);\n else{GaloisPresolve source{full,galois};k=staged_kernel::solve(source,opt,{},&stats);}','require(mode=="masked","affine TR requires masked mode");k=staged_kernel::solve(full,opt,{},&stats);')
(D/'sew_affine_TR.cpp').write_text(s);print(D)
