# Executed by build.py, using O, edit and parent.
q=O/'src/reduced_d3/irr_projection.hpp';s=q.read_text()
s=s.replace('#include <set>','#include <set>\n#include "shared/workload_kernels.hpp"\n#include "shared/exact_operator_compression.hpp"')
s=s.replace('bool shoup=true;','size_t projection_raw_rows=0;bool shoup=true;',1)
s=s.replace('size_t n,ulong v)const{','size_t n,ulong v,bool dense=false)const{',1)
s=s.replace('for(size_t i=0;i<n;++i)if(src[i])dst[i]=nmod_add(dst[i],n_mulmod_shoup(v,src[i],pre,F.mod.n),F.mod);','''if(dense)for(size_t i=0;i<n;++i)dst[i]=nmod_add(dst[i],n_mulmod_shoup(v,src[i],pre,F.mod.n),F.mod);
            else for(size_t i=0;i<n;++i)if(src[i])dst[i]=nmod_add(dst[i],n_mulmod_shoup(v,src[i],pre,F.mod.n),F.mod);''')
# Density is a performance hint: the dense kernel also processes exact zeros.
s=s.replace('std::vector<ulong> candidates(count*stride);','''std::vector<ulong> candidates(count*stride);std::vector<unsigned char> dense(count);
                    std::vector<std::vector<ulong>> scratch;if(workload::kernels())scratch.assign(c.opt->pool.get_thread_count(),std::vector<ulong>(stride));''')
s=s.replace('auto candidate=candidates.data()+col*stride;std::vector<ulong> term_image(stride);','''auto candidate=candidates.data()+col*stride;std::vector<ulong> local;
                        if(!workload::kernels())local.resize(stride);auto& term_image=workload::kernels()?scratch[thread_id()]:local;''')
s=s.replace('c.axpy(candidate,term_image.data(),stride,v);\n                        }','''c.axpy(candidate,term_image.data(),stride,v);
                        }
                        if(workload::kernels())dense[col]=workload::prefer_dense(candidate,stride);''')
old='''                    irr_parallel::rows(M.nrow,&c.opt->pool,[&](size_t r){
                        auto dst=current.data()+(next.off[b]+mu*next.m[b]+r)*stride;size_t lo=0,hi=M[r].nnz();'''
new='''                    auto gather=[&](size_t r){
                        auto dst=current.data()+(next.off[b]+mu*next.m[b]+r)*stride;size_t lo=0,hi=M[r].nnz();'''
assert old in s;s=s.replace(old,new)
s=s.replace('''c.axpy(dst,candidate,stride,v);
                        }
                    });''','''c.axpy(dst,candidate,stride,v,dense[M[r](p)-first]);
                        }
                    };
                    if(workload::kernels())workload::balanced(M.nrow,&c.opt->pool,[&](size_t r){return M[r].nnz()*stride;},gather);
                    else irr_parallel::rows(M.nrow,&c.opt->pool,gather);''')
# Reuse per-worker projection scratch across all row tasks in a tile.
s=s.replace('''            irr_parallel::rows(S.nrow,&c.opt->pool,[&](size_t r){
                std::vector<ulong> sums(stride);size_t lo=0,hi=S[r].nnz();''','''            std::vector<std::vector<ulong>> scratch;if(workload::kernels())scratch.assign(c.opt->pool.get_thread_count(),std::vector<ulong>(stride));
            auto gather=[&](size_t r){
                std::vector<ulong> local;if(!workload::kernels())local.resize(stride);
                auto& sums=workload::kernels()?scratch[thread_id()]:local;std::fill(sums.begin(),sums.end(),0);size_t lo=0,hi=S[r].nnz();''')
s=s.replace('''rows[s].coefficients[r]=nmod_add(rows[s].coefficients[r],sums[s],c.F.mod);
            });''','''rows[s].coefficients[r]=nmod_add(rows[s].coefficients[r],sums[s],c.F.mod);
            };
            if(workload::kernels())workload::balanced(S.nrow,&c.opt->pool,[&](size_t r){return S[r].nnz()*stride;},gather);
            else irr_parallel::rows(S.nrow,&c.opt->pool,gather);''')
# A conservative feature space for convolution operators, independent of the
# unknown solution. Support omission is proved over EVERY supplied head copy.
where=s.index('inline equations slots(')
s=s[:where]+'''inline std::vector<size_t> convolution_targets(context& c,const std::vector<ulong>& head,
    const std::vector<ulong>& tail,size_t d,const equations& rows,bool enabled){
    std::vector<size_t> selected;size_t stride=c.samples*d;
    if(!enabled||!workload::compression()||d==1){selected.resize(rows.size());std::iota(selected.begin(),selected.end(),0);return selected;}
    size_t heads=head.size()/stride,tails=tail.size()/stride,audited=0;
    for(size_t sample=0;sample<c.samples;++sample){
        std::vector<size_t> powers;
        for(size_t j=0;j<d;++j){bool active=false;for(size_t i=0;i<heads&&!active;++i)active=head[i*stride+sample*d+j]!=0;if(active)powers.push_back(j);}
        exact_operator::Dense source(d,std::vector<ulong>(tails*powers.size()+1));
        for(size_t q=0;q<d;++q){for(size_t j=0;j<powers.size();++j)if(powers[j]<=q)
            for(size_t l=0;l<tails;++l)source[q][j*tails+l]=tail[l*stride+sample*d+q-powers[j]];
            source[q].back()=rows[sample*d+q].rhs;}
        auto basis=exact_operator::compress(source,c.F);audited+=basis.certified_entries;
        for(auto q:basis.selected)selected.push_back(sample*d+q);
    }
    // One exact zero equation preserves the existing stop accounting.
    if(selected.empty()&&!rows.empty())selected.push_back(0);
    std::cout<<"OPERATOR_COMPRESSION family=polynomial_convolution raw="<<rows.size()<<" selected="<<selected.size()<<" certified_entries="<<audited<<" affine_rhs_included=1"<<std::endl;
    return selected;
}
inline void select_equations(equations& rows,const std::vector<size_t>& selected){
    if(rows.size()==selected.size())return;equations compressed;compressed.reserve(selected.size());
    for(auto r:selected)compressed.push_back(std::move(rows[r]));rows=std::move(compressed);
}

''' +s[where:]
s=s.replace('inline equations slots(context& c,bool mrk,int metric=3){','inline equations slots(context& c,bool mrk,int metric=3,bool compress=false){')
start=s.index('inline equations slots(');end=s.index('// Keep the last k boundary',start)
part=s[start:end]
part=part.replace('    // Form each pair', '    c.projection_raw_rows=rows.size();auto targets=convolution_targets(c,head,tail,d,rows,compress);select_equations(rows,targets);\n    // Form each pair')
part=part.replace('for(size_t s=0;s<c.samples;++s)for(size_t q=0;q<d;++q){ulong x=0;', 'for(size_t z=0;z<targets.size();++z){size_t s=targets[z]/d,q=targets[z]%d;ulong x=0;')
part=part.replace('image[s*d+q]=nmod_add(image[s*d+q],','image[z]=nmod_add(image[z],')
part=part.replace('apply_rows(c,rows,stride,','apply_rows(c,rows,targets.size(),')
s=s[:start]+part+s[end:]
# Generic vector modular multiply-add also replaces repeated modular divisions
# in the tail recurrence; the germ and context use the same finite field.
s=s.replace('''std::fill(candidate.begin(),candidate.end(),0);for(auto [r,v]:M[col])for(size_t z=0;z<stride;++z){auto x=values[(current.off[b]+mu*current.m[b]+r)*stride+z];if(x)candidate[z]=germ.add(candidate[z],germ.mul(v,x));}''','''std::fill(candidate.begin(),candidate.end(),0);for(auto [r,v]:M[col]){
                            const ulong* src=values.data()+(current.off[b]+mu*current.m[b]+r)*stride;
                            if(workload::kernels())c.axpy(candidate.data(),src,stride,v);
                            else for(size_t z=0;z<stride;++z){auto x=src[z];if(x)candidate[z]=germ.add(candidate[z],germ.mul(v,x));}}''')
s=s.replace('''for(size_t z=0;z<stride;++z)if(term_image[z])next[r*stride+z]=germ.add(next[r*stride+z],germ.mul(v,term_image[z]));''','''if(workload::kernels())c.axpy(next.data()+r*stride,term_image.data(),stride,v);
                            else for(size_t z=0;z<stride;++z)if(term_image[z])next[r*stride+z]=germ.add(next[r*stride+z],germ.mul(v,term_image[z]));''')
s=s.replace('scratch[thread_id()]','scratch[workload::slot()]')
q.write_text(s)
# Use generic certified compression in the old collinear adapter; keep the
# original implementation available for ablations.
q=O/'src/shared/hexagon_collinear_tail_compression.hpp';s=q.read_text();pos=s.index('    struct Pivot')
s=s[:pos]+'''    if(workload::compression()){
        auto result=exact_operator::compress(source,F);auto selected=std::move(result.selected);size_t rank=selected.size();
        if(selected.empty())selected.push_back(0);
        std::cout<<"COL_TAIL_COMPRESSION raw_words="<<count<<" selected_words="<<selected.size()<<" affine_rank="<<rank<<" certificate_entries="<<result.certified_entries<<" rhs_columns="<<samples<<" complete_operator_preserved=1 generic=1"<<std::endl;
        return selected;
    }
''' +s[pos:];q.write_text(s)
# The direct original-coordinate discovery adapter formerly projected all
# suffix rows. Now the same certified feature-space compression applies.
q=O/'src/shared/direct_physical_rows.hpp';s=q.read_text();start=s.index('inline equations direct_collinear_tail_exact(')
part=s[start:];point=part.index('    apply_direct_invariant_columns(c,rows,stride,')
part=part[:point]+'''    if(workload::compression()){
        std::vector<unsigned char> active(la.dim);
        for(size_t b=0;b<6;++b)for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){bool used=false;
            for(size_t f=0;f<fa.m[b]&&!used;++f)for(size_t sample=0;sample<n;++sample)if(head[(fa.off[b]+mu*fa.m[b]+f)*n+sample]){used=true;break;}
            if(used)for(size_t l=0;l<la.m[b];++l)active[la.off[b]+mu*la.m[b]+l]=1;
        }
        auto selected=independent_collinear_suffix(c,la,rows,n,&active);std::vector<size_t> targets;
        for(size_t sample=0;sample<n;++sample)for(auto beta:selected)targets.push_back(sample*t+beta);
        select_equations(rows,targets);std::vector<size_t> selected_words;for(auto beta:selected)selected_words.push_back(c.tail_words[beta]);
        // The shared tail word layout stays immutable across adapters/batches.
        t=selected.size();stride=n*t;
        apply_direct_invariant_columns(c,rows,stride,[&](size_t b,size_t j,ulong* image){
            for(size_t mu=0;mu<irr::specs()[b].irrep_dimension;++mu){size_t i=fa.off[b]+mu*fa.m[b]+j/la.m[b],l=la.off[b]+mu*la.m[b]+j%la.m[b];
                for(size_t sample=0;sample<n;++sample)for(size_t beta=0;beta<t;++beta){auto v=c.exact_tail[l*q+selected_words[beta]];
                    auto x=nmod_mul(head[i*n+sample],v,c.F.mod);if(mu)x=nmod_mul(x,3,c.F.mod);image[sample*t+beta]=nmod_add(image[sample*t+beta],x,c.F.mod);}}
        });
        size_t nnz=0;for(const auto& row:rows)for(auto x:row.coefficients)nnz+=bool(x);
        c.record("collinear_direct_tail_exact",rows.size(),c.rows,nnz,(head.size()+c.exact_tail.size()+rows.size()*c.rows)*8,start);return rows;
    }
''' +part[point:];s=s[:start]+part
start=s.index('inline equations direct_slots(');end=s.index('// Keep the last k boundary',start);part=s[start:end]
part=part.replace('    // Form each pair','    auto targets=convolution_targets(c,head,tail,d,rows,true);select_equations(rows,targets);\n    // Form each pair')
part=part.replace('for(size_t s=0;s<c.samples;++s)for(size_t q=0;q<d;++q){ulong x=0;','for(size_t z=0;z<targets.size();++z){size_t s=targets[z]/d,q=targets[z]%d;ulong x=0;')
part=part.replace('image[s*d+q]=nmod_add(image[s*d+q],','image[z]=nmod_add(image[z],').replace('apply_direct_invariant_columns(c,rows,stride,','apply_direct_invariant_columns(c,rows,targets.size(),')
s=s[:start]+part+s[end:];q.write_text(s)
q=O/'src/reduced_d3/irr_physics_pipeline.hpp';s=q.read_text()
# Only training paths opt into operator compression. Independent scalar
# validation continues to generate the complete original target list.
s=s.replace('irr_physics::slots(c,r.at("operator")=="mrk"));c.cache_summary();','irr_physics::slots(c,r.at("operator")=="mrk",3,!validation));c.cache_summary();')
s=s.replace('irr_physics::slots(c,r.at("operator")=="mrk");\n','irr_physics::slots(c,r.at("operator")=="mrk",3,true);\n',1)
s=s.replace('            auto eq=r.at("operator")','            c.projection_raw_rows=0;auto eq=r.at("operator")')
s=s.replace('stagnant=rank==last_rank?stagnant+eq.size():0;','stagnant=rank==last_rank?stagnant+std::max(eq.size(),c.projection_raw_rows):0;')
q.write_text(s)
