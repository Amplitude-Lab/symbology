#include "partitioned_tensor_view.hpp"
using namespace symrep;
int main(int argc,char** argv){try{
    require(argc==2,"partitioned_tensor_probe OUTPUT_DIR");fs::path root=argv[1];fresh(root);
    for(size_t trial=0;trial<12;++trial){
        size_t inner=trial%3?20:0;Mat expected(0,inner);std::vector<PartitionedTensorView::Part> parts;
        for(size_t s=0;s<4;++s){Mat m((trial+s)%5,inner);
            for(size_t r=0;r<m.nrow;++r)for(size_t c=0;c<inner;++c)if((r+c+s+trial)%4==0){Q x(long(c)-7,long(r)+1);if(trial%2)x*=Q(int_t(2).pow(95ul)+int_t(3));if(x!=0)m[r].push_back(index(c),x);}
            auto file=root/(std::to_string(trial)+"-"+std::to_string(s)+".wxf");write_tensor_view(file,MatrixTensorView(m,{m.nrow,inner/4,4}));
            PartitionedTensorView::Part part{file,{}};for(const auto& row:m.rows)part.row_nonzeros.push_back(row.nnz());parts.push_back(std::move(part));
            for(auto& row:m.rows)expected.rows.push_back(std::move(row));expected.nrow=expected.rows.size();}
        PartitionedTensorView view(std::move(parts),{inner/4,4});auto file=root/(std::to_string(trial)+"-joined.wxf");write_tensor_view(file,view);
        verify_streamed_encoding(file,MatrixTensorView(expected,{expected.nrow,inner/4,4}));require(equal(flatten(read(file)),expected),"partitioned readback differs from unpartitioned matrix");
    }
    std::cout<<"PASS 12 partitioned tensor encodings, empty partitions/axes/rows, rationals and big integers"<<std::endl;
}catch(const std::exception& e){std::cerr<<e.what()<<std::endl;return 1;}}
