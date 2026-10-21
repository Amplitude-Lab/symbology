# SparseRREF 依赖版本

当前版本：`23ceb2e`，无本地补丁。相对 `ebb98f2` 只更新 `wxf_parser.h`：
`Encoder::push_array` 按转换回调的返回类型判断整数数组及元素宽度，修复
`rat_t` 张量以 `mma_layout=true` 导出时丢失 packed integer array 的问题。
发布导出现在逐个读回 WXF，对尺寸、索引和系数进行精确校验。
下面保留此前同步过程记录。

## 2026-09-21 更新

当前基线为本机上游干净工作树 `572e71c`，包含新的 `transpose_replace()`、
线程池排序、WXF 校验和张量移动/reshape 修复。同步原有11个依赖文件，未同步二进制。
另有一项本地兼容补丁：将 token 版 `sparse_tensor_read_wxf` 的默认参数移至首次
声明，避免 GCC 拒绝函数模板在后续定义中追加默认参数；不改变运行语义。
同步前依赖保存于工作区 `temp/sparserref_before_572e71c/`。
上述“仅本地”限制已由用户解除；正式压缩包纳入 ebb98f2，云端同步并重新编译。
详细影响见 handoff 的 `docs/SPARSERREF_572E71C_CN.md`。

## 历史基线

2026-09-20 将 `src/SparseRREF/` 已 vendored 的11个文件逐字节同步自
`D:/Documents/Git/Li-Xiang-Ideal/SparseRREF` 的干净工作树，HEAD 为 `eb0c7fe`
（2026-09-19）。未引入其本地 debug、可执行文件或 LibraryLink 二进制。
同步前副本保存在原开发 workspace 的 `temp/sparserref_before_irr_convert_20260920/`，
不属于本包的运行依赖。

上游已修复空 CSR 转 COO 时维度状态未转换的问题；不在工作流中保留该问题的 workaround。
当前 writer 接受 nnz=0 且各维度为正的稀疏张量，但明确拒绝包含零维度的张量。
直接重数工作流对此采用 `.wxf.shape` 文本记录，内容为 rank 后接各维度，
无 WXF 文件；读回时按这些维度构造零张量。旧的零维 WXF 仍可读，不修改历史数据。
见本包 `docs/IRR_CONVERSION_CN.md`。
