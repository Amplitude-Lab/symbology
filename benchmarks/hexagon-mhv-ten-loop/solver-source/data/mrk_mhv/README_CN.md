# MHV MRK 输入

`bfkl_f_L<loop>_k<degree>.tsv` 是 MHV Fourier--Mellin 输出的最高权 binary-word symbol，
不是 NMHV 输入。`k` 是 MRK 对数次数，word 长度为 `2L-1-k`。
`manifest.tsv` 记录目前随包的 L1--L10 共 55 个 sector 的 weight 和 nnz。

生成步骤：

```bash
math -script MRK/prepare_mhv_input.wls /path/to/Ct.m /path/to/expansion_g10.txt /path/to/mhv_monomials.tsv 5
make -C MRK SPARSERREF_ROOT=../src/SparseRREF
./MRK/mrk_top_symbol --mhv -i /path/to/mhv_monomials.tsv -o /path/to/new_mhv_data
```

输入来自本项目已有的 BFKL eigenvalue (`Ct.m`) 和
[1407.3766](https://arxiv.org/abs/1407.3766) 的 impact-factor expansion。
MHV 路径不使用 helicity flip；原有生成器不加 `--mhv` 时仍保持 NMHV 行为。
L1--L5 结果已与 [1903.10890](https://arxiv.org/abs/1903.10890) 的
`hexMRKL1-7.m` 中 `FMres[L,k]/Pi` 逐项精确比较，15 项全部一致。

工作流读取全部低圈 sector，以构造 BDS-like 归一化带来的 shuffle 乘积。
因此计算 L 圈需要 `1..L` 的数据，不仅是 L 圈文件。

## L6--L10 扩展

高阶输入来自本地 all-order BFKL 系数级数，不截用仅到 g^10 的 ancillary。
所用 `Omega/Phi/Iminus` association 已保存在 `output/MHV/MRK/bfkl_series_L10.wxf`；
生成公式与原型来源见 `MRK/SELF_CONTAINED_L1_L9.md`。MHV 只读取 Omega 和 Phi。

```bash
math -script MRK/prepare_mhv_series_input.wls output/MHV/MRK/bfkl_series_L10.wxf /path/to/new_monomials.tsv 10
./MRK/mrk_top_symbol --mhv -i /path/to/new_monomials.tsv -o /path/to/new_output
```

本次单项式表保存在 `output/MHV/MRK/bfkl_monomials_L10.tsv`，共37,023项。
series WXF 的 SHA256 为 `D4105ED7FC1BED03C3AE90BFE214FAAB5E7846A27B68C3C06EE70F9C31740FFA`；
单项式表 SHA256 为 `ED30053366EFA36264DCB08BCE8A10A3A7E0A87ECE80712E49CD5F74FE3FA237`。
已有 L1--L7 TSV 逐字节保持不变。L6--L7 共13个 sector 已与 Lance 的同一份 MRK ancillary
逐项一致；L8--L10 尚未进行独立物理数据对照，不能把生成成功写成外部验证通过。
