# 最大 helicity 的 framed-W OPE 数据

本目录不是完整的 T 展开，也不是可直接交给现有 sewing solver 的 BDS-like RHS。

## 内容

`framed_gluon_L1_T5.tsv` 至 `framed_gluon_L9_T5.tsv` 保存
`g^(2L) T^a (F^a+F^-a) log(T)^j z^p log(z)^q` 的系数：

- a=1,...,5，F=exp(i phi)，z=1/S。
- j=L-1,L-2,L-3，略去负数；无 zeta 部分。
- p=a,a+2，即每个 a 的前两个留数。没有保存完整的 S 依赖。
- 正负 helicity 系数相同，只保存正 helicity 一份；不能再给单份系数乘2。
- 表内保留声明范围中的零系数。九圈有360行，一至九圈合计2010行。
- 只要 a<=4 的内容即为本范围内的 T4 数据，a=5 是额外的 T5 数据。

TSV 列为 `t, helicity, log_t, z, log_z, coefficient`，后三个指数中的 z 表示 z 的幂。
`scope.json` 记录归一化与适用范围，`generation.tsv` 记录逐圈、逐 twist 的生成耗时。

## 理论范围

依据 [1407.1736 §4.2、附录 A](https://arxiv.org/html/1407.1736)，同 helicity 的
N 粒子项从 g^(2N^2) 开始。因此双粒子在 L 圈至多有 log(T)^(L-4)，
不会影响这里的最高三层；三粒子虽从九圈开始，但此时只影响更低 log 阶。
最大 helicity 扇区没有 scalar/fermion 污染，这一简化不能套用到其他 helicity。

若需要九圈的 log(T)^5 或更低阶，必须补充双粒子贡献；
若需要完整 T4/T5，还要包含非最大 helicity 的 gluon、scalar、fermion 等状态。
把当前单粒子结果称为完整 T4/T5 会漏项。

## 生成与验证

在包根目录运行 `OPE_MAX_LOOP=9 OPE_MAX_TWIST=5 math -script OPE/prepare_framed.wls`。
默认参数即9圈、T5、两个留数；计算脚本有时间和内存上限。
该生成器不读取任何物理 amplitude 或作者 ancillary。

2026-09-14 本机所有系数计算耗时合计2.90秒，Mathematica峰值约77.7 MiB。
九圈至T4耗时0.423秒，增加T5耗时0.190秒，全部九张表合计39293 bytes。
这些数字是生成理论系数的成本，不是将 OPE 施加到高圈 sewing 的成本。

- 一圈闭式 Wilson loop 公式检查 a=1,...,5，精确相等。
- 已有独立 bootstrap 的一至三圈物理 symbol，在
  p=1152921504606847009 上直接投影到 T5，再恢复 framed W 并提取最大 helicity；
  a=1,...,5 的所有可用 log 层、各两个留数全部相等，共60个 log(z) 多项式检查。
- 此前 a=1,2 的 NNLL 公式还通过作者三至六圈数据的精确检验。

## 与 BDS-like RHS 的区别

现有 `../calE_ope_L1..L8.tsv` 仍然只针对 T1/T2。
新表不能直接改名成 `calE_ope_L9.tsv`：归一化
`W=calE exp[g^2 (X-calE_1)]` 需要较低圈、相应展开阶数的数据。
现有低圈 germ 只到 T2，且只到七圈；九圈 T4/T5 的数值 BDS-like RHS 尚未准备。
现有 C++ OPE projector 也仍按 T2 截断，未在本次生成数据时静默改动。
