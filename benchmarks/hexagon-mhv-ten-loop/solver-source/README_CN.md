# Hexagon 直接重数 Bootstrap

发布版 **v2，2026-09-28**：`hexagon_irr_bootstrap_v2_20260928.tar.gz`。
已完成的优化、兼容边界和后续工作见[版本说明](docs/RELEASE_V2_CN.md)。

2026-09-29 本地增补：共享 QQ/Fp 递推缓存、逐阶段 data/projection、十圈低 twist OPE 主库。
范围、生成器、阶段计划及未提交的集群配方见[十圈准备](docs/L10_PREPARATION_CN.md)。
这些本地修改尚未同步到上述归档或云端；扩展库完成状态以验收记录为准。

2026-09-28本地新增FEC逐块清缓存：`irr_fec ... --clean-block-cache`，或工作流设置
`FEC_CLEAN_BLOCK_CACHE=1`。默认关闭；块输出核验后持锁清理，保留重构/校验凭据并支持中断续算。
同时改用流式SHA256，避免核验大文件时额外占用整个文件两倍内存。
详见[逐块安全清理](docs/FEC_BLOCK_CACHE_CN.md)；版本标识`fec_block_cache_revision=20260928`。

2026-09-27 本地新增显式有限域贯通入口 `irr_pipeline` / `FINITE_FIELD_PIPELINE=1`，
支持 bulk、条件、左移交错执行，默认只重构最终 sewing。旧入口保留。
用法及缓存兼容边界见 [有限域贯通](docs/MODULAR_PIPELINE_CN.md)；未自动同步压缩包或云端。

2026-09-24本地新增MHV `ope_small_s`：默认使用完整log的 `T^2 F^(+/-2)` 至S^8定解，
T1及更高S幂留作检查；RHS备至L9。原 `ope` 保持不变。用法和验证见
[MHV small-S OPE](docs/MHV_SMALL_S_OPE_CN.md)。本版包含条件表和生成器，未测试九圈定解。

2026-09-23物理条件更新：`irr_conditions` 默认采用 adaptive sampling、scalar projection 与 Shoup 算术。可显式选择 `--projection tail-exact` 加速小k的collinear训练；MRK验证仍走独立scalar采样。`--resume-legacy` 用于保留旧阶段身份和缓存的续算，先逐项核验一枚旧素数再复用其余缓存。详见 [兼容说明](docs/LEGACY_RESUME_20260923_CN.md) 和 [七圈benchmark](docs/PHYSICS_V2_BENCHMARK_20260923_CN.md)。`make bin/irr_conditions_v2` 提供相同程序的备用入口。

本包在 parity × D3 的不可约表示重数空间直接构造 irrFEC、irrLEC 和 sewing，支持 MHV 和 joint E/Et NMHV。可按指定顺序施加 collinear、MRK，以及 MHV 的既有 T 展开 OPE。条件不足时输出完整仿射解空间，不预设物理解唯一。

已包含固定输入、C++ 源码和 SparseRREF `ebb98f2`（无本地补丁），不依赖相邻仓库、Mathematica 或 Python。系统需提供 C++20 编译器、GNU Make、Bash、flock、GMP、MPFR、FLINT、TBB、mimalloc。不要单独更新库头文件而不重新验证。

## 目录与约定

- `src/`：FEC、LEC、bulk sewing、转换、物理条件和缓存清理程序。
- `data/`：letter actions、NMHV Qbar seed、已经适配 D3 的 1-term/full-CA condition tensor 及第三维 metadata。
- `bin/`：可执行文件与 `.d`；`log/`：运行日志；`log/slurm/`：调度日志。
- `output/FEC/wXX/`：直接重数递推矩阵 `irrFEC_w_<block>.wxf`，不是展开后的 rank-3 FEC；零维块以同名 `.wxf.shape` 保存。2026-09-21 naming revision 2统一采用此命名，张量、基和CG约定不变。
- `output/LEC/MHV/wXX/`、`output/LEC/NMHV/wXX/`：`irrLEC_w_<block>.wxf`。
- `output/SEW_hpk/<MHV|NMHV>/bulk/`：六个共享 solution 轴的 `irrSEW_hpk_<block>.wxf`。
- `output/SEW_hpk/<MHV|NMHV>/collinear_mrk/` 等：按实际施加顺序累计的条件结果，每阶段独立拥有完整 origin、directions 和缓存。
- `release_<helicity>_wW_hpK/`：单独、自包含的发布目录，存可连乘的 FCE/LCE、sewing 和 WL 读取脚本。

六块固定顺序为 `pe_A1,pe_A2,pe_E,po_A1,po_A2,po_E`。递推矩阵是 `{new_copy,old_copy/channel}`，列按 input irrep → input copy → channel 排列。表示分量不进入主线大矩阵。CG 与 compact actions 可在 LEC 的 `common/` 中查看；action 作用于小表示分量，不作用于重数。

FEC 默认 integrability + ES + 1-term CA；LEC/sewing 默认 full CA。条件第三维的 irrep metadata 必须与张量配套，且采用包内固定 CG convention。递推的约束模板按完整方程清分母和公因子；kernel 逐行 primitive。LEC1 只在物理终端基中按整个 copy primitive，再转为重数坐标，不再次缩放 irrLEC1。

零维块使用 `.wxf.shape` 文本，`columns.tsv` 仍记录逻辑 WXF 文件名；非零尺寸、nnz=0 的块仍为 WXF。复制链时必须同时复制这些 sidecar 和 metadata。

## 编译与本地小测试

### Sewing 单侧预 RREF

`irr_sew ... --prereduce left|right|none` 控制对 `irrFEC.dlogmatIRR.irrLEC` 的单侧预约简。
省略时，last weight `k=1` 默认 `right`，其余默认 `left`，MHV/NMHV 使用相同规则。
shell 工作流可用 `SEW_PREREDUCE=left|right|none` 覆盖；不设置时由 C++ 选择默认值。

- `left/right`：每枚素数内分别约简 `irrFEC.dlogmatIRR` / `dlogmatIRR.irrLEC`，再接另一侧。
  不有理重构部分矩阵，只重构最终 sewing kernel；不常驻完整 QQ 接缝矩阵。
- `none`：原来的 QQ 一次装配，保留完整 QQ 接缝矩阵，每枚素数模约化后求 kernel。
- 方向加入运行身份，禁止在同一输出目录混用不同方向的缓存；旧目录可显式用 `none` 续算。
  不删除旧数据，换方向性能比较使用新目录。输出张量布局和共享 solution 轴不变。
- `--verify` 会额外装配原始 QQ 接缝以验证，只用于小规模检查；不要将其峰值当生产路径峰值。
- `manifest.json` 记录 `prereduce` 和 `constraint_arithmetic`。`summary.tsv` 的 rows/nnz
  在预约简模式下是本次求解所装配的有限域最终接缝的观测最大值；QQ constraint height 不适用，记为 `NA`。
  `rows.tsv` / `row_chart.tsv` 仍描述原始 CG 方程模板，不是素数相关的约简后行坐标。

本地验证和部署边界见 [预约简记录](docs/PREREDUCE_20260922_CN.md)。

2026-09-21 本地改用上游 `sparse_mat::transpose_replace()`，删除独立 helper。
正式压缩包已纳入此更新；更新范围、兼容性及验证记录见
[依赖更新说明](docs/SPARSERREF_572E71C_CN.md)。可单独运行 `make test-transpose`。

2026-09-22本地解压包新增低内存并行装配、分批CRT、接缝移动和物理投影优化，
正式发布版保留这些优化，但 sewing 装配仍用旧版静态调度。范围、排序内存注意事项及本地验证见
[并行改造记录](docs/PARALLEL_20260922_CN.md)。计算约定和缓存格式不变。

需要自定义库位置时，将 `config.mk.example` 复制为 `config.mk`。默认不指定 CPU 架构，避免在登录节点编译出计算节点不能运行的指令。编译并发默认仅 2，与求解线程数独立。

```bash
make -j2
THREADS=2 make test
```

测试从 seed 开始构造低阶链，检查 NMHV L3 的 4 维和 MHV L2 的 3 维空间、逐项参考数据、续算、sewing 移动往返及 release 校验和。测试输出保存在独立的 `output/smoke_*`，不会删除。

## 阶段与参数

`STAGE` 不区分大小写，可选 `FEC,LEC,SEW,CONDITIONS,SHIFT,RELEASE`。`SEW` 会先补齐所需 FEC/LEC 和 QQ bulk；`CONDITIONS` 只读取已有 sewing，不重新构造 bulk。`LOOP` 默认 3，`TAIL_WEIGHT` 默认 2，`HEAD_WEIGHT` 默认 `2*LOOP-TAIL_WEIGHT`；当前 **h,k 均须大于零，尚不支持 k=0**。

```bash
STAGE=FEC MAX_WEIGHT=16 bash run_workflow.sh
STAGE=LEC HELICITY=NMHV LEC_WEIGHT=6 bash run_workflow.sh
STAGE=SEW HELICITY=NMHV LOOP=8 TAIL_WEIGHT=2 bash run_workflow.sh
STAGE=SHIFT HELICITY=NMHV LOOP=8 TAIL_WEIGHT=2 TARGET_HEAD=12 bash run_workflow.sh
STAGE=RELEASE HELICITY=NMHV LOOP=8 TAIL_WEIGHT=2 TARGET_HEAD=12 DESTINATION="$PWD/release_NMHV_w16_12p4" bash run_workflow.sh
```

发布导出依赖 SparseRREF `23ceb2e` 或更新版本。`irr_convert export` 写出每个
`mma_layout=true` WXF 后都会立即读回并逐项核对；旧编码器生成的 WXF 即使通过
SHA256 与压缩包检查，也可能无法导入，不能作为有效发布数据。

SHIFT/RELEASE 读取原 split 的 sewing，然后在重数空间精确移动；自动补齐所需链，不构造完整高 weight symbol。`INPUT_SEW` 可指定输入，`DESTINATION` 指定输出。不会覆盖已有导出目录。

`OUTPUT_ROOT`、`FEC_ROOT`、`LEC_ROOT`、`SEW_ROOT`、`LOG_ROOT` 可分别覆盖。共享链会加文件锁，独立 sewing 任务可以共用 FEC/LEC。`BUILD=0` 跳过编译；同时提交多个任务时应先独立完成一次编译，再都设 `BUILD=0`。

条件可通过 `FEC_DLOGMAT/FEC_BLOCKS`、`LEC_DLOGMAT/LEC_BLOCKS`、`SEW_DLOGMAT/SEW_BLOCKS` 成对设置。更换条件必须使用新的 output 根，不能将新条件接在旧链上。兼容范围是固定 alphabet、群表示和 monomial codomain chart，不是任意用户 alphabet。

## Slurm

先在包目录中检查 `slurm/run_workflow.slurm` 的平台配置块，按集群启用 modules 和库路径。该文件故意允许编辑，不以 SHA256 限制 Makefile/Slurm 修改。程序运行位置取 `PACKAGE_ROOT` 或 `SLURM_SUBMIT_DIR`，不取 Slurm spool 路径。

从包目录提交，`log/slurm` 必须在 sbatch 前存在：

```bash
mkdir -p log/slurm
sbatch --chdir="$PWD" --export=ALL,STAGE=FEC,MAX_WEIGHT=16 slurm/run_workflow.slurm
sbatch --chdir="$PWD" --export=ALL,STAGE=SEW,HELICITY=NMHV,LOOP=8,TAIL_WEIGHT=2 slurm/run_workflow.slurm
sbatch --chdir="$PWD" --partition=amd_a8_768 --cpus-per-task=128 --mem=768G --export=ALL,STAGE=FEC,MAX_WEIGHT=18 slurm/run_workflow.slurm
sbatch --chdir="$PWD" --partition=amd_a8_384 --cpus-per-task=128 --mem=384G --export=ALL,STAGE=SEW,HELICITY=NMHV,LOOP=8,TAIL_WEIGHT=2 slurm/run_workflow.slurm
job=$(sbatch --parsable --chdir="$PWD" --export=ALL,STAGE=SEW,LOOP=8,TAIL_WEIGHT=2 slurm/run_workflow.slurm); sbatch --chdir="$PWD" --dependency=afterok:"$job" --export=ALL,STAGE=RELEASE,LOOP=8,TAIL_WEIGHT=2,TARGET_HEAD=12,BUILD=0 slurm/run_workflow.slurm
```

`amd_a8_768` 为 128 核、每核最多 6 GB；`amd_a8_384` 为 128 核、每核最多 3 GB。上例只是资源申请示例，不保证高 weight 可以在该内存内完成。没有设置 walltime；需要时自行用 `sbatch --time=...` 覆盖。日志含实际线程数、约束 rows/columns/nnz、重构模数、最终 height 和内存统计；Linux 的 peak RSS 与某些平台的 private commit 不是同一口径。

## 续算与重构

重复同一 FEC/LEC 命令会逐块校验完成标记并复用；也可提高最高 weight。SEW 重复运行会校验完整输出后返回，未完成时复用已保存的素数。缓存自动增加素数，不要求手工猜个数。FEC/LEC 的 RREF 采用矩阵高度界重构；sewing 沿用已验证的最终 kernel 重构及两个额外素数验证，二者不要混称为同一种证明标准。

完整块以 SHA256 `.ready` 标记确认；输入条件身份不符立即停止。中断留下无标记的文件会保留并报错，不会静默覆盖。有限域 cache 是昂贵数据，不自动删除。原始二进制 cache 不承诺跨架构/库版本兼容；WXF 与配套文本是跨平台交换格式。

## 可组合物理条件

```bash
STAGE=SEW HELICITY=MHV LOOP=6 TAIL_WEIGHT=2 CONDITIONS=collinear+mrk+ope bash run_workflow.sh
STAGE=CONDITIONS HELICITY=NMHV LOOP=6 TAIL_WEIGHT=2 INPUT_SEW="$PWD/output/SEW_10p2/NMHV/collinear" CONDITIONS=mrk BUILD=0 bash run_workflow.sh
STAGE=CONDITIONS HELICITY=MHV LOOP=6 TAIL_WEIGHT=2 CONDITIONS=mrk+collinear EXPORT_INTERMEDIATE=1 bash run_workflow.sh
```

`mrk` 不暗含 collinear，`ope` 不暗含前两者；目录名是累计顺序，例如 `mrk_collinear`。
从 `collinear` 续接 `mrk` 保存到 `collinear_mrk`。历史来自 `stage.tsv`，而不是猜目录名。
同一目录下配置或输入不同会报 identity 冲突，改用另一个 `OUTPUT_ROOT`，不覆盖旧结果。

默认只导出最终阶段 QQ；`EXPORT_INTERMEDIATE=1` 同时导出各阶段。中间 cache-only 阶段也含完整有限域仿射对象，可直接作为后续输入；缺少某枚素数时按父链补算。父路径相对当前阶段保存，整体搬运目录不需改变它；若单独搬运，应保持父对象可访问。

`MODULAR_PRIMES=auto` 默认每增加一枚重构素数尝试重构，另取两枚验证；也可指定正整数作为固定预算。缓存内部 origin 是最后一行，不能读取 benchmark 中 origin-first 的二进制文件。跨素数核对小系统的 pivot profile；遇到不同 profile 或不相容方程会停止并保留缓存，不静默混合。

`SAMPLES` 默认至少 384，并至少为输入方向数的两倍加 32；续接沿用父阶段样本数。`SAMPLE_BATCH` 默认 32，控制投影缓冲，而非物理方程覆盖域。`VALIDATION_SAMPLES` 默认 32，使用独立确定性 seed 检查完整 origin 和所有剩余方向。collinear/MRK 这里使用有限域 slot 采样，**不是符号恒等式证明**；失败时增加样本并用新输出目录。两枚独立素数验证有理重构，不等于矩阵高度界的严格证明。

支持 `PHYSICS_CACHE_MIB`（默认 1024，设为 0 禁用；32768 为32 GiB）：缓存当前素数下的
collinear/MRK 递推矩阵，避免每个采样批次重读 WXF 和取模。缓存按分配容量计费，
超过预算的层仍流式读取，换素数即清空；这不是进程总内存上限，读取中的 QQ 张量、
投影缓冲和 sewing 另占内存。输入文件在一次进程运行中必须保持不变。
日志输出采样批次进度、各阶段耗时、缓存字节数及进程峰值 RSS。
缓存不足明确打印 `NOT_STORED`，每批汇总命中数和剩余容量。移动 sewing 默认限频输出详细阶段，
`SHIFT_VERBOSE_DETAIL=1` 恢复完整日志，详见 `docs/LOGGING_20260923_CN.md`。

对当前 `X=X0+z S`，求解 `P(S) z=b-P(X0)`，并在模 p 下立即生成 `X0new,Snew`。
最后重构的是这些完整 sewing 块，不是限制矩阵的各层因子。不会在乘回后作 RREF、排序或清公因子；尤其不能缩放非齐次代表元。
流式 CG/channel 投影只生成必要小表示分量，slot 按批次处理，OPE 用末字母积分，不是 head/tail germ 相乘。

输出方向为 `irrSEW_hpk_<block>.wxf`，维度 `{solution,first_copy,last_copy}`；origin 为 `irrSEW_hpk_<block>_origin.wxf`，维度 `{first_copy,last_copy}`。
没有齐次自由度时只发布 origin，六块共享同一 solution 轴。每阶段 `manifest.json` 标明 QQ 完成或 cache-only，`summary.tsv` 记录维数、nnz、height、模数和内存。

可选 `CONDITION_PLAN=conditions.tsv`，每行 `label operator config_path`，空配置用 `-`。
配置是两列、tab 分隔，支持 `samples`、`validation_samples`、`seed`、`coverage`；OPE 的 coverage 对 MHV 为 `T1T2-top3`，对 NMHV 为 `NMHV-T1T2-all`。
重复算子不会被折叠，目录例如 `collinear_ope_ope`，label 只记录在 metadata 中。

物理数据：NMHV/MHV MRK 包含至 L10；原 MHV `ope` 的 BDS-like RHS 至 L8，仅 T1/T2、最高三个 log(T) 层和既有 1/S 留数。新增 MHV `ope_small_s` 提供至 L9 的全 log 条件，详见上述 small-S 文档。NMHV 的 (1111) T1/T2 全 log RHS 在 `data/ope_nmhv/`，备至 L9，本地验证至 L6；使用方式仍为 `CONDITIONS=ope`，见 [接入与验证](docs/NMHV_OPE_20260924_CN.md)。更高 T framed W 是研究材料，**不作为本程序 RHS**；不支持的范围报错，不静默截断。

从物理阶段导出：

```bash
STAGE=RELEASE HELICITY=MHV LOOP=6 TAIL_WEIGHT=2 TARGET_HEAD=9 INPUT_SEW="$PWD/output/SEW_10p2/MHV/collinear_mrk_ope" DESTINATION="$PWD/release_MHV_w12_9p3" bash run_workflow.sh
sbatch --chdir="$PWD" --export=ALL,STAGE=SEW,HELICITY=NMHV,LOOP=8,TAIL_WEIGHT=2,CONDITIONS=collinear+mrk,MODULAR_PRIMES=auto slurm/run_workflow.slurm
```

## 有限域缓存清理

```bash
bash clean_cache.sh
bash clean_cache.sh --apply
OUTPUT_ROOT=/path/to/output bash clean_cache.sh --apply
```

第一行仅预览；第二行一键清理所有**凭据完整且校验通过**的已导出 QQ 对象的有限域缓存，包括本包新生成的 FEC/LEC/bulk/物理阶段。
只删缓存二进制及其标记，保留张量、日志、重构记录和输入身份。未完成、cache-only、校验失败、正在使用的目录均跳过；符号链接/junction/共享 hardlink 不跟随、不清理。
凭据为 `qq_complete.tsv`（完整仿射物理结果）或 `qq_cache_safe.tsv`（六块全部完成的链/bulk）。导入的旧数据没有凭据时宁可保留缓存，不凭同名 WXF 猜完成。
凭据同时密封 `cache_files.tsv` 的精确文件清单；不按通配文件名前缀删除。混入别的计算文件、缺失标记或截断文件会导致该阶段保留，锁位置由生产者记录，不靠目录名猜。
清理过的完整 QQ 仍可续接；后续阶段需要它时直接模 p 读取 QQ，不依赖已删除的 cache。清理脚本不会自动清理未知历史工作流。
`KEEP_MODULAR_CACHE=0` 可在 CONDITIONS 成功结束后对当前 split/helicity 调用同一清理器；默认 `1` 保留。失败或仅 cache-only 的阶段无论如何都保留，以免丢失尚无 QQ 替代品的计算。

每次调用都有独立日志。续算仍会重建少量模板/装配矩阵；完成的 RREF 不必重复。请保留第一次运行日志，续算汇总中的当前装配统计不能代替原始峰值。

## 尚未重构的缓存左移

可用 `bin/irr_convert shift-cache NMHV H K INPUT OUTPUT --head NEW_H --reconstruct`，
直接逐素数左移，再尝试目标坐标下的有理重构；默认 recurrence 路径仍是 `output/FEC`
和 `output/LEC/NMHV`，可以通过 `--fec-root` / `--lec-root` 覆盖。
不加 `--reconstruct` 时只转换缓存。原缓存保持不变，重复命令可续算。
完整约定与验证见 [SHIFT_CACHE_CN.md](docs/SHIFT_CACHE_CN.md)。

## Release 输出

默认文件名为 `hexagon_coproduct_entries_first_w.wxf`、`hexagon_<mhv|nmhv>_coproduct_entries_last_w.wxf`、`hexagon_<mhv|nmhv>_sewing_hpk.wxf`。`FIRST_NAME/LAST_NAME/SEW_NAME/ORIGIN_NAME` 可使用 `{w},{h},{k},{weight},{loop},{helicity}` 模板，只接受目录内 `.wxf` 文件名。

```bash
STAGE=RELEASE LOOP=3 TAIL_WEIGHT=2 TARGET_HEAD=3 SEW_NAME='SEW_{h}p{k}.wxf' DESTINATION="$PWD/release_custom" bash run_workflow.sh
```

齐次空间保留 `{solution,first_basis,last_basis}`，即使只有一个 solution 也不自动删轴。通过 `DIRECTIONS`/`ORIGIN` 导入外部仿射系数时，代表元单独输出 `{first_basis,last_basis}`，零维齐次部分不发布。只有显式 `UNIQUE=1` 才将单个方向作为无 solution 轴的结果；该开关不意味着程序替你证明它是物理解。

### WXF 格式重导出

`irr_convert export` 使用 `mma_layout=true`。若要与旧本地发布文件的二进制格式一致，
可以用 `make bin/reencode_wxf` 编译转换程序，再执行
`bin/reencode_wxf input.wxf output.wxf`。它使用 `mma_layout=false`，写出后读回并精确
比较所有索引和有理数值；拒绝覆盖已有文件。此转换不做换基或归一化。

已合并的 L2--L9、40个WXF发布目录可用以下命令整体转换，同时更新文件大小和校验清单：

```bash
python3 tools/reencode_release.py /path/to/hexagon_symbols_weight_18_mma /path/to/hexagon_symbols_weight_18 /path/to/local_reference bin/reencode_wxf
```

第三个参数是本地参照发布目录，其中同名WXF须逐字节一致；对比结果保存在新目录的
上级 `local_comparison.json`。来源计算记录保留，新的 `manifest.json` 另记编码转换来源。

## 接入已有数据与算法说明

见 `docs/IMPORT_CN.md`。`docs/IRR_*_CN.md` 是本地实现的算法说明快照，里面的旧路径仅用于说明来源；本包运行入口和路径以本 README 为准。

验证记录见 `docs/VALIDATION_PHYSICS_CN.md`。尚未在真实 Slurm 集群执行本包；提交高阶任务前先跑小测试。
