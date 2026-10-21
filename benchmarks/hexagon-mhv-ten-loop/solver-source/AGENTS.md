# Agent 交接

- 2026-09-29本地增加OPE逐端点幂秩报告，适用于irr_conditions与irr_pipeline；投影一次、
  分批增量消元，不是降低jet截断的自适应生成器。--ope-rank-only可用两枚素数研究已有
  pre-OPE有理仿射空间，不重构新的sewing。九圈MHV/NMHV实验脚本已准备，未提交/部署。
  六圈真实数据与历史sewing逐项一致；详情docs/OPE_RANK_CN.md。

- 2026-09-29本地十圈准备：pipeline共享只读QQ/Fp递推，`--qq-cache-mib`默认0；每阶段可独立
  配置data/projection，库身份及完整性在计算前检查。六圈新旧OPE QQ六块逐字节相同。
  `data/ope_L10_main`已完整，不覆盖旧库；扩展库以生成manifest为准，partial严禁作为完整RHS。
  本地程序已更新，归档/云端未同步，未提交任务；见docs/L10_PREPARATION_CN.md与工作区实施记录。

- 2026-09-28发布v2：正式归档名hexagon_irr_bootstrap_v2_20260928.tar.gz。
  已批准优化均已部署；下文各日期的“本地/未同步”描述是实施当时的历史状态，不代表当前状态。
  新版本仅整理版本/日期/说明，数值源码与当日逐块清理部署版相同，输出及缓存schema不因包版本改变。
  后续研究性优化和未支持的功能见docs/RELEASE_V2_CN.md，不得声称不存在任何改进空间。

- 2026-09-28本地新增FEC显式`--clean-block-cache`/`FEC_CLEAN_BLOCK_CACHE=1`，默认保留缓存。
  每块有理输出核验后原子保存授权凭据并逐文件清理；旧整阶清理兼容，未完成块保留。
  仅此显式开关允许生产者持锁执行块级清理，不得在外部绕开活动锁。SHA256已改成流式。
  见docs/FEC_BLOCK_CACHE_CN.md。云端部署状态以工作区remote/TASK_BRIEFINGS.md为准。

- 2026-09-27：本地新增 irr_pipeline，以原 bulk/physics/shift 算法贯通有限域，默认仅最终 QQ。
  旧入口保留；新 bulk 缓存为 kernel-free-coordinate-v1，禁止与旧 QQ primitive/sorted 基混用。
  首枚素数复用仅为明确无状态的 bulk evaluator 启用；其他调用默认保持原行为。
  用法与验证见 docs/MODULAR_PIPELINE_CN.md、docs/LOCAL_OPTIMIZATION_RESULTS_20260927_CN.md。
  本轮未同步压缩包/云端，不改正在运行的 irrFEC18。

- 2026-09-26：发布导出依赖 SparseRREF `23ceb2e` 的 `push_array` 修复；
  `irr_convert export` 对每个 WXF 写出后立即读回，核对尺寸、索引和系数。
  旧版生成的 weight-18 发布目录和归档含无效 WXF，不能沿用或仅重新压缩。

- 2026-09-24本地接入MHV `ope_small_s`，默认T2全log至S^8、T1及S更高幂留出验证。
  旧 `ope` 和旧缓存身份不变；RHS至L9，详见docs/MHV_SMALL_S_OPE_CN.md。
  本版包含条件表和生成器，部署凭据见工作区remote/TASK_BRIEFINGS.md；不得把L6回归或RHS生成完成称为九圈定解。

- 2026-09-24：本地 `irr_convert shift-cache` 可直接左移已完成的有限域缓存，再可选重构；
  保持 solution/origin 约定、不做 presentation RREF。`--resume-shifted-cache` 可将唯一 OPE
  缓存接入新 split 续算，首枚原生重算逐项比较。用法、限制与六圈验证见 `docs/SHIFT_CACHE_CN.md`；
  部署状态见工作区 `remote/TASK_BRIEFINGS.md`。

- 2026-09-24：本地源码与 `bin/irr_conditions` 已接入 NMHV (1111) T1/T2 全 log OPE；
  `data/ope_nmhv/` RHS 备至 L9，L1--L6 已测试。详见 `docs/NMHV_OPE_20260924_CN.md`。
  后续压缩包/云端部署状态见工作区 `remote/TASK_BRIEFINGS.md`。高 twist 的部分 framed W 数据尚未接入，不得当作 BDS-like RHS。

- 2026-09-23新版物理投影与旧缓存兼容已纳入正式源码；irr_conditions_v2是irr_conditions的等价可执行入口。默认adaptive+scalar+Shoup，tail-exact显式选择且仅用于collinear训练(k<=3)。旧缓存必须显式--resume-legacy，首枚逐项核验，保留旧身份及全部数据。部署/作业状态见remote/TASK_BRIEFINGS.md；说明及七圈对照见docs/LEGACY_RESUME_20260923_CN.md和docs/PHYSICS_V2_BENCHMARK_20260923_CN.md。下方记录均为历史状态。

- 2026-09-23用户要求将下述本地物理条件优化和日志监控一并同步正式压缩包、云端源码及irr_conditions/irr_convert。实际部署与编译完成状态见工作区remote/TASK_BRIEFINGS.md；不改变历史输出、有限域缓存或任务脚本，不自动提交作业。下列“仅本地”是部署前历史记录。

- 2026-09-23 本地新增缓存逐阶stored/NOT_STORED与逐批命中/剩余容量监控；移动sewing日志默认限频，SHIFT_VERBOSE_DETAIL=1恢复详细日志。见docs/LOGGING_20260923_CN.md。未自动同步压缩包或云端。

- 2026-09-23 加物理条件的优化仅在本地：有界有限域递推缓存、分批并行前向消元/单次回代、MRK RHS shuffle 优化及进度日志。见 `docs/PHYSICS_OPTIMIZATION_20260923_CN.md`。不要误认为已同步正式 tar.gz 或云端；不得影响正在运行的九圈任务。

- 2026-09-22 shift_one改成逐行按目标坐标即时累加，不存全部重复乘积再排序，并输出阶段/solution/内存日志。保留精确回乘和solution顺序。见docs/SHIFT_ACCUMULATION_20260922_CN.md；云端部署状态以工作区remote/TASK_BRIEFINGS.md为准。

这是 parity × D3 直接重数包，含 bulk 与可组合的 collinear/MRK/OPE，不支持 k=0。
先读 README_CN.md 和 docs/IMPORT_CN.md。不要把旧 R-sew、展开 FEC 或旧 last-0 基接进来。

- 固定 CG/channel convention 和 physical-copy-primitive-v1 seed；不得逐分量单独缩放 E irrep。
- 六块共享 solution 轴；移动 sewing 不重选 solution basis。
- 保留所有完成块、有限域缓存和日志。不得为了重跑删除昂贵数据。
- 计算目录统一 output/，二进制 bin/，日志 log/，release 独立目录。
- first递推文件统一为 `irrFEC_<weight>_<irrep>.wxf`，不再使用原型前缀；零维块、完成标记及所有读写接口遵循相同命名。更名不改变张量或CG convention。
- sewing 路径为 `output/SEW_hpk/<MHV|NMHV>/<ordered_conditions>/`。
- 每阶段必须拥有完整仿射对象；origin 不能拖到后续 MRK 阶段才生成。
- 每枚素数内先求限制再乘回 sewing，重构最终乘积；禁止附加 presentation RREF。
- 有限域缓存内部 origin 是最后一行，不是 benchmark packed 文件的第一行。
- `clean_cache.sh` 默认预览，只清理有核验通过的完整 QQ 凭据的缓存；不要绕开锁和链接检查。
- MHV OPE 使用 BDS-like T1/T2 top-three-log；NMHV 使用 (1111) T1/T2 全 log 的有理 symbol 部分。
  不把 framed W 表直接当 RHS，不把缺失系数当零。
- 修改 Makefile、Slurm 不需要重新生成数据 hash；但更改输入/算法需要新的输出链或明确的回归验证。
- 2026-09-22低内存并行发布版仅将sewing装配恢复为旧版静态detach_loop，其余优化保留。见docs/PARALLEL_20260922_CN.md；不改变数据/缓存格式。不要给张量transpose_replace盲目传pool：上游并行sort_indices会另分配索引和值缓冲；默认无pool的转置已经排序，勿重复排序。并行任务不得嵌套等待同一个线程池，也不得并发多个大型素数消元以换取速度。
- 不沿符号链接清理任何目标。优先使用独立新目录作验证。
- 2026-09-22新增sewing `--prereduce left|right|none`，k=1默认right，其余left；shell用SEW_PREREDUCE覆盖。预约简在每枚素数内进行，仅重构最终kernel；none保留QQ一次装配。此修订纳入正式压缩包，云端部署记录见工作区remote/TASK_BRIEFINGS.md。不能把--verify额外QQ装配的RSS当生产路径峰值；变更方向必须使用独立输出目录，旧目录用none续算。
