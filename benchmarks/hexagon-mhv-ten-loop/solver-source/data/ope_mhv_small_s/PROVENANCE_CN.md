# MHV small-S RHS 来源

2026-09-24由本地 `hexagon/OPE/flux_tube/generate_small_s.wls` 从论文公式生成，
L1--L9、twist1/2全log、S幂至10，归一化为BDS-like calE的普通symbol无周期部分。
不读取振幅结果或合作者的冻结RHS。复现副本在包内 `tools/ope_small_s/`。

`generator_manifest.json` 保留生成范围、核心源码SHA256、Mathematica版本及时间。
核心生成耗时326.82秒，WL内部峰值189280792 bytes。原W/calE和日志保留于工作区
`temp/mhv_ope_small_s_20260924/generated/`，不是求解的外部依赖。

表中role按端点幂划分：<=8为solve，>8为check。程序默认将全部T1行进一步转作check，
仅T2参与求解；可通过coverage配置改成T1/T2均参与。实际覆盖与测试见
`../../docs/MHV_SMALL_S_OPE_CN.md`。
