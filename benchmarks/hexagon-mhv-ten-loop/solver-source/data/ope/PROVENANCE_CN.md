# OPE 输入来源

`kernels_F2/F3.tsv`、`X_F2/F3.tsv` 由包内 `OPE/prepare_kernels.wls` 生成。
`calE_ope_L1..L8.tsv` 由 `OPE/prepare_rhs.wls` 生成，无作者 ancillary 依赖。
L1--L6 的已有表取两层 log，L7--L8 表取三层；后者新增公式见 `OPE/gluon_ope_nnll.wl`。

低圈 calE 展开来自本地独立 bootstrap。下列是来源文件 SHA256，仅作为来源标识，
不是运行时依赖的路径；不需要把原 hexagon 目录一起复制。

| L | 来源条件 | 来源物理张量 SHA256 |
| --- | --- | --- |
| 1 | collinear+MRK | C20A9A11BCD75D6173BAC85EC99CF9266C47CBFDC648659552BB425523235D37 |
| 2 | collinear+MRK | D270CA6943A5838B3EBF5945EF5DB438B045D6334474083B5E903338F470333E |
| 3 | collinear+MRK | 016107CB28E881232E5162286A5020AA8EFFB8AE207F47393476745D177966A9 |
| 4 | collinear+MRK | 10C916CECCDFF98B832B6C8AB117AD9658DF784DC6466F285025FF940E44A2C9 |
| 5 | collinear+MRK | A481BB2EBAB4E3062CE0437EA162651C5DA66858154A4F982F64EB392D9E2706 |
| 6 | collinear+MRK+OPE | BF83EAF580696EC164AE5761DF6930E0F304CADE2619BCD5FDD08A2C41640A23 |
| 7 | collinear+MRK+三层 OPE，FEC_12 x LEC_2 | 34C75BA61EF3F260176BD18E4DC62EAAA25882B0EA2221E1EF41C36DC206E713 |

L1--L6 每个 F=2,3 的有理展开由六个31-bit素数重构，另两个素数验证。
保存格式：Association，key 为2、3，value 为 T,z,lt,lz 的有理 Laurent polynomial。
L1–L6 每个 F 的非零单项数分别为 9、35、82、149、238、346。
用新包的 `OPE/export_germ.cpp` 与 `OPE/import_lower_germ.wls` 可以从重新计算的
物理 SEW 再生；不同 sewing 的输入坐标不同，但正确展开应相同。

L6 原型独立求解后，已和 Lance 的 symbol 在两个新素数上逐项对照；作者数据仅用于
事后检验，未用作 OPE 约束输入。低圈归一化依赖不等于使用本圈物理答案。

七圈物理张量使用3枚约60-bit素数重构，另2枚验证；其 F=2,3 展开各473项，
由包内 `OPE/export_germ.cpp` 使用2枚约60-bit素数重构、另2枚验证。
七圈物理解在确定自身时只使用 L1--L6 归一化输入；其展开仅用于八圈 RHS。
`lower_calE_L7.wxf` SHA256：
`1E952AA8C47CDB84211BD5886C4C0AD50FF91953207B372C633917954859D350`。
`calE_ope_L8.tsv` SHA256：
`D1B235E8065E386C56A832C86F193094ED7A1D28B5CB2EC13ABBEC5C1BF8CCF0`。
