#set document(title: "Component.R 组件文档", author: "Dongdong Kong")
#set page(paper: "a4", margin: (x: 2.3cm, y: 2.5cm), numbering: "1")
#set text(
  font: ("New Computer Modern", "Source Han Sans SC"),
  size: 10.5pt,
  lang: "zh",
  region: "cn",
)
#set par(justify: true, leading: 0.85em, spacing: 1.15em)
#set heading(numbering: "1.1")
#show raw.where(block: true): set text(size: 9pt, lang: "en")
#show raw: set text(font: ("DejaVu Sans Mono", "Source Han Sans SC"))
#set math.equation(numbering: "(1)")
#let hline = table.hline

#align(center)[
  #text(size: 19pt, weight: "bold")[`Component.R` 组件文档] \
  #v(0.3em)
  #text(size: 9.5pt, fill: gray)[
    分段回归 · 去季节化 · `rpkg.template` v0.1.0.9000
  ]
]

#v(1em)

= 概述

本文档描述 `R/Component.R` 中四个函数的**实际行为**（而非设计意图）。
文档与 `tests/testthat/test-Component.R` 一一对应：每一条数学性质都有对应断言，
每一个已知缺陷都有编号（B1–B8）并在 §7 汇总。

#table(
  columns: (2.6cm, 1.5cm, 1fr),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*函数*], [*对外可见*], [*职责*]),
  hline(stroke: 0.4pt),
  [`RemoveSeasonal`], [否], [按日历月日分组求气候均值，返回去季节序列],
  [`lm_solve`], [否], [QR 最小二乘内核，返回系数、拟合、残差、SSE],
  [`LS_OneBreak`], [是（`@export`）], [单断点同斜率分段回归],
  [`LS_MultiBreaks`], [否], [多断点同斜率分段回归],
  hline(stroke: 0.8pt),
)

注：`NAMESPACE` 尚未重新生成，`LS_OneBreak` 在当前版本实际也未导出；
测试统一用 `rpkg.template:::` 取用四个函数，不受导出状态影响。

== 记号与约定

约定如下，后续各节不再重复：

- $n$：观测长度；$p$：设计阵列数。
- $X$：自变量，通常是时间索引 $1, dots, n$。
- $Y$：响应向量，长度 $n$。
- $t_i$：月日整数编码（MMDD），`101` 表示 1 月 1 日，`1231` 表示 12 月 31 日。
- $I_j$：第 $j$ 段末尾位置，$j = 1, dots, m$，要求 $I_m = n$ 且严格递增。
- $s_i$：观测 $i$ 所属段号，$s_i = op("min"){ j : i <= I_j }$。
- $I(A)$：示性函数，$A$ 成立取 1，否则取 0。

一致性约定：所有函数对输入**不做校验**（无 `stopifnot`、无维度检查），
错误输入一律以“静默错误结果”或 `.lm.fit` 的底层报错形式暴露。

= `RemoveSeasonal`

== 数学定义

设 $S_k = { i : t_i = k, bold(t)_i " 非缺失" }$，对 $t$ 中出现的每个互异取值 $k$（升序）：

$ mu_k = 1 / |S_k| sum_(i in S_k) y_i, quad a_i = y_i - mu_(t_i) $

即“按月日的气候态”去季节：同月日跨年取均值，减掉后得到异常序列 $a$。

== 接口

```r
RemoveSeasonal(y, t_mmdd)
```

返回值 `list(mu, anomaly)`：

- `mu`：长度等于 `length(sort(unique(t_mmdd)))`，按 MMDD 升序排列——月尺度最多 12 项，日尺度最多 366 项。
- `anomaly`：与 `y` 等长，逐元素减去所属组均值。

`mu` 与分组键的对应关系由**位置**隐式承载（无 `names`），见 B1。

== 实现要点

```r
inds <- sort(unique(t_mmdd))
mu <- numeric(length(inds))
anomaly <- numeric(length(y))
for (i in seq_along(inds)) {
  ind <- which(t_mmdd == inds[i])
  mu[i] <- mean(y[ind], na.rm = TRUE)
  anomaly[ind] <- y[ind] - mu[i]
}
```

- 时间复杂度 $O(n + k log k)$，$k$ 为互异月日数；单层循环，无矩阵运算。
- `anomaly` 预分配为 `numeric()`，故**未被任何分组命中的位置恒为 0**（B2）。

== 边界情形

#table(
  columns: (auto, 1fr),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*输入*], [*行为*]),
  hline(stroke: 0.4pt),
  [`y` 含 `NA`], [该点不参与均值（`na.rm = TRUE`），`anomaly` 对应位置为 `NA`],
  [某组 `y` 全为 `NA`], [`mean` 返回 `NaN`，该组 `anomaly` 全为 `NaN`],
  [某组仅 1 个观测], [`mu` 即该观测值，`anomaly` 为 0],
  [`t_mmdd` 含 `NA`], [该位置不进入任何组，`anomaly` 静默为 0，见 B2],
  [`t_mmdd` 无序], [输出 `mu` 仍按升序排列，不依赖输入顺序],
  hline(stroke: 0.8pt),
)

= `lm_solve`

== 数学定义

最小二乘：$hat(bold(beta)) = op("argmin")_bold(beta) norm(Y - D bold(beta))_2^2$，
返回 $hat(bold(beta))$、$hat(Y) = D hat(bold(beta))$、$e = Y - hat(Y)$、$"SSE" = sum_i e_i^2$。

== 接口

```r
lm_solve(design, response, only.SSE = FALSE, ...)
```

`only.SSE = TRUE` 时只返回 `list(SSE)`，供断点搜索此类“只比较 SSE”的场景省掉拟合向量。

== 实现要点

```r
l <- .lm.fit(design, as.matrix(response))
fitted <- design %*% l$coefficients
SSE <- sum(l$residuals^2)
```

- `.lm.fit` 走 LINPACK `dqrls`，复杂度 $O(n p^2)$；刻意绕开 `lm()` 的公式与模型框架开销。
- **不检查 `NA`/`NaN`/`Inf`**：`Y` 含缺失立刻报错，与 `RemoveSeasonal` 的 `na.rm` 语义不一致（B7）。
- 秩亏时 `.lm.fit` 在 `coefficients` 中留 `NA`，随后 `design %*% coef` 把 `NA` 扩散进 `fitted`（B4）。
- `...` 原样接收但**从未使用**（B5）。

= `LS_OneBreak`

== 数学定义

单断点、同斜率、截距跳变（断点 $I_c$ 为第一段末位置）：

$ Y_i = beta_0 + beta_1 X_i + delta I(i > I_c) + epsilon_i $ <eq:onebreak>

设计阵 $D = [bold(1), X, I(i > I_c)]$，$p = 3$；$delta$ 即两段截距之差。

== 接口

```r
LS_OneBreak(Y, X, Ic, ...)
```

`Ic = NA` 表示“无断点”，此时退回普通线性回归（$p = 2$）。`...` 转交 `lm_solve()`。

== 边界情形

#table(
  columns: (auto, 1fr),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*输入*], [*行为*]),
  hline(stroke: 0.4pt),
  [`Ic = NA`], [单段线性回归，`coef` 长度 2],
  [`Ic = 0`], [示性列全为 1，与截距列共线，`coef` 含 `NA`，见 B4],
  [`Ic = n`], [示性列全为 0，秩亏，`coef` 含 `NA`，见 B4],
  [`Ic` 非整数], [`rep(0, Ic)` 截断导致两段长度之和 $!= n$，`cbind` 报错],
  [`Y` 含 `NA`], [`.lm.fit` 报错],
  hline(stroke: 0.8pt),
)

参数文档写作 `T_anorm`，形参却叫 `X`；且此处 $X$ 未中心化，与 `LS_MultiBreaks` 不一致（B3、B6）。

= `LS_MultiBreaks`

== 数学定义

$m$ 段、同斜率、各段截距不同，段 1 为基线：

$ Y_i = beta_0 + beta_1 (X_i - overline(X)) + sum_(j=2)^m delta_j I(s_i = j) + epsilon_i $ <eq:multibreak>

设计阵 $D = [bold(1), X - overline(X), I(s = 2), dots, I(s = m)]$，$p = m + 1$。

系数含义：$beta_1$ 为全段共用斜率；$beta_0$ 为**段 1 在 $X = overline(X)$ 处的水平**；
$delta_j$ 为段 $j$ 相对段 1 的截距差，故 $delta_1 eq.not 0$ 不出现，段 1 由 $beta_0$ 承载。

== 实现要点

```r
n_breaks <- length(Ips) - 1
design <- cbind(1, X - mean(X))
for (i in seq_len(n_breaks)) {
  segment <- numeric(length(Y))
  segment[(Ips[i] + 1):Ips[i + 1]] <- 1
  design <- cbind(design, segment)
}
lm_solve(design, Y, ...)
```

- 时空复杂度 $O(n m^2)$；循环内反复 `cbind` 增长矩阵，$m$ 大时为 $O(n m^2)$ 次拷贝，可接受但不优雅。
- `Ips` 长度 1（即 $m = 1$）时 `seq_len(0)` 为空，退化为普通线性回归。

== 边界情形

#table(
  columns: (auto, 1fr),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*输入*], [*行为*]),
  hline(stroke: 0.4pt),
  [`Ips[m] = n`], [预期用法，全部观测被覆盖],
  [`Ips[m] < n`], [尾部示性列全 0，**静默归入段 1**，见 B8],
  [`Ips` 有重复值], [出现空段，对应列全 0，秩亏，`coef` 含 `NA`],
  [`Ips` 非递增], [示性列重叠或反向赋值，结果无意义，不报错],
  hline(stroke: 0.8pt),
)

= 参数化对比

两个分段函数解的是同一类模型，但**参数化不同**，系数不可跨函数直接比较。

#table(
  columns: (2.8cm, 1fr, 1fr),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*项目*], [`LS_OneBreak`], [`LS_MultiBreaks`]),
  hline(stroke: 0.4pt),
  [回归列], [$X$ 原值], [$X - overline(X)$],
  [断点表示], [标量 `Ic`], [向量 `Ips`],
  [段数], [至多 2], [$m$ 段，$m >= 1$],
  [系数], [$beta_0, beta_1, delta$], [$beta_0, beta_1, delta_2, dots, delta_m$],
  [拟合列空间], [同左], [同左（平移不改变列空间）],
  [导出的 `coef`], [不可直接对比], [不可直接对比],
  [拟合值 / SSE], [一致], [一致],
  hline(stroke: 0.8pt),
)

最后两行是测试的关键性质：$X$ 的平移不改变列空间，故单断点情形两函数
`fitted`、`SSE` 必须相等，而 `coef` 只差常数平移。

= 已知缺陷

按严重度排序，供后续修复时逐条销账。测试中对应位置均有注释。

#table(
  columns: (1.1cm, 1.9cm, 1fr, 1.9cm),
  stroke: none,
  inset: (x: 5pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*编号*], [*位置*], [*现象*], [*建议*]),
  hline(stroke: 0.4pt),
  [B1], [`RemoveSeasonal`], [`mu` 无 `names`，与 `sort(unique(t_mmdd))` 的对应关系只能靠位置猜], [给 `mu` 命名],
  [B2], [`RemoveSeasonal`], [`t_mmdd` 为 `NA` 处 `anomaly` 静默等于 0，看起来像“已去季节”], [置 `NA` 或报错],
  [B3], [`LS_*`], [一处用 $X$、一处用 $X - overline(X)$，同包内两套参数化], [统一中心化],
  [B4], [`LS_*`], [断点越界/退化 → 秩亏 → `coef` 为 `NA` → `fitted` 为 `NA`，静默], [`stopifnot(Ic > 0, Ic < n)`],
  [B5], [`lm_solve`], [`...` 传入后从未使用，`LS_*` 的 `...` 实为死参数], [删除 `...`],
  [B6], [`LS_OneBreak`], [`@param T_anorm` 与形参 `X` 不符；`LS_MultiBreaks` 缺 `@export`], [对齐文档与导出],
  [B7], [`lm_solve`], [`.lm.fit` 不接受 `NA`，与 `RemoveSeasonal` 的 `na.rm` 语义冲突], [入口统一 `na.omit` 或显式报错],
  [B8], [`LS_MultiBreaks`], [`Ips[m] < n` 时尾部观测静默归入段 1，断点定位错误不可见], [断言 `Ips[m] == length(Y)`],
  hline(stroke: 0.8pt),
)

= 测试计划

测试文件 `tests/testthat/test-Component.R`，共 15 个用例，只断言可观测的输出性质，不复刻实现细节。

#table(
  columns: (auto, 1.2cm, 1fr),
  stroke: none,
  inset: (x: 6pt, y: 4pt),
  hline(stroke: 0.8pt),
  table.header([*被测函数*], [*用例数*], [*断言要点*]),
  hline(stroke: 0.4pt),
  [`RemoveSeasonal`], [6],
  [分组均值数值正确；输出按 MMDD 升序对齐；各组异常值均值为 0；组内只有 1 个观测时异常值为 0；
   `y` 的 `NA` 不污染均值；`t_mmdd` 的 `NA` 不报错且长度不变],
  [`lm_solve`], [2], [`SSE = sum(e^2)`、`fitted + resi = Y`；`only.SSE` 只返回 `SSE`],
  [`LS_OneBreak`], [4],
  [无噪两段数据精确还原系数；断点模型 SSE 不大于无断点模型；
   `Ic = NA` 退化为线性回归；`Ic in {0, n}` 秩亏暴露为 `NA` 系数],
  [`LS_MultiBreaks`], [3],
  [三段数据还原斜率与截距跳变；单断点时与 `LS_OneBreak` 的 `fitted`/`SSE` 相等；
   `Ips` 只有一个元素时退化为线性回归],
  hline(stroke: 0.8pt),
)

运行方式：

```r
testthat::test_local()          # 或 devtools::test()
R CMD check --no-manual .       # CI 中由 R-CMD-check.yaml 执行
```

测试对四个函数的取用方式统一为 `rpkg.template:::f`，与导出状态解耦：

```r
RemoveSeasonal <- rpkg.template:::RemoveSeasonal
lm_solve       <- rpkg.template:::lm_solve
```

= 编译

```bash
typst compile docs/component.typ docs/component.pdf
```

文档刻意不 `#import` 外部模板，任意 Typst 0.11 及以上版本均可独立编译。
