# Mathematical definition of undirected doublet-pair enrichment

# 无方向 Doublet Pair 富集的数学定义

Last updated: 2026-06-23

最后更新：2026-06-23

This note defines the observed-to-expected ratio used by
`stg_pair_enrichment()` and separates the scientific estimator from the
finite value used only for heatmap display.

本说明定义 `stg_pair_enrichment()` 使用的观察/期望比，并严格区分科学估计量与
仅用于热图显示的有限数值。

---

## 1. Statistical unit and notation

## 1. 统计单位与符号

Let \(r\) index an independent ROI or biological replicate, \(c\) index a
condition, and \(i,j,k\) index cell types. Within ROI \(r\), define:

令 \(r\) 表示独立 ROI 或生物学重复，\(c\) 表示疾病组，\(i,j,k\) 表示细胞类型。
在 ROI \(r\) 内定义：

| Symbol | Definition | 中文定义 |
|---|---|---|
| \(S_{ir}\) | number of singlets assigned to type \(i\) | 类型 \(i\) 的 singlet 数 |
| \(S_r=\sum_kS_{kr}\) | total singlet count | singlet 总数 |
| \(p_{ir}=S_{ir}/S_r\) | empirical singlet proportion | singlet 经验比例 |
| \(D_r\) | total number of certain heterotypic doublets | 确定异型 doublet 总数 |
| \(O_{ijr}\) | observed count of undirected pair \(\{i,j\}\) | 无方向 pair \(\{i,j\}\) 的观察数 |

The braces \(\{i,j\}\) indicate an unordered pair:
\(\{i,j\}=\{j,i\}\). Therefore, each biological pair appears once in the
triangular heatmap.

花括号 \(\{i,j\}\) 表示无序 pair：\(\{i,j\}=\{j,i\}\)。因此每个生物学 pair
在三角热图中只出现一次。

---

## 2. Conditional heterotypic null model

## 2. 条件化异型配对零模型

Consider two independent draws \(X,Y\) from the ROI singlet composition.
Before conditioning,

考虑从 ROI singlet 组成中独立抽取 \(X,Y\) 两次。在条件化之前，

$$
P(X=i,Y=j)=p_{ir}p_{jr}.
$$

For \(i\ne j\), the unordered event \(\{i,j\}\) contains two mutually
exclusive ordered outcomes, \((i,j)\) and \((j,i)\). Thus,

当 \(i\ne j\) 时，无序事件 \(\{i,j\}\) 包含两个互斥的有序结果
\((i,j)\) 与 \((j,i)\)，因此

$$
P(\{i,j\})=2p_{ir}p_{jr}.
$$

The probability that two draws have different types is

两次抽取属于不同类型的概率为

$$
H_r=P(X\ne Y)
=1-P(X=Y)
=1-\sum_kp_{kr}^2.
$$

Because the analyzed RCTD class consists of heterotypic doublets, the relevant
sample space is \(X\ne Y\). The conditional pair probability is therefore

由于分析对象是 RCTD 异型 doublet，相关样本空间为 \(X\ne Y\)。因此条件化
pair 概率为

$$
q_{ijr}
=P(\{i,j\}\mid X\ne Y)
=\frac{2p_{ir}p_{jr}}{1-\sum_kp_{kr}^2},
\qquad i<j.
$$

The probabilities sum to one over all unordered heterotypic pairs:

全部无方向异型 pair 的概率和为 1：

$$
\sum_{i<j}q_{ijr}
=
\frac{\sum_{i<j}2p_{ir}p_{jr}}
{1-\sum_kp_{kr}^2}
=1.
$$

The expected count in ROI \(r\) is

ROI \(r\) 中的期望数为

$$
E_{ijr}=D_rq_{ijr}.
$$

This null model asks whether the doublet-pair composition differs from random
heterotypic pairing under the local singlet composition. It does not assume
equal cell-type abundance.

该零模型检验 doublet pair 组成是否偏离由局部 singlet 组成决定的随机异型配对；
它不假设各细胞类型等丰度。

---

## 3. Pair-specific ROI eligibility

## 3. Pair 特异的 ROI 可纳入条件

Rare singlet populations yield unstable \(p_{ir}\) and very small expected
counts. For minimum count \(m=100\), define

稀有 singlet 群体会产生不稳定的 \(p_{ir}\) 和极小期望数。令最小计数
\(m=100\)，定义

$$
I_{ijr}
=\mathbf 1(S_{ir}\ge m,\;S_{jr}\ge m).
$$

ROI \(r\) contributes to pair \(\{i,j\}\) only when \(I_{ijr}=1\).
The eligible set for condition \(c\) is

仅当 \(I_{ijr}=1\) 时，ROI \(r\) 才纳入 pair \(\{i,j\}\) 的估计。疾病组
\(c\) 的合格集合为

$$
\mathcal R_{ijc}
=\{r:\operatorname{condition}(r)=c,\ I_{ijr}=1\}.
$$

Crucially, \(H_r=1-\sum_kp_{kr}^2\) is computed from every singlet type in
the ROI **before** applying \(I_{ijr}\). Recalculating the denominator after
removing rare types would change the null sample space and over-normalize the
remaining pairs.

关键点是：\(H_r=1-\sum_kp_{kr}^2\) 使用 ROI 内的全部 singlet 类型，并在应用
\(I_{ijr}\) **之前**计算。若删除稀有类型后重新计算分母，会改变零模型样本空间，
并对保留 pair 造成再次归一化。

---

## 4. Condition-level aggregation

## 4. 疾病组层面的汇总

Observed and expected counts are summed before their ratio is calculated:

应先汇总观察数和期望数，再计算比值：

$$
O_{ijc}
=\sum_{r\in\mathcal R_{ijc}}O_{ijr},
\qquad
E_{ijc}
=\sum_{r\in\mathcal R_{ijc}}E_{ijr}.
$$

This is a ratio of pooled counts, not the arithmetic mean of ROI-level
\(O/E\) values. It weights each ROI through its expected information content
and avoids undefined ROI-level ratios when \(E\) is tiny.

这是汇总计数之比，而不是 ROI 层面 \(O/E\) 的算术平均。该方式通过期望信息量
自然赋予 ROI 权重，并避免在单个 ROI 的 \(E\) 极小时产生不稳定比值。

Let

定义

$$
N_{ijc}
=\sum_{r\in\mathcal R_{ijc}}D_r.
$$

The adaptive expected-count threshold is

动态期望数门槛为

$$
E_{\min,ijc}
=\max(E_0,\alpha N_{ijc})
=\max(5,10^{-4}N_{ijc}).
$$

The estimate is considered reportable when

当满足下式时，估计量可报告：

$$
|\mathcal R_{ijc}|\ge K
\quad\text{and}\quad
E_{ijc}\ge E_{\min,ijc},
$$

where the package default is \(K=1\). Increasing \(K\) is appropriate when
the study design requires evidence from multiple biological replicates.

package 默认 \(K=1\)。当研究设计要求多个生物学重复共同提供证据时，可提高
\(K\)。

---

## 5. Primary enrichment effect size

## 5. 主富集效应量

For an estimable pair, define

对于可估计 pair，定义

$$
R_{O/E,ijc}
=\frac{O_{ijc}}{E_{ijc}},
$$

and its log-scale representation

及其对数尺度表示

$$
L_{ijc}
=\log_2R_{O/E,ijc}
=\log_2\left(\frac{O_{ijc}}{E_{ijc}}\right).
$$

Interpretation:

解释：

- \(R_{O/E}>1\), or \(L>0\): enriched relative to the null.
- \(R_{O/E}=1\), or \(L=0\): consistent with the null expectation.
- \(0<R_{O/E}<1\), or \(L<0\): depleted relative to the null.

- \(R_{O/E}>1\) 或 \(L>0\)：相对零模型富集。
- \(R_{O/E}=1\) 或 \(L=0\)：与零模型期望一致。
- \(0<R_{O/E}<1\) 或 \(L<0\)：相对零模型缺失。

No pseudocount is added to positive \(O\), and \(E\) is never modified.
This preserves the exact effect size and prevents a fixed constant from
dominating pairs with small \(E\).

正观察数 \(O\) 不添加 pseudocount，且 \(E\) 永不修改。这保留了精确效应量，
并避免固定常数主导小 \(E\) pair 的排序。

---

## 6. Observed zeros and heatmap display

## 6. 零观察与热图显示

If an estimable pair has \(O_{ijc}=0\), then

若某个可估计 pair 的 \(O_{ijc}=0\)，则

$$
R_{O/E,ijc}=0,
\qquad
L_{ijc}=-\infty.
$$

An infinite value cannot be mapped to a finite color scale. For visualization
only, stGrads defines

无穷值无法映射到有限色阶。仅为可视化，stGrads 定义

$$
L^{\mathrm{display}}_{ijc}
=\log_2\left(\frac{\delta}{E_{ijc}}\right),
\qquad \delta=0.5.
$$

The output simultaneously records `zero_observed = TRUE`, and the triangular
heatmap labels the tile “Zero observed”. This display value must not be
reported as the exact biological effect size.

输出同时记录 `zero_observed = TRUE`，三角热图在格中标记“Zero observed”。
该显示值不应作为精确生物学效应量报告。

Non-estimable pairs have `log2_oe = NA` and are not passed to
`geom_tile()`, so their locations remain blank rather than appearing as a
third biological state.

不可估计 pair 的 `log2_oe = NA`，且不会传入 `geom_tile()`；因此这些位置保持
空白，不会被误解为第三种生物学状态。

---

## 7. Sensitivity estimate

## 7. 敏感性估计

For comparison with earlier analyses, the package retains

为与早期分析比较，package 保留

$$
R^{\mathrm{sensitivity}}_{ijc}
=\frac{O_{ijc}+0.5}{E_{ijc}+0.5}.
$$

This quantity is stored in `oe_sensitivity` and
`log2_oe_sensitivity`. It is not used by the default heatmap and should be
described explicitly as a sensitivity analysis if reported.

该数值存储于 `oe_sensitivity` 和 `log2_oe_sensitivity`，默认热图不使用它；
若报告，应明确说明其为敏感性分析。

---

## 8. Reproducible implementation

## 8. 可复现实现

```r
enrichment <- stg_pair_enrichment(
  metadata,
  sample_col = "MatchKey",
  condition_col = "Condition",
  min_singlet_cells = 100,
  min_expected = 5,
  min_expected_fraction = 1e-4,
  min_eligible_roi = 1,
  zero_display_pseudocount = 0.5
)

condition_oe <- enrichment$condition
filter_parameters <- enrichment$parameters

stg_plot_pair_triangle(
  condition_oe,
  value_col = "log2_oe",
  show_counts = TRUE,
  zero_label = "Zero\nobserved"
)
```

Important output fields:

重要输出字段：

| Field | Meaning | 中文含义 |
|---|---|---|
| `observed` | aggregated \(O_{ijc}\) | 汇总观察数 |
| `expected` | aggregated \(E_{ijc}\) | 汇总期望数 |
| `n_total_doublet_eligible` | \(N_{ijc}\) | 合格 ROI 的 doublet 总数 |
| `min_expected_required` | adaptive \(E_{\min,ijc}\) | 动态期望数门槛 |
| `estimable` | whether all reporting criteria are met | 是否满足可报告条件 |
| `exclusion_reason` | reason for exclusion | 排除原因 |
| `oe_raw` | exact \(O/E\), including zero | 精确 \(O/E\)，可为零 |
| `log2_oe_raw` | exact log ratio, possibly \(-\infty\) | 精确对数比，可为负无穷 |
| `log2_oe` | estimable heatmap value | 可估计结果的热图数值 |
| `zero_observed` | explicit observed-zero flag | 零观察标记 |
| `log2_oe_sensitivity` | former symmetric-pseudocount result | 旧对称 pseudocount 敏感性结果 |

---

## 9. Suggested Methods wording

## 9. 推荐的 Methods 表述

> Within each ROI, expected frequencies of undirected heterotypic doublet
> pairs were derived from the local singlet composition. For cell types
> \(i\) and \(j\), the conditional null probability was
> \(2p_ip_j/(1-\sum_kp_k^2)\), and the expected count equaled this probability
> multiplied by the total number of certain heterotypic doublets. A
> condition-level pair estimate included only ROIs containing at least 100
> singlets of each constituent type. Observed and expected counts were summed
> across eligible ROIs before calculating the unregularized \(O/E\) ratio.
> Results were reported only when the aggregated expected count was at least
> \(\max(5,10^{-4}N_{\mathrm{total}})\). Observed zeros were represented by
> \(\log_2(0.5/E)\) for heatmap coloring only and were explicitly labelled;
> non-estimable combinations were left blank.

> 在每个 ROI 中，无方向异型 doublet pair 的期望频率由局部 singlet 组成计算。
> 对细胞类型 \(i\) 和 \(j\)，条件化零模型概率为
> \(2p_ip_j/(1-\sum_kp_k^2)\)，期望数为该概率乘以确定异型 doublet 总数。
> 疾病组层面的 pair 估计仅纳入两类组成 singlet 均至少为 100 的 ROI。先在合格
> ROI 间分别汇总观察数和期望数，再计算未经正则化的 \(O/E\)。仅当汇总期望数达到
> \(\max(5,10^{-4}N_{\mathrm{total}})\) 时报告结果。零观察仅在热图着色时使用
> \(\log_2(0.5/E)\) 并明确标记；不可估计组合保持空白。
