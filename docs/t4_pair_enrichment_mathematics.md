# (Advanced)Mathematical definition


Last updated: 2026-07-11

This note defines the observed-to-expected ratio used by
`stg_pair_enrichment()` and separates the scientific estimator from the
finite value used only for heatmap display.


---

## 1. Statistical unit and notation

Let $r$ index an independent ROI or biological replicate, $c$ index a
condition, and $i,j,k$ index cell types. Within ROI $r$, define:


| Symbol | Definition |
|---|---|
| $S_{ir}$ | number of singlets assigned to type (i) |
| $S_r=\sum_kS_{kr}$ | total singlet count |
| $p_{ir}=S_{ir}/S_r$ | empirical singlet proportion |
| $D_r$ | total number of certain heterotypic mixed units |
| $O_{ijr}$ | observed count of the undirected pair $\{i,j\}$ |

The braces $\{i,j\}$ indicate an unordered pair:
$\{i,j\}=\{j,i\}$. Therefore, each biological pair appears once in the
triangular heatmap.


---

## 2. Conditional heterotypic null model

Consider two independent draws $X,Y$ from the ROI singlet composition.
Before conditioning,

```{math}
P(X=i,Y=j)=p_{ir}p_{jr}.
```

For $i\ne j$, the unordered event $\{i,j\}$ contains two mutually
exclusive ordered outcomes, $(i,j)$ and $(j,i)$. Thus,


```{math}
P(\{i,j\})=2p_{ir}p_{jr}.
```

The probability that two draws have different types is


```{math}
H_r=P(X\ne Y)
=1-P(X=Y)
=1-\sum_kp_{kr}^2.
```

Because the analyzed RCTD class consists of heterotypic doublets, the relevant
sample space is $X\ne Y$. The conditional pair probability is therefore

```{math}
q_{ijr}
=P(\{i,j\}\mid X\ne Y)
=\frac{2p_{ir}p_{jr}}{1-\sum_kp_{kr}^2},
\qquad i<j.
```

The probabilities sum to one over all unordered heterotypic pairs:


```{math}
\sum_{i<j}q_{ijr}
=
\frac{\sum_{i<j}2p_{ir}p_{jr}}
{1-\sum_kp_{kr}^2}
=1.
```

The expected count in ROI $r$ is


```{math}
E_{ijr}=D_rq_{ijr}.
```

This null model asks whether the doublet-pair composition differs from random
heterotypic pairing under the local singlet composition. It does not assume
equal cell-type abundance.


---

## 3. Pair-specific ROI eligibility


Rare singlet populations yield unstable $p_{ir}$ and very small expected
counts. For minimum count $m=100$, define


```{math}
I_{ijr}
=\mathbf 1(S_{ir}\ge m,\;S_{jr}\ge m).
```

ROI $r$ contributes to pair $\{i,j\}$ only when $I_{ijr}=1$.
The eligible set for condition $c$ is


```{math}
\mathcal R_{ijc}
=\{r:\operatorname{condition}(r)=c,\ I_{ijr}=1\}.
```

Crucially, $H_r=1-\sum_kp_{kr}^2$ is computed from every singlet type in
the ROI **before** applying $I_{ijr}$. Recalculating the denominator after
removing rare types would change the null sample space and over-normalize the
remaining pairs.


---

## 4. Condition-level aggregation


Observed and expected counts are summed before their ratio is calculated:


```{math}
O_{ijc}
=\sum_{r\in\mathcal R_{ijc}}O_{ijr},
\qquad
E_{ijc}
=\sum_{r\in\mathcal R_{ijc}}E_{ijr}.
```

This is a ratio of pooled counts, not the arithmetic mean of ROI-level
$O/E$ values. It weights each ROI through its expected information content
and avoids undefined ROI-level ratios when $E$ is tiny.


Let


```{math}
N_{ijc}
=\sum_{r\in\mathcal R_{ijc}}D_r.
```

The adaptive expected-count threshold is


```{math}
E_{\min,ijc}
=\max(E_0,\alpha N_{ijc})
=\max(5,10^{-4}N_{ijc}).
```

The estimate is considered reportable when


```{math}
|\mathcal R_{ijc}|\ge K
\quad\text{and}\quad
E_{ijc}\ge E_{\min,ijc},
```

where the package default is $K=1$. Increasing $K$ is appropriate when
the study design requires evidence from multiple biological replicates.


---

## 5. Primary enrichment effect size


For an estimable pair, define


```{math}
R_{O/E,ijc}
=\frac{O_{ijc}}{E_{ijc}},
```

and its log-scale representation


```{math}
L_{ijc}
=\log_2R_{O/E,ijc}
=\log_2\left(\frac{O_{ijc}}{E_{ijc}}\right).
```

Interpretation:


- $R_{O/E}>1$, or $L>0$: enriched relative to the null.
- $R_{O/E}=1$, or $L=0$: consistent with the null expectation.
- $0<R_{O/E}<1$, or $L<0$: depleted relative to the null.

No pseudocount is added to positive $O$, and $E$ is never modified.
This preserves the exact effect size and prevents a fixed constant from
dominating pairs with small $E$.


---

## 6. Observed zeros and heatmap display


If an estimable pair has $O_{ijc}=0$, then


```{math}
R_{O/E,ijc}=0,
\qquad
L_{ijc}=-\infty.
```

An infinite value cannot be mapped to a finite color scale. For visualization
only, stGrads defines


```{math}
L^{\mathrm{display}}_{ijc}
=\log_2\left(\frac{\delta}{E_{ijc}}\right),
\qquad \delta=0.5.
```

The output simultaneously records `zero_observed = TRUE`, and the triangular heatmap labels the tile “Zero observed”. This display value must not be reported as the exact biological effect size.

Non-estimable pairs have `log2_oe = NA` and are not passed to
`geom_tile()`, so their locations remain blank rather than appearing as a
third biological state.


---

## 7. Sensitivity estimate


For comparison with earlier analyses, the package retains


```{math}
R^{\mathrm{sensitivity}}_{ijc}
=\frac{O_{ijc}+0.5}{E_{ijc}+0.5}.
```

This quantity is stored in `oe_sensitivity` and
`log2_oe_sensitivity`. It is not used by the default heatmap and should be
described explicitly as a sensitivity analysis if reported.


---

## 8. Reproducible implementation


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


| Column | Meaning |
|---|---|
| `observed` | aggregated ($O_{ijc}$) |
| `expected` | aggregated ($E_{ijc}$) |
| `oe_raw` | primary unregularized ($O/E$) effect size |
| `log2_oe` | exact $log_2(O/E)$, including $-\infty$ for observed zero |
| `log2_oe_display` | finite value used only for heatmap coloring |
| `zero_observed` | indicator that the exact observed count is zero |
| `estimable` | whether abundance, replicate, and expected-count rules pass |
| `exclusion_reason` | reason a pair is not estimable |
