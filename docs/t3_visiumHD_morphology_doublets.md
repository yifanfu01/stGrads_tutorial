# Visium HD Morphology and RCTD Doublet Refinement

# Visium HD 形态学与 RCTD Doublet 增强

This tutorial introduces the pathology-aware workflow added in stGrads 0.2.
It extends the original distance–strength–expression framework from the
stGrads article to segmented Visium HD cells.

本教程介绍 stGrads 0.2 新增的病理感知分析流程，将原 stGrads 文章中的
“距离—强度—表达”框架扩展至分割后的 Visium HD 细胞。

The complete workflow is:

完整流程为：

> Visium HD → cell/nucleus segmentation → pathology morphology → RCTD
> singlets and doublets → local spatial support → marker overlap →
> evidence-based composite-segment refinement.

> Visium HD → 细胞/细胞核分割 → 病理形态 → RCTD singlet 与 doublet →
> 局部空间支持 → marker 重叠 → 基于多证据的复合分割单元增强。

---

## 1. Install and load stGrads 0.2

## 1. 安装并载入 stGrads 0.2

```r
remotes::install_github("yifanfu01/stGrads")

library(stGrads)
library(ggplot2)
library(patchwork)
```

For local development:

本地开发版本：

```r
remotes::install_local("/Users/yifanfu/LabCode_2026/stGrads2")
```

---

## 2. Load the package demo

## 2. 读取 package demo

The package contains an anonymized crop of 694 Visium HD spatial units and
19 marker genes. It does not contain a tissue image, patient identifier, or
full transcriptome.

package 内置一个匿名化裁剪数据，包含 694 个 Visium HD 空间单元和 19 个 marker
基因，不包含组织图像、患者标识或完整转录组。

```r
data("stgrads_hd_demo")

metadata <- stgrads_hd_demo$metadata
counts <- stgrads_hd_demo$counts
log_expression <- stgrads_hd_demo$log_expression

dim(metadata)
table(metadata$spot_class)
```

---

## 3. Upgrade a previously saved RDS object

## 3. 升级既往保存的 RDS 对象

The compatibility layer supports Seurat v4/v5 objects, legacy Visium image
coordinates, VisiumV2 centroid boundaries, metadata coordinates, and both
slot- and layer-based assay access.

兼容层支持 Seurat v4/v5 对象、旧版 Visium image 坐标、VisiumV2 centroid
boundary、metadata 坐标以及 slot/layer 两种 assay 读取方式。

```r
old_object <- readRDS("previous_stGrads_or_Visium_HD_object.rds")

# The operation adds standardized fields but does not remove existing content.
# 该操作只增加标准化字段，不会删除原有内容。
old_object <- stg_upgrade_object(old_object)

head(old_object[[]][, c("stg_x", "stg_y")])
```

If the old object already contains micron-scale morphology coordinates:

如果旧对象已经包含微米尺度的形态学质心坐标：

```r
old_object <- stg_upgrade_object(
  old_object,
  coordinate_scale = "micron"
)
```

---

## 4. Calculate cell and nuclear morphology

## 4. 计算细胞与细胞核形态

For a Space Ranger `outs/` directory:

对于 Space Ranger `outs/` 目录：

```r
morphology <- stg_prepare_hd_morphology(
  outs_dir = "/path/to/spaceranger/outs",
  chip = "PC1",
  run_qc = TRUE
)

head(
  morphology[, c(
    "area_um2_cell",
    "area_um2_nuc",
    "nucleus_fraction",
    "nc_ratio",
    "nucleus_centroid_displacement_norm",
    "morph_qc_pass"
  )]
)
```

The main morphology fields include:

主要形态字段包括：

| Field | English definition | 中文含义 |
|---|---|---|
| `area_um2_cell` | segmented cell area | 细胞分割面积 |
| `area_um2_nuc` | nuclear area | 细胞核面积 |
| `nucleus_fraction` | nuclear area / cell area | 核面积占细胞面积比例 |
| `nc_ratio` | nucleus-to-cytoplasm area ratio | 核质面积比 |
| `circularity_cell` | cell circularity | 细胞圆形度 |
| `solidity_cell` | cell solidity | 细胞实度 |
| `boundary_complexity_cell` | perimeter relative to equal-area circle | 相对于等面积圆的边界复杂度 |
| `nucleus_centroid_displacement_norm` | normalized nuclear displacement | 标准化核—细胞质心偏移 |

Attach the table without rebuilding the Seurat object:

无需重建 Seurat 对象即可合并形态表：

```r
old_object <- stg_add_morphology(old_object, morphology)
```

---

## 5. Standardize RCTD doublets

## 5. 标准化 RCTD doublet

```r
metadata <- stg_annotate_rctd(metadata)

table(metadata$doublet_pair, useNA = "ifany")
```

The pair is canonical and undirected:

pair 经过规范化且无方向：

```r
stg_canonical_pair(
  c("Duct2", "Fibroblast"),
  c("Fibroblast", "Duct2")
)

# Both values are "Duct2__Fibroblast".
# 两个结果均为 "Duct2__Fibroblast"。
```

---

## 6. Test local spatial support

## 6. 检验局部空间支持

```r
support <- stg_doublet_spatial_support(
  metadata = metadata,
  coordinates = metadata,
  radius = 30,
  sample_col = "sample",
  cell_id_col = "CellID",
  x_col = "centroid_x_um_cell",
  y_col = "centroid_y_um_cell",
  coordinate_cell_col = "CellID"
)

table(support$contact_evidence)

metadata <- stg_add_spatial_support(
  metadata,
  support,
  cell_id_col = "CellID"
)
```

A supported doublet has nearby singlets matching both `first_type` and
`second_type` within 30 μm.

支持性 doublet 要求在 30 μm 内分别存在匹配 `first_type` 和 `second_type` 的
singlet。

This is local spatial evidence, not direct proof of membrane contact.

这是局部空间证据，并非细胞膜直接接触的证明。

---

## 7. Display Duct2, fibroblasts, and mixed units

## 7. 展示 Duct2、成纤维细胞与混合单元

```r
metadata <- stg_classify_pair_display(
  metadata,
  pair = "Duct2__Fibroblast",
  support_col = "contact_evidence"
)

p_pair <- stg_plot_doublet_map(
  metadata,
  x_col = "stg_x",
  y_col = "stg_y",
  point_size = 0.8
) +
  labs(
    title = "Duct2-Fibroblast mixed units",
    subtitle = "Red: spatially supported RCTD doublets"
  )

p_pair
```

---

## 8. Calculate Duct2 and fibroblast marker programs

## 8. 计算 Duct2 与成纤维细胞 marker program

Single genes are intuitive but sensitive to dropout. Marker programs are
recommended for the main figure.

单基因容易理解但受 dropout 影响较大，主图推荐使用 marker program。

```r
programs <- stg_marker_programs(
  log_expression,
  gene_sets = list(
    Duct2 = stgrads_hd_demo$duct2_markers,
    Fibroblast = stgrads_hd_demo$fibroblast_markers
  )
)

metadata$Duct2_program <- programs[metadata$CellID, "Duct2"]
metadata$Fibroblast_program <- programs[metadata$CellID, "Fibroblast"]

metadata$program_overlap <- stg_program_overlap(
  metadata$Duct2_program,
  metadata$Fibroblast_program,
  method = "minimum"
)
```

The minimum-based overlap is high only when both programs are high.

基于最小值的重叠评分仅在两个 program 均高时升高。

```r
p_duct2 <- stg_plot_spatial_feature(
  metadata,
  feature = "Duct2_program",
  x_col = "stg_x",
  y_col = "stg_y"
) +
  labs(title = "Duct2 marker program")

p_fibroblast <- stg_plot_spatial_feature(
  metadata,
  feature = "Fibroblast_program",
  x_col = "stg_x",
  y_col = "stg_y",
  low = "#F7FCF5",
  mid = "#74C476",
  high = "#00441B"
) +
  labs(title = "Fibroblast marker program")

p_overlap <- stg_plot_spatial_feature(
  metadata,
  feature = "program_overlap",
  x_col = "stg_x",
  y_col = "stg_y",
  low = "#F7F7F7",
  mid = "#FDBB84",
  high = "#7F0000"
) +
  labs(title = "Joint Duct2-Fibroblast signal")

p_duct2 | p_fibroblast | p_overlap
```

---

## 9. Integrate evidence for segmentation refinement

## 9. 整合证据进行分割增强

```r
refined <- stg_refine_doublets(
  metadata,
  support_col = "both_types_supported",
  overlap_col = "program_overlap",
  sample_col = "sample"
)

table(
  refined$stg_refinement_class,
  useNA = "ifany"
)
```

The score integrates:

评分整合以下证据：

- spatial support;
- cell area;
- nuclear area;
- UMI count;
- detected genes;
- dual marker-program overlap.

- 空间支持；
- 细胞面积；
- 细胞核面积；
- UMI 数；
- 检测基因数；
- 双 marker program 重叠。

The function prioritizes candidate composite segments but does not
automatically modify polygon geometry.

该函数用于筛选候选复合分割单元，但不会自动修改 polygon 几何边界。

---

## 10. Calculate conditional observed-to-expected enrichment

## 10. 计算条件化观察/期望富集

```r
enrichment <- stg_pair_enrichment(
  refined,
  sample_col = "sample",
  condition_col = "condition"
)

head(
  enrichment$condition[
    order(enrichment$condition$log2_oe, decreasing = TRUE),
  ]
)
```

For an undirected heterotypic pair:

对于无方向的异型 pair：

$$
P(A,B\mid A\ne B)
=
\frac{2p_Ap_B}{1-\sum_i p_i^2}.
$$

$$
\mathrm{Enrichment}_{AB}
=
\log_2\left(
\frac{O_{AB}+0.5}{E_{AB}+0.5}
\right).
$$

The conditioning avoids a systematic upward enrichment bias when the RCTD
output contains only heterotypic doublets.

当 RCTD 输出仅包含异型 doublet 时，该条件化可避免系统性的富集上偏。

```r
stg_plot_pair_triangle(
  enrichment$condition,
  show_counts = TRUE
)
```

---

## 11. Connect morphology to spatial gradients

## 11. 将形态特征连接到空间梯度

The original stGrads logic can now be applied to morphology as well as
expression.

原 stGrads 的梯度逻辑现在既可以用于表达，也可以用于形态。

```r
gradient <- stg_spatial_gradient(
  object = refined,
  label_col = "first_type",
  reference = "Duct2",
  query = "Fibroblast",
  x_col = "centroid_x_um_cell",
  y_col = "centroid_y_um_cell",
  cell_id_col = "CellID",
  sample_col = "sample",
  max_distance = 100,
  decay = "exponential"
)

morphology_matrix <- t(as.matrix(
  refined[, c(
    "area_um2_cell",
    "area_um2_nuc",
    "nucleus_fraction",
    "boundary_complexity_cell"
  )]
))
colnames(morphology_matrix) <- refined$CellID

morphology_gradients <- stg_gradient_test(
  morphology_matrix,
  gradient,
  metric = "distance",
  method = "spearman"
)

morphology_gradients
```

This extends the original expression-gradient question:

这将原始表达梯度问题扩展为：

> Do fibroblast morphology and transcriptional programs change continuously
> with distance from Duct2 tumor cells?

> 成纤维细胞的形态与转录 program 是否随其到 Duct2 肿瘤细胞的距离连续变化？

---

## Interpretation boundary

## 解释边界

Use the following wording:

推荐使用以下表述：

> Spatially supported mixed cell-type units with concordant morphometric and
> marker-program evidence.

> 具有一致形态学和 marker program 证据的空间支持性混合细胞类型单元。

Avoid:

避免使用：

> Confirmed physical cell-cell contacts.

> 已确认的物理细胞接触。

Direct contact requires image-level validation, multi-nucleus detection,
membrane segmentation, or orthogonal staining.

直接接触仍需图像层验证、多细胞核检测、细胞膜分割或正交染色。

