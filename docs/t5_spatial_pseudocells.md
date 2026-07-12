# Spatial pseudo-cells and segment-level back-mapping

This tutorial adds a spatially constrained pseudo-cell layer for sparse
segmented Visium HD data. It is intended for fibroblast/CAF subtyping,
module scoring, differential expression, and ligand-receptor or NicheNet
workflows.

The key principle is that aggregation improves transcriptomic signal without
discarding the original spatial units. Every pseudo-cell retains a complete
segment-level map, so downstream labels can be projected back to the original
segmentation.

## 1. Define spatial pseudo-cells

For an original segment $s$, let $mathbf{x}_s=(x_s,y_s)$ denote its spatial
centroid. With grid width $h$, the spatial grid index is

```{math}
g(s) = (floor(x_s / h), floor(y_s / h)).
```

Here, $h$ is the `grid_size` parameter. Segments are first restricted to the
same sample and, optionally, the same cell type. Within each sample-grid
cell, segments are randomly partitioned into groups with target size $m$:

```{math}
P_{g,k} = { s in g : k m < rank_g(s) <= (k + 1) m }.
```

Only groups satisfying the following size constraint are retained:

```{math}
n_min <= |P_{g,k}| <= n_max.
```

Counts are then aggregated gene by gene:

```{math}
Y_{gk,j} = sum_{s in P_{g,k}} X_{s,j}.
```

This produces an expression matrix with pseudo-cells as columns and genes as
rows. The aggregation is spatially constrained, but it does not claim that
the original segments form a single biological cell.

## 2. Run the package function

The function accepts a Seurat object containing segmentation-level counts and
coordinates. Here, `MatchKey` prevents pseudo-cells from crossing samples and
`first_type` restricts the input to Fibroblast singlets.

```r
library(stGrads)

fibroblast_pseudocells <- stg_spatial_pseudocells(
  object = hd_object,
  group_cols = "MatchKey",
  cell_type_col = "first_type",
  cell_types = "Fibroblast",
  assay = "Spatial.Polygons",
  layer = "counts",
  grid_size = 200,
  target_cells = 80,
  min_cells = 8,
  max_cells = 220,
  seed = 20260624,
  as_seurat = TRUE
)

fibroblast_pc <- fibroblast_pseudocells$object
dim(fibroblast_pseudocells$counts)
head(fibroblast_pseudocells$metadata)
```

The returned object is ready for a standard Seurat workflow:

```r
fibroblast_pc <- NormalizeData(fibroblast_pc, verbose = FALSE)
fibroblast_pc <- FindVariableFeatures(fibroblast_pc, nfeatures = 3000, verbose = FALSE)
fibroblast_pc <- ScaleData(fibroblast_pc, verbose = FALSE)
fibroblast_pc <- RunPCA(fibroblast_pc, npcs = 30, verbose = FALSE)
fibroblast_pc <- FindNeighbors(fibroblast_pc, dims = 1:20, verbose = FALSE)
fibroblast_pc <- FindClusters(fibroblast_pc, resolution = 0.2, verbose = FALSE)
fibroblast_pc <- RunUMAP(fibroblast_pc, dims = 1:20, verbose = FALSE)
DimPlot(fibroblast_pc, group.by = "seurat_clusters")
```

## 3. Preserve the spatial mapping

The central output is `cell_map`. It contains one row per retained original
segment:

```r
mapping <- fibroblast_pseudocells$cell_map

mapping[, c(
  "CellID", "pseudocell_id", "grid_id", "grid_x", "grid_y",
  "x", "y", "MatchKey", "first_type"
)] |> head()
```

The mapping is an explicit many-to-one relation:

```{math}
f: S_retained -> P,  f(s) = p(s).
```

Every retained segment $s$ is assigned to exactly one pseudo-cell $p(s)$.
A pseudo-cell-level label can therefore be projected back to the original
segments by a join:

```r
pc_labels <- data.frame(
  pseudocell_id = colnames(fibroblast_pc),
  CAF_state = fibroblast_pc$CAF_state,
  stringsAsFactors = FALSE
)

segment_labels <- merge(
  mapping,
  pc_labels,
  by = "pseudocell_id",
  all.x = TRUE,
  sort = FALSE
)
```

This table can be joined to the original Seurat metadata or used directly for
spatial plotting. The pseudo-cell centroid is retained in
`fibroblast_pseudocells$metadata$x_centroid` and
`fibroblast_pseudocells$metadata$y_centroid`; the original segment centroids
remain in `cell_map$x` and `cell_map$y`.

## 4. Plot pseudo-cell locations and projected labels

```r
library(ggplot2)

ggplot(segment_labels, aes(x = x, y = y, color = CAF_state)) +
  geom_point(size = 0.25, alpha = 0.8) +
  scale_y_reverse() +
  coord_equal() +
  theme_minimal() +
  labs(
    title = "CAF labels projected from spatial pseudo-cells",
    subtitle = "Each original segment inherits the label of its pseudo-cell",
    x = "X coordinate",
    y = "Y coordinate"
  )
```

The projected map is a visualization of the aggregation result, not a new
segmentation. For morphology-aware review, the original polygon boundaries
should remain the primary geometric representation.

## 5. CAF subtyping and NicheNet input

After pseudo-cell construction, canonical CAF programs can be scored on the
pseudo-cell object. For example:

```r
caf_programs <- list(
  iCAF = c("IL6", "CXCL12", "CXCL14", "LIF", "CFD", "PDGFRA"),
  myCAF = c("ACTA2", "TAGLN", "MYL9", "POSTN", "COL1A1", "COL3A1"),
  apCAF = c("HLA-DRA", "HLA-DRB1", "CD74", "CIITA")
)

program_scores <- stg_marker_programs(
  fibroblast_pseudocells$counts,
  gene_sets = caf_programs
)
fibroblast_pc$CAF_state <- colnames(program_scores)[
  max.col(as.matrix(program_scores), ties.method = "first")
]
```

The pseudo-cell object can then be used as the sender population for a
NicheNet workflow, while `cell_map` provides the link back to the original
spatial segments. NicheNet resource files and model choices remain external
to stGrads; the package provides the spatial aggregation and the
provenance-preserving mapping layer.
