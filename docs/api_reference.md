# API Reference

This page lists the public stGrads 2.0 functions. Arguments shown here are
the most important user-facing arguments; run `?function_name` in R for the
complete help page and defaults.

## Spatial pseudo-cells and transcriptomic aggregation

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_spatial_pseudocells()` | `object`, `group_cols`, `cell_type_col`, `cell_types`, `assay`, `layer`, `grid_size`, `target_cells`, `min_cells`, `max_cells`, `seed`, `as_seurat` | List with `counts`, `metadata`, `cell_map`, `n_cells`, and optionally `object` | Aggregates spatially adjacent segmented units into pseudo-cells while preserving segment-level provenance. |
| `stg_pseudobulk()` | `object`, `group_cols`, `assay`, `layer`, `min_cells` | List with sparse counts and grouped metadata | Aggregates counts by metadata groups without using spatial coordinates. |
| `stg_normalize_counts()` | `counts`, `scale_factor`, `log1p` | Normalized matrix | Performs sparse-friendly library-size normalization. |
| `stg_marker_programs()` | `expression`, `gene_sets`, `min_genes` | Cell-by-program score matrix | Scores deterministic marker programs from one or more gene sets. |
| `stg_program_overlap()` | `program_a`, `program_b`, `method` | Numeric overlap score | Combines two marker-program scores, including the minimum fuzzy intersection. |

### `stg_spatial_pseudocells()` mapping output

The `cell_map` table has one row per retained original segment and includes:

| Column | Meaning |
|---|---|
| `CellID` | Original segmented-unit identifier. |
| `pseudocell_id` | Pseudo-cell receiving the segment. |
| `grid_id` | Sample-specific spatial grid identifier. |
| `grid_x`, `grid_y` | Integer spatial-grid coordinates. |
| `x`, `y` | Original segment centroid. |
| `group_cols` | Original sample or replicate identifiers. |
| `cell_type_col` | Original cell-type annotation, when supplied. |

The pseudo-cell-level `metadata` table contains `n_segments`,
`x_centroid`, `y_centroid`, `grid_x`, and `grid_y`. This makes it possible to
assign a CAF or NicheNet-related label to a pseudo-cell and project that label
back to every original segment.

## Spatial coordinates and gradients

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_get_coordinates()` | `object`, `image`, `x_col`, `y_col`, `scale` | `CellID`, `x`, `y` data frame | Extracts coordinates from metadata, Seurat images, or segmentation boundaries. |
| `stg_coordinates()` | `coordinates`, `object`, `x_col`, `y_col` | Coordinate matrix | Standardizes coordinate input for spatial calculations. |
| `stg_spatial_gradient()` | `object`, `query`, `reference`, `radius`, `decay` | Cell-level gradient table | Computes spatial distance or decayed reference strength. |
| `stg_gradient_test()` | `gradient`, `expression`, `group`, `method` | Feature-level test table | Tests transcriptomic association with a spatial gradient. |
| `stg_decay()` | `distance`, `scale`, `model` | Numeric weights | Applies an exponential or user-defined spatial decay. |
| `stg_spatial_neighborhood()` | `object`, `cell_type_col`, `radius`, `group_cols` | Neighborhood counts and proportions | Quantifies local cell-type composition. |
| `stg_tumor_interface()` | `object`, `tumor_col`, `radius` | Interface distance table | Quantifies distance to a tumor or epithelial interface. |
| `stg_stromal_barrier()` | `object`, `stromal_col`, `radius` | Barrier proxy table | Summarizes stromal obstruction between spatial populations. |

## Morphology and segmentation

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_morphology()` | `polygons`, `centroids`, `object` | Morphology table | Computes cell-area and shape descriptors from segmentation polygons. |
| `stg_add_morphology()` | `object`, `polygons`, `centroids` | Updated object or metadata | Adds morphology features to a spatial object. |
| `stg_combine_cell_nucleus()` | `cell_features`, `nucleus_features` | Combined morphology table | Combines cell and nuclear measurements. |
| `stg_morphology_qc()` | `morphology`, `thresholds` | QC table | Flags implausible or low-quality morphology. |
| `stg_morphology_atypia()` | `morphology`, `group_cols` | Atypia score table | Summarizes morphology deviations within a reference group. |
| `stg_prepare_hd_morphology()` | `object`, `polygons`, `nuclei` | Prepared morphology object | Standardizes polygon and nucleus inputs for Visium HD analysis. |

## RCTD mixed units and refinement

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_annotate_rctd()` | `metadata`, `spot_class`, `first_type`, `second_type` | Normalized annotations | Standardizes RCTD singlet and mixed-unit labels. |
| `stg_canonical_pair()` | `type_a`, `type_b` | Canonical pair labels | Treats mixed pairs as undirected. |
| `stg_doublet_spatial_support()` | `metadata`, `coordinates`, `radius` | Support table | Tests whether both constituent singlet types occur nearby. |
| `stg_add_spatial_support()` | `object`, `metadata`, `coordinates`, `radius` | Updated metadata/object | Adds local support evidence to an object. |
| `stg_doublet_expression_residual()` | `counts`, `metadata`, `gene_sets` | Residual score table | Measures expression beyond expected constituent signals. |
| `stg_refine_doublets()` | `metadata`, `evidence`, `weights` | Ranked refinement table | Integrates spatial, morphology, and transcriptomic evidence. |
| `stg_classify_pair_display()` | `pair`, `observed`, `expected` | Display class | Labels enriched, depleted, zero-observed, or non-estimable pairs. |

## Pair enrichment

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_pair_enrichment()` | `metadata`, `sample_col`, `condition_col`, `min_singlet_cells`, `min_expected` | ROI-level and condition-level tables | Computes conditional heterotypic observed-to-expected pair enrichment. |
| `stg_plot_pair_triangle()` | `enrichment`, `condition`, `value`, `show_counts` | ggplot object | Displays undirected pair enrichment as a triangular heatmap. |

The primary effect size is the unregularized ratio (O/E). The display-only
value for an observed zero may use a finite pseudocount so that the heatmap
can be colored, but the exact zero remains recorded in the output.

## Publication-oriented plotting

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_theme_publication()` | `base_size`, `base_family` | ggplot theme | Provides a consistent publication-oriented theme. |
| `stg_plot_hd_he()` | `demo`, `image_path`, `alpha` | ggplot object | Displays cropped H&E context. |
| `stg_plot_polygon_class()` | `metadata`, `polygons`, `class_col` | ggplot object | Plots segmentation polygons by class. |
| `stg_plot_polygon_feature()` | `metadata`, `polygons`, `feature` | ggplot object | Maps a continuous feature onto polygons. |
| `stg_plot_doublet_map()` | `metadata`, `coordinates`, `pair` | ggplot object | Displays mixed units and spatial support. |
| `stg_plot_spatial_feature()` | `object`, `feature`, `coordinates` | ggplot object | Maps a feature across spatial coordinates. |

## Compatibility and data helpers

| Function | Main arguments | Returns | Description |
|---|---|---|---|
| `stg_layer_data()` | `object`, `assay`, `layer` | Matrix | Reads Seurat v5 layers or Seurat v4 slots. |
| `stg_upgrade_object()` | `object` | Updated object | Applies compatibility updates to older stGrads/Seurat objects. |
| `stg_format_hd_cell_id()` | `cell_ids`, `sample` | Character vector | Standardizes HD cell identifiers. |
| `stg_demo_image_path()` | `demo` | File path | Returns the packaged demo H&E image path. |
| `stg_read_hd_scalefactors()` | `path` | List | Reads HD scale-factor metadata. |

## Legacy API

The original stGrads 0.1.x functions remain available for old scripts and RDS
workflows. New analyses should use the `stg_*` functions above whenever an
equivalent function exists.

`FindPrimSpot()`, `CalcNearDis()`, `CalcNearDis_HD()`, `PlotNearDis()`,
`PlotStrengthDis()`, `PlotStrengthExpr()`, `PlotDisExpr()`, `PlotDisProp()`,
`CalcDRG()`, `PlotDRGs()`, `SeuLocation_StConvert()`, and
`calculate_distance()` are retained for backward compatibility.
