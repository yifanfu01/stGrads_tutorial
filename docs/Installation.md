# Installation

`stGrads` 2.0 supports spatial gradients, Visium HD pathology morphology,
RCTD mixed-unit support, and morphology–transcriptome integration.


This page shows how to install it and verify the setup.


---

## Requirements

- **R**: 4.1 or newer is recommended  
- **OS**: macOS / Linux / Windows  
- **Optional tutorial packages**: `Seurat (v5+)`, `ggplot2`, `patchwork`,
  `data.table`, and `jpeg`

> If you install from source on Windows/macOS, you may need basic build tools:
> - **Windows**: Rtools (matching your R version)
> - **macOS**: Xcode Command Line Tools (`xcode-select --install`)
> - **Linux**: common build essentials (e.g., `build-essential`)

---

## Install the development version

```r
# 1) If you don't have remotes yet:
install.packages("remotes")

# 2) Install stGrads from GitHub:
remotes::install_github("yifanfu01/stGrads")

# 3) Load the package:
library(stGrads)

# 4) Quick check:
packageVersion("stGrads")
# Expected: 2.0.0 or newer
```

> Tip: To update later, just run the same `install_github()` command again.

------

## (Optional) Install common dependencies for the tutorials

```R
install.packages(c(
  "Seurat",      # main analysis framework (v5+)
  "ggplot2",     # plotting
  "patchwork",   # plot layout
  "data.table"   # fast data handling
))

```

## Install from a local clone (advanced)

If you have cloned the repository locally:

```R
# Using remotes:
remotes::install_local("yifanfu01/stGrads")

# Or from terminal:
# R CMD build yifanfu01/stGrads
# R CMD INSTALL stGrads_<version>.tar.gz

```
