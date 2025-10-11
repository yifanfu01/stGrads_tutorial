# Installation

`stGrads` is an R package for analyzing and visualizing spatial transcriptomics data (10x Visium / 10x HD / STO).  
This page shows how to install it and verify the setup.

---

## Requirements

- **R**: 4.1 or newer is recommended  
- **OS**: macOS / Linux / Windows  
- **(Optional) R packages for tutorials**: `Seurat (v5+)`, `ggplot2`, `patchwork`, `data.table`

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