# CiteFuse


CiteFuse is a streamlined package consisting of a suite of tools for pre-processing, modality integration, clustering, differential RNA and ADT expression analysis, ADT evaluation, ligand-receptor interaction analysis, and interactive web-based visualization of CITE-seq data


## Installation

Install CiteFuse required Bioconductor 3.10

```r
library(BiocManager)
install(version = "3.10")
```


```r
BiocManager::install(c("SummarizedExperiment", "S4Vectors", "SingleCellExperiment", "scater", "scran"))
```

```r
library(devtools)
devtools::install_github("SydneyBioX/CiteFuse")
```

## Vignette

You can find the vignette at our website: https://sydneybiox.github.io/CiteFuse/articles/CiteFuse.html.


## CiteFuse overview


<img src="https://raw.githubusercontent.com/SydneyBioX/CiteFuse/master/inst/figures/CiteFuse_overview.png" align="center"/>

