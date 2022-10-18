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


## Citation

<div class="oxford-citation-text">

Hani Jieun Kim, Yingxin Lin, Thomas A Geddes, Jean Yee Hwa Yang, Pengyi Yang, CiteFuse enables multi-modal analysis of CITE-seq data, _Bioinformatics_, Volume 36, Issue 14, 15 July 2020, Pages 4137–4143, [https://doi.org/10.1093/bioinformatics/btaa282](https://doi.org/10.1093/bioinformatics/btaa282)

</div>

