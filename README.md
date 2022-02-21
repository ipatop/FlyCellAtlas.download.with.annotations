
# FlyCellAtlas.download.with.annotations

<!-- badges: start -->
<!-- badges: end -->

The goal of FlyCellAtlas.download.with.annotations is to put the data from the [FlyCellAtlas](https://flycellatlas.org) into a Seurat object. This proved to be complicated using functions from Seurat and impossible only using the loom files.

The functions of the package will requiere both 10X loom AND h5ad files from [FlyCellAtlas](https://flycellatlas.org) 

The function `loom_too_Seurat` will only read the objects from FlyCellAtlas and transform it in a Seurat object. This object only has the most variable gene expression.
The function `loom_too_Seurat_brandNew` will read the objects from FlyCellAtlas and transform it in a NEW Seurat object and perform new dimensional reduction. This is useful to get all the genes in a propper format. It is also useful to just have a fresh start but with the celltype annotation from the FlyCellAtlas paper.
The function `read_loom_and_analyze` will perform multiple analyses. It will load and save the Seurat object from the paper, reanalyze it into a new object and also calculate marker genes for both analyzes. This will be saved in an RData object.
All functions have the option to generate and save general plots.


## Installation

You can install the development version of FlyCellAtlas.download.with.annotations from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("ipatop/FlyCellAtlas.download.with.annotations")
```

## Example

Load library

``` r
library(FlyCellAtlas.download.with.annotations)
``` 

Analyze bodywall example

``` r
#Files are in https://flycellatlas.org
# Here we download the body wall 10X loom and h5ad
download.file(url = "https://cloud.flycellatlas.org/index.php/s/egxzns8NgtCjonB/download",destfile = "./bodywall.loom")
download.file(url = "https://cloud.flycellatlas.org/index.php/s/sZZMSbcNk4SWtHE/download",destfile = "./bodywall.h5ad")

loom_file<-"./bodywall.loom"
h5ad_file<-"./bodywall.h5ad"
## basic example code to create a new Seurat object with annotations from FlyCellAtlas
bodywall_new<-loom_too_Seurat_brandNew(loom_file =loom_file,h5ad_file=h5ad_file)
```


Now `bodywall_new` is a Seurat object we can not plot etc

``` r

DimPlot(bodywall_new,repel = T,label = T,reduction = "umap", group.by = "annotation")+ NoLegend()

```

This is the general structure of the flyCellAtlas data:

``` r
head(bodywall_new@meta.data)
```

You see that you can then color by sex or batch
``` r

DimPlot(bodywall_new,repel = T,label = T,reduction = "umap", group.by = "sex")+ NoLegend()

```
