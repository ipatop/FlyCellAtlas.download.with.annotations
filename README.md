
# FlyCellAtlas.download.with.annotations

<!-- badges: start -->
<!-- badges: end -->

The goal of FlyCellAtlas.download.with.annotations is to ...

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
