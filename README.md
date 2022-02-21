
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

Here a few options

``` r
library(FlyCellAtlas.download.with.annotations)

#try files
system.file("loom", "s_fca_biohub_body_wall_10x.loom", package = "FlyCellAtlas.download.with.annotations",mustWork = T)

loom_bodywall<-system.file("loom", "s_fca_biohub_body_wall_10x.loom", package = "FlyCellAtlas.download.with.annotations",mustWork = T)

## basic example code to create a new Seurat object with annotations from FlyCellAtlas
bopdywall_new<-loom_too_Seurat_brandNew(loom_file =loom_bodywall)
```

