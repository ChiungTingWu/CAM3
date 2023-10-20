# CAM3
A new and fast version of debCAM ( https://github.com/Lululuella/debCAM )

# Introduction

Complex tissues are dynamic eco-systems consisting of molecularly distinct yet interacting cell types. Computational deconvolution aims to dissect bulk tissue data into cell type compositions and cell-specific expressions. With few exceptions, most existing deconvolution tools exploit supervised approaches requiring various types of references that may be unreliable or even unavailable for specific tissue microenvironments. We have developed an efficient and fully unsupervised deconvolution tool – Convex Analysis of Mixtures (CAM3.0) – specifically designed for studying biologically diverse conditions. We improve and extend our previous CAM framework [@CAM2016] to introduce three new and effective algorithms, namely, radius-fixed clustering to identify reliable markers, linear programming to detect an initial scatter simplex, and a smart floating search for the optimum latent variable model.

`CAM3.0` is an R package for fully unsupervised deconvolution of complex tissues. It provides basic functions to perform unsupervised deconvolution on mixture expression profiles by Convex Analysis of Mixtures (CAM) and some auxiliary functions to help understand the cell type-specific results. It also implements functions to perform supervised deconvolution based on prior knowledge of molecular markers, S matrix or A matrix. Combining molecular markers from CAM and from prior knowledge can achieve semi-supervised deconvolution of mixtures.

You can install the latest version of CAM3.0 from GitHub by
```{r, eval = FALSE}
devtools::install_github("ChiungTingWu/CAM3")
```
or from Bioconductor by
```{r, eval = FALSE}
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("CAM3")
```

# Quick Start

Similar to the previous versions of CAM [@debCAM], the function `CAM3()` includes all necessary steps to decompose a matrix of mixture expression profiles. Optional steps include upstream procedures like norm-based filtering, dimension deduction, perspective projection, local outlier removal and aggregation of gene expression vectors by clustering. Each step in `CAM3()` can also be performed separately if you prefer a more flexible workflow. More details will be introduced in the section 'CAM3.0 Workflow'.

Starting your analysis by `CAM3()` where you need to specify the percentage of low/high-expressed molecules to be removed, in order to reduce the effect of noise and also to speed up the computation. `dim.rdc` refers to the reduced data dimension, which should also be specified and be no less than the maximum number of cell type `K`:
```{r, eval = FALSE}
rCAM3 <- CAM3Run(data, K = 3, dim.rdc = 3, thres.low = 0.30, thres.high = 0.95)
```

Or, if the range of possible cell type is given instead of an exact number:
```{r, eval = FALSE}
rCAM3 <- CAM3Run(data, K = 2:5, dim.rdc = 10, thres.low = 0.30, thres.high = 0.95)
```
# User Manual

https://htmlpreview.github.io/?https://github.com/ChiungTingWu/CAM3/blob/main/inst/doc/CAM3.html
