# StemnessScore: An R package to estimate the stemness score of each glioma cell using scATAC-seq data

Cancer stem cells are one of these populations and the study of cancer stem cells, which are thought to play an active role in the growth and recurrence of tumor cells, is critical in understanding the initiation, development, and resistance to treatment of cancer. This R package was developed to estimate the stemness score specific to glioblastoma using scATAC-seq data. In this way, researchers will be able to investigate the cancer stemness of GBM cells using only accessible regions in chromatin (scATAC-seq data). 

# Dependencies for the StemnessScore package

1. While using R programming, we suggest you use Rstudio which is the R statistical computing environment to use. 
2. devtools is required to install StemnessScore.
3. StemnessScore uses the following R packages which have to be installed before installing the StemnessScore package:
* Iranges
* GenomicRanges
* readr
* readxl
* dplyr
* Matrix
* tidyverse

# Installation
```
library(devtools)

devtools::install_github("Necla/Stemness-Score")
```
For complete list of functions and instructions:
```
library(help = "StemnessScore") 
```

Note that devtools does not build vignettes by default. To view the vignette:

```
devtools::install_github("Necla/Stemness-Score", build_vignettes = TRUE)

library(StemnessScore)

vignette("StemnessScore")
```

 
