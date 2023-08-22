# StemnessScore: An R package to estimate the stemness score of each glioma cell using scATAC-seq data

Cancer stem cells are one of these populations and the study of cancer stem cells, which are thought to play an active role in the growth and recurrence of tumor cells, is critical in understanding the initiation, development, and resistance to treatment of cancer. This R package was developed to estimate the stemness score specific to glioblastoma using scATAC-seq data. In this way, researchers will be able to investigate the cancer stemness of GBM cells using only accessible regions in chromatin (scATAC-seq data). 

# Dependencies for the StemnessScore package

1. R version should be version 3.5+
2. While using R programming, we suggest you use Rstudio which is the R statistical computing environment to use. 
3. devtools is required to install StemnessScore.
4. StemnessScore uses these R packages so you have to install all of them. You may visit the following websites to install them easily:
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

 
