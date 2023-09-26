# StemnesScoRe: An R package to estimate the stemness of glioma cancer cells at single-cell resolution

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

devtools::install_github("Necla/StemnessScore")
```
For complete list of functions and instructions:
```
library(help = "StemnessScore") 
```

Note that devtools does not build vignettes by default. To view the vignette:

```
devtools::install_github("Necla/StemnessScore", build_vignettes = FALSE)

library(StemnessScore)

vignette("StemnessScore")
```

# Tutorial
# Downloading GBM scATAC-seq data  
To execute the stemness.score() function within the package, you can download the GBM scATAC-seq data from either of the following sources:

1. Download the data from the National Center for Biotechnology Information (NCBI) Gene Expression Omnibus (GEO) using the accession number GSE139136.

2. Access the data directly from the "DATA" folder included within the package.


# Obtaining/Reading gene list, gene importance, and genomic region data 

Download genes obtained from machine learning algorithms along with their respective importance values from the "DATA" folder within the package. One can use different gene lists and their associated importance scores for downstream analyses and applications. 

Obtain promoter and enhancer regions associated with the genes identified through ML algorithms using the GeneAlaCart database (https://genealacart.genecards.org/). 

#  Set the File Path 
Before running the stemness.score() function, you need to specify the file path to the folder containing your data. Let's call this folder "GBM_Data." Now, you can proceed to run the stemness.score() function on the acquired dataset. 

```
stemness.score("GBM_Data","GSM4131776_4218_matrix.mtx.gz", "GSM4131776_4218_peaks.bed.gz", "GSM4131776_4218_barcodes.tsv.gz", "genelist.csv", "Classification_feature_Importance_RF.tsv", genomic.regions)

```
 
