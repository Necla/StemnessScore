# TODO: Add comment
#
# Author: Necla Kochan and Gokhan Karakulah
###############################################################################

#' @title Stemness Score of Glioma Cells Using scATAC-seq Data
#'
#' @description The \code{stemness.score} function estimates the stemness index of
#' each glioma cell using scATAC-seq data.
#'
#' @importFrom GenomicRanges GRanges 
#'
#' @importFrom IRanges IRanges mergeByOverlaps
#'
#' @importFrom readr read_tsv
#'
#' @importFrom readxl read_excel
#'
#' @import dplyr
#' 
#' @import tidyverse
#'
#' @importFrom Matrix readMM
#' 
#' @param main_folder_name the name of the folder where three files for a patient stored
#'
#' @param matrix_file_name matrix file name in MEX (Market Exchange Format) format
#'
#' @param peaks_file_name BED file name containing peaks
#'
#' @param barcodes_file_name tsv file name containing barcodes
#'
#' @param gene_list_file_name csv file name containing gene list
#'
#' @param importance a data frame composed of genes and their importance
#'
#' @param genomic.regions a data frame obtained from GeneAlaCart
#'
#' @return a data frame 
#'
#' @author Necla Kochan and Gokhan Karakulah
#' 
#' @keywords stemness glioblastoma
#'
#' @export
#
#
# # read the genes obtained from machine learning algorithms
# gene.list <- c("ATAD2", "BANF1", "BIRC5", "CCND2", "CDCA8", "CDK1", "CKS2", "CLU", "DEK", "DUT",
#        "EEF1A1", "EZH2", "GJA1", "GLUL", "HMGB2", "HTRA1", "JAG1", "KIF15", "MALAT1",
#        "MARCKS", "PRC1", "PTMA", "RPL32", "RPS11", "RPS26", "RRM2", "SET", "SMC2",
#        "SMC4", "SOX11", "TMSB4X", "TOP2A", "TPX2", "TUBB5", "PBK")
#
# # read the file where importance of genes are sorted.
#
# path1 <- system.file("Data", "Classification_feature_Importance_RF.tsv", package = "StemnessScore")
# importance <- read.delim(path1)

# importance <- read.delim("Classification_feature_Importance_RF.tsv")
#
# genomic.regions <- readxl::read_excel("GeneALaCart-59598-220421-134734.xlsx", sheet = "GeneHancers")
#
#
# path2 <- system.file("Data", "GeneALaCart-59598-220421-134734.xlsx", package = "StemnessScore")
# # read the files which include promoter and enhancer regions of the genes identified by machine learning algorithms.
# genomic.regions <- read_excel(path2, sheet = "GeneHancers")
#
# StemnessScore("GSE139136_RAW","GSM4131776_4218_matrix.mtx.gz", "GSM4131776_4218_peaks.bed.gz", "GSM4131776_4218_barcodes.tsv.gz", "genelist.csv", "Classification_feature_Importance_RF.tsv", genomic.regions)

stemness.score  <- function(main_folder_name, matrix_file_name, peaks_file_name, barcodes_file_name, gene_list_file_name, importance_file_name, genomic_regions) {
  
 
  # print(gene_list_file_name)
  gene.list <- read.csv(gene_list_file_name, header = F)
  # print(gene.list)
  gene.list <- as.vector(as.matrix(gene.list))
  # print(gene.list)
  importance <- read.delim(importance_file_name, header = T)


  # convert gene symbols into capital letters
  importance$Symbol <- toupper(importance$Symbol)
  importance <- importance[which(importance$Symbol %in% gene.list),]


  #select promoters with TSS distance between -5000 and 5000
  genomic.regions$TSSdistance <- as.numeric(genomic.regions$TSSdistance)
  genomic.regions <- genomic.regions#[genomic.regions$TSSdistance>=-5000 & genomic.regions$TSSdistance<=5000,]

  #Enhancer <- genomic.regions[genomic.regions$GHtype=="Enhancer",]
  Promoter <- genomic.regions[genomic.regions$GHtype!="Enhancer",]

  # create Granges objects for promoters and enhancers
  NEFTEL_promoter <- GRanges(seqnames = Promoter$Chromosome, IRanges(as.numeric(Promoter$StartHg19), end = as.numeric(Promoter$EndHg19)), Score=Promoter$Score, Symbol = Promoter$Symbol)
  #gr_enhancer <- GRanges(seqnames = Enhancer$Chromosome, IRanges(as.numeric(Enhancer$StartHg19), end = as.numeric(Enhancer$EndHg19)), Score=Enhancer$Score, Symbol = Enhancer$Symbol)


  # peak-bc matrix
  mex_dir_path <- main_folder_name
  print(mex_dir_path)
  mtx_path <- paste(mex_dir_path, matrix_file_name, sep = '/') #
  print(mtx_path)
  feature_path <- paste(mex_dir_path, peaks_file_name, sep = '/')
  #print(feature_path)
  barcode_path <- paste(mex_dir_path, barcodes_file_name, sep = '/')



  features <- readr::read_tsv(feature_path, col_names = F, show_col_types = FALSE) %>% tidyr::unite(feature)
  barcodes <- readr::read_tsv(barcode_path, col_names = F, show_col_types = FALSE) %>% tidyr::unite(barcode)

  mtx <- Matrix::readMM(mtx_path) %>%
    magrittr::set_rownames(features$feature) %>%
    magrittr::set_colnames(barcodes$barcode)


  peaks = read.delim(feature_path, header = FALSE, stringsAsFactors = FALSE)
  peaks_modified <- peaks
  peaks_modified$V1 <- peaks$V1
  peaks_modified$V2 <- peaks_modified$V3 <- floor(peaks$V2+(peaks$V3-peaks$V2)/2)
  # Peaks_all <- cbind(peaks,peaks_modified)

  d_mat_1 <- as.matrix(mtx)
  colnames <- paste0("`",colnames(d_mat_1),"`")
  d_mat <- data.frame(chr = peaks_modified$V1, start=peaks_modified$V2, end=peaks_modified$V3, d_mat_1);
  colnames(d_mat) <- c("chr", "start","end",colnames)

  # peaks from scATACseq data
  gr_peaks_old <- GRanges(seqnames = peaks_modified$V1, IRanges(peaks_modified$V2, end = peaks_modified$V3))
  GR <- gr_peaks_old

  prom_inter <- mergeByOverlaps(gr_peaks_old,NEFTEL_promoter)
  prom_inter_data <- as.data.frame(prom_inter)
  dim(prom_inter_data)
  prom_inter_data <- prom_inter_data[!duplicated(prom_inter_data[,1:3]), ]
  prom_inter_data <- prom_inter_data[order(prom_inter_data$gr_peaks_old.start, decreasing = FALSE),]
  dim(prom_inter_data)

  D_prom <- d_mat[d_mat$chr %in% prom_inter_data$gr_peaks_old.seqnames & d_mat$start %in% prom_inter_data$gr_peaks_old.start & d_mat$end %in% prom_inter_data$gr_peaks_old.end,]
  dim(D_prom)
  D_prom <- D_prom[order(D_prom$start, decreasing = FALSE),]
  D_prom <- D_prom[!duplicated(D_prom[,2:3]), ]
  DATA_prom <- cbind(prom_inter_data, D_prom)
  dim(DATA_prom)

  # gene.list: gene set obtained by machine learning algorithms. TUBB5=TUBB+TUBB4A. Genealacart verisinde bu sekilde verilmis
  indis <- c(which(DATA_prom$Symbol=="TUBB"),which(DATA_prom$Symbol=="TUBB4A"))
  DATA_prom$Symbol[indis] <- "TUBB5"
  DATA_prom$NEFTEL_promoter.Score <- as.numeric(DATA_prom$NEFTEL_promoter.Score)

  sum_counts1 = DATA_prom %>% group_by(Symbol) %>% summarise(across(.cols = starts_with("`"),.fns = sum))
  rownames(sum_counts1) <- sum_counts1$Symbol
  sum_counts <- sum_counts1[,-1]
  rownames(sum_counts)

  mean_scores =  DATA_prom %>% group_by(Symbol) %>% summarise(across(.cols = ends_with(".Score"),.fns = mean))


  raw_counts <- as.data.frame(sum_counts1)
  rownames(raw_counts) <- raw_counts[,1]
  raw_counts <- raw_counts[,-1]
  raw_counts <- log2(raw_counts+1)

  #normalize <- function(x){(x-min(x))/(max(x)-min(x))}

  # add importance as a new column to DATA_prom
  raw_counts$importance <- rownames(raw_counts)
  for (i in 1:length(gene.list)){
    raw_counts$importance <- ifelse(rownames(raw_counts) == gene.list[i],importance[which(importance$Symbol==gene.list[i]),"importance"], raw_counts$importance)
  }


  raw_counts$WeightedScore <- normalize(mean_scores$NEFTEL_promoter.Score) * normalize(as.numeric(raw_counts$importance))
  raw_counts[is.na(raw_counts$WeightedScore),"WeightedScore"] <- 0


  score <- NULL
  for (j in 1:(ncol(raw_counts)-2)){
    score[j] <- sum(as.numeric(raw_counts[,j]) *as.numeric(raw_counts$WeightedScore))
  }

  stemness_score_prom <- normalize(score)

  P1_score <- data.frame(cellID = colnames(raw_counts[,1:(ncol(raw_counts)-2)]), StemnessScore = stemness_score_prom)

  return(P1_score)

}







