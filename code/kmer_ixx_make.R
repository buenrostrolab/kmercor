
library(GenomicRanges)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
source("get_kmer_indices.R")

if (basename(getwd()) != "code") setwd("code")

bed <- "../input/ImmGenATAC1219.peak.bed"   # Point to new bed file with peaks
peaks <- get_peaks(bed, sort_peaks = FALSE)

kmer_idx <- get_kmer_indices(peaks, genome = BSgenome.Mmusculus.UCSC.mm10, k = 5)

sm2 <- reshape::melt.list(kmer_idx,  level=1)
fac <- as.factor(sm2$L1)
sm <- sparseMatrix(i = sm2$value, j = as.integer(fac), x = 1)
sm2 <- sparseMatrix(summary(sm)$i, summary(sm)$j, x = 1)

ixx <- SummarizedExperiment(assays = list(match = sm2), rowRanges = peaks,  colData = DataFrame(kmer = names(kmer_idx)))
saveRDS(ixx, "../output/kmer.ixx.rds")

