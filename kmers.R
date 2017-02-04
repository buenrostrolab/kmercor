
library(GenomicRanges)
library(chromVAR)
library(BSgenome.Mmusculus.UCSC.mm10)
source("get_kmer_indices.R")

if (basename(getwd()) != "code") setwd("code")

bed <- "../input/ImmGenATAC1219.peak.bed"   # Point to new bed file with peaks
peaks <- get_peaks(bed, sort_peaks = FALSE)

kmer_idx <- get_kmer_indices(peaks, genome = BSgenome.Mmusculus.UCSC.mm10, k = 5)

# Convert to sparse Matrix
sm <- sparseMatrix(i = length(peaks), j = length(kmer_idx), x = 1)
lapply(1:length(kmer_idx), function(i){
  sm[kmer_idx[[i]],i] <- 1
})
saveRDS(sm, "../output/kmerhits.idx.rds")


