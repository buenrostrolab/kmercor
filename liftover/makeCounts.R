# Load Libraries and CL source code
library(readr)
library(Matrix)
library(GenomicRanges)
library(chromVAR)

if (basename(getwd()) != "code") setwd("code")

# Register Paralleization; change 2 to 1 if multiple cores are not available
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = TRUE))  # Update this with more cores if appropriate

bed <- "ImmGen_mm9.bed"   # Point to new bed file with peaks
peaks <- get_peaks(bed, sort_peaks = FALSE)

counts <- get_counts("single-SCdata.mergeAll.bam", peaks, paired = TRUE, by_rg = TRUE)  # Takes a while to execute
saveRDS(counts, "other_LPS_Counts_mm9.rds")


r <- readRDS("LPS_Counts.rds")

str(r)
counts@colData@rownames

