library(chromVAR)

if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

bg <- readRDS("../output/backgroundPeaks.rds")
ixx <- readRDS("../output/kmer.ixx.rds")
immgen <- readRDS("../output/immgenSamples.rds")

deviations <- compute_deviations(immgen, background_peaks = bg, annotations = ixx)
