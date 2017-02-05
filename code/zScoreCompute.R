library(chromVAR)

if (basename(getwd()) != "code") setwd("code")

# For laptop
BiocParallel::register(BiocParallel::MulticoreParam(2, progressbar = FALSE))

#bg <- readRDS("../output/backgroundPeaks.rds")
ixx <- readRDS("../output/kmer.ixx.rds")
zh <- readRDS("../input/LPS_Counts.rds")

boo <- Matrix::rowSums(assays(zh, "counts")[[1]]) >= 1
obj <- add_gc_bias(zh[boo,], genome = BSgenome.Mmusculus.UCSC.mm10)
bg <- get_background_peaks(obj)

deviations <- compute_deviations(obj, background_peaks = bg, annotations = ixx[boo,])
saveRDS(deviations, file = "../output/other_LPS_Deviations.rds")


