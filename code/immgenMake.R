library(chromVAR)
library(data.table)

# TSS Peaks
d <- read.table("../../immgen_dat/ImmGenATAC1219.peak.bed", header = FALSE)
peaks <- GenomicRanges::makeGRangesFromDataFrame(setNames(d, c("chr", "start", "stop")))

# Counts 
csv <- "../../immgen_dat/brian/ImmGen_ATACseq_17Jan17.csv.gz"
counts <- data.matrix(data.frame(fread(input = paste0('zcat < ', csv))))

sample_annotation <- DataFrame(sampleName = colnames(counts))
immgen <- SummarizedExperiment(assays = list(counts = Matrix(counts)), rowRanges = peaks, colData = sample_annotation)
saveRDS(immgen, file = "../output/immgenSamples.rds")

gcbias <- add_gc_bias(immgen, genome = BSgenome.Mmusculus.UCSC.mm10)
bg <- get_background_peaks(gcbias)
