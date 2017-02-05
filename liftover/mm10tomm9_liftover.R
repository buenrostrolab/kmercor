library(rtracklayer)
ch <- import.chain("mm10ToMm9.over.chain")
original <- import("ImmGenATAC1219.peak.bed")
b <- liftOver(x=original, chain=ch)
b <- unlist(b)
write.table(data.frame(b)[,1:3], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE,
            file = "ImmGen_mm9.bed")

original <- import("/Volumes/dat/Research/BuenrostroResearch/00_test_run_essentials/01_additionalData/BlacklistFiles/mm9_blacklist.JDB.bed")
blacklist <- liftOver(x=original, chain=ch)
blacklist <- unlist(blacklist)
write.table(data.frame(blacklist)[,1:3], sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE,
            file = "/Volumes/dat/Research/BuenrostroResearch/00_test_run_essentials/01_additionalData/BlacklistFiles/mm10_blacklist.JDB.bed")



