# Function to merge lists, either by name or order -----------------------------
source("utils.a.R")

get_kmer_indices <- function(peaks, genome = BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19, k = 6){
  
  if (k > 8){
    stop("k must be less than 8")
  }
  if (k < 5){
    stop("k must be greater than or equal to 5")
  }
  seqs <- Biostrings::getSeq(genome, peaks)
  kmers <- Biostrings::DNAStringSet(Biostrings::mkAllStrings(c("A","C","G","T"),width = k))
  pd <- Biostrings::PDict(kmers)
  indices  <- Biostrings::vwhichPDict(pd,seqs)
  
  tmp <- data.frame(peak_ix = unlist(lapply(1:length(indices), function(x) rep(x, length(indices[[x]])))), kmer_ix = factor(unlist(indices), levels = 1:length(kmers), ordered=T))
  out <- split(tmp$peak_ix,tmp$kmer_ix)
  names(out) <- as.character(kmers)
  
  #remove reverse complements
  tmp <- cbind(names(out), as.character(Biostrings::reverseComplement(kmers)))
  names(out)[tmp[,1] > tmp[,2]] <- tmp[tmp[,1] > tmp[,2],2]
  out <- merge_lists(out, by = "name")
  
  return(out)
}