metadata <- read.csv("../data/raw/FLu_8445_segments_METADATA.csv", header = TRUE)
lastf <- function(x) x[length(x)]
metadata$Collection.Year <- as.numeric(unlist(lapply(strsplit(as.character(metadata$Collection.Date), "/"), lastf)))
Freqs.full <- table(metadata$Collection.Year)
Freqs.full
metadata$Sequence.Accession  <- gsub("\\*", "", metadata$Sequence.Accession)
########
# Renaming sequences, remove commenting to run
# Accessions <- metadata$Sequence.Accession
# Accessions <- gsub("\\*", "", Accessions)
# oriAln <- ape::read.dna("../data/raw/flu_8445_segments_seqs_aligned.fasta", format = "fasta")
# Names <- sapply(1:nrow(metadata),
#                 function(i) paste("seq", i, Accessions[i], metadata$Collection.Year[i], sep = "_")
# )
# Positions <- sapply(Accessions, function(x) grep(x, names(oriAln)))
# table(unlist(lapply(Positions, length))) ## check one
# length(unique(Positions)) == nrow(metadata) ## should be true, i.e., every sequence was found
# newAln <- oriAln
# names(newAln) <- Names[Positions]
# ape::write.dna(newAln, "../data/raw/flu_8445_segments_seqs_aligned_renamed.fasta", format = "fasta")
################################
howManysamples <- function(time_span, freqs){
  if(time_span < 2) stop("Time span less than 2 years makes no sense")
  Years <- as.numeric(names(freqs))
  windows <- seq(min(Years), max(Years), by = time_span)
  K <- length(windows) + 1
  Sizes <- vector(K, mode = "list")
  for(k in 2:K) Sizes[[k]] <- freqs[which(Years > windows[k-1] & Years <= windows[k])]
  Sizes <- Sizes[!sapply(Sizes, is.null)]
  return(Sizes)
}
howManysamples(time_span = 13, freqs = Freqs.full)
#################################
Master <- metadata
nSeqs <- 50
Size <- floor(nrow(Master)/nSeqs)
Tentative.samples <- vector(Size, mode = "list")
for(i in 1:Size){
  pos <- sample(1:nrow(Master), nSeqs, replace = FALSE)
  Tentative.samples[[i]] <- Master[pos, ]
  Master <- Master[-pos, ]
}
getSpan <- function(dt) max(dt$Collection.Year) - min(dt$Collection.Year)
Spans <- unlist(lapply(Tentative.samples, getSpan))
table(Spans)
####
Aln <- ape::read.dna("../data/raw/flu_8445_segments_seqs_aligned_renamed.fasta", format = "fasta")
getSeq <- function(genbank_id, aln){
  pos <- grep(genbank_id, names(aln))
  if(length(pos) > 1) stop("Something is wrong, ID maps to multiple sequences")
  return(aln[pos])
}
getAln <- function(dt){
  aln <- unlist(lapply(dt$Sequence.Accession, function(x) getSeq(genbank_id = x, aln = Aln)), recursive = FALSE)
  class(aln) <- "DNAbin"
  return(aln)
}
Sequences_samples <- lapply(Tentative.samples, getAln)
export_seqs <- function(index){
  ape::write.dna(Sequences_samples[[index]],
                 file = paste("../data/processed/flu_50seqs_dataset_", index, ".fasta", sep = "") ,
                 format = "fasta")
}
sapply(seq_along(Sequences_samples), export_seqs)