metadata <- read.csv("../data/alignments/flu/experiment_1/raw/flu_8445_segments_METADATA.csv",
                     header = TRUE)
lastf <- function(x) x[length(x)]
metadata$Collection.Year <- as.numeric(
  unlist(lapply(strsplit(as.character(metadata$Collection.Date), "/"), lastf))
  )
Freqs.full <- table(metadata$Collection.Year)
Freqs.full
metadata$Sequence.Accession  <- gsub("\\*", "", metadata$Sequence.Accession)
#######
Accessions <- metadata$Sequence.Accession
Accessions <- gsub("\\*", "", Accessions)
oriAln <- ape::read.dna("../data/alignments/flu/experiment_1/raw/flu_8445_segments_seqs.fasta",
                        format = "fasta")
Names <- sapply(1:nrow(metadata),
                function(i) paste("seq", i, Accessions[i], metadata$Collection.Year[i], sep = "_")
)
Positions <- sapply(Accessions, function(x) grep(x, names(oriAln)))
table(unlist(lapply(Positions, length))) ## check one
length(unique(Positions)) == nrow(metadata) ## should be true, i.e., every sequence was found
newAln <- oriAln[Positions]
names(newAln) <- Names
ape::write.dna(newAln, "../data/alignments/flu/experiment_1/flu_8445_segments_seqs_renamed.fasta",
               format = "fasta")