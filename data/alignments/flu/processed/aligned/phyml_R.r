######## This script takes a .fasta alignment and estimates a maximum likelihood tree using PhyML
########
args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1 | args[1] == "help" | args[1] == "usage"){
  cat("usage:  <name_of_file> \n"); q(save = "no")
}
File <- args[1]
Stem <- gsub(".fasta", "", File)
fs <- ape::read.dna(paste(File), format = "fasta")
ape::write.dna(fs, file = paste(Stem, ".phy", sep = ""), format = "sequential") ## export .phy file
####
system(paste("phyml --quiet -i ", Stem, ".phy -q -t e -a e -o tlr -b 0", sep = ""))
system(paste("mv ", Stem, ".phy_phyml_tree.txt ",  Stem, "_MLTree.nwk", sep = ""))
system(paste("rm ", Stem, ".phy ", Stem, ".phy_phyml_stats.txt", sep = ""))
