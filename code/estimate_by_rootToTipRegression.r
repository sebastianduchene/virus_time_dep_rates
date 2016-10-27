library(ape)
devtools::source_url("https://raw.githubusercontent.com/maxbiostat/CODE/master/R/PHYLO/find_best_rooting.R")
rm(get.ages)
lastf <- function(x) x[length(x)]
get.ages <- function(tree){
  dates <- unlist(lapply(strsplit(tree$tip.label, "_"), lastf))
  return(as.numeric(dates))
}
estimate_rate <- function(aln){
  Stem <- gsub(".fasta", "", aln)
  MLtree <- read.tree(paste(Stem, "_MLTree.nwk", sep = ""))
  res <- find_best_rooting(MLtree)$lm
  return(res)
}
setwd("../data/alignments/flu/experiment_1/processed/aligned/")
alns <- system("ls *.fasta", intern = TRUE)
Alns <- lapply(alns, function(x) ape::read.dna(x, format = "fasta"))
Diversities <- unlist(lapply(Alns, function(x) pegas::nuc.div(x)))
Names <- gsub(".fasta", "", alns)
RDVS <- parallel::mclapply(alns, estimate_rate, mc.cores = 8)
setwd("../../../../../../code/")
getRate <- function(lm_obj){
  m <-  coef(lm_obj)[2]
  ci <- confint(lm_obj)[2, ]
  return(c(ci[1], m, ci[2]))
}
getSpan <- function(rdv){
  r <-range(rdv$model$ages)
  r[2]-r[1]
}
getMinYear <- function(rdv){
  min(rdv$model$ages)
}
spans <- unlist(lapply(RDVS, getSpan))
initYear <- unlist(lapply(RDVS, getMinYear))
Rates <- matrix(NA, ncol = 3, nrow = length(RDVS))
colnames(Rates) <- c("lwr", "mean", "upr")
for (i in 1:nrow(Rates)) Rates[i, ] <- getRate(RDVS[[i]])
##############
alns[which(Rates[, "mean"] < 0)] ## problematic data sets
alns[which.min(Rates[, "mean"])]
alns[which.max(Rates[, "mean"])]
##
plot(Diversities ~ spans)
LM <- lm(Diversities ~ spans)
summary(LM)
abline(LM)
##
boxplot(Diversities ~ ifelse(Rates[, "mean"] < 0, 1, 0)) ## differences in diversity between rates < 0 and rates > 0
##############
library(ggplot2)
forPlot <- data.frame(Rates, time_span = spans, start_year = initYear) 
number_ticks <- function(n) { function(limits) pretty(limits, n) }
p <- ggplot(forPlot, aes(x = time_span, y = mean)) +
geom_smooth(method = 'lm') + 
geom_pointrange(aes(ymin = lwr, ymax = upr, col = start_year), position = position_dodge(0.5)) +
geom_abline(intercept = 0, slope = 0, linetype = "longdash", size = .5, color = "black") + 
scale_y_continuous("Evolutionary Rate (s/s/y) [regression slope]",
                   breaks = number_ticks(10), expand = c(0, 0)) +
scale_x_continuous("Sampling span (years)", breaks = number_ticks(10), expand = c(0, 0)) +
ggtitle("Rate estimates by RDV regression -- Influenza H3N2") +
theme_bw()
p
pdf("../plots/preliminary_influenza.pdf")
p
dev.off()