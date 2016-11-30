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
pdf("../plots/preliminary_influenza_rdv.pdf")
p
dev.off()
############
Fnames <- gsub(".fasta", "", alns)
BeastData <- read.csv("../data/Summaries_meanRate_338_logFiles_flu.csv")
Clocks <- unlist(lapply(strsplit(as.character(BeastData$file), "_"), function(x) x[1]))
BeastData$clock <- Clocks
Fnames2 <- gsub(".log", "", BeastData$file)
Fnames2 <- gsub("ucln_", "", Fnames2)
Fnames2 <- gsub("strict_", "", Fnames2)
BeastData$file <- Fnames2
forPlot2 <- forPlot
forPlot2$file <- Fnames
forPlot2$clock <- rep("rdv", nrow(forPlot2))
BeastData$time_span <- c(forPlot2$time_span, forPlot2$time_span)
BeastData$start_year <- c(forPlot2$start_year, forPlot2$start_year)
forPlot3 <- rbind(forPlot2[, c("lwr", "mean", "upr", "clock", "time_span", "start_year")],
                  BeastData[, c("lwr", "mean", "upr", "clock", "time_span", "start_year")])
forPlot3$width <- forPlot3$upr-forPlot3$lwr
#########
q0 <- ggplot(forPlot3, aes(x = time_span, y = width, colour = clock)) +
  geom_point() + geom_smooth(method = 'lm') +
  scale_y_continuous("BCI width",
                        breaks = number_ticks(10), expand = c(0, 0)) +
  scale_x_continuous("Sampling span (years)", breaks = number_ticks(10), expand = c(0, 0)) +
  ggtitle("Credibility/confidence interval width -- Influenza H3N2") +
  theme_bw() 
pdf("../plots/preliminary_influenza_widths.pdf")
q0
dev.off()
#########
q1 <- ggplot(subset(forPlot3, clock == "strict"), aes(x = time_span, y = mean)) +
  geom_smooth(method = 'lm') + 
  geom_pointrange(aes(ymin = lwr, ymax = upr, col = start_year), position = position_dodge(0.5)) +
  geom_abline(intercept = 0, slope = 0, linetype = "longdash", size = .5, color = "black") + 
  scale_y_continuous("Evolutionary Rate (s/s/y)",
                     breaks = number_ticks(10), expand = c(0, 0)) +
  scale_x_continuous("Sampling span (years)", breaks = number_ticks(10), expand = c(0, 0)) +
  ggtitle("Rate estimates (strict clock) -- Influenza H3N2") +
  theme_bw()
q1
pdf("../plots/preliminary_influenza_strict.pdf")
q1
dev.off()
#########
q2 <- ggplot(subset(forPlot3, clock == "ucln"), aes(x = time_span, y = mean)) +
  geom_smooth(method = 'lm') + 
  geom_pointrange(aes(ymin = lwr, ymax = upr, col = start_year), position = position_dodge(0.5)) +
  geom_abline(intercept = 0, slope = 0, linetype = "longdash", size = .5, color = "black") + 
  scale_y_continuous("Evolutionary Rate (s/s/y)",
                     breaks = number_ticks(10), expand = c(0, 0)) +
  scale_x_continuous("Sampling span (years)", breaks = number_ticks(10), expand = c(0, 0)) +
  ggtitle("Rate estimates (UCLN) -- Influenza H3N2") +
  theme_bw()
q2
pdf("../plots/preliminary_influenza_ucln.pdf")
q2
dev.off()
#########
q3 <- ggplot(forPlot3, aes(x = time_span, y = mean, fill = clock, col = clock)) +
  geom_smooth(method = 'lm') + 
  geom_pointrange(aes(ymin = lwr, ymax = upr, col = clock), position = position_dodge(0.5)) +
  geom_abline(intercept = 0, slope = 0, linetype = "longdash", size = .5, color = "black") + 
  scale_y_continuous("Evolutionary Rate (s/s/y)",
                     breaks = number_ticks(10), expand = c(0, 0)) +
  scale_x_continuous("Sampling span (years)", breaks = number_ticks(10), expand = c(0, 0)) +
  ggtitle("Rate estimates -- Influenza H3N2") +
  theme_bw()
q3
pdf("../plots/preliminary_influenza.pdf")
q3
dev.off()
#################
#################
Cols <- rev(heat.colors(5))[cut(spans, breaks = seq(min(spans), max(spans), length.out = 5), labels = c(1:4))]
pdf("../plots/influenza_strict_vs_ucln_MeanRate.pdf")
plot(subset(forPlot3, clock == "ucln")$mean ~ subset(forPlot3, clock == "strict")$mean,
     col = Cols, pch  = 16,
     xlab = "Mean rate strict (s/s/y)", ylab = "Mean rate UCLN (s/s/y)")
abline(a = 0, b = 1, lwd = 2)
dev.off()
##
pdf("../plots/influenza_rdv_vs_ucln_MeanRate.pdf")
plot(subset(forPlot3, clock == "ucln")$mean ~ subset(forPlot3, clock == "rdv")$mean,
     col = Cols, pch  = 16,
     xlab = "Mean rate RDV (s/s/y)", ylab = "Mean rate UCLN (s/s/y)")
abline(a = 0, b = 1, lwd = 2)
dev.off()
##
pdf("../plots/influenza_rdv_vs_strict_MeanRate.pdf")
plot(subset(forPlot3, clock == "strict")$mean ~ subset(forPlot3, clock == "rdv")$mean,
     col = Cols, pch  = 16,
     xlab = "Mean rate RDV (s/s/y)", ylab = "Mean rate strict (s/s/y)")
abline(a = 0, b = 1, lwd = 2)
dev.off()
##
plot(subset(forPlot3, clock == "ucln")$lwr ~ subset(forPlot3, clock == "strict")$lwr,
     col = Cols, pch  = 16,
     xlab = "95% lwr rate strict (s/s/y)", ylab = "95% lwr rate UCLN (s/s/y)")
abline(a = 0, b = 1, lwd = 2)
##
plot(subset(forPlot3, clock == "ucln")$upr ~ subset(forPlot3, clock == "strict")$upr,
     col = Cols, pch  = 16,
     xlab = "95% upr rate strict (s/s/y)", ylab = "95% upr rate UCLN (s/s/y)")
abline(a = 0, b = 1, lwd = 2)