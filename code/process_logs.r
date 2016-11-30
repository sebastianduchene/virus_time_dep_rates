## This script will focus solely on MeanRate
findOptimalBurnIn <- function(log){
  ESS.burnin <- function(trace, p){
    return(-coda::effectiveSize(trace[round(p*nrow(trace)):nrow(trace), "meanRate"]))
  }
  Opt <- optimise(ESS.burnin, interval = c(0, .5), trace = log)
  return(list(p.burnin = Opt$minimum, min.ESS = -Opt$objective))
}
getSummary <- function(x, alpha = .95){
  return(data.frame(
    lwr = quantile(x, probs = (1 - alpha)/2),
    mean = mean(x),
    upr = quantile(x, probs = (1 + alpha)/2),
    ess = coda::effectiveSize(x)
  ))
}
returnSummaries <- function(logfile){
  Log <- read.table(paste(logfile), header = TRUE)
  opt.b <- findOptimalBurnIn(log = Log)$p.burnin
  MeanRate <- Log[round(opt.b*nrow(Log)):nrow(Log), ]$meanRate
  return(getSummary(MeanRate))
}
Files <- system("ls *.log", intern = TRUE)
Res <- parallel::mclapply(Files, returnSummaries, mc.cores = 14)
DtRes <-  as.data.frame(data.table::rbindlist(Res, idcol = TRUE))
DtRes$file <- Files[DtRes$.id]
write.csv(DtRes,
            file = paste("Summaries_meanRate_", length(Files), "_logFiles.csv", sep = ""),
            row.names = FALSE)