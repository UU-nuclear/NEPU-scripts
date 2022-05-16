##############################################################
# TODO: These functions should go somewhere else in the future
##############################################################

getThresEn <- function(reac, expDt, defaultThresEn = NA) {
  
  if (expDt[REAC == reac, all(DATA == 0)]) {
    if (is.na(defaultThresEn))
      stop("automatic energy threshold specification failed")
    else
      defaultThresEn
  }
  else
    expDt[REAC == reac, L1[which(DATA[order(L1)] != 0)[1]]]
}


createSubentStub <- function(exforReacStr, en, xs=NULL) {                                                                                                                                                                                                                                                                                                                             
  require(digest)
  list(
    ID = paste0("MOD_", digest(exforReacStr, algo="crc32", serialize=FALSE)),
    BIB = list(REACTION = exforReacStr),
       DATA = list(DESCR = c("EN", "DATA"),
                   UNIT = c("MEV", "MB"),
                   TABLE = {
                     tmpTable <- data.table(EN = en)
                     if (!is.null(xs)) tmpTable[, DATA:=xs]
                     tmpTable[]
                   }))
}


# sample from a multivariate normal distribution
sample_mvn <- function(num, meanvec, covmat) {
  L <- chol(as.matrix(covmat))
  p <- length(meanvec)
  N <- num
  meanvec <- as.vector(meanvec)
  stopifnot(dim(covmat)[1] == p)
  samples <- meanvec + t(L) %*% matrix(rnorm(p*N), nrow=p, ncol=N)
  return(samples)
}

