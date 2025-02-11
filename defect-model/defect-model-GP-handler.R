
create_GP_handler <- function(model) {

  gpDt <- NULL

  # this function should acctually belong in utility functions
  # remove before integrating in NucDataBaynet
  returnSparseMatrix <- function(dt, dims, ret.mat = TRUE) {

    if (!ret.mat) dt[] else
      dt[, sparseMatrix(i = IDX1, j = IDX2, x = X, dims = dims)]
  }

  addGP <- function(expid, sigma, len, nugget) {

    newDt <- data.table(EXPID = expid,
                        GPTYPE = "sparse",
                        PARNAME = c("sigma", "len", "nugget"),
                        PARVAL = c(sigma, len, nugget))
    proposedDt <- rbind(gpDt, newDt)
    stopifnot(anyDuplicated(proposedDt) == 0)
    gpDt <<- proposedDt
  }
	
	
}