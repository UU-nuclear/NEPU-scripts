# Generator function to create an object containing
# two functions fun and jac.
# The function 'fun' calculates up to 
# a constant the logarithmized posterior density
# for supplied parameter sets. The parameter sets
# are given as a matrix, each column corresponds to
# a specific parameter set.
# The function 'jac' calculates the gradient of 
# the logarithmic posterior density for a specific
# parameter set pref.
# Meaning of parameters:
#   p0 ... best prior parameter vector
#   P0 ... prior covariance matrix of parameters
#   D  ... Diagonal matrix with statistical errors of experiments
#   S  ... matrix S mapping systematic errors to the appropriate indices of experimental data points
#   X  ... covariance matrix associated with systematic error components
#   talys ... talys wrapper as created by function 'createTalysFun' in 'talys_wrapper.R'

setupLogPost <- function(p0, P0, yexp, D, S, X, talys) {
  
  lastPars <- NULL
  
  getLastPars <- function() {
    lastPars
  }
  
  fun <- function(pref, fref = NULL) {
    print("fun")
    lastPars <<- pref
    if (is.null(fref))
      fref <- talys$fun(pref)
    
    pref <- as.matrix(pref)
    fref <- as.matrix(fref)
    stopifnot(ncol(pref) == ncol(fref))
    invP0 <- solve(P0)
    res <- rep(NA_real_, ncol(pref))
    for (i in seq_len(ncol(pref))) {
      prefCur <- pref[,i]
      frefCur <- fref[,i]
      dpriorRef <- prefCur - p0
      LpriorRef <- as.vector(crossprod(dpriorRef, invP0 %*% dpriorRef))
      Lref <- as.vector(mult_xt_invCov_x(yexp - frefCur, D, S, X)) + LpriorRef
      Lref <- (-0.5) * Lref
      res[i] <- Lref
    }
    res
  }
  
  
  jac <- function(pref) {
    invP0 <- solve(P0)
    J <- talys$jac(pref)
    fref <- talys$fun(pref)
    dpriorRef <- pref - p0
    res <- as.matrix(2 * invP0 %*% dpriorRef - 2 * t(J) %*% mult_invCov_x(yexp - fref, D, S, X))
    res <- (-0.5) * res
    res 
  }
  
  list(fun = fun, jac = jac, getLastPars = getLastPars)
}
