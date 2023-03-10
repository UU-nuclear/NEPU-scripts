#  nucdataBaynet - Nuclear Data Evaluation Using a Bayesian Network 
#  Copyright (C) 2019  Georg Schnabel
#  
#  nucdataBaynet is free software: you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation, either version 3 of the License, or
#  (at your option) any later version.
#  
#  nucdataBaynet is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#  
#  You should have received a copy of the GNU General Public License
#  along with this program.  If not, see <https://www.gnu.org/licenses/>

#' Create Functions for ML Optimization
#'
#' @return List with functions \code{logLike} and \code{gradLogLike}
#' @export
#'
createMLOptimFuns <- function() {

  print("MLOwithPrior")

  # private variables
  this <- list(

    sysCompHandler = NULL,

    d = NULL,
    D = NULL,
    S = NULL,
    P = NULL,

    expDt = NULL,
    sysDt = NULL,
    gpDt = NULL,

    statUncIdx = NULL,
    sysUncIdx = NULL,
    gpIdx = NULL,

    idxSetStat = NULL,
    idxSetSys = NULL,
    idxSetGp = NULL,

    priorExpectation = NULL,
    priorCovariance = NULL,
    priorLogDetCov = NULL,
    priorCov_inv = NULL
  )

  # setter/getter functions

  setDts <- function(expDt = NULL, sysDt = NULL, gpDt = NULL, sysCompHandler) {
    stopifnot(all(c("IDX","ADJUSTABLE","DATA","REFDATA") %in% names(expDt)),
              all(c("IDX","ADJUSTABLE","DATA","REFDATA") %in% names(sysDt)),
              is.null(gpDt) || all(c("IDX","ADJUSTABLE") %in% names(gpDt)))
    this$expDt <<- copy(expDt)
    this$sysDt <<- copy(sysDt)
    this$gpDt <<- copy(gpDt)
    setkey(this$expDt, IDX)
    setkey(this$sysDt, IDX)
    if (!is.null(this$gpDt))
      setkey(this$gpDt, IDX)
    this$statUncIdx <<- which(this$expDt$ADJUSTABLE)
    this$sysUncIdx <<- which(this$sysDt$ADJUSTABLE)
    if (!is.null(this$gpDt))
      this$gpIdx <<- which(this$gpDt$ADJUSTABLE)

    this$idxSetStat <<- seq_along(this$statUncIdx)
    this$idxSetSys <<- length(this$idxSetStat) + seq_along(this$sysUncIdx)
    this$idxSetGp <<- length(this$idxSetStat) + length(this$idxSetSys) + seq_along(this$gpIdx)

    this$S <<- sysCompHandler$map(this$expDt, this$sysDt, ret.mat = TRUE)
    this$P <<- sysCompHandler$cov(this$sysDt, this$gpDt, ret.mat = TRUE)
    this$d <<- getDt_DATA(this$expDt) - getDt_REFDATA(this$expDt) -
      this$S %*% (getDt_DATA(this$sysDt) - getDt_REFDATA(this$sysDt))
    this$D <<- Diagonal(x = getDt_UNC(this$expDt)^2)

    this$sysCompHandler <<- sysCompHandler
  }

  setPrior <- function(expectation, covariance) {
    # exectation = prior mean (or expectation value)
    # of the parameters to optimize
    # covariance = prior covariance matrix of the
    # parameters to optimize

    # check that there is an expectation value for all ADJUSTABLEs
    stopifnot(length(expectation) == length(this$idxSetStat) + length(this$idxSetSys) + length(this$idxSetGp))

    # check that covariance matrix has the right dimensions
    stopifnot(all(dim(covariance) == length(this$idxSetStat) + length(this$idxSetSys) + length(this$idxSetGp)))

    this$priorExpectation <<- expectation
    this$priorCovariance <<- covariance
    this$priorCov_inv <<- solve(covariance)
    detCov <- determinant(covariance)
    stopifnot(detCov$sign == 1)
    this$priorLogDetCov <<- as.vector(detCov$modulus)
  }

  # functions

  updateDts <- function(x) {

    stopifnot(length(x) == length(this$idxSetStat) + length(this$idxSetSys) + length(this$idxSetGp))
    setkey(this$expDt, IDX)
    this$expDt[J(this$statUncIdx), UNC := x[this$idxSetStat]]
    setkey(this$sysDt, IDX)
    this$sysDt[J(this$sysUncIdx), UNC := x[this$idxSetSys]]
    if (!is.null(this$gpDt)) {
      setkey(this$gpDt, IDX)
      this$gpDt[J(this$gpIdx), PARVAL := x[this$idxSetGp]]
    }
    this$P <<- this$sysCompHandler$cov(this$sysDt, this$gpDt, ret.mat = TRUE)
    this$D <<- Diagonal(x = getDt_UNC(this$expDt)^2)
  }

  getUpdateDt <- function(x) {
    # does basically the same as updateDts, but avoids data race for parallel execution.
    # This is slightly slower (some %), since it requires copying of the data
    stopifnot(length(x) == length(this$idxSetStat) + length(this$idxSetSys) + length(this$idxSetGp))
    
    D <- this$D
    if(length(this$idxSetStat)) { # we need to update exp. stat unc
      expDt <- copy(this$expDt)
      setkey(expDt, IDX)
      expDt[J(this$statUncIdx), UNC := x[this$idxSetStat]]
      D <- Diagonal(x = getDt_UNC(expDt)^2)
    }
    
    sysDt <- this$sysDt
    if(length(this$idxSetSys)) { # we need to update sys. uncs.
      sysDt <- copy(this$sysDt)
      setkey(sysDt, IDX)
      sysDt[J(this$sysUncIdx), UNC := x[this$idxSetSys]]
    }

    gpDt <- NULL
    if (!is.null(this$gpDt)) { # gpDt exists
      if(length(this$idxSetGp)) { # we need to update GP hyperpars
        gpDt <- copy(this$gpDt)
        setkey(gpDt, IDX)
        gpDt[J(this$gpIdx), PARVAL := x[this$idxSetGp]]
      } else {
        gpDt <- this$gpDt
      }
    }

    P <- this$sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)

    list(
      d = this$d, # doesn't change
      S = this$S, #doesn't change
      P = P, #changes
      D = D, #changes
      sysDt = sysDt,
      gpDt = gpDt
    )
  }

  getModifiedDts <- function(x) {

    updateDts(x)
    list(statUncIdx = this$statUncIdx,
         sysUncIdx = this$sysUncIdx,
         gpIdx = this$gpIdx,
         expDt = copy(this$expDt),
         sysDt = copy(this$sysDt),
         gpDt = copy(this$gpDt))
  }

  locLogLike <- function(x) {

    #updateDts(x)
    #logLike(this$d, this$D, this$S, this$P)

    update <- getUpdateDt(x)
    logLike(update$d, update$D, update$S, update$P)
  }

  locGradLogLike <- function(x) {

    #updateDts(x)
    #gradStruc <- gradLogLike(this$d, this$D, this$S,
    #                         sysDt = this$sysDt,
    #                         gpDt = this$gpDt,
    #                         statUncIdx = this$statUncIdx,
    #                         sysUncIdx = this$sysUncIdx,
    #                         gpIdx = this$gpIdx,
    #                         sysCompHandler = this$sysCompHandler)
    #with(gradStruc, c(
    #  gradLogLike_wrtStatUnc[this$statUncIdx],
    #  gradLogLike_wrtSysUnc[this$sysUncIdx],
    #  gradLogLike_wrtHyp[this$gpIdx]
    #))

    update <- getUpdateDt(x)
    gradStruc <- gradLogLike(update$d, update$D, update$S,
                             sysDt = update$sysDt,
                             gpDt = update$gpDt,
                             statUncIdx = this$statUncIdx,
                             sysUncIdx = this$sysUncIdx,
                             gpIdx = this$gpIdx,
                             sysCompHandler = this$sysCompHandler)
    with(gradStruc, c(
      gradLogLike_wrtStatUnc[this$statUncIdx],
      gradLogLike_wrtSysUnc[this$sysUncIdx],
      gradLogLike_wrtHyp[this$gpIdx]
    ))
  }

  locLogPost <- function(x) {
    # evaluates the posterior prob. (up to a normalisation)
    L <- locLogLike(x)

    P <- locLogPrior(x)

    P + L
  }

  locGradLogPost <- function(x) {
    dL <- locGradLogLike(x)

    dP <- locGradLogPrior(x)

    as.vector(dL + dP)
  }

  locLogPrior <- function(x) {

    d <- x - this$priorExpectation
    
    -0.5*( as.vector( t(d) %*% this$priorCov_inv %*% d )
          + this$priorLogDetCov
          + length(d)*log(2*pi) 
          )
  }

  locGradLogPrior <- function(x) {

    this$priorCov_inv %*% (this$priorExpectation - x)

  }

  list(setDts = setDts,
       getModifiedDts = getModifiedDts,
       setPrior = setPrior,
       # computation functions
       logLike = locLogLike,
       gradLogLike = locGradLogLike,
       logPost = locLogPost,
       gradLogPost = locGradLogPost,
       logPrior = locLogPrior,
       gradLogPrior = locGradLogPrior
       )
}







