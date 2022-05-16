#  This is a generator function to create wrapper for running TALYS.
#  It provides two important functions:
#
#    fun(x)    runs a calculation with the parameter specifications given
#              in the vector x and returns a vector with the predicted observables.
#    jac(x)    calculates the Jacobian matrix for the parameter specifications
#              given in the vector x and returns a matrix with the results
#
#  Default parameters and which parameters are part of x and what observables should
#  be returned must be configured. This can be done with the set... functions,
#  in particular setPars and setNeedsDt,  
#  Scroll to the bottom of the file to see a complete list of functions provided.
#
#  This function is used in config.R to create the wrapper object.
#  The wrapper object is used in file 07_tune_talyspars.R 

createTalysFun <- function(talysClust, print.info = TRUE) {

  # closure variables

  thisParamDt <- NULL
  thisNeedsDt <- NULL
  thisSexp <- NULL
  thisMask <- NULL
  thisEps <- 1e-3
  
  thisParTrafoFun <- NULL
  thisParTrafoJac <- NULL

  numPars <- NULL
  numNeeds<- NULL

  cacheSize <- 50
  cacheIdx <- 0
  cache <- replicate(cacheSize, list(), simplify=FALSE)
  
  # debugging vars

  refCalcId <- NULL
  refCalcJobs <- NULL

  jacCalcId <- NULL
  jacCalcJobs <- NULL

  # print
  infoprint <- function(str) {
    if (print.info)
      cat(str)
  }

  # setter functions

  setPars <- function(paramDt) {
    stopifnot(c("IDX","ADJUSTABLE","PARNAME") %in%
                names(paramDt))
    thisParamDt <<- copy(paramDt)
    setkey(thisParamDt, IDX)
    numPars <<- nrow(thisParamDt)
  }
  
  setParTrafo <- function(fun, jac) {
    stopifnot(is.function(fun))
    stopifnot(is.function(jac))
    thisParTrafoFun <<- fun
    thisParTrafoJac <<- jac
  }

  setNeeds <- function(needsDt) {
    stopifnot(all(needsDt$PROJECTILE[1] == needsDt$PROJECTILE))
    stopifnot(all(needsDt$ELEMENT[1] == needsDt$ELEMENT))
    stopifnot(all(needsDt$MASS[1] == needsDt$MASS))
    thisNeedsDt <<- copy(needsDt)
    setkey(thisNeedsDt, IDX)
    numNeeds <<- nrow(thisNeedsDt)
  }

  setSexp <- function(Sexp) {
    thisSexp <<- copy(Sexp)
  }

  setMask <- function(mask) {
    thisMask <<- mask
  }
  
  setEps <- function(eps) {
    thisEps <<- eps
  }
  
  setCacheSize <- function(n) {
    
    if (cacheSize < n) {
      cache <<- c(cache, replicate(n-cacheSize, list(), simplify=FALSE))
    } 
    else if (cacheSize > n) {
      cache <- cache[1:n]
      if (cacheIdx >= n) cacheIdx <<- 0
    }
    cacheSize <<- n
  }

  # helper functions

  getCachedVariable <- function(x, varname) {
    idx <- which(sapply(cache, function(el) !is.null(el$x) && all(el$x == x)))
    stopifnot(length(idx) <= 1)
    if (length(idx)==0) NULL else cache[[idx]]$data[[varname]]
  }

  updateCache <- function(x, ...) {
    idx <- which(sapply(cache, function(el) !is.null(el$x) && all(el$x == x)))
    stopifnot(length(idx) <= 1)
    if (length(idx) == 0) {
      cacheIdx <<- (cacheIdx %% cacheSize) + 1
      cache[[cacheIdx]] <<- list(x = x, data = list())
      idx <- cacheIdx
    }
    cache[[idx]]$data <<- modifyList(cache[[idx]]$data, list(...))
  }

  getCache <- function() {
    cache
  }

  getRefCalcJobs <- function() {
    refCalcJobs
  }

  getJacCalcJobs <- function() {
    jacCalcJobs
  }
  
  getEps <- function() {
    thisEps
  }

  # functions

  fun <- function(x, applySexp = TRUE, ret.dt = FALSE, saveDir = NULL) {

    stopifnot(is.numeric(x))
    x <- as.matrix(x)
    stopifnot(! (ncol(x) > 1 & ret.dt))
    trafo_x <- x
    if (!is.null(thisParTrafoFun))
      x <- thisParTrafoFun(x)
    
    y <- matrix(NA_real_, nrow = nrow(thisNeedsDt), ncol = ncol(x))
    isCached <- rep(FALSE, ncol(x))
    inpList <- list()
    for (j in seq_len(ncol(x))) {
      curx <- x[,j]
      funval <- getCachedVariable(curx, "funval")
      if (!is.null(funval))
      {
        isCached[j] <- TRUE
        y[,j] <- funval
      }
      else
      {
        setkey(thisParamDt, IDX)
        setkey(thisNeedsDt, IDX)
        thisParamDt[ADJUSTABLE == TRUE, PARVAL := as.list(curx)]
        inputDt <- convertToInput(thisParamDt, thisNeedsDt)
        inpList <- c(inpList, inputDt$inputs)
      }
    }
    
    if (!all(isCached)) {
      
      redNeedsDt <- copy(thisNeedsDt)
      redNeedsDt[, c("PROJECTILE","MASS","ELEMENT","V1"):=NULL]
      refCalcJobs <<- talysClust$run(inpList, redNeedsDt, saveDir = saveDir,
                                     calcsPerJob = 1000, runOpts=list(TMPDIR="/dev/shm/talysTemp"))
      
      while(talysClust$isRunning(refCalcJobs, combine=TRUE)) Sys.sleep(30)
      res <- talysClust$result(refCalcJobs)
      calcIdcs <- which(!isCached)
      for (k in seq_along(calcIdcs)) {
        curCalcIdx <- calcIdcs[k]
        curx <- x[,curCalcIdx]
        curRes <- res[[k]]$result
        stopifnot(all(c("IDX","V1") %in% names(curRes)))
        # stopifnot(all(c("IDX","V1") %in% names(curRes)))
        setkey(curRes, IDX)
        funval <- curRes[,V1]
        y[,curCalcIdx] <- funval
        updateCache(curx, funval = funval)
      }
    }
    if (ret.dt) {
      setkey(thisNeedsDt, IDX)
      thisNeedsDt[, V1 := y[,1]]
      copy(thisNeedsDt)
    } else {
      if (is.null(thisSexp) || !applySexp)
        y
      else
        as.matrix(thisSexp %*% y)
    }
  }

  jac <- function(x, applySexp = TRUE) {

    stopifnot(is.numeric(x))
    trafo_x <- x
    if (!is.null(thisParTrafoJac))
      x <- thisParTrafoFun(x)
    
    jacmat <- getCachedVariable(x, "jacmat")
    if (is.null(jacmat)) {
      setkey(thisParamDt, IDX)
      thisParamDt[ADJUSTABLE == TRUE, PARVAL := as.list(x)]
      jacInputsDt <- createInputsForJacobian(thisParamDt, thisNeedsDt, mask = thisMask, eps = thisEps)

      newJacCalcId <<- digest(jacInputsDt)
      if (!identical(jacCalcId, newJacCalcId) || !is.null(refCalcJobs)) {
        jacCalcJobs <<- talysClust$run(jacInputsDt$inputs, jacInputsDt$outspecs,
                                      calcsPerJob = 1000, runOpts=list(TMPDIR="/dev/shm/talysTemp"))
        jacCalcId <<- newJacCalcId
      }
      
      while(talysClust$isRunning(jacCalcJobs, combine=TRUE)) Sys.sleep(30)
      jacRes <- talysClust$result(jacCalcJobs)
      jacInputsDt$outspecs <- lapply(jacRes, function(x) x$result)
      SparDt <- computeJacobian(jacInputsDt, drop0 = TRUE)
      jacmat <- SparDt[, sparseMatrix(i = IDX1, j = IDX2, x = X, dims = c(numNeeds, numPars))]
      setkey(thisParamDt, IDX)
      jacmat <- jacmat[, thisParamDt$ADJUSTABLE==TRUE, drop=FALSE]
      refNeedsDt <- jacInputsDt[ADJIDX==0, outspecs[[1]]]
      setkey(refNeedsDt, IDX)
      funval <- refNeedsDt[["V1"]]
      updateCache(x, jacmat = jacmat, funval = funval)
      jacCalcId <<- NULL
    }
    res <- jacmat
    if (!is.null(thisParTrafoJac)) res <- res %*% thisParTrafoJac(trafo_x)
    if (is.null(thisSexp) || !applySexp) res else drop0(thisSexp %*% res)
  }

  list(
    # initialization
    setPars = setPars,
    setParTrafo = setParTrafo,
    setNeeds = setNeeds,
    setSexp = setSexp,
    setMask = setMask,
    setEps = setEps,
    # calculation
    fun = fun,
    jac = jac,
    # debugging
    getCache = getCache,
    getRefCalcJobs = getRefCalcJobs,
    getJacCalcJobs = getJacCalcJobs,
    getEps = getEps
  )
}

