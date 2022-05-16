mult_invCov <- function(D, S, P = NULL, cholZ = NULL, Di = FALSE) {
  
  stopifnot(grepl("Matrix$", class(D)),
            is.null(P) || grepl("Matrix$", class(P)),
            grepl("Matrix$", class(S)))
  
  invD <- NULL
  invCov <- NULL
  
  if(Di) {
    invD <- D
    
      if(is.null(cholZ)){
        Z <- forceSymmetric(crossprod(S, crossprod(invD, S)) + solve(P))
        cholZ <- makeCholesky(Z)
        invCov <-invD - crossprod(invD, S %*% solve(cholZ, crossprod(S, invD)))
    } 
  
  } else {
    invD <- solve(D)
    Z <- forceSymmetric(crossprod(S, crossprod(invD, S)) + solve(P))
    cholZ <- makeCholesky(Z)
    invCov <- invD - crossprod(invD, S %*% solve(cholZ, crossprod(S, invD)))
    }
  
  invCov

}
woodburry <- function(D, S, P, invD=FALSE){
  Di <- NULL
  
  if(invD){
    Di <- D
  } else {
     Di <- solve(D) 
    }
  
  Pi <- solve(P)
  central_term <- solve(Pi + t(S) %*% Di %*% S)
  Di - Di %*% S %*% central_term %*% t(S) %*% Di
}

predict_lin <- function(Xi, S, Kzz, K, y, yhat){
  

  
  d <- as.vector(y - yhat)  
  Sx <- tcrossprod(S %*% Kzz, S)
  D2i <- diag(1/diag(K-Sx))
  SD2iS <- t(S) %*% D2i %*% S
  covi <- woodburry(Xi, Sx, Kzz+SD2iS, invD = TRUE)
  
    
  Sk <- crossprod(Sx, covi)
  
  mu_xp <- as.vector(yhat + Sk %*% d)
  var_fxp <- K - tcrossprod(Sk, Sx)
  
  
  list("mu"=mu_xp,
       "cov"=var_fxp)
}

gp_post_est <- function(Xi, S, Kzz, K, y, yhat) {
  
  stopifnot(grepl("Matrix$", class(Xi)),
            is.null(K) || grepl("Matrix$", class(K)),
            grepl("Matrix$", class(S)))
  
    r <- y - yhat
    
    SK <- S %*% Kzz
    XiS <- crossprod(Xi, S)
    
    StXiS <- crossprod(S, XiS)
    SKStXiS <- SK %*% StXiS
    StXiSKSt <- t(SKStXiS)

    Z <- forceSymmetric(StXiS + solve(Kzz))
    cholZ <- makeCholesky(Z)
    ZiStXiSKS <- solve(cholZ, StXiSKSt)
    StKSXiStZi <- t(ZiStXiSKS)
    invCov <- tcrossprod(SKStXiS, SK) - SKStXiS %*% ZiStXiSKS

    mu <- as.vector(yhat+crossprod(tcrossprod(SK - StKSXiStZi, XiS), r))   
  
    cov <-tcrossprod(SK,S) - invCov
  
  list("mu"=mu, "cov"=cov)
  
}

gpPostExact <- function(expDt, gpDt, mapAssignment){
  gpHandeler <- createSysCompGPHandler()
  
  expData <- getDt_DATA(expDt)
  refExpData <- getDt_REFDATA(expDt)
  
  Kxx <- gpHandeler$mapExact(expDt, gpDt, mapAssignment, ret.mat = TRUE)
  K <- gpHandeler$covExact(expDt, gpDt, mapAssignment, ret.mat = TRUE)
  invK <- gpHandeler$invCovExact(expDt, gpDt, mapAssignment, ret.mat = TRUE)
  Sxx <- Kxx %*% invK
  
  mu <- as.vector(refExpData + Sxx %*% (expData - refExpData))
  cov <- K - Sxx %*% Kxx
  
  list("mu"=mu, "cov"=cov)
}

getPosteriorGp <- function(expDt, sysDt, gpDt = NULL, sysCompHandler) {
  
  setkey(sysDt, IDX)
  setkey(expDt, IDX)
  
  expSel <- optSysDt[, (!grepl("TALYS-", EXPID) & !grepl("REACEXP-", EXPID))]
  gpOBsSel <- optSysDt[, grepl("REACEXP-", EXPID)]
  
  S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)
  cov_y <- Diagonal(x = getDt_UNC(expDt)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)
  

  X <- P[expSel, expSel] 
  SX <- S[, expSel]
  
  y <- getDt_DATA(expDt)
  mu <- getDt_REFDATA(expDt)
  mup <- getDt_REFDATA(expDt)
  d <- as.vector(y - mu)  
  
  S_xk <- S[, gpOBsSel]
  S_xpk <- S[, gpOBsSel]
  
  K_kk <- P[gpOBsSel, gpOBsSel] 
  k <- S_xk %*% K_kk %*% t(S_xk)
  k_xpxp <- kK
  print("TEST 1")
  
  SKS_xkk <- k
  Di <- mult_invCov(D=Diagonal(x = 1/getDt_UNC(expDt)^2), S=SX, P=X, Di=TRUE)
  #Di <- Diagonal(x = 1/getDt_UNC(expDt)^2)
  covi <- mult_invCov(D=Di, S=S_xk, P=K_kk, Di=TRUE)
  
  print("TEST 2")
  
  #S_xpx <- S_xpk %*% K_kk %*% t(S_xk)
  S_xpx <- k
  #mu_xp <- as.vector(mup + covi %*% d)
  mu_xp <- as.vector(mup + S_xpx %*% covi %*% d)
  var_fxp <- k_xpxp - S_xpx %*% covi %*% t(S_xpx)
  
  print("TEST 3")
  
  #expDtPost[, GPOBSPRIOR := gpObsPrior]
  
  newExpDt <- copy(expDt)
  newExpDt[, GPMEAN := mu_xp]
  newExpDt[, GPVAR := sqrt(diag(var_fxp))]  
  list(expDt = newExpDt, cov = var_fxp)
}

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


#' Get Posterior of Systematic Components
#'
#' @param expDt datatable with experimental data
#' @param sysDt datatable with systematic components
#' @param gpDt datatable with Gaussian process components
#' @param sysCompHandler systematic component handler, see \link{createSysCompHandler}
#'
#' @return list with two elements \code{sysDt} and \code{cov} which contain
#' the posterior of the systematic components and the covariance matrix, respectively.
#' and covariance matrix
#' @export
#'
getPosteriorGpSys <- function(expDt, sysDt, gpDt = NULL, sysCompHandler) {
  
  setkey(expDt, IDX)
  S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)
  D <- Diagonal(x = getDt_UNC(expDt)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)
  
  expSel <- sysDt[, (!grepl("TALYS-", EXPID) & !grepl("pw", ERRTYPE))]
  talysSel <- !expSel
  
  P0 <- P[talysSel, talysSel]
  X <- P[expSel, expSel] 
  SX <- S[, expSel]
  
  optSysDtGpObs <- sysDt[grepl("REACEXP-", EXPID), ][, IDX := seq_len(.N)]
  Sk <- sysCompHandler$map(expDt, optSysDtGpObs, ret.mat = TRUE)
  Kzz <- sysCompHandler$cov(optSysDtGpObs, gpDt, ret.mat = TRUE)
  
  invD <- solve(D)
  invP <- solve(forceSymmetric(P))
  p0 <-getDt_DATA(sysDt)
  expData <- getDt_DATA(expDt)
  refExpData <- getDt_REFDATA(expDt)
  refSysData <- getDt_REFDATA(sysDt)
  effData <- expData - refExpData + S %*% refSysData
  gpHandeler <- createSysCompGPHandler()
  K <- gpHandeler$covExact(expDt, gpDt, mapAssignment, ret.mat = TRUE)
  Dk <- diag(diag(K) - diag(Sk %*% Kzz %*% t(Sk)))
  DXi <- mult_invCov(D+Dk, SX, X)

  P1 <- solve(forceSymmetric(t(S) %*% invD %*% S) + invP)
  p1 <- as.vector(P1 %*% (invP %*% p0 + t(S) %*% invD %*% effData))
  
  newSysDt <- copy(sysDt)
  setDt_DATA(newSysDt, p1)
  setDt_UNC(newSysDt, sqrt(diag(P1)))
  list(sysDt = newSysDt, cov = P1, Kzz = Kzz, K=K, Sk = Sk, SX=SX, DXi = DXi)
}




#' Map Systematic Components to Experimental Grid
#'
#' @param expDt  datatable with experimental data
#' @param sysList list with components \code{sysDt} and potentially \code{cov}
#' @param sysCompHandler handler for systematic components as created by \link{createSysCompHandler}
#' @param ret.cov should the posterior matrix be computed?
#' @param idx vector of indices that specifies which systematic components should be mapped.
#' Default \code{NULL} means to map all of them.
#'
#' @return list with two components \code{expDt} and \code{cov}
#' @export
#'
mapGpToExp <- function(expDt, gpDt, mapAssignment, sysList, sysCompHandler, ret.cov = FALSE, idx = NULL) {
  
  if (is.null(idx)) idx <- seq_len(nrow(sysList$sysDt))
  sysDt <- sysList$sysDt
  SX <- sysList$SX
  S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)[,idx]
  refExpData <- getDt_REFDATA(expDt)
  expData <- getDt_DATA(expDt)
  refSysData <- getDt_REFDATA(sysDt, idx = idx)
  sysData <- getDt_DATA(sysDt, idx = idx)
  effData <- expData - refExpData 
  newExpDt <- copy(expDt)
  xsCovmat <- NULL
  
  if (!is.null(sysList$cov)) {
    curCov <- sysList$cov[idx,idx]
    setDt_UNC(newExpDt, sqrt(rowSums((S %*% curCov) * S)))
    if (isTRUE(ret.cov)) {
      xsCovmat <- S %*% curCov %*% t(S)
    }
    
  K <- sysList$K
  Kzz <- sysList$Kzz
  Sk <- sysList$Sk
  DXi <- sysList$DXi
  #D <- Diagonal(x = getDt_UNC(expDt)^2)
  #DXi <- Diagonal(x = (1/expDt$UPDUNC)^2)
  
  gp_post <- predict_lin(Xi = DXi, S=Sk, Kzz=Kzz, K=K, y=expData, yhat=refExpData)
  #gp_post <- gp_post_est(Xi = DXi, S=Sk, Kzz=Kzz, K=K, y=expData, yhat=refExpData)
  gp_post_mu <- gp_post$mu
  gp_post_cov <- gp_post$cov
  
  #gp_post <- gpPostExact(expDt = expDt, gpDt = gpDt, mapAssignment = mapAssignment)
  #gp_post_mu <- gp_post$mu
  #gp_post_cov <- gp_post$cov
  newExpDt[, GPMU := gp_post_mu]  
  newExpDt[, GPSIG := sqrt(diag(gp_post_cov))]

  #invCov <- mult_invCov(D=DXi, S=Sk, P=K, Di=TRUE)
  #SKS <- Sk %*% K %*% t(Sk)  
  #SKSinvCov <- SKS %*% invCov
  #newExpDt[, GPMU := as.vector(refExpData + crossprod(SKSinvCov, effData))]  
  #newExpDt[, GPSIG := sqrt(diag(SKS - SKSinvCov %*% SKS))]
  }
  list(expDt = newExpDt[], cov = xsCovmat)
}


plot_gp_post <- function(expDt, optGpDt){
  expDt[,{  
    
    
    curGpDt <-  optGpDt[EXPID == mapAssignment[REAC==.BY, EXPID]]
    sig <- curGpDt[PARNAME=="sigma", PARVAL]
    len <- curGpDt[PARNAME=="len", PARVAL]
    nug <- curGpDt[PARNAME=="nugget", PARVAL]
    
    lableString <- sprintf("sigma == %3.2f~~lambda == %3.2f~~tau ==  %3.2f  ", sig, len, nug)
    
    ggr <- ggplot(.SD) + theme_bw() + ggtitle(.BY)
    
    #r <- ORIGDATA-REFDATA 
    r <- DATA 
    
    #rpred <- DATA-DATAREF
    rpred <- GPMEAN
    
    ggr <- ggr + geom_point(aes(x = L1, y = r), size = 0.5)
    
    alpha <- 1 - .N /nrow(optExpDt)
    
    ggr <- ggr + geom_errorbar(aes(x = L1, ymin = r - UNC, ymax= r + UNC), color="gray", size = 0.8, alpha=alpha)  
    
    
    ggr <- ggr + geom_line(aes(x = L1, y = DATAREF), color="black", size = 0.5)
    ggr <- ggr + geom_ribbon(aes(x = L1, ymin = rpred-GPVAR, ymax = rpred+GPVAR), fill="orange", alpha=0.3)  
    ggr <- ggr + geom_line(aes(x = L1, y = rpred), color="orange", size = 0.5) 
    
    ggr <- ggr + annotate("text", x = Inf, y = Irnf, label = lableString, hjust =1.1, vjust = 3.1, size=5, parse = TRUE)   
    
    print(ggr)
    filepath <- file.path(rootpath, outdataDir, "plots","06", "OPT_OBS_GP") 
    filename <- paste0(curReac, "_GP.pdf")
    if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
    
    ggsave(file.path(filepath, filename), ggr, width = 36, height = 24, units = "cm", dpi = 300)
  }
  , by=REAC]
}

plot_bayes_update <- function(expDt, optGpDt){
  expDt[,{  
    
    curReac <- .BY
    
    curGpDt <-  optGpDt[EXPID %in% mapAssignment[REAC==.BY, EXPID]]
    sig <- curGpDt[PARNAME=="sigma", PARVAL]
    len <- curGpDt[PARNAME=="len", PARVAL]
    nug <- curGpDt[PARNAME=="nugget", PARVAL]
    
    lableString <- sprintf("sigma == %3.2f~~lambda == %3.2f~~tau ==  %3.2f  ", sig, len, nug)
    
    ggr <- ggplot(.SD) + theme_bw() + ggtitle(.BY)
    
    #r <- ORIGDATA-REFDATA 
    r <- ORIGDATA
    
    rpred <- DATA
    #rpred <- GPMEAN
    
    ggr <- ggr + geom_point(aes(x = L1, y = r), size = 0.5)
    
    alpha <- 1 - .N /nrow(optExpDt)
    
    
    
    ggr <- ggr + geom_line(aes(x = L1, y = DATAREF), color="black", size = 0.5)
    ggr <- ggr + geom_ribbon(aes(x = L1, ymin = GPMU-GPSIG, ymax = GPMU+GPSIG), fill="orange", alpha=0.3)  
    ggr <- ggr + geom_line(aes(x = L1, y = GPMU), color="orange", size = 0.5) 
    
    #ggr <- ggr + geom_errorbar(aes(x = L1, ymin = r - ORIGUNC, ymax= r + ORIGUNC), color="gray", size = 0.8, alpha=alpha)  
    ggr <- ggr + geom_errorbar(aes(x = L1, ymin = r - UPDUNC, ymax= r + UPDUNC), color="gray", size = 0.8, alpha=alpha)
    
        
    ggr <- ggr + annotate("text", x = Inf, y = Inf, label = lableString, hjust =1.1, vjust = 3.1, size=5, parse = TRUE)   
    
    print(ggr)
    filepath <- file.path(plotPath,"06", "OPT_OBS_GP") 
    filename <- paste0(curReac, "_GP.pdf")
    if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
    
    path_to_file <- file.path(filepath, filename)
    print(path_to_file)
    ggsave(path_to_file, ggr, width = 36, height = 24, units = "cm", dpi = 300)
  }
  , by=REAC]
}

plot_gp_prior <- function(expDt, optGpDt){
  expDt[,{  
    
    curReac <- .BY
    
    curGpDt <-  optGpDt[EXPID %in% mapAssignment[REAC==.BY, EXPID]]
    sig <- curGpDt[PARNAME=="sigma", PARVAL]
    len <- curGpDt[PARNAME=="len", PARVAL]
    nug <- curGpDt[PARNAME=="nugget", PARVAL]
    
    #lableString <- sprintf("sigma == %3.2f~~lambda == %3.2f~~tau ==  %3.2f  ", sig, len, nug)
    
    ggr <- ggplot(.SD) + theme_bw() + ggtitle(.BY)
    
    r <- ORIGDATA-REFDATA 
    #r <- ORIGDATA
    
    rpred <- DATA
    #rpred <- GPMEAN
    

    
    alpha <- 1 - .N /nrow(optExpDt)
    
    
    
    #ggr <- ggr + geom_line(aes(x = L1, y = DATAREF), color="black", size = 0.5)
    ggr <- ggr + geom_ribbon(aes(x = L1, ymin = 0-2*GPPRIOR, ymax = 0+2*GPPRIOR), fill="olivedrab", alpha=0.4) 
    ggr <- ggr + geom_ribbon(aes(x = L1, ymin = 0-GPPRIOR, ymax = 0+GPPRIOR), fill="olivedrab", alpha=0.4)  
    ggr <- ggr + geom_line(aes(x = L1, y = 0), color="black", linetype="dashed", size = 0.5) 
    
  
    ggr <- ggr + geom_errorbar(aes(x = L1, ymin = r - ORIGUNC, ymax= r + ORIGUNC), color="gray", size = 0.8, alpha=alpha)  
    ggr <- ggr + geom_point(aes(x = L1, y = r), size = 0.5)      
    #ggr <- ggr + annotate("text", x = Inf, y = Irnf, label = lableString, hjust =1.1, vjust = 3.1, size=5, parse = TRUE)   
    
    print(ggr)
    filepath <- file.path(plotPath,"06", "OPT_OBS_GP_PRIOR") 
    filename <- paste0(curReac, "_GP_PRIOR.pdf")
    if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
    
    path_to_file <- file.path(filepath, filename)
    print(path_to_file)
    ggsave(path_to_file, ggr, width = 36, height = 24, units = "cm", dpi = 300)

  }
  , by=REAC]
}

source("config.R")
subents <- read_object(1, "subents")
modList <- read_object(3, "modList")
optSysDt <- read_object(6, "optSysDt")
optGpDt <- read_object(6, "optGpDt")
optExpDt <- read_object(6, "optExpDt")
mapAssignment <- read_object(6, "mapAssignment")
sysCompHandler <- read_object(6, "sysCompHandler")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")

normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler2 <- createSysCompHandler()
sysCompHandler2$addHandler(normHandler)

S <- sysCompHandler2$map(expDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler2$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler2$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(expDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S))) 
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))
#expDtPost <- gpPostList$expDt

optExpDt[, UPDUNC := updUnc]
#gpPostList <-getPosteriorGp(expDt=optExpDt, sysDt=optSysDt, gpDt = optGpDt, sysCompHandler=sysCompHandler)

#sysList <- getPosteriorSys(expDt=optExpDt, sysDt=optSysDt, gpDt = optGpDt, sysCompHandler=sysCompHandler)
#expDtList <- mapSysToExp(expDt=optExpDt, sysList=sysList, sysCompHandler=sysCompHandler, ret.cov = TRUE, idx = NULL)

sysList <- getPosteriorGpSys(expDt=optExpDt, sysDt=optSysDt, gpDt = optGpDt, sysCompHandler=sysCompHandler)
expDtList <- mapGpToExp(expDt=optExpDt, gpDt = optGpDt, mapAssignment = mapAssignment, sysList=sysList, sysCompHandler=sysCompHandler, ret.cov = TRUE, idx = NULL)

expDtBay <- expDtList$expDt


optSysDtGpObs <- optSysDt[grepl("REACEXP-", EXPID), ][, IDX := seq_len(.N)]
Sk <- sysCompHandler$map(optExpDt, optSysDtGpObs, ret.mat = TRUE)
K <- sysCompHandler$cov(optSysDtGpObs, optGpDt, ret.mat = TRUE)
statUnc <- getDt_UNC(optExpDt)

gpObsPrior <- sqrt(diag(Sk %*% K %*% t(Sk)))


expDtBay[, GPPRIOR := gpObsPrior]
expDtBay[, ORIGUNC := optExpDt$UNC]
expDtBay[, ORIGDATA := optExpDt$DATA]

library(ggplot2)

plot_gp_prior(expDtBay, optGpDt)
plot_bayes_update(expDtBay, optGpDt)


#gpHandeler <- createSysCompGPHandler()
#Kxz <- gpHandeler$map(optExpDt, optSysDt, optGpDt, mapAssignment, ret.mat = TRUE)
#Kxx <- gpHandeler$covExact(optExpDt, optGpDt, mapAssignment, ret.mat = TRUE)
#invKxx <- gpHandeler$invCovExact(optExpDt, optGpDt, mapAssignment, ret.mat = TRUE)

#optExpDt[, {curGPSpec <- optGpDt[EXPID %in% mapAssignment[REAC == .BY, EXPID]]}, by="REAC"]

