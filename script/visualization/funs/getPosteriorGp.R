getPosteriorGp <- function(expDt, sysDt, gpDt = NULL, sysCompHandler) {
  
  setkey(sysDt, IDX)
  setkey(expDt, IDX)
  
  expSel <- optSysDt[, (!grepl("TALYS-", EXPID) & !grepl("REACEXP-", EXPID))]
  gpOBsSel <- optSysDt[, grepl("REACEXP-", EXPID)]
  
  
  setkey(expDt, IDX)
  S <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE)
  D <- Diagonal(x = getDt_UNC(expDt)^2)
  P <- sysCompHandler$cov(sysDt, gpDt, ret.mat = TRUE)
  
  K <- P[gpOBsSel, gpOBsSel] 
  X <- P[expSel, expSel] 
  SX <- S[, expSel]
  SK <- S[, gpOBsSel]
  
  DXi <- mult_invCov(D=D, S=SX, P=X)
  DX <- D + SX %*% X %*% t(SX)
  invK <- solve(forceSymmetric(K))
  KX <- SK %*% invK
  SKS <- KX %*% t(SK)
  
  expData <- getDt_DATA(expDt)
  refExpData <- getDt_REFDATA(expDt)
  refSysData <- getDt_REFDATA(sysDt)
  
  r <- expData - refExpData    
  
  Zi <- solve(forceSymmetric(invK + t(SK) %*% DXi %*% SK))
  KZiK <- KX %*% Zi %*% t(KX)
  SZK <- SKS %*% KZiK 
  Kpost <- SKS - SZK %*% SKS
  muPost <- refExpData + SZK %*% r
  #Kpost <- DX + KZiK
  #muPost <- as.vector(refExpData + KZiK %*% DXi %*% r)
  
  newExpDt <- copy(expDt)
  setDt_DATA(newExpDt, muPost)
  setDt_UNC(newExpDt, sqrt(diag(Kpost)))
  list(expDt = newExpDt, cov = Kpost)
}