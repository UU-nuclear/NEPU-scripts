##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

library(matrixStats)

extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
expDt <- read_object(3, "expDt")
modDt <- read_object(3,"modDt")
allResults <- read_object(9, 'allResults')

normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

S <- sysCompHandler$map(expDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(expDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S))) 
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))

setkey(expDt, IDX)
expDt[, ORIGUNC := origUnc]
expDt[, UPDUNC := updUnc]

# --------------------------------------

# create model grid for experimental data
exforReacStrs <- unique(expDt$REAC)
modSubents <- lapply(exforReacStrs, createSubentStub, en=energyGrid)
#modSubents <- lapply(exforReacStrs, createSubentStub, en=energyGridrandomFiles)
modDt_post <- exforHandler$extractData(modSubents, ret.values=FALSE)
Smod <- exforHandler$getJac(modDt_post, extNeedsDt, modSubents)
setkey(modDt_post, IDX, DIDX)

# for plotting uncertainty bands
tmpDt <- copy(extNeedsDt)
setkey(tmpDt, IDX)
uncinfo <- cov.wt(t(allResults))
tmpDt[, V1:=uncinfo$center]
tmpDt[, UNC:=sqrt(diag(uncinfo$cov))]
postDt <- copy(modDt_post)

postDt[, DATA:=as.vector(Smod %*% tmpDt$V1)]
postDt[, UNC:=as.vector(diag(Smod %*% diag(tmpDt$UNC) %*% t(Smod)))]
# --------------------------------------

reactions <- expDt[,unique(REAC)]
curReac <- "(26-FE-56(N,TOT),,SIG)"

curExpDt <- expDt[REAC == curReac]
curModDt <- modDt[REAC == curReac]
curPostDt <- postDt[REAC == curReac]
curTmpDt <- tmpDt[REAC == "CS/TOT"]


#============================================================

#Apparently this is wrong, the data stored in the matrix is not directly the results
# of the calculation but the result mulitpied with some sensitivity matrix
all_total_xs <- allResults[1:length(energyGrid),] # matrix[length(energyGrid),n_calc]
mean_total_xs <- rowMeans(all_total_xs)
sds_total_xs <- rowSds(all_total_xs)

# I believe that the best paramter estimate is calc 001
best_total_xs <- allResults[1:length(energyGrid),1]

total_xs_df <- data.frame(energyGrid,mean_total_xs,sds_total_xs,best_total_xs)

ggp <- ggplot(total_xs_df,aes(x=energyGrid,y=mean_total_xs))
ggp <- ggp + geom_line()
ggp <- ggp + geom_ribbon(aes(ymin=mean_total_xs-sds_total_xs,ymax=mean_total_xs+sds_total_xs,alpha=0.3))
ggp <- ggp + geom_line(aes(y=best_total_xs),col='red')

ggp <- ggp + geom_line(data=curPostDt,aes(x=L1, y=DATA),col="green")
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3, data=curPostDt,fill="green")