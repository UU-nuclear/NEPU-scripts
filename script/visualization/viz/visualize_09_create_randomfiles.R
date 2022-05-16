
source('config.R')
library(ggplot2)

extNeedsDt <- read_object(2, "extNeedsDt")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
sysDt <- updSysDt[grepl('^EXPID-', EXPID)]
sysDt[, IDX:=seq_len(.N)]
needsDt <- read_object(1, "needsDt")
retrieve_option <- 2L

# create data table for extrapolation
extNeedsDt <- needsDt[,{
  stopifnot(all(L2 == 0) & all(L3 == 0))
  list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
       L2 = 0, L3 = 0, V1 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]

extNeedsDt[, IDX := seq_len(.N)]
# option 1: retrieve the information from resulting tarfiles
# NOTE: savePathTalys must be accessible from within the container
#       this may not be the case if calculations have not been run
#       inside the Docker container
if (retrieve_option == 1L) {
  
  tarfiles <- list.files(savePathTalys, full.names=TRUE, pattern='\\.tar.gz')
  mytempdir <- file.path(tempdir(), 'visualize_09_create_randomfiles_tmpdir')
  dir.create(mytempdir)  
  
  talysHnd <- createModelTALYS()
  
  allResults <- NULL
  for (idx in seq_along(tarfiles[1:20])) {
    tarfile <- tarfiles[idx]
    cat('processing number',idx,'---',tarfile,'\n')
    # unpack the tarfile
    untarcmd <- paste0('tar -xC "',mytempdir,'" -f "',tarfile,'"')
    retcode <- system(untarcmd, intern=FALSE, wait=TRUE)
    if (retcode != 0) {
      stop(paste0('Problem unpacking file ', tarfile, ' - retcode is ',retcode))
    }
    # retrieve result
    resultDt <- talysHnd$read(mytempdir, extNeedsDt, packed=FALSE)
    setkey(resultDt, IDX)
    allResults <- cbind(allResults, resultDt$V1)
    # clean up
    file.remove(list.files(mytempdir, full.names=TRUE))
  }
  
} else if (retrieve_option==2) {
  
  # option 2: retrieve the results from object allResults created in step 09
  allResults <- read_object(9, 'allResults')
}

# plot the results

# reconstruct experimental covariance matrix
sysCompHandler <- createSysCompHandler()
normHandler <- createSysCompNormHandler(dataref='DATAREF')
normHandler$addSysUnc('EXPID', unique(expDt$EXPID), 0, 0, TRUE)
sysCompHandler$addHandler(normHandler)
S_syserr <- sysCompHandler$map(expDt, sysDt, ret.mat=TRUE)
Cov_syserr <- sysCompHandler$cov(sysDt, ret.mat=TRUE)
totunc_exp <- sqrt(expDt$UNC^2 + rowSums((S_syserr %*% Cov_syserr) * S_syserr))
tmpExpDt <- copy(expDt)
setkey(tmpExpDt, IDX)
tmpExpDt[, UNC := totunc_exp]


# create model grid for experimental data
exforReacStrs <- unique(expDt$REAC)
modSubents <- lapply(exforReacStrs, createSubentStub, en=energyGrid)
modDt <- exforHandler$extractData(modSubents, ret.values=FALSE)
Smod <- exforHandler$getJac(modDt, extNeedsDt, modSubents)
setkey(modDt, IDX, DIDX)

# for plotting uncertainty bands
tmpDt <- copy(extNeedsDt)
setkey(tmpDt, IDX)
uncinfo <- cov.wt(t(allResults))
tmpDt[, V1:=uncinfo$center]
tmpDt[, UNC:=sqrt(diag(uncinfo$cov))]
plotDt <- copy(modDt)

plotDt[, DATA:=as.vector(Smod %*% tmpDt$V1)]
plotDt[, UNC:=as.vector(sqrt(diag(Smod %*% diag(tmpDt$UNC) %*% t(Smod))))]

ggp <- ggplot(data=plotDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('cross section [mbarn]')
# overlay experimental data
ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC, col=EXPID), data=tmpExpDt)
ggp <- ggp + geom_point(aes(x=L1, y=DATA, col=EXPID), data=tmpExpDt, size=0.2)
# plot the model
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC), alpha=0.3)
ggp <- ggp + facet_wrap(~REAC, scales='free_y') + theme(strip.background = element_blank())

ggp

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(plotPath, 'plot_posterior_xs.png'), ggp,
       units='cm', width=17.8, height=10)



