library(ggplot2)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################
outdataPathRun <- outdataPath
plotPath <- paste0(outdataPathRun, "/plots")


finalPars <- read_object(08, "finalPars", outdata_path = outdataPathRun)
finalParCovmat <- read_object(08, "finalParCovmat", outdata_path = outdataPathRun)
optSysDt_allpars <- read_object(07, "optSysDt_allpars", outdata_path = outdataPathRun) # all parameters
optSysDt_optpars <- read_object(07, "optSysDt_optpars", outdata_path = outdataPathRun) # parameters considered for optimization
optParamDt <- read_object(07,"optParamDt")
optGpDt <- read_object(06,"optGpDt")

plotDt <- copy(optSysDt_allpars)
setkey(plotDt, IDX)
plotDt$DATA <- paramTrafo$fun(finalPars)
plotDt$UNC <- sqrt(diag(finalParCovmat))
plotDt$DATAMIN <- paramTrafo$fun(finalPars - sqrt(diag(finalParCovmat)))
plotDt$DATAMAX <- paramTrafo$fun(finalPars + sqrt(diag(finalParCovmat)))
plotDt$EXPID <- sub('TALYS-','', plotDt$EXPID)

# sort the EXPIDs according to uncertainty
expid_order <- order(plotDt[,abs(DATAMAX-DATAMIN)])
plotDt$EXPID <- factor(plotDt$EXPID, levels=unique(plotDt$EXPID[expid_order]), ordered=TRUE)

###########################################
#    HYPER-PARAMETER VALUES
###########################################
GPsigmas <- optGpDt[PARNAME=="sigma"]
GPlenths <- optGpDt[PARNAME=="len"]

###########################################
#    ENERGY-DEPENDENT PARAMETERS
###########################################

#varnames <- plotDt[ERRTYPE=='talyspar_endep', unique(EXPID)]
# sensitive energy-dependent parameters
varnames <- c('TALYS-v1adjust a', 'TALYS-vso1adjust p',
              'TALYS-d1adjust p', 'TALYS-v1adjust p', 'TALYS-vso1adjust n',
              'TALYS-w1adjust n', 'TALYS-d1adjust n', 'TALYS-v1adjust n')
              #'TALYS-rcadjust p')
varnames <- sub('TALYS-','',varnames)

varnames_all <- unique(plotDt[ERRTYPE=="talyspar_endep"]$EXPID)

# varname <- 'TALYS-v1adjust n'
#ggp <- ggplot(data=plotDt[ERRTYPE=='talyspar_endep' & EXPID %in% varnames_all[1:8],])
#ggp <- ggplot(data=plotDt[ERRTYPE=='talyspar_endep' & EXPID %in% varnames_all[9:18],])
#ggp <- ggplot(data=plotDt[ERRTYPE=='talyspar_endep' & EXPID %in% varnames_all[19:28],])
ggp <- ggplot(data=plotDt[ERRTYPE=='talyspar_endep' & EXPID %in% varnames_all[29:36],])
ggp <- ggp + theme_bw()
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   plot.title=element_text(size=12))
ggp <- ggp + xlab('energy') + ylab('parameter value relative to default')
ggp <- ggp + geom_ribbon(aes(x=EN, ymin=DATAMIN, ymax=DATAMAX), alpha=0.3)
ggp <- ggp + geom_line(aes(x=EN, y=DATA))
ggp <- ggp + facet_wrap(~ EXPID, ncol=2)
ggp <- ggp + ylim(0.5,1.5)
#ggp <- ggp + annotate("label", x = 10, y = 1, label = "avg rate")
ggp
