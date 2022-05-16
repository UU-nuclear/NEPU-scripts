# plot the posterior distribution of energy-dependent TALYS parameters

source("config.R")
library(ggplot2)

finalPars <- read_object(8, "finalPars")
finalParCovmat <- read_object(8, "finalParCovmat")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")

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
#    ENERGY-DEPENDENT PARAMETERS
###########################################

#varnames <- plotDt[ERRTYPE=='talyspar_endep', unique(EXPID)]
# sensitive energy-dependent parameters
varnames <- c('TALYS-v1adjust a', 'TALYS-vso1adjust p',
              'TALYS-d1adjust p', 'TALYS-v1adjust p', 'TALYS-vso1adjust n',
              'TALYS-w1adjust n', 'TALYS-d1adjust n', 'TALYS-v1adjust n')
              #'TALYS-rcadjust p')
varnames <- sub('TALYS-','',varnames)

# varname <- 'TALYS-v1adjust n'
ggp <- ggplot(data=plotDt[ERRTYPE=='talyspar_endep' & EXPID %in% varnames,])
ggp <- ggp + theme_bw()
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   plot.title=element_text(size=12))
ggp <- ggp + xlab('energy') + ylab('parameter value relative to default')
ggp <- ggp + geom_ribbon(aes(x=EN, ymin=DATAMIN, ymax=DATAMAX), alpha=0.3)
ggp <- ggp + geom_line(aes(x=EN, y=DATA))
ggp <- ggp + facet_wrap(~ EXPID, ncol=2)
ggp


dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(plotPath, 'plot_posterior_endep_pars.png'), ggp,
       units='cm', width=8.66, height=12)


###########################################
#    ENERGY-INDEPENDENT PARAMETERS
###########################################

plotDt2 <- plotDt[ERRTYPE=='talyspar']
plotDt2 <- plotDt2[order(abs(DATAMAX-DATAMIN))]
plotDt2[, IDX:=seq_len(.N)]
setkey(plotDt2, IDX)
plotDt2$EXPID <- factor(plotDt2$EXPID, levels=rev(plotDt2$EXPID), ordered=TRUE)

ggp <- ggplot(data=plotDt2[1:30])
ggp <- ggp + theme_bw() #+ theme(axis.text.x = element_text(angle=90))
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   axis.title.x=element_text(hjust=1),
                   plot.title=element_text(size=12))
ggp <- ggp + xlab('parameter value relative to default') + ylab('')
ggp <- ggp + geom_errorbarh(aes(y=EXPID, xmin=DATAMIN, xmax=DATAMAX), size=0.5, height=0.3)
ggp <- ggp + geom_point(aes(y=EXPID, x=DATA), col='red', size=0.75)
ggp

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(plotPath, 'plot_posterior_pars.png'), ggp,
       units='cm', width=8.65, height=12, dpi=300)
