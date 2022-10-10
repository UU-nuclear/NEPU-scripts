#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}


library(ggplot2)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################
outdataPathRun <- outdataPath
plotPath <- paste0(outdataPathRun, "/plots")


finalPars <- read_object(08, "finalPars", outdata_path = outdataPathRun)
finalParCovmat <- read_object(08, "finalParCovmat", outdata_path = outdataPathRun)
optSysDt_allpars <- read_object(07, "optSysDt_allpars", outdata_path = outdataPathRun) # all parameters
optParamDt <- read_object(07,"optParamDt")
optGpDt <- read_object(06,"optGpDt")

plotDt <- copy(optSysDt_allpars)
setkey(plotDt, IDX)
plotDt$DATA <- paramTrafo$fun(finalPars)
plotDt$UNC <- sqrt(diag(finalParCovmat))
plotDt$DATAMIN <- paramTrafo$fun(finalPars - sqrt(diag(finalParCovmat)))
plotDt$DATAMAX <- paramTrafo$fun(finalPars + sqrt(diag(finalParCovmat)))
plotDt$EXPID <- sub('TALYS-','', plotDt$EXPID)


# select only the parameters that are adjusted in the LM algorithm
adjustable_par_names <- optParamDt[ADJUSTABLE==TRUE]$PARNAME
plotDt_adjustables <- plotDt[PARNAME %in% adjustable_par_names]
adjustable_EXPID <- unique(plotDt_adjustables$EXPID)
plotDt_adjustables <- plotDt[EXPID %in% adjustable_EXPID]

###########################################
#    HYPER-PARAMETER VALUES
###########################################
GPsigmas <- optGpDt[PARNAME=="sigma"]
GPsigmas$EXPID <- sub('TALYS-','',GPsigmas$EXPID)
GPlenths <- optGpDt[PARNAME=="len"]
GPlenths$EXPID <- sub('TALYS-','',GPlenths$EXPID)

for(exp_id in GPsigmas$EXPID) { # add the values to the EXPIDs so they will be printed in the facet wrap
    plotDt_adjustables[EXPID==exp_id]$EXPID <- paste0(exp_id," : \u03B4 = ",format(GPsigmas[EXPID==exp_id]$PARVAL,digits=3,nsmall=3)," : \u03BB = ",format(GPlenths[EXPID==exp_id]$PARVAL,digits=3,nsmall=3))
}


# sort the EXPIDs according to uncertainty
expid_order <- order(plotDt[,abs(DATAMAX-DATAMIN)])
plotDt$EXPID <- factor(plotDt$EXPID, levels=unique(plotDt$EXPID[expid_order]), ordered=TRUE)

###########################################
#    ENERGY-DEPENDENT PARAMETERS
###########################################

ggp1 <- ggplot(data=plotDt_adjustables[ERRTYPE=='talyspar_endep'])
ggp1 <- ggp1 + theme_bw()
ggp1 <- ggp1 + scale_x_continuous(breaks=seq(0,50,5))
ggp1 <- ggp1 + theme(axis.text=element_text(size=10),
                   axis.title=element_text(size=10),
                   plot.title=element_text(size=10))
ggp1 <- ggp1 + xlab('energy') + ylab('parameter value relative to default')
ggp1 <- ggp1 + geom_ribbon(aes(x=EN, ymin=DATAMIN, ymax=DATAMAX), alpha=0.3)
ggp1 <- ggp1 + geom_line(aes(x=EN, y=DATA))
ggp1 <- ggp1 + geom_point(data=plotDt_adjustables[ERRTYPE=='talyspar_endep' & PARNAME %in% adjustable_par_names],aes(x=EN, y=DATA))
ggp1 <- ggp1 + facet_wrap(~ EXPID, ncol=4)
ggp1 <- ggp1 + ylim(0.5,1.5)

#print(ggp1)
dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, 'endep_parameters.png')
#ggsave(filepath, ggp1, width = 8.65, height = 5.6, units = "cm", dpi = 300)
#ggsave(filepath, ggp1, width = 20, height = 21.0, units = "cm", dpi = 300)
ggsave(filepath, ggp1, width = 29.7, height = 21.0, units = "cm", dpi = 300)
#ggsave(filepath, ggp1)

# IT APPEARS THAT SOME OF THE ENERGY DEPENDENT PARAMTERS ARE NOT ADJUST OVER THE ENTIRE ENERGY INTERVALL
# DURING THE LM FIT. IN PARTICULAR THE POINT AT E=0 IS NOT PART OF THE ADJUSTABLE PARAMETERS. THIS HAS 
# ONE OBVIOUS EXPLANATION, THERE IS NO DATA INCLUDED AT BELOW 2 MeV, SINCE PARAMETERS ARE SELECTED BASED
# ON THEIR SENSITIVITY TO EXPERIMENTAL DATA THE POINT AT E=0 WILL ALWAYS BE EXCLUDED. CURIOSLY, ALTHOUGH
# THE PARAMETER (AT E=0) IS NOT INCLUDED IN THE OPTIMIZATION IT GETS A VALUE != PRIOR. I GUESS THIS HAS TO DO
# WITH THE GAUSSIAN PROCESS LENGTH SCALE. 

###########################################
#    ENERGY-INDEPENDENT PARAMETERS
###########################################

plotDt2 <- plotDt[ERRTYPE=='talyspar']
plotDt2 <- plotDt2[order(abs(DATAMAX-DATAMIN))]
plotDt2[, IDX:=seq_len(.N)]
setkey(plotDt2, IDX)
plotDt2$EXPID <- factor(plotDt2$EXPID, levels=rev(plotDt2$EXPID), ordered=TRUE)

# Add column showing the relative change = fitted val. / starting val.
plotDt2[,RELATIVE_CHANGE := DATA/REFDATA]
plotDt2[,RELATIVE_CHANGE_MIN := DATAMIN/REFDATA]
plotDt2[,RELATIVE_CHANGE_MAX := DATAMAX/REFDATA]

#ggp <- ggplot(data=plotDt2[1:30])
ggp <- ggplot(data=plotDt2[order(UNC)][RELATIVE_CHANGE!=1.0])
#ggp <- ggplot(data=plotDt2[UNC!=0.1 & DATA!=1.0]) # all parameters that have been adjusted
ggp <- ggp + theme_bw() #+ theme(axis.text.x = element_text(angle=90))
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   axis.title.x=element_text(hjust=1),
                   plot.title=element_text(size=12))
ggp <- ggp + xlab('parameter value relative to initial') + ylab('')
#ggp <- ggp + geom_errorbarh(aes(y=EXPID, xmin=DATAMIN, xmax=DATAMAX), size=0.5, height=0.3)
#ggp <- ggp + geom_point(aes(y=EXPID, x=DATA), col='red', size=0.75)
ggp <- ggp + geom_errorbarh(aes(y=EXPID, xmin=RELATIVE_CHANGE_MIN, xmax=RELATIVE_CHANGE_MAX), size=0.5, height=0.3)
ggp <- ggp + geom_point(aes(y=EXPID, x=RELATIVE_CHANGE), col='red', size=0.75)

dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
ggsave(file.path(plotPath, 'posterior_pars.png'), ggp,
       units='cm', width=8.65, height=12, dpi=300)