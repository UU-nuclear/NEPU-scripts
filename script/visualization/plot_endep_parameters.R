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

ggp <- ggplot(data=plotDt_adjustables[ERRTYPE=='talyspar_endep'])
ggp <- ggp + theme_bw()
ggp <- ggp + theme(axis.text=element_text(size=10),
                   axis.title=element_text(size=10),
                   plot.title=element_text(size=10))
ggp <- ggp + xlab('energy') + ylab('parameter value relative to default')
ggp <- ggp + geom_ribbon(aes(x=EN, ymin=DATAMIN, ymax=DATAMAX), alpha=0.3)
ggp <- ggp + geom_line(aes(x=EN, y=DATA))
ggp <- ggp + facet_wrap(~ EXPID, ncol=2)
ggp <- ggp + ylim(0.5,1.5)
ggp

print(ggp)
dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
filepath <- file.path(plotPath, 'endep_parameters.png')
#ggsave(filepath, ggp, width = 8.65, height = 5.6, units = "cm", dpi = 300)
ggsave(filepath, ggp)

# IT APPEARS THAT SOME OF THE ENERGY DEPENDENT PARAMTERS ARE NOT ADJUST OVER THE ENTIRE ENERGY INTERVALL
# DURING THE LM FIT. IN PARTICULAR THE POINT AT E=0 IS NOT PART OF THE ADJUSTABLE PARAMETERS. THIS HAS 
# ONE OBVIOUS EXPLANATION, THERE IS NO DATA INCLUDED AT BELOW 2 MeV, SINCE PARAMETERS ARE SELECTED BASED
# ON THEIR SENSITIVITY TO EXPERIMENTAL DATA THE POINT AT E=0 WILL ALWAYS BE EXCLUDED. CURIOSLY, ALTHOUGH
# THE PARAMETER (AT E=0) IS NOT INCLUDED IN THE OPTIMIZATION IT GETS A VALUE != PRIOR. I GUESS THIS HAS TO DO
# WITH THE GAUSSIAN PROCESS LENGTH SCALE. 