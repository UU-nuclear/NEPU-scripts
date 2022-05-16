#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
  source("./config/config.R")
  stop("No config file supplied, using default file config.R", call.=FALSE)
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(ggplot2)
library(stringr)
library(data.table)
library(Matrix)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################
dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)

optPars <- read_object(8, "finalPars")
finalParamDt <- read_object(8, "finalParamDt")
optParCovmat <- read_object(8, "finalParCovmat")
P0 <- read_object(8, "P0_all")
refParamDt <- read_object(2, "refParamDt")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")
optParamDt <- read_object(7, "optParamDt")

par_selection <- optParamDt[, ADJUSTABLE==TRUE] 

#paramDt <- finalParamDt[par_selection]
paramDt <- optParamDt[par_selection]

paramDt[, OPTPARVAL_ref := as.numeric(paramDt$PARVAL)]
paramDt[, OPTPARUNC_ref := as.numeric(as.vector(sqrt(diag(P0))))]
paramDt[, OPTPARVAL := optPars[par_selection]]
paramDt[, OPTPARUNC := sqrt(diag(optParCovmat))[par_selection]]


paramDt[, ENDEPPARTYPE := str_extract(paramDt$PARNAME, ".*(?=\\(.+?\\))")]
paramDt[!is.na(ENDEPPARTYPE), ENDEPPARTYPE := paste0(ENDEPPARTYPE," ",str_sub(PARNAME, start= -1))]
paramDt[ !is.na(ENDEPPARTYPE), ENERGY := as.numeric(sapply(PARNAME, function(x) gsub("[\\(\\)]", "", regmatches(x , gregexpr("\\(.*?\\)",x))[[1]])))]
paramDt[,IDX:=seq(1,.N)]


edep_partypes <- unique(paramDt[!is.na(ENDEPPARTYPE)]$ENDEPPARTYPE)
eindep_partypes <- unique(paramDt[is.na(ENDEPPARTYPE)]$PARNAME)

for (edep_par in edep_partypes) {
  print(edep_par)
  print(paramDt[ENDEPPARTYPE==edep_par, OPTPARUNC_ref])
}

for (edep_par in edep_partypes) {
  curparamDt <- paramDt[ENDEPPARTYPE == edep_par]
  print(curparamDt)
  setkey(curparamDt, ENERGY) 
  ggp_epar <- ggplot(curparamDt) + theme_bw() + ggtitle(edep_par)
  #ggp_epar <- ggplot(curparamDt) + theme_bw() + ylim(0.2, 1.8) + ggtitle(edep_par)
  ggp_epar <- ggp_epar + geom_point(aes(x = ENERGY, y = OPTPARVAL, shape ="opt"), size = 2)
  ggp_epar <- ggp_epar + geom_point(aes(x = ENERGY, y = OPTPARVAL_ref, shape ="prior"), size = 2)
  ggp_epar <- ggp_epar + geom_errorbar(aes(x = ENERGY, ymin = OPTPARVAL - OPTPARUNC, ymax = OPTPARVAL + OPTPARUNC, col = "opt",linetype = "opt"), size = 0.4)
  ggp_epar <- ggp_epar + geom_errorbar(aes(x = ENERGY, ymin = OPTPARVAL_ref - OPTPARUNC_ref, ymax = OPTPARVAL_ref + OPTPARUNC_ref, col = "prior",linetype = "prior"), size = 0.4)
  ggp_epar <- ggp_epar + annotate("text", x = Inf, y = Inf, label = "Preliminary", 
                                  hjust =1.2, vjust = 0.5, size=20, fontface="bold.italic", color="gray", alpha= 0.3, angle = 45)
  ggp_epar <- ggp_epar + xlab(expression("Incident Energy"~italic(E)~"(MeV)")) + 
    ylab(expression("Value"))
  
  ggp_epar <- ggp_epar + scale_color_manual(name="", values=c("prior"="gray20", "opt"="black")) 
  ggp_epar <- ggp_epar + scale_linetype_manual(name="", values=c("prior"="dashed", "opt"="solid")) 
  ggp_epar <- ggp_epar + scale_shape_manual(name="", values=c("prior"=8, "opt"=1))    
  print(ggp_epar)
  
  filepath <- file.path(plotPath, "fit_par_prior_comp") 
  filename <- paste0(edep_par, "_fit_par_prior_comp.pdf")
  if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(filepath, filename), ggp_epar, width = 36, height = 24, units = "cm", dpi = 300)  
  
  filepath <- file.path(plotPath, "fit_par_prior_comp_png") 
  filename <- paste0(edep_par, "_fit_par_prior_comp.png")
  if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(filepath, filename), ggp_epar, width = 36, height = 24, units = "cm", dpi = 300)  
}


curparamDt <- paramDt[is.na(ENDEPPARTYPE)]
ggp_par <- ggplot(curparamDt) + theme_bw()  + ggtitle("Energy independent parameters") 
#ggp_par <- ggp_par + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggp_par <- ggp_par + geom_point(aes(x = PARNAME, y = OPTPARVAL, shape ="opt"))
ggp_par <- ggp_par + geom_point(aes(x = PARNAME, y = OPTPARVAL_ref, shape ="prior"), size = 2)
ggp_par <- ggp_par + geom_errorbar(aes(x = PARNAME, ymin = OPTPARVAL - OPTPARUNC, ymax = OPTPARVAL + OPTPARUNC, col = "opt", linetype = "opt"),  size = 0.4)
ggp_par <- ggp_par + geom_errorbar(aes(x = PARNAME, ymin = OPTPARVAL_ref - OPTPARUNC_ref, ymax = OPTPARVAL_ref + OPTPARUNC_ref, col = "prior", linetype = "prior"), size = 0.4)
ggp_par <- ggp_par + annotate("text", x = Inf, y = Inf, label = "Preliminary", 
                              hjust =1.2, vjust = 0.5, size=20, fontface="bold.italic", color="gray", alpha= 0.3, angle = 45)
ggp_par <- ggp_par + xlab(expression("Parameter")) + 
  ylab(expression("Value"))
ggp_par <- ggp_par + coord_flip()

ggp_par <- ggp_par + scale_color_manual(name="", values=c("prior"="gray20", "opt"="black")) 
ggp_par <- ggp_par + scale_linetype_manual(name="", values=c("prior"="dashed", "opt"="solid")) 
ggp_par <- ggp_par + scale_shape_manual(name="", values=c("prior"=8, "opt"=1))    

filepath <- file.path(plotPath, "fit_par_prior_comp") 
filename <- paste0("e_indep_par", "_fit_par_prior_comp.pdf")
if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(filepath, filename), ggp_par, width = 36, height = 24, units = "cm", dpi = 300)  

filepath <- file.path(plotPath, "fit_par_prior_comp_png") 
filename <- paste0("e_indep_par", "_fit_par_prior_comp.png")
if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(filepath, filename), ggp_par, width = 36, height = 24, units = "cm", dpi = 300)  
print(ggp_par)


##############################################################################################
# Untransformed parameters paramTrafo$fun()
########################################################################

for (edep_par in edep_partypes) {
  curparamDt <- paramDt[ENDEPPARTYPE == edep_par]
  setkey(curparamDt, ENERGY) 
  ggp_epar <- ggplot(curparamDt) + theme_bw() + ggtitle(edep_par)
  #ggp_epar <- ggplot(curparamDt) + theme_bw() + ylim(0.2, 1.8) + ggtitle(edep_par)
  ggp_epar <- ggp_epar + geom_point(aes(x = ENERGY, y = paramTrafo$fun(OPTPARVAL), shape ="opt"), size = 2)
  ggp_epar <- ggp_epar + geom_point(aes(x = ENERGY, y = paramTrafo$fun(OPTPARVAL_ref), shape ="prior"), size = 2)
  ggp_epar <- ggp_epar + geom_errorbar(aes(x = ENERGY, ymin = paramTrafo$fun(OPTPARVAL - OPTPARUNC), ymax = paramTrafo$fun(OPTPARVAL + OPTPARUNC), col = "opt",linetype = "opt"), size = 0.4)
  ggp_epar <- ggp_epar + geom_errorbar(aes(x = ENERGY, ymin = paramTrafo$fun(OPTPARVAL_ref - OPTPARUNC_ref), ymax = paramTrafo$fun(OPTPARVAL_ref + OPTPARUNC_ref), col = "prior",linetype = "prior"), size = 0.4)
  ggp_epar <- ggp_epar + annotate("text", x = Inf, y = Inf, label = "Preliminary", 
                                  hjust =1.2, vjust = 0.5, size=20, fontface="bold.italic", color="gray", alpha= 0.3, angle = 45)
  ggp_epar <- ggp_epar + xlab(expression("Incident Energy"~italic(E)~"(MeV)")) + 
    ylab(expression("Value"))
  
  ggp_epar <- ggp_epar + scale_color_manual(name="", values=c("prior"="gray20", "opt"="black")) 
  ggp_epar <- ggp_epar + scale_linetype_manual(name="", values=c("prior"="dashed", "opt"="solid")) 
  ggp_epar <- ggp_epar + scale_shape_manual(name="", values=c("prior"=8, "opt"=1))    
  print(ggp_epar)
  
  filepath <- file.path(plotPath, "fit_par_prior_comp_trafo") 
  filename <- paste0(edep_par, "_fit_par_prior_comp_trafo.pdf")
  if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(filepath, filename), ggp_epar, width = 36, height = 24, units = "cm", dpi = 300)  
  
  filepath <- file.path(plotPath, "fit_par_prior_comp_trafo_png") 
  filename <- paste0(edep_par, "_fit_par_prior_comp_trafo.png")
  if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(filepath, filename), ggp_epar, width = 36, height = 24, units = "cm", dpi = 300)  
}


curparamDt <- paramDt[is.na(ENDEPPARTYPE)]
ggp_par <- ggplot(curparamDt) + theme_bw()  + ggtitle("Energy independent parameters") 

ggp_par <- ggp_par + geom_point(aes(x = PARNAME, y = paramTrafo$fun(OPTPARVAL), shape ="opt"))
#ggp_par <- ggp_par+ theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))
ggp_par <- ggp_par + geom_point(aes(x = PARNAME, y = paramTrafo$fun(OPTPARVAL_ref), shape ="prior"), size = 2)
ggp_par <- ggp_par + geom_errorbar(aes(x = PARNAME, ymin = paramTrafo$fun(OPTPARVAL - OPTPARUNC), ymax = paramTrafo$fun(OPTPARVAL + OPTPARUNC), col = "opt", linetype = "opt"),  size = 0.4)
ggp_par <- ggp_par + geom_errorbar(aes(x = PARNAME, ymin = paramTrafo$fun(OPTPARVAL_ref - OPTPARUNC_ref), ymax = paramTrafo$fun(OPTPARVAL_ref + OPTPARUNC_ref), col = "prior", linetype = "prior"), size = 0.4)
ggp_par <- ggp_par + annotate("text", x = Inf, y = Inf, label = "Preliminary", 
                              hjust =1.2, vjust = 0.5, size=20, fontface="bold.italic", color="gray", alpha= 0.3, angle = 45)
ggp_par <- ggp_par + xlab(expression("Parameter")) + 
  ylab(expression("Value"))
ggp_par <- ggp_par + coord_flip()

ggp_par <- ggp_par + scale_color_manual(name="", values=c("prior"="gray20", "opt"="black")) 
ggp_par <- ggp_par + scale_linetype_manual(name="", values=c("prior"="dashed", "opt"="solid")) 
ggp_par <- ggp_par + scale_shape_manual(name="", values=c("prior"=8, "opt"=1))    

filepath <- file.path(plotPath, "fit_par_prior_comp_trafo") 
filename <- paste0("e_indep_par", "_fit_par_prior_comp_trafo.pdf")
if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(filepath, filename), ggp_par, width = 36, height = 24, units = "cm", dpi = 300)  

filepath <- file.path(plotPath, "fit_par_prior_comp_trafo_png") 
filename <- paste0("e_indep_par", "_fit_par_prior_comp_trafo.png")
if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(filepath, filename), ggp_par, width = 72, height = 24, units = "cm", dpi = 300)  
print(ggp_par)