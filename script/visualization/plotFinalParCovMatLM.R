#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}


library(reshape2)
library(gplots)
library(RColorBrewer)


##################################################
#       READ IN RESULTS NEEDED
##################################################
optRes <- read_object(10,"optRes")
optParamDt <- read_object(10, "optParamDt")
FinalParCovMatLM <- as.matrix(optRes$parCovLM)

##################################################
#       SEPPARATE ENERGY DEPENTENT AND OTHER PARS
##################################################

#set the names
ParNames <- optParamDt[ADJUSTABLE==TRUE]$PARNAME


colnames(FinalParCovMatLM) <- ParNames
rownames(FinalParCovMatLM) <- ParNames

endep_pars <- grepl("\\(.+\\)",ParNames)
other_pars <- !grepl("\\(.+\\)",ParNames)

Cov_endep <- FinalParCovMatLM[endep_pars,endep_pars]
Cov_other <- FinalParCovMatLM[other_pars,other_pars]


##################################################
#       PLOT ENERGY DEPENTENT PARS
##################################################

# sort by:
# (1) the last character (projectile)
# (2) first 3 characters in the name
# (3) energies of the parameters
energies <- as.numeric(substr(rownames(Cov_endep),regexpr("\\(",rownames(Cov_endep))+1,regexpr("\\)",rownames(Cov_endep))-1))
ordering <- order(substr(rownames(Cov_endep),nchar(rownames(Cov_endep))-1,nchar(rownames(Cov_endep))),substr(rownames(Cov_endep),1,3),energies,decreasing=TRUE)
Cov_endep <- Cov_endep[ordering,ordering]

# set parameter names only for the central energy
central_energy <- energyGridForParams[0.5*length(energyGridForParams)]

endep_par_names <- colnames(Cov_endep)
central_par_names <- grepl(paste0("\\(",central_energy,"\\)"),endep_par_names)
endep_par_names[central_par_names] <- str_remove(endep_par_names[central_par_names],"\\(.+\\)")
#endep_par_names[!central_par_names] <- ""
endep_par_names <- str_remove(endep_par_names,"\\(.+\\)")

colnames(Cov_endep) <- endep_par_names
rownames(Cov_endep) <- endep_par_names

#color_palette <- colorRampPalette(brewer.pal(8, "Blues"))(25)
color_palette<-colorRampPalette(c("red","white","blue"))
heatmap.2(Cov_endep,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

png(file=file.path(plotPath, 'FinalCov_endep_pars_LM.png'),
width=1200, height=700)
heatmap.2(Cov_endep,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
dev.off()

heatmap.2(Cov_endep,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

png(file=file.path(plotPath, 'FinalCorr_endep_pars_LM.png'),
width=1200, height=700)
heatmap.2(cov2cor(Cov_endep),Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),
  trace = "none",dendrogram="none",density.info = "none",
  cexRow=1,cexCol=1)
dev.off()

##################################################
#       PLOT ENERGY INDEPENTENT PARS
##################################################

heatmap.2(Cov_other,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

png(file=file.path(plotPath, 'FinalCov_energy_indep_pars_LM.png'),
width=1200, height=700)
heatmap.2(Cov_other,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
dev.off()

png(file=file.path(plotPath, 'FinalCorr_energy_indep_pars_LM.png'),
width=1200, height=700)
heatmap.2(cov2cor(Cov_other),Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
dev.off()


##################################################
#       PLOT THE FULL PAR COV MATRIX
##################################################

# sort by:
# (0) wether or not the parameter is energy dependent
# (1) the last character (projectile)
# (2) first 3 characters in the name
# (3) energies of the parameters
energies <- as.numeric(substr(rownames(FinalParCovMatLM),regexpr("\\(",rownames(FinalParCovMatLM))+1,regexpr("\\)",rownames(FinalParCovMatLM))-1))

energy_dependence <- !is.na(energies)
projectile_order <- substr(rownames(FinalParCovMatLM),nchar(rownames(FinalParCovMatLM))-1,nchar(rownames(FinalParCovMatLM)))
alphabetic_order <- substr(rownames(FinalParCovMatLM),1,3)

ordering <- order(energy_dependence,projectile_order,alphabetic_order,energies,decreasing=TRUE)
FinalParCovMatLM <- FinalParCovMatLM[ordering,ordering]

heatmap.2(FinalParCovMatLM,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

png(file=file.path(plotPath, 'FinalParCovLM_full.png'),
width=1200, height=700)
heatmap.2(FinalParCovMatLM,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
dev.off()

png(file=file.path(plotPath, 'FinalParCorrLM_full.png'),
width=1200, height=700)
heatmap.2(cov2cor(FinalParCovMatLM),Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
dev.off()