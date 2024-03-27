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
#       OUTPUT FROM PREVIOUS STEPS
##################################################

allParsets <- read_object(12, "allParsets")
allParamDt <- read_object(12, "allParamDt")
optParamDt <- read_object(10, "optParamDt")

par_energies <- optParamDt[,str_extract(PARNAME,"\\(.+\\)")]
par_energies <- str_remove(par_energies,"\\(")
par_energies <- str_remove(par_energies,"\\)")
par_energies <- as.numeric(par_energies)
optParamDt[,PAR_EN:=par_energies]
optParamDt[,FULLPARNAME:=PARNAME]
optParamDt[,PARNAME:=str_remove(PARNAME,"\\(.+\\)")]

par_energies <- allParamDt[,str_extract(PARNAME,"\\(.+\\)")]
par_energies <- str_remove(par_energies,"\\(")
par_energies <- str_remove(par_energies,"\\)")
par_energies <- as.numeric(par_energies)
allParamDt[,PAR_EN:=par_energies]
allParamDt[,FULLPARNAME:=PARNAME]
allParamDt[,PARNAME:=str_remove(PARNAME,"\\(.+\\)")]

allParamDt_adj <- allParamDt[ADJUSTABLE==TRUE]
allParamDt_adj[,IDX:=seq_len(.N)]

#construct the covariance matrix from the result of the sampling in step09
uncinfo <- cov.wt(t(allParsets),cor=TRUE)

par_covariance <- uncinfo$cov
par_correlation <- uncinfo$cor

colnames(par_correlation) <- allParamDt_adj[,FULLPARNAME]
rownames(par_correlation) <- allParamDt_adj[,FULLPARNAME]

# first order the indices in the data.table
setorder(allParamDt_adj,PARNAME,PAR_EN)

# select by parameter name and get the indices in the correct order
# pars_to_plot <- c("rwdadjust n","rvadjust n")
#pars_to_plot <- c("avadjust d","rvadjust d")
#pars_to_plot <- c("rwsoadjust n")

# all energy dependent parameter : this matrix is too large to plot
# pars_to_plot <- unique(allParamDt_adj[grepl("\\(.+\\)",FULLPARNAME),PARNAME]) 

# all energy dependent parameters that were adjusted in the LM
pars_to_plot <- unique(optParamDt[ADJUSTABLE==TRUE &!is.na(PAR_EN),PARNAME])

# all energy dependent parameters that were NOT adjusted in the LM
#pars_to_plot <- unique(optParamDt[ADJUSTABLE==FALSE &!is.na(PAR_EN),PARNAME])

pars_to_plot <- pars_to_plot[1:14] # only a subset

# plot the correlation matrix
#color_palette <- colorRampPalette(brewer.pal(8, "Blues"))(25)
idx_to_plot <- allParamDt_adj[PARNAME %in% pars_to_plot,IDX]

# select only those that are outside the original range (that go trough the LM algorithm)
idx_to_plot <- allParamDt_adj[PARNAME %in% pars_to_plot & PAR_EN>max(energyGridForParams),IDX]

ncolors <- 256
color_palette<-colorRampPalette(c("red","white","blue"))(ncolors)
heatmap.2(par_correlation[idx_to_plot,idx_to_plot],
  Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
  margins=c(8,8),trace = "none",dendrogram="none",density.info = "none",
  breaks=seq(-1.,1.,length.out = ncolors+1))

