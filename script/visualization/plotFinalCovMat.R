#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(ggplot2)
library(latex2exp)

library(reshape2)
library(gplots)
library(RColorBrewer)


##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

allResults <- read_object(12, 'allResults')
needsDt <- read_object(1, "needsDt")
expDt <- read_object(3, "expDt")

# create a data.table with the model energy-grid
# and the reaction channel specification

reactions <- expDt[,unique(REAC)]

extNeedsDt <- needsDt[,{
    stopifnot(all(L2 == 0) & all(L3 == 0))
    list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
         L2 = 0, L3 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]
extNeedsDt[, IDX := seq_len(.N)]

# create the interpolation matrix using the exforHandler
modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles) # fake exfor sub-entries
modDt <- exforHandler$extractData(modSubents, ret.values=FALSE) # fake experimental data
Smod <- exforHandler$getJac(modDt, extNeedsDt, modSubents)
setkey(modDt, IDX, DIDX)

modDt[,L2:=NULL]
modDt[,L3:=NULL]

# reoder the data so that the REAC keyword is correct
reorder <- function(x)
{
  as.vector(Smod %*% x)
}
allResults <- apply(allResults,2,reorder)


#construct the covariance matrix from the result of the sampling in step09
uncinfo <- cov.wt(t(allResults),cor=TRUE)

modDt$FIT_MEAN <- uncinfo$center

model_covariance <- uncinfo$cov
model_correlation <- uncinfo$cor

#color_palette <- colorRampPalette(brewer.pal(8, "Blues"))(25)
#heatmap.2(model_covariance,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

covariance_data <- melt(model_covariance)
names(covariance_data)[3] <- "cov"
covariance_data$E1 <- modDt[covariance_data$Var1]$L1
covariance_data$E2 <- modDt[covariance_data$Var2]$L1
covariance_data$REAC1 <- modDt[covariance_data$Var1]$REAC
covariance_data$REAC2 <- modDt[covariance_data$Var2]$REAC

correlation_data <- melt(model_correlation)
covariance_data$corr <- correlation_data$value

#replace NaNs in the correlation with NAs
covariance_data[is.nan(covariance_data$corr),]$corr <- NA

# create column holding the bin widths for plotting tiles
# refParamDt <- read_object(2,"refParamDt")
# energyGrid <- refParamDt[PARNAME=="energy",unlist(PARVAL)]
energyGrid <- energyGridrandomFiles
binC <- energyGrid
binLow <- (c(0,binC[1:length(binC)-1])+binC)*0.5
binUp <- (c(binC[2:length(binC)],2*binC[length(binC)]-binC[length(binC)-1]) + binC)*0.5

covariance_data$binLow1 <- binLow[ifelse(covariance_data$Var1%%length(binC)==0,length(binC),covariance_data$Var1%%length(binC))]
covariance_data$binLow2 <- binLow[ifelse(covariance_data$Var2%%length(binC)==0,length(binC),covariance_data$Var2%%length(binC))]

covariance_data$binUp1 <- binUp[ifelse(covariance_data$Var1%%length(binC)==0,length(binC),covariance_data$Var1%%length(binC))]
covariance_data$binUp2 <- binUp[ifelse(covariance_data$Var2%%length(binC)==0,length(binC),covariance_data$Var2%%length(binC))]

#cov_n_inl <- covariance_data[covariance_data$REAC1=="(26-FE-56(N,INL)26-FE-56,,SIG)" & covariance_data$REAC2=="(26-FE-56(N,INL)26-FE-56,,SIG)",]
#cov_n_p <- covariance_data[covariance_data$REAC1=="(26-FE-56(N,P)25-MN-56,,SIG)" & covariance_data$REAC2=="(26-FE-56(N,P)25-MN-56,,SIG)",]

#cov_np_plot <- ggplot(cov_n_p, aes(x=E2, y=E1, xmin=binLow2, xmax=binUp2, ymin=binLow1, ymax=binUp1))
#cov_np_plot <- cov_np_plot + geom_rect(aes(fill=value)) 
#cov_np_plot <- cov_np_plot + scale_fill_gradient2(midpoint=0.5*(max(cov_n_p$value)+min(cov_n_p$value)))
#
#corr_np_plot <- ggplot(cov_n_p, aes(x=E2, y=E1, xmin=binLow2, xmax=binUp2, ymin=binLow1, ymax=binUp1))
#corr_np_plot <- corr_np_plot + coord_cartesian(xlim = c(2.65, max(binUp)),ylim = c(2.65, max(binUp)),expand=FALSE)
#corr_np_plot <- corr_np_plot + geom_rect(aes(fill=corr)) 
#corr_np_plot <- corr_np_plot + scale_fill_gradient2(midpoint=0.5)
#corr_np_plot

reactions <- unique(covariance_data$REAC1)

path <- file.path(plotPath, 'correlation_matrices')
dir.create(path, recursive=TRUE, showWarnings=FALSE)
for(reaction in reactions) {
  cov_reac <- covariance_data[covariance_data$REAC1==reaction & covariance_data$REAC2==reaction,]
  cov_reac <- cov_reac[!is.na(cov_reac$corr),]
  corr_plot <- ggplot(cov_reac, aes(x=E2, y=E1, xmin=binLow2, xmax=binUp2, ymin=binLow1, ymax=binUp1))
  corr_plot <- corr_plot + coord_cartesian(xlim = c(min(cov_reac$binLow1), max(cov_reac$binUp1)),ylim = c(min(cov_reac$binLow1), max(cov_reac$binUp1)),expand=FALSE)
  corr_plot <- corr_plot + geom_rect(aes(fill=corr))
  #corr_plot <- corr_plot + scale_fill_gradient2(midpoint=0.5*(min(cov_reac$corr)+max(cov_reac$corr)))
  corr_plot <- corr_plot + scale_fill_gradient2(breaks=seq(from=-1,to=1,by=0.5),limits=c(-1.00001,1.000001))

  file_name <- file.path(path, paste0('correlation_', reaction,'.png'))
  ggsave(file_name, corr_plot, width = 8.65, height = 5.6, units = "cm", dpi = 300)

  # saving a cross section plot just to check that the reaction mapping is correct
  xs_plot <- ggplot(modDt[REAC==reaction],aes(x=L1,y=FIT_MEAN)) +
  geom_line()

  file_name <- file.path(path, paste0('fit_mean_', reaction,'.png'))
  #ggsave(file_name, xs_plot, width = 8.65, height = 5.6, units = "cm", dpi = 300)
}