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

allParsets <- read_object(12, "allParsets")
allParamDt <- read_object(12, "allParamDt")


par_energies <- allParamDt[,str_extract(PARNAME,"\\(.+\\)")]
par_energies <- str_remove(par_energies,"\\(")
par_energies <- str_remove(par_energies,"\\)")
par_energies <- as.numeric(par_energies)
allParamDt[,PAR_EN:=par_energies]
allParamDt[,PARNAME:=str_remove(PARNAME,"\\(.+\\)")]

#construct the covariance matrix from the result of the sampling in step09
uncinfo <- cov.wt(t(allParsets),cor=TRUE)

par_covariance <- uncinfo$cov
par_correlation <- uncinfo$cor

#color_palette <- colorRampPalette(brewer.pal(8, "Blues"))(25)
#heatmap.2(par_covariance,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

covariance_data <- data.table(melt(par_covariance))
names(covariance_data)[3] <- "cov"
covariance_data[,PARNAME1:=allParamDt[ADJUSTABLE==TRUE,PARNAME][Var1]]
covariance_data[,PARNAME2:=allParamDt[ADJUSTABLE==TRUE,PARNAME][Var2]]
covariance_data[,E1:=allParamDt[ADJUSTABLE==TRUE,PAR_EN][Var1]]
covariance_data[,E2:=allParamDt[ADJUSTABLE==TRUE,PAR_EN][Var2]]

correlation_data <- melt(par_correlation)
covariance_data$corr <- correlation_data$value

###################
energyGrid <- covariance_data[!is.na(E1),unique(E1)]
binC <- energyGrid
bins <- c(0,(binC[2:length(binC)] + binC[1:(length(binC)-1)])*0.5,2*binC[length(binC)]-binC[length(binC)-1])

covariance_data[,binLow1:= bins[.bincode(E1,breaks=bins)]]
covariance_data[,binLow2:= bins[.bincode(E2,breaks=bins)]]
covariance_data[,binUp1:= bins[.bincode(E1,breaks=bins)+1]]
covariance_data[,binUp2:= bins[.bincode(E2,breaks=bins)+1]]
#plot_dt <- covariance_data[!is.na(E1) & !is.na(E2)]
plot_dt <- covariance_data[(PARNAME1=="rvadjust n" | PARNAME1=="v1adjust n") & (PARNAME2=="rvadjust n" | PARNAME2=="v1adjust n") ]


corr_plot <- ggplot(plot_dt, aes(x=E2, y=E1, xmin=binLow2, xmax=binUp2, ymin=binLow1, ymax=binUp1))
corr_plot <- corr_plot + coord_cartesian(xlim = c(min(plot_dt$binLow1), max(plot_dt$binUp1)),ylim = c(min(plot_dt$binLow1), max(plot_dt$binUp1)),expand=FALSE)
corr_plot <- corr_plot + geom_rect(aes(fill=corr))
corr_plot <- corr_plot + scale_fill_gradient2(breaks=seq(from=0,to=1,by=0.25))
corr_plot
