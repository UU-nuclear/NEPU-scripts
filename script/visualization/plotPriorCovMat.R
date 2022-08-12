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
P0 <- read_object(7, "P0")
optParamDt <- read_object(7, "optParamDt")

# convert from sparse matrix to regular matrix format
P0mat <- as.matrix(P0)

#set the names
ParNames <- optParamDt[ADJUSTABLE==TRUE]$PARNAME
colnames(P0mat) <- ParNames
rownames(P0mat) <- ParNames

endep_pars <- grepl("\\(.+\\)",ParNames)
other_pars <- !grepl("\\(.+\\)",ParNames)

P0_endep <- P0mat[endep_pars,endep_pars]
P0_other <- P0mat[other_pars,other_pars]

# sort by:
# (1) the last character (projectile)
# (2) first 3 characters in the name
# (3) energies of the parameters
energies <- as.numeric(substr(rownames(P0_endep),regexpr("\\(",rownames(P0_endep))+1,regexpr("\\)",rownames(P0_endep))-1))
ordering <- order(substr(rownames(P0_endep),nchar(rownames(P0_endep))-1,nchar(rownames(P0_endep))),substr(rownames(P0_endep),1,3),energies,decreasing=TRUE)
P0_endep <- P0_endep[ordering,ordering]

color_palette <- colorRampPalette(brewer.pal(8, "Blues"))(25)
#heatmap(P0_endep,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none")
heatmap.2(P0_endep,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

png(file=file.path(plotPath, 'PriorCov_endep_pars.png'),
width=1200, height=700)
heatmap.2(P0_endep,Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
dev.off()

# plot a specific paramter
par2plot <- "d1adjust\\(.+\\) n"

selection <- grepl(par2plot,rownames(P0_endep))
#heatmap(P0_endep[selection,selection],Rowv=NA,Colv=NA,symm=TRUE,col=color_palette)
#heatmap.2(P0_endep[selection,selection],Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")
