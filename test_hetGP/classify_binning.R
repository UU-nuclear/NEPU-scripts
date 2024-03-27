library(class)

expDt <- read_object(3,"expDt")

reac <- "(24-CR-52(N,TOT),,SIG)"
experiment <- "13840002"

curExpDt <- expDt[EXPID==experiment]

# calculate the binning
setorder(curExpDt,L1)
energies <- curExpDt$L1
nn <- length(energies)
dE <- energies[2:nn]-energies[1:nn-1]
dE <- c(dE[1],dE) # the distance for the first point is to the next point
curExpDt[,dE := dE]

################################

km_res <- kmeans(curExpDt[,c("L1","dE")],3)
#km_res <- kmeans(curExpDt[,c("L1","dE")]
curExpDt[,kmclust:=km_res$cluster]

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_point(aes(x=L1, y=DATA, col=kmclust))
ggp

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_point(aes(x=L1, y=dE, col=kmclust))
ggp


###################################

library("bcp",lib.loc="/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/R-user-libs/")


bcp.1a <- bcp(curExpDt[,dE],w0=0.9,p0=1.e-3)
plot(bcp.1a, main="Univariate Change Point Example")

bcp.1b <- bcp(curExpDt[,L1])
bcp.1a <- bcp(curExpDt[,dE])
plot(bcp.1b, main="Univariate Change Point Example")

###################################

library(zoo,lib.loc="/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/R-user-libs/")
library(changepoint,lib.loc="/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/R-user-libs/")
# change in mean
ansmean=cpt.mean(curExpDt[,dE], method = 'BinSeg')
plot(ansmean,cpt.col='blue')
print(ansmean)

#####################################

library(EnvCpt,lib.loc="/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/R-user-libs/")

fit_envcpt = envcpt(curExpDt[,dE])  # Fit all models at once
fit_envcpt$summary  # Show log-likelihoods
plot(fit_envcpt)

fit_envcpt$meancpt