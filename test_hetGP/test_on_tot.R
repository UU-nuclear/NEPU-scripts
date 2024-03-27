
library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)


expDt <- read_object(3,"expDt")

# take one of the data sets
#curExpDt <- expDt[EXPID=="10047029"] # 240
curExpDt <- expDt[EXPID=="13840002"] # 7130 points
#curExpDt <- expDt[EXPID=="22131005"] # 510 points
#curExpDt <- expDt[EXPID=="22549005"] # single point
#curExpDt <- expDt[EXPID=="40149008"] # single point

#curExpDt <- curExpDt[L1>3] # just to limit the number of points a bit for testing

Ec <- 0.01 # 10 keV
# data with energy spacing below Ec are considered as resonance data
# data with energy spacing larger than Ec are considered as average cross-section
setorder(curExpDt,L1)
energies <- curExpDt$L1
nn <- length(energies)
dE <- energies[2:nn]-energies[1:nn-1]
dE <- c(dE[1],dE) # the distance for the first point is to the next point
curExpDt[,dE := dE]

curExpDt <- curExpDt[dE<Ec]

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
#ggp <- ggp + geom_point(aes(x=L1, y=DATA), size=1.0)
#ggp <- ggp + geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC))

# construct a grid on which I want the cross sections
points_to_use <- which(energyGridrandomFiles<max(curExpDt$L1))
points_to_use <- c(points_to_use,max(points_to_use)+1)
Xgrid <- matrix(energyGridrandomFiles[points_to_use], ncol=1)

# fit the homoschedastic model
#hom <- mleHomGP(curExpDt$L1,curExpDt$DATA, covtype = "Matern5_2")
#p_hom <- as.data.table(predict(x = Xgrid, object = hom))
#p_hom[,L1:=Xgrid]

# plot the homoschedastic model
#ggp <- ggp + geom_line(aes(x=L1,y=mean),data=p_hom, col='blue')
#ggp <- ggp + geom_ribbon(aes(x=L1, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=p_hom, col='blue', alpha=0., linetype = "dashed")

# fit the heteroschedastic model
settings <- list(return.hom = TRUE)
het <- mleHetGP(curExpDt$L1,curExpDt$DATA, covtype = "Matern5_2", settings = settings)
p_het <- as.data.table(predict(x = Xgrid, object = het))
p_het[,L1:=Xgrid]

# plot the heteroschedastic model
ggp <- ggp + geom_line(aes(x=L1,y=mean),data=p_het, col='red')
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=p_het, col='red', alpha=0., linetype = "dashed")
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=qnorm(0.05,mean,sqrt(nugs)), ymax=qnorm(0.95,mean,sqrt(nugs))),data=p_het, col='green', alpha=0., linetype = "dashed")
ggp

#ggp <- ggp + geom_line(aes(x=L1,y=V1),data=modDt[REAC=="(24-CR-52(N,INL)24-CR-52,,SIG)" & L1<30], col='cyan')

ggp + xlim(c(1.5,NA))

Xgrid <- matrix(seq(from=min(curExpDt$L1)-0.1,to=max(curExpDt$L1) + 0.1, by=0.01), ncol=1)
p_het <- as.data.table(predict(x = Xgrid, object = het))
p_het[,L1:=Xgrid]

Xgrid <- matrix(seq(from=2,to=3, by=0.001), ncol=1)
p_het <- as.data.table(predict(x = Xgrid, object = het))
p_het[,L1:=Xgrid]