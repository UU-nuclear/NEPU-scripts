# This is a script doing what is done on page 10 of https://cran.r-project.org/web/packages/hetGP/vignettes/hetGP_vignette.pdf

library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)

hom <- mleHomGP(mcycle$times, mcycle$accel, covtype = "Matern5_2")
het <- mleHetGP(mcycle$times, mcycle$accel, covtype = "Matern5_2")

Xgrid <- matrix(seq(0, 60, length = 301), ncol = 1)

p_hom <- as.data.table(predict(x = Xgrid, object = hom))
p_hom[,times:=Xgrid]

p_het <- as.data.table(predict(x = Xgrid, object = het))
p_het[,times:=Xgrid]


ggp <- ggplot(data=mcycle)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_point(aes(x=times, y=accel), size=1.0)
# plot the models

# homoschedastic
ggp <- ggp + geom_line(aes(x=times,y=mean),data=p_hom, col='red')
ggp <- ggp + geom_ribbon(aes(x=times, ymin=mean-sqrt(sd2), ymax=mean+sqrt(sd2)),data=p_hom, col='red', alpha=0., linetype = "dashed")

# heteroschedastic
ggp <- ggp + geom_line(aes(x=times,y=mean),data=p_het, col='blue')
ggp <- ggp + geom_ribbon(aes(x=times, ymin=mean-sqrt(sd2), ymax=mean+sqrt(sd2)),data=p_het, col='blue', alpha=0., linetype = "dashed")
ggp

################3

ggp <- ggplot(data=mcycle)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_point(aes(x=times, y=accel), size=1.0)
# plot the models

# homoschedastic
ggp <- ggp + geom_line(aes(x=times,y=mean),data=p_hom, col='red')
ggp <- ggp + geom_ribbon(aes(x=times, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=p_hom, col='red', alpha=0., linetype = "dashed")

# heteroschedastic
ggp <- ggp + geom_line(aes(x=times,y=mean),data=p_het, col='blue')
ggp <- ggp + geom_ribbon(aes(x=times, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=p_het, col='blue', alpha=0., linetype = "dashed")
ggp <- ggp + geom_ribbon(aes(x=times, ymin=qnorm(0.05,mean,sqrt(nugs)), ymax=qnorm(0.95,mean,sqrt(nugs))),data=p_het, col='green', alpha=0., linetype = "dashed")
ggp