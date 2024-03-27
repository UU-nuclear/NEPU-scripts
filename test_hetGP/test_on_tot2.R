
library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)

# a script to test the idea of inflating the random uncertainty from resonance like data
# 
# We are trying to fit the average (smooth) cross section. The data in the low energy region
# is however not representative of this mean. They follow a different distribution, due to the resonance-
# like structure. Therefore, their weight in the fit should be representative of that these are samples
# from this distribution. If we do not do this the propagated uncertainty will be an underestimation of
# the true uncertainty.

# The simple way of doing this would be to calculate the standard deviation
# from a model (GLS fit) in each energy interval and inflate the uncertainty of the data points to this,
# this is what Helgesson did in his thesis.
# Here is a different way of achieving this by using a heteroschedastic gaussian process
# 1) First we "train" a GP to fit the default talys calculation, giving us a representative length scale of
#    the mean cross section.
# 2) We then fit a het. GP to the data using this length scale. The het. GP estimates the variance of the 
#    data around a mean (with the length scale extracted from talys)
# 3) This can then be used to replace the uncertainty of the experimental data points, to be represtative
#    of the distribution around the mean.
# So what we are doing is to estimate, using the GP, the distribution of data around a representative mean.
# We do not move the data points, but only inflate their uncertainty by considering them as samples from
# the determined distribution. This will then reflect in a larger uncertainty in the mean cross section once
# we fit the talys model


expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

# take one of the data sets
#curExpDt <- expDt[EXPID=="10047029"] # 240
curExpDt <- expDt[EXPID=="13840002"] # 7130 points
#curExpDt <- expDt[EXPID=="22131005"] # 510 points
#curExpDt <- expDt[EXPID=="22549005"] # single point
#curExpDt <- expDt[EXPID=="40149008"] # single point

curModDt <- modDt[REAC=="(24-CR-52(N,TOT),,SIG)"]

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + geom_line(aes(x=L1, y=DATA),data=curModDt,col='red')

# train a GP on the talys-model data

X <- matrix(curModDt$L1, ncol = 1)
Z <- curModDt$DATA
nvar <- 1

settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training

model_hetGP <- mleHetGP(X = X, Z = Z,
covtype = "Matern5_2", settings = settings)

# and draw it
Xgrid <- matrix(seq(from=min(curModDt$L1),to=max(curModDt$L1),by=0.01), ncol=1)
pred_model_het <- as.data.table(predict(x = Xgrid, object = model_hetGP))
pred_model_het[,L1:=Xgrid]
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=pred_model_het, col='green', alpha=0., linetype = "dashed")


# now try to model the data, with the length scale fixed to the value obtained from the fit to the talys-model

# limit to the region of resonance structure
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

# limit the energy range further for the test to not take too long
#curExpDt <- curExpDt[L1<3]

X <- matrix(curExpDt$L1, ncol = 1)
Z <- curExpDt$DATA
nvar <- 1

settings <- list(return.hom = TRUE, ) # To keep homoskedastic model used for training

model_hetGP_exp <- mleHetGP(X = X, Z = Z,
covtype = "Matern5_2", settings = settings, known = list(theta=model_hetGP$theta))

Xgrid <- matrix(seq(from=min(curExpDt$L1),to=max(curExpDt$L1),by=0.01), ncol=1)
pred_exp <- as.data.table(predict(x = Xgrid, object = model_hetGP_exp))
pred_exp[,L1:=Xgrid]

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp <- ggp + geom_line(aes(x=L1, y=DATA),data=curModDt[L1<=5.1],col='red')
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=pred_exp, col='green', alpha=0., linetype = "dashed")
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=mean-sqrt(nugs), ymax=mean+sqrt(nugs)),data=pred_exp, col='blue', alpha=0., linetype = "dashed")
ggp <- ggp + geom_line(aes(x=L1, y=mean),data=pred_exp, col='blue')

# for comparison, estimate the standard deviation in each energy bin

energyGrid <- curModDt[L1>=min(curExpDt$L1) & L1<max(curExpDt$L1),L1]

x_grid <- c()
stddev <- c()
y_mean <- c()
mean_unc <- c()
for(i in seq_len(length(energyGrid)-1)) {
    E1 <- energyGrid[i]
    E2 <- energyGrid[i+1]

    x_grid <- c(x_grid,0.5*(E1+E2))
    stddev <- c(stddev,sd(curExpDt[L1>=E1 & L1<E2,DATA]))
    y_mean <- c(y_mean,mean(curExpDt[L1>=E1 & L1<E2,DATA]))
    mean_unc <- c(mean_unc,mean(curExpDt[L1>=E1 & L1<E2,DATA]))
}

simple <- data.table(x_grid,y_mean,stddev)

ggp <- ggp + geom_ribbon(aes(x=x_grid, ymin=y_mean-stddev, ymax=y_mean+stddev),data=simple, col='cyan', alpha=0., linetype = "dashed")
ggp <- ggp + geom_line(aes(x=x_grid, y=y_mean),data=simple, col='cyan')

# It seems clear that the nugget parameter in the heteroschedasticGP
# models the variance of the data with respect to the mean prediction
# we can use it to inflate the uncertainty on the data points as estimates
# for the mean cross section.

# To account for the statistical uncertainty in the data we can take gaussian noise
# samples for each data point before fitting the gp, this should be equvilant in 
# performance if we first group the data, since the samples from a given point
# constitutes what is refered to as repetitions

######################33
# X <- matrix(curExpDt$L1, ncol = 1)
# Z <- curExpDt$DATA
# nvar <- 1
# 
# # this step groups the data when there are repetitions
# # we could use this to account for the stat. error
# # in the data taking several samples for each data point
# # then group them like this to account for the distribution
# data_m <- find_reps(X, Z, rescale = TRUE)
# 
# plot(rep(data_m$X0, data_m$mult), data_m$Z, #ylim = c(-160, 90),
# ylab = 'cross section (mb)', xlab = "En (MeV)")
# 
# 
# ggp <- ggplot(data=curExpDt)
# ggp <- ggp + theme_bw() + theme(legend.position="none")
# ggp <- ggp + theme(axis.text=element_text(size=4),
                   # axis.title=element_text(size=4),
                   # strip.text=element_text(size=3))
# ggp <- ggp + geom_point(aes(x=L1, y=DATA), size=1.0)
# 
# # construct a grid on which I want the cross sections
# points_to_use <- which(energyGridrandomFiles<max(curExpDt$L1))
# points_to_use <- c(points_to_use,max(points_to_use)+1)
# Xgrid <- matrix(energyGridrandomFiles[points_to_use], ncol=1)
# 
# ## Model fitting
# settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training
# 
# model <- mleHetGP(X = matrix(curExpDt$L1, ncol = 1), Z = curExpDt$DATA, lower = rep(0.1, nvar), upper = rep(50, nvar),
# covtype = "Matern5_2", settings = settings)
# ## A quick view of the fit
# summary(model)