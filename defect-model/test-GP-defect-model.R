
rm(list=ls())
gc()

#library(microbenchmark, lib.loc="R-user-libs/")

source("config/config-Fe56.R")
source("defect-model/defect-model.R")
source("defect-model/GP-defect-model.R")

expDt <- read_object(3, "expDt")
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc-Fe56"

# test on reduced energy range
setorder(expDt,REAC,L1)
#expDt <- expDt[L1<1.3]
#expDt <- expDt[L1<3]
expDt <- expDt[L1<5.5 & L1>5]
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
# energies <- expDt[,sort(L1)]

# get energies only from the total xs
energies <- expDt[grepl("\\(N,TOT\\)",REAC),sort(L1)]

# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
energies <- unique(round(energies*1e4))*1e-04

model <- defect_model(energies, talys_calc_dir, expDt)

pars <- model$parsDt[,V1]

# and map the experimental energies
predictions <- model$jac %*% pars # this will be the default TALYS prediction

expDt[,default_prediction:=as.vector(predictions)]
expDt[,RESIDUAL:=DATA-default_prediction]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA, col=EXPID)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction), col='red') +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp
################################################

guessSigma <- 0.25404736
guessLen <- 0.01271152
min_hyper_pars <- c(0,1e-4)
max_hyper_pars <- c(1,1)

hyper_pars_init <- c(guessSigma,guessLen)

#myGPmodel <- GP_model(hyper_pars_init, model, expDt[,RESIDUAL], Diagonal(x=expDt[,UNC^2]))
GP_prior_mean <-  expDt[,RESIDUAL]
myGPmodel <- GP_model(hyper_pars_init, model, GP_prior_mean, Diagonal(x=expDt[,UNC^2]))

#####################

optim_control <- list(fnscale=-1)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, gr=myGPmodel$grad_ML, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)

#####################
library(optimParallel)
nCores <- detectCores(all.tests = FALSE, logical = TRUE)
cl <- makeCluster(6)
setDefaultCluster(cl=cl)

dummy <- clusterEvalQ(cl, c(library(data.table)))
print(dummy)
print("------")
print(cl)
print("------")
clusterExport(cl, c("myGPmodel","model","GP_prior_mean"), 
              envir=environment())

optim_res <- optimParallel(par = hyper_pars_init, 
                          fn = myGPmodel$marginal_likelihood, 
                          gr = myGPmodel$grad_ML, 
                          method = "L-BFGS-B",
                          lower = min_hyper_pars, 
                          upper = max_hyper_pars, 
                          control = list(fnscale = -1)
                          )

stopCluster(cl)
# for energy range L1<3>=: optim_pars <- c(0.06474224, 0.01260597)
####################

optim_pars <- optim_res$par
pars_cov_mat <- myGPmodel$get_pars_cov_mat(optim_pars)
exp_cov_mat <- myGPmodel$get_exp_cov_mat(optim_pars)

Sigma_12 <- pars_cov_mat %*% t(model$jac)

pars_mean <- Sigma_12 %*% solve(exp_cov_mat,expDt[,RESIDUAL])
model$parsDt[,CONDITIONAL:=as.vector(pars_mean)]

#pars_cov <- pars_cov_mat - Sigma_12 %*% solve(exp_cov_mat,t(Sigma_12))
pars_cov <- as.matrix(pars_cov_mat) - as.matrix(Sigma_12 %*% solve(exp_cov_mat,t(Sigma_12))) # smaller in memory and faster to compute

model$parsDt[,CONDITIONAL_UNC:=sqrt(diag(pars_cov))]

residual_predictions <- model$jac %*% pars_mean
predictions_cov <- model$jac %*% pars_cov %*% t(model$jac)

# this is a way to compute the covariance matrix at the experiments without first computing the
# the very large model covariance matrix, this is much faster and saves memory (if it is not already computed)!
# prior_exp <- model$jac %*% pars_cov_mat %*% t(model$jac)
# predictions_cov <- prior_exp - mult_xt_invCov_x(prior_exp, D=Diagonal(x=expDt[,UNC^2]), S=model$jac, P=pars_cov_mat)

expDt[,RESIDUAL_GP:=as.vector(residual_predictions)]
expDt[,RESIDUAL_GP_UNC:=sqrt(diag(predictions_cov))]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA, col=EXPID)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
		geom_line(aes(x=L1,y=default_prediction),col='green') +
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
		#xlim(1,1.5) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

#################################################################################

new_energy_grid = energies
fake_expDt <- data.table(
	REAC=rep(expDt[,unique(REAC)], each=length(new_energy_grid)),
	L1=rep(new_energy_grid, times=length(expDt[,unique(REAC)]))
	)

fake_expDt[,IDX:=seq_len(.N)]

smooth_model <- defect_model(energies, talys_calc_dir, fake_expDt)

# get the default model parameters (talys prediction)
pars <- smooth_model$parsDt[,V1]
# and map the experimental energies
predictions <- smooth_model$jac %*% pars
fake_expDt[,default_prediction:=as.vector(predictions)]

smooth_predictions <- smooth_model$jac %*% pars_mean
smooth_predictions_cov <- smooth_model$jac %*% pars_cov %*% t(smooth_model$jac)

fake_expDt[,RESIDUAL_GP:=as.vector(smooth_predictions)]
fake_expDt[,RESIDUAL_GP_UNC:=sqrt(diag(smooth_predictions_cov))]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='green') +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt,col='red') +
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

####################
L <- myGPmodel$fun(hyper_pars_init)
# test the gradient
delta <- 1e-08
hyper_pars_mod <- c(guessSigma+delta,guessLen)
deltaL <- myGPmodel$fun(hyper_pars_mod) - L
deltaL/delta

hyper_pars_mod <- c(guessSigma,guessLen+delta)
deltaL <- myGPmodel$fun(hyper_pars_mod) - L
deltaL/delta

grad_finite_diff <- function(hyper_pars, delta=1e-08) {
	L <- myGPmodel$fun(hyper_pars)

	hyper_pars_mod <- c(guessSigma+delta,guessLen)
	deltaL <- myGPmodel$fun(hyper_pars_mod) - L
	dL_dsigma <- deltaL / delta

	hyper_pars_mod <- c(guessSigma,guessLen+delta)
	deltaL <- myGPmodel$fun(hyper_pars_mod) - L
	dL_dl <- deltaL / delta

	matrix(c(dL_dsigma,dL_dl), ncol=1)
}

grad_finite_diff(hyper_pars_init)

#######################################################
microbenchmark(grad_finite_diff(hyper_pars_init), times = 10)
microbenchmark(myGPmodel$grad(hyper_pars_init))

# the analytical gradient is only slightly faster ~23%, may be due to that I'm evaluating it at the same parameter set each time
sigmas <- rnorm(10, mean=guessSigma, sd=0.01*guessSigma)
lengths <- rnorm(10, mean=guessLen, sd=0.01*guessLen)

tot_time <- 0
for(i in seq_along(sigmas)) {
	tot_time <- tot_time + system.time(grad_finite_diff(c(sigmas[i],lengths[i])))
}
print(tot_time)

tot_time <- 0
for(i in seq_along(sigmas)) {
	tot_time <- tot_time + system.time(myGPmodel$grad(c(sigmas[i],lengths[i])))
}
print(tot_time)
####################

hyper_pars <- hyper_pars_init
observed <- expDt[,RESIDUAL]
observed_cov_mat <- Diagonal(x=expDt[,UNC^2])
nugget <- 1e-04

##############################

# here I go parameter by parameter and propagate a unit change to the experimental
# data. If a change in the parameter has no efect on the prediction at the experimental
# energies sum(as.vector(model$jac %*% par_vector)) will be zero
insensitive <- c()
for(i in seq_len(nrow(model$parsDt))) {
	par_vector <- rep(0, nrow(model$parsDt))
	par_vector[i] <- 1
	sensitive <- sum(as.vector(model$jac %*% par_vector))
	if(sensitive==0) {
		cat("parameter nbr ", i, "is not sensitive\n")
		insensitive <- c(insensitive,i)
	}
}

model$parsDt[insensitive]
# these parameters won't be affected by the conditioning on the experimental data,
# and could in principle be removed from the calculations. Their final value may be
# affected by the GP covariance though.

# Maybe this is not correct!!!
# In principle all parameters will contribute to at least the total cross section

# What is happening is that I have an energy grid point in the model for each
# experimental data point. This means that there may be an experimental data point
# for example in (n,inl) whcih falls in between energy grid points for the total xs.
# then changing for (n,a) at this energy won't affect directly the prediction of the
# total xs (at the experimental energies) since they are determined by interpolation
# of other energy points only, it will still affect the total du to the correlation
# introduced by the GP covariance function.

# I guess I need to be more careful in choosing the energy grid points, I don't need to
# have at every experimental energy for the exclusive channels 