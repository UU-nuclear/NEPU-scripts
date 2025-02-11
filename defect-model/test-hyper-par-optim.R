
rm(list=ls())
gc()

source("config/config-Fe56.R")
source("defect-model/defect-model.R")

expDt <- read_object(3, "expDt")
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc-Fe56"

# test on reduced energy range
setorder(expDt,REAC,L1)
expDt <- expDt[L1<3]
# expDt <- expDt[L1<1.3]
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
energies <- expDt[,sort(L1)]
# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
energies <- unique(round(energies*1e4))*1e-04

model <- defect_model(energies, talys_calc_dir, expDt)

pars <- model$parsDt[,V1]

# and map the experimental energies
predictions <- model$jac %*% pars # this will be the default TALYS prediction

expDt[,default_prediction:=as.vector(predictions)]
expDt[,RESIDUAL:=DATA-default_prediction]

# now try to create a GP prior covariance matrix on the parameters
# of the model

opt_fun <- function(hyper_pars) {

	cov_func <- function(x1,x2) {
		sigma2*exp(-0.5*(x1 - x2)^2/l2)
	}

	log_likelihood <- function(observed, sigma_obs) {
		# observed = vector of experimental data points, or residual
		# sigma_obs = prior mapped to experiments + experimental covariance
		chi_square <- observed %*% solve(sigma_obs, observed)

		tmp <- determinant(sigma_obs)
		stopifnot(tmp$sign == 1)
		log_det_cov <- tmp$modulus

		as.vector(-0.5*(length(observed)*log(2*pi) + log_det_cov + chi_square))
	}

	sigma2 <- (hyper_pars[1])^2
	l2 <- (hyper_pars[2])^2

	priors <- list()
	for(reac in model$parsDt[,unique(REAC)]) {
		model_energies <- model$parsDt[REAC==reac,L1]
		model_default <- model$parsDt[REAC==reac,V1]
		cov_mat <- outer(model_energies,model_energies,cov_func)

		scale <- outer(model_default,model_default)

		cov_mat <- Matrix(cov_mat*scale, sparse = TRUE)

		priors <- append(priors,list(cov_mat))
	}

	full_prior_cov_mat <- bdiag(priors)
	prior_exp <- model$jac %*% full_prior_cov_mat %*% t(model$jac)

	log_likelihood(expDt[,RESIDUAL], prior_exp + Diagonal(x=expDt[,UNC^2]))
}

guess_hyper_pars <- c(0.1,3e-1)
min_hyper_pars <- c(-1,1e-3)
max_hyper_pars <- c(1,1e-2)

opt_fun(guess_hyper_pars)

optim_control <- list(fnscale=-1, ndeps=c(1e-03,1e-05))
optim_res <- optim(guess_hyper_pars, fn=opt_fun, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)

sigma2 <- (optim_res$par[1])^2
	l2 <- (optim_res$par[2])^2

cov_func <- function(x1,x2) {
	sigma2*exp(-0.5*(x1 - x2)^2/l2)
}

priors <- list()
for(reac in model$parsDt[,unique(REAC)]) {
	model_energies <- model$parsDt[REAC==reac,L1]
	model_default <- model$parsDt[REAC==reac,V1]
	cov_mat <- outer(model_energies,model_energies,cov_func)

	scale <- outer(model_default,model_default)

	cov_mat <- Matrix(cov_mat*scale, sparse = TRUE)

	priors <- append(priors,list(cov_mat))
}

full_prior_cov_mat <- bdiag(priors)
prior_exp <- model$jac %*% full_prior_cov_mat %*% t(model$jac)

Sigma_12 <- full_prior_cov_mat %*% t(model$jac)

pars_mean <- Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC^2]),expDt[,RESIDUAL])
model$parsDt[,CONDITIONAL:=as.vector(pars_mean)]

pars_cov <- full_prior_cov_mat - Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC^2]),t(Sigma_12))
model$parsDt[,CONDITIONAL_UNC:=sqrt(diag(pars_cov))]

residual_predictions <- model$jac %*% pars_mean
predictions_cov <- model$jac %*% pars_cov %*% t(model$jac)

expDt[,RESIDUAL_GP:=as.vector(residual_predictions)]
expDt[,RESIDUAL_GP_UNC:=sqrt(diag(predictions_cov))]

#################################################################################

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
		geom_line(aes(x=L1,y=default_prediction),col='green') +
		#geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp


filepath <- file.path(getwd(), 'defect-model', 'gp-hyper-par-optim-square-exp.png')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

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
		#geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

filepath <- file.path(getwd(), 'defect-model', 'gp-hyper-par-optim-sqr-exp.png')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)
################################################################################

ggp_model <- ggplot(data=model$parsDt) +
		geom_line(aes(x=L1,y=V1 + CONDITIONAL),col='red') +
		geom_line(aes(x=L1,y=V1),col='green', linetype=2) +
		geom_ribbon(aes(x=L1,ymin=V1+CONDITIONAL-CONDITIONAL_UNC, ymax=V1+CONDITIONAL+CONDITIONAL_UNC), fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp_model

################################################

# prior_exp <- model$jac %*% full_prior_cov_mat %*% t(model$jac)
# pars_cov <- full_prior_cov_mat - Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC^2]),t(Sigma_12))
# should be the same as mult_xt_invCov_x(x=      D=      S=      P=      cholZ=)
pars_cov <- full_prior_cov_mat - mult_xt_invCov_x(x=t(Sigma_12), D=Diagonal(x=expDt[,UNC^2]), S=model$jac, P=full_prior_cov_mat)
# but gives
# Error in h(simpleError(msg, call)) : 
#  error in evaluating the argument 'x' in selecting a method for function 'forceSymmetric':
#  A = LU is computationally singular: min(d)/max(d) = 4.809e-38, d = abs(dia​‌​
#  g(U))

# do it in steps

D <- Diagonal(x=expDt[,UNC^2])
invDx <- solve(D, t(Sigma_12))

# cholZ <- makeCholZ(D, S, P)
# Z <- forceSymmetric(crossprod(S, solve(D, S)) + solve(P))
Z1 <- crossprod(model$jac, solve(D, model$jac))
#Z2 <- solve(full_prior_cov_mat) # computationally singular
Z2 <- solve(full_prior_cov_mat + Diagonal(n=nrow(full_prior_cov_mat), x=1e-04)) # adding a small nugget parameter
# Z <- forceSymmetric(crossprod(model$jac, solve(D, model$jac)) + solve(full_prior_cov_mat))
Z <- forceSymmetric(Z1 + Z2)
# makeCholesky(Z)
cholZ <- Cholesky(Z, LDL = FALSE, perm = TRUE)

pars_cov <- full_prior_cov_mat - crossprod(t(Sigma_12), (invDx - solve(D, S %*% solve(cholZ, crossprod(model$jac, invDx)))))