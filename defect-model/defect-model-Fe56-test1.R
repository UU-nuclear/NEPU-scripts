
rm(list=ls())
gc()

library(mvnfast)
#library(microbenchmark, lib.loc="R-user-libs/")

source("config/config-Fe56.R")
source("defect-model/defect-model.R")
#source("defect-model/GP-defect-model.R")

expDt <- read_object(3, "expDt")
expDt[,UNC:=ORIG_UNC] # UNC was replaced using hetGP
#talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc-Fe56"
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/talys-fit-fe56"

# test on reduced energy range
setorder(expDt,REAC,L1)

minE <- 3
maxE <- 5
#expDt <- expDt[L1<1.3]
#expDt <- expDt[L1<3]
expDt <- expDt[L1<maxE & L1>minE]
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
# energies <- expDt[,sort(L1)]

# get energies only from the total xs
# energies <- expDt[grepl("\\(N,TOT\\)",REAC),sort(L1)]

energies <- expDt[EXPID==22316003,sort(L1)]
# energies <- c(energies,expDt[grepl("\\(N,TOT\\)",REAC)][L1>max(energies),L1])
# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
energies <- unique(round(energies*1e4))*1e-04
#energies <- c(3, energies)

# take every second energy
last_energy <- energies[length(energies)] 
energies <- energies[rep(c(TRUE,FALSE),length.out=length(energies))]

energies[1] <- energies[1] - 1e-5
if(last_energy>energies[length(energies)]) energies <- c(energies,last_energy)
##########################################

##########################################

model <- defect_model(energies, talys_calc_dir, expDt)

pars <- model$parsDt[,V1]

# and map the experimental energies
predictions <- model$jac %*% pars # this will be the default TALYS prediction

expDt[,default_prediction:=as.vector(predictions)]
expDt[,RESIDUAL:=DATA-default_prediction]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction), col='red') +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

minE <- min(energies) - 1e-04
maxE <- max(energies) + 1e-04
#hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out =30)
hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out =8)
#sigmas <- rep(0.1, length(hyper_par_energy_grid))
#len <- 3*min(diff(energies))

# ===============================================

hetGP_model <- function(n_sigmas, hyper_pars, model, observed, observed_cov_mat, nugget=1e-04) {

	cov_func <- function(x1, x2, length_scale) {
		d_over_l <- abs(x1-x2)/length_scale
		ifelse(d_over_l < 1, ( (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)), 0)
		# OBS: The sigma hyper-parameter is not included here but outside the function
	}

	cov_func_alt <- function(x, length_scale) {
		d_over_l <- abs(outer(x,x,"-"))/length_scale

		idx <- which(d_over_l<1, arr.ind=TRUE)
		d_over_l <- d_over_l[idx]
		data <- (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)

		sparseMatrix(i=idx[,1], j=idx[,2], x=data)
	}

	cov_func_partial_sigma <- function(x, length_scale) {
		# note that this is independent of the value of sigma!!!
		cov_func_alt(x, length_scale)
	}

	matern_1_2 <- function(x1, x2, length_scale) {
		exp(-abs(x1-x2)/length_scale)
	}

	matern_3_2 <- function(x1, x2, length_scale) {
		(1+sqrt(3)*abs(x1-x2)/length_scale)*exp(-abs(x1-x2)/length_scale)
	}

	make_cov_mat <- function(hyper_pars, model, nugget=1e-04) {

		sigmas <- hyper_pars[1:n_sigmas]
		len <- hyper_pars[n_sigmas + 1]

		priors <- list()
		for(reac in model$parsDt[,unique(REAC)]) {

#			model_energies <- model$parsDt[REAC==reac,L1]
#			model_default <- model$parsDt[REAC==reac,V1]
#			#cov_mat <- cov_func_alt(model_energies, length_scale=len)
#			cov_mat <- outer(model_energies,model_energies,cov_func, length_scale=len)
#			#cov_mat <- outer(model_energies,model_energies,cov_func, length_scale=len)*outer(model_energies,model_energies,matern_1_2, length_scale=len/3.121378/log(2))
#			#cov_mat <- outer(model_energies,model_energies,cov_func, length_scale=len)*outer(model_energies,model_energies,matern_3_2, length_scale=len/3.121378/log(2))
#
#
#			# sigma_scale <- outer(sigmas[match(model_energies,hyper_par_energy_grid)],sigmas[match(model_energies,hyper_par_energy_grid)])
#			sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
#			sigma_scale <- outer(sigmas_tmp,sigmas_tmp)
#
#			scale <- outer(model_default,model_default)
#
#			cov_mat <- Matrix(cov_mat*sigma_scale*scale)
#			cov_mat <- cov_mat + Diagonal(n=nrow(cov_mat), x=nugget)


			###########
			model_energies <- model$parsDt[REAC==reac,L1]
			model_default <- model$parsDt[REAC==reac,V1]
			cov_mat <- cov_func_alt(model_energies, length_scale=len)
			# the matrix due to the correlation function is the same for each channel (at least each channel is a sub-set of the total),
			# so the construction could be sped up by not recalculating it for each channel

			# sigma_scale <- outer(sigmas[match(model_energies,hyper_par_energy_grid)],sigmas[match(model_energies,hyper_par_energy_grid)])
			sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
			sigma_scale <- outer(sigmas_tmp,sigmas_tmp)

			scale <- outer(model_default,model_default)

			cov_mat <- cov_mat*sigma_scale*scale
			cov_mat <- cov_mat + Diagonal(n=nrow(cov_mat), x=nugget)

			###########

			priors <- append(priors,list(cov_mat))
		}

		bdiag(priors)
	}

	update_cov_mat <- function(hyper_pars) {

		print("update_cov_mat!")
		# if the hyper-parameters have changed from the last call
		# re-build covariance matrix
		if(!all(hyper_pars==cur_hyper_pars)) {
			print("hyper-pars have changed!")
			cur_hyper_pars <<- hyper_pars
			cur_cov_mat <<- make_cov_mat(cur_hyper_pars, model)
			cur_cov_mat_exp <<-  model$jac %*% cur_cov_mat %*% t(model$jac)
		}

	}

	log_likelihood <- function(hyper_pars) {

		n_calls <<- n_calls + 1
		print(paste("------- log_likelihood!",n_calls))
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)
		print("cov mat has been updated!")

		chi_square <- observed %*% solve(cur_cov_mat_exp + observed_cov_mat, observed)

		tmp <- determinant(cur_cov_mat_exp + observed_cov_mat)
		stopifnot(tmp$sign == 1)
		log_det_cov <- tmp$modulus

		L <- as.vector(-0.5*(length(observed)*log(2*pi) + log_det_cov + chi_square))

		# GP prior on the sigmas vector
		
		# the amplitude of the GP on sigmas covariance function
		# sigma0 <- 1 # 0.01 correspond to a 10% relative error on the talys prediction
		# the GP prior on sigmas covariance matrix:
		sigmas_prior_cov_mat <- latent_cov_mat(cur_hyper_pars[length(cur_hyper_pars)-1], cur_hyper_pars[length(cur_hyper_pars)])
		
		tmp <- determinant(sigmas_prior_cov_mat)
		stopifnot(tmp$sign == 1)
		log_det_sigmas_prior_cov_mat <- tmp$modulus
		
		sigmas <- hyper_pars[1:n_sigmas]
		# right now the mean of teh latent GP is sqrt(0.1): should be set to some estimate based on the data
		chi_square <- (sigmas - sqrt(0.1)) %*% solve(sigmas_prior_cov_mat, (sigmas-sqrt(0.1)))

		P <- as.vector(-0.5*(length(sigmas)*log(2*pi) + log_det_sigmas_prior_cov_mat + chi_square))

		L + P
	}

	cov_mat_partial_sigma <- function(hyper_pars) {
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		model_energies <- model$parsDt[,L1]
		sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
		sigma_scale <- outer(sigmas_tmp,sigmas_tmp)

		cur_cov_mat / sigma_scale
	}

	cov_mat_partial_len <- function(hyper_pars) {
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		cov_func_d_len <- function(x, length_scale) {
			d_over_l <- abs(outer(x,x,"-"))/length_scale

			idx <- which(d_over_l<1, arr.ind=TRUE)
			d_over_l <- d_over_l[idx]
			data <- (4/3)*(pi*(1-d_over_l)*cos(pi*d_over_l) + sin(pi*d_over_l))*sin(pi*d_over_l)*(d_over_l/len)

			sparseMatrix(i=idx[,1], j=idx[,2], x=data)
		}

		sigmas <- hyper_pars[1:length(hyper_pars)-1]
		len <- hyper_pars[length(hyper_pars)]

		priors <- list()
		for(reac in model$parsDt[,unique(REAC)]) {


			###########
			model_energies <- model$parsDt[REAC==reac,L1]
			model_default <- model$parsDt[REAC==reac,V1]
			cov_mat_dl <- cov_func_d_len(model_energies, length_scale=len)
			# the matrix due to the correlation function is the same for each channel (at least each channel is a sub-set of the total),
			# so the construction could be sped up by not recalculating it for each channel

			# sigma_scale <- outer(sigmas[match(model_energies,hyper_par_energy_grid)],sigmas[match(model_energies,hyper_par_energy_grid)])
			sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
			sigma_scale <- outer(sigmas_tmp,sigmas_tmp)

			scale <- outer(model_default,model_default)

			cov_mat_dl <- cov_mat_dl*sigma_scale*scale

			###########

			priors <- append(priors,list(cov_mat))
		}

		bdiag(priors)
	}

	get_pars_cov_mat <- function(hyper_pars) {
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		return(cur_cov_mat)
	}

	get_exp_cov_mat <- function(hyper_pars) {
		# update the covariance matrix if hyper-parameters have changed
		update_cov_mat(hyper_pars)

		return(cur_cov_mat_exp)
	}
	
	latent_cov_mat <- function(sigma, len) {
	  Matrix(outer(hyper_par_energy_grid,hyper_par_energy_grid,cov_func, length_scale=len))*sigma
	}

	# member variables used by both log_likelihood() and grad_log_likelihood()
	cur_hyper_pars <- hyper_pars
	cur_cov_mat <- make_cov_mat(cur_hyper_pars, model)
	cur_cov_mat_exp <-  model$jac %*% cur_cov_mat %*% t(model$jac)

	n_calls <- 0

	list(marginal_likelihood=log_likelihood, get_pars_cov_mat=get_pars_cov_mat, get_exp_cov_mat=get_exp_cov_mat)
}
# ===============================================

#hyper_pars <- c(sigmas,len)

#guessSigma <- rep(0.2,length(sigmas))
guessSigma <- seq(from=0.2, to=0.1, length.out=length(hyper_par_energy_grid))
# guessSigma <- c(0.9,0.671,0.3,0.1)
guessLen <- min(diff(energies))*3

guess_latent_sigma <- 1
guess_latent_len <- 1.25
hyper_pars_init <- c(guessSigma, guessLen, guess_latent_sigma, guess_latent_len)
# hyper_pars_init <- c(0.88333830 0.65310938 0.29088514 0.09821135 0.01513433)
# hyper_pars_init <- c(0.11656545, 0.15463725, 0.09103543, 0.11805153, 0.11115000, 0.08886865, 0.03450786, 0.02291360)

min_hyper_pars <- c(rep(1e-03,length(guessSigma)),min(diff(energies)), 0.1, 0.1)
max_hyper_pars <- c(rep(1,length(guessSigma)),10, 10, 10)


GP_prior_mean <-  expDt[,RESIDUAL]
myGPmodel <- hetGP_model(length(guessSigma),hyper_pars_init, model, GP_prior_mean, Diagonal(x=expDt[,UNC^2]))

# graphical representation of the guess parameters
expDt[,PRIOR_GP_UNC:=sqrt(diag(myGPmodel$get_exp_cov_mat(hyper_pars_init)))]

sample <- rmvn(1, rep(0,nrow(expDt)), myGPmodel$get_exp_cov_mat(hyper_pars_init) + Diagonal(n=nrow(expDt), x=1e-04), ncores = 1, isChol = FALSE)
#sample <- rmvn(1, rep(0,nrow(expDt)), myGPmodel$get_exp_cov_mat(optim_pars) + Diagonal(n=nrow(expDt), x=1e-04), ncores = 1, isChol = FALSE)
expDt[,GP_SAMPLE:=sample[1,]]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction),col='green') +
		geom_ribbon(aes(x=L1,ymin=default_prediction-PRIOR_GP_UNC, ymax=default_prediction+PRIOR_GP_UNC), fill='green', alpha=0.3) +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE),col='red') +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp


# filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'GP-prior-MLO-result.png')
# ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

# if you restart
# hyper_pars_init <- optim_res$par
# optim_control <- list(fnscale=-1)
# optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-3,length(hyper_pars_init)-1),1e-05))
optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-3,length(hyper_pars_init)-1),1e-05), maxit=40)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)


optim_pars <- optim_res$par
expDt[,PRIOR_GP_UNC:=sqrt(diag(myGPmodel$get_exp_cov_mat(hyper_pars_init)))]
ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction),col='red') +
		geom_ribbon(aes(x=L1,ymin=default_prediction-PRIOR_GP_UNC, ymax=default_prediction+PRIOR_GP_UNC), fill='red', alpha=0.3) +
#		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE),col='red') +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp
#######################################################
# > optim_res
# $par
# [1] 0.11656545 0.15463725 0.09103543 0.11805153 0.11115000 0.08886865 0.03450786
# [8] 0.02291360
# 
# $value
# [1] -7930.037
# 
# $counts
# function gradient 
     # 104      104 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
########################################3
# > optim_res
# $par
# [1] 0.40572822 0.25303662 0.20922609 0.12341743 0.01016092
# 
# $value
# [1] -29833.91
# 
# $counts
# function gradient 
#       43       43 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"
#######################################################
# $par
#  [1] 0.06275961 0.02626605 0.05899257 0.12096815 0.04492133 0.64699198
#  [7] 0.43339615 0.04221537 0.25814369 0.03159313 0.04311991 0.06749023
# [13] 0.05530061 0.06056973 0.30075665 0.03826405 0.17396169 0.14634146
# [19] 0.17709988 0.26899723 0.00100000 0.41238771 0.54677206 0.26798097
# [25] 0.06845227 0.31307383 0.08757651 0.25251574 0.31307897 0.00100000
# [31] 0.64646463
# 
# $value
# [1] -20988.93
# 
# $counts
# function gradient 
#      114      114 
# 
# $convergence
# [1] 0
# 
# $message
# [1] "CONVERGENCE: REL_REDUCTION_OF_F <= FACTR*EPSMCH"


optim_pars <- optim_res$par
pars_cov_mat <- myGPmodel$get_pars_cov_mat(optim_pars)
exp_cov_mat <- myGPmodel$get_exp_cov_mat(optim_pars) + Diagonal(x=expDt[,UNC^2])

Sigma_12 <- pars_cov_mat %*% t(model$jac)

pars_mean <- Sigma_12 %*% solve(exp_cov_mat,expDt[,RESIDUAL])
model$parsDt[,CONDITIONAL:=as.vector(pars_mean)]

#pars_cov <- pars_cov_mat - Sigma_12 %*% solve(exp_cov_mat,t(Sigma_12))
pars_cov <- as.matrix(pars_cov_mat) - as.matrix(Sigma_12 %*% solve(exp_cov_mat,t(Sigma_12))) # smaller in memory and faster to compute
# tmp_mat <- mult_xt_invCov_x(x=t(Sigma_12), D=Diagonal(x=expDt[,UNC^2]), S=model$jac, P=pars_cov_mat)
# pars_cov <- pars_cov_mat - tmp_mat

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
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'exp-with-unc-matern-kernel.png')
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
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp


filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'mod-without-unc-matern-kernel-test.pdf')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

#########################

# Try to draw some samples from the GP posterior

# GPsamples <- rmvn(10, rep(0,nrow(fake_expDt)), smooth_predictions_cov, ncores = 1, isChol = FALSE)

parGPsamples <- rmvn(10, rep(0,nrow(model$parsDt)), pars_cov, ncores = 1, isChol = FALSE)

fake_expDt[,GP_SAMPLE1:=as.vector(smooth_model$jac %*% parGPsamples[1,])]
fake_expDt[,GP_SAMPLE2:=as.vector(smooth_model$jac %*% parGPsamples[2,])]
fake_expDt[,GP_SAMPLE3:=as.vector(smooth_model$jac %*% parGPsamples[3,])]
fake_expDt[,GP_SAMPLE4:=as.vector(smooth_model$jac %*% parGPsamples[4,])]
fake_expDt[,GP_SAMPLE5:=as.vector(smooth_model$jac %*% parGPsamples[5,])]
fake_expDt[,GP_SAMPLE6:=as.vector(smooth_model$jac %*% parGPsamples[6,])]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='green') +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE1+RESIDUAL_GP),data=fake_expDt,col='red') +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE2+RESIDUAL_GP),data=fake_expDt,col='red') +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE3+RESIDUAL_GP),data=fake_expDt,col='red') +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE4+RESIDUAL_GP),data=fake_expDt,col='red') +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE5+RESIDUAL_GP),data=fake_expDt,col='red') +
		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE6+RESIDUAL_GP),data=fake_expDt,col='red') +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp