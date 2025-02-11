
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
maxE <- 20
#expDt <- expDt[L1<1.3]
#expDt <- expDt[L1<3]
expDt <- expDt[L1<maxE & L1>minE]
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
# energies <- expDt[,sort(L1)]

# get energies only from the total xs
# energies <- expDt[grepl("\\(N,TOT\\)",REAC),sort(L1)]

energies <- expDt[EXPID==22316003,sort(L1)]
energies <- c(energies,expDt[grepl("\\(N,TOT\\)",REAC)][L1>max(energies),L1])
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
hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out =30)
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
    # right now the mean of the latent GP is sqrt(0.1): should be set to some estimate based on the data
    # chi_square <- (sigmas - sqrt(0.1)) %*% solve(sigmas_prior_cov_mat, (sigmas-sqrt(0.1)))
    chi_square <- (sigmas - sqrt(sigmas_mean)) %*% solve(sigmas_prior_cov_mat, (sigmas-sqrt(sigmas_mean)))
    
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
  
  # the mean of the GP on the sigma hyper-paramter is set to the estimate of the
  # standard deviation of the residual/default prediction, we use the more robust
  # estimate of the standard deviation obtained from 1.4826*MAD(x)
  #sigmas_mean <- 1.4826*mad(observed/(model$jac %*% model$parsDt[,V1]))
  #sigmas_mean <- 0.1
  sigmas_mean <- 0
  
  n_calls <- 0
  
  list(marginal_likelihood=log_likelihood, get_pars_cov_mat=get_pars_cov_mat, get_exp_cov_mat=get_exp_cov_mat)
}
# ===============================================

#hyper_pars <- c(sigmas,len)

# $par
# [1] 0.288531580 0.140337542 0.137897181 0.097659496 0.083594204 0.088440044 0.088136669 0.103576047 0.256055112 0.266043531 0.257830776
# [12] 0.107828862 0.111414733 0.143444990 0.126338756 0.123767139 0.121427454 0.106233422 0.115467110 0.234023670 0.107638283 0.098285070
# [23] 0.108452172 0.248376020 0.235658072 0.229228686 0.241255261 0.231557943 0.229092749 0.224189138 0.009627443 0.871492305 1.929426168

# X = 0.2828, 0.143628, 0.153552, 0.113152, 0.112824, 0.0926551, 0.081706, 0.0810937, 0.239879, 0.252941, 0.246246, 0.0836372, 0.0959295, 0.141797, 0.119223, 0.118403, 0.115396, 0.0970989, 0.109205, 0.226819, 0.100474, 0.0935237, 0.1056, 0.264707, 0.230831, 0.225586, 0.25267, 0.236932, 0.231585, 0.222977, 0.00969007, 0.867319, 1.9446 

#guessSigma <- rep(0.2,length(sigmas))
guessSigma <- seq(from=0.2, to=0.1, length.out=length(hyper_par_energy_grid))
# guessSigma <- c(0.9,0.671,0.3,0.1)
guessLen <- min(diff(energies))*3

guess_latent_sigma <- 1
guess_latent_len <- 1.25
hyper_pars_init <- c(guessSigma, guessLen, guess_latent_sigma, guess_latent_len)
# hyper_pars_init <- c(0.88333830 0.65310938 0.29088514 0.09821135 0.01513433)
# hyper_pars_init <- c(0.11656545, 0.15463725, 0.09103543, 0.11805153, 0.11115000, 0.08886865, 0.03450786, 0.02291360)
hyper_pars_init <- c(0.288531580, 0.140337542, 0.137897181, 0.097659496, 0.083594204, 0.088440044, 0.088136669, 0.103576047, 0.256055112, 0.266043531, 0.257830776,
                     0.107828862, 0.111414733, 0.143444990, 0.126338756, 0.123767139, 0.121427454, 0.106233422, 0.115467110, 0.234023670, 0.107638283, 0.098285070,
                     0.108452172, 0.248376020, 0.235658072, 0.229228686, 0.241255261, 0.231557943, 0.229092749, 0.224189138, 0.009627443, 0.871492305, 1.929426168)
# It may be important to set proper box constraints on the length hyper-parameter
# This is because this affects the sparseness of the covariance matrix, and thereby
# the computation time (and memory use), if the optimizer starts to move out to 
# try very long length scales the calculation of the log-likelihood will slow down
# and the program may even crash if it runs out of memory.
min_hyper_pars <- c(rep(1e-03,length(guessSigma)),min(diff(energies)), 0.1, 0.1)
max_hyper_pars <- c(rep(1,length(guessSigma)),0.1, 10, 10)


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
# hyper_pars_init <- c(0.2828, 0.143628, 0.153552, 0.113152, 0.112824, 0.0926551, 0.081706, 0.0810937, 0.239879, 0.252941, 0.246246, 0.0836372, 0.0959295, 0.141797, 0.119223, 0.118403, 0.115396, 0.0970989, 0.109205, 0.226819, 0.100474, 0.0935237, 0.1056, 0.264707, 0.230831, 0.225586, 0.25267, 0.236932, 0.231585, 0.222977, 0.00969007, 0.867319, 1.9446)
# optim_control <- list(fnscale=-1)
# optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-3,length(hyper_pars_init)-1),1e-05))
optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-3,length(hyper_pars_init)-1),1e-05), maxit=30)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)

#################
# iterations 12
# function evaluations 50
# segments explored during Cauchy searches 44
# BFGS updates skipped 0
# active bounds at final generalized Cauchy point 0
# norm of the final projected gradient 7.25208
# final function value 23011.2
# 
# X = 0.250866 0.143847 0.144044 0.0973012 0.0939648 0.0881319 0.0918502 0.0875745 0.180621 0.204302 0.197009 0.0393718 0.0722271 0.134241 0.0995773 0.0998962 0.0991354 0.0643457 0.0922449 0.184842 0.087376 0.0691168 0.0942496 0.35707 0.199633 0.189411 0.347219 0.315791 0.310211 0.19764 0.0094844 0.783042 2.74792 
# F = 23011.2
# Warning:  more than 10 function and gradient evaluations
# in the last line search
# final  value 23011.237236 
# converged
# > optim_res$par
# [1] 0.250866313 0.143846708 0.144043533 0.097301171 0.093964803 0.088131928 0.091850239 0.087574525 0.180620814 0.204301888 0.197009000
# [12] 0.039371802 0.072227068 0.134241200 0.099577307 0.099896151 0.099135390 0.064345678 0.092244861 0.184841912 0.087376017 0.069116822
# [23] 0.094249576 0.357069985 0.199632663 0.189411409 0.347219040 0.315790735 0.310210948 0.197639724 0.009484396 0.783042054 2.747916364
#################

optim_pars <- optim_res$par
expDt[,PRIOR_GP_UNC:=sqrt(diag(myGPmodel$get_exp_cov_mat(hyper_pars_init)))]
ggp <- ggplot(data=expDt[L1<=max(energies)]) +
  geom_point(aes(x=L1,y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),col='red') +
  geom_ribbon(aes(x=L1,ymin=default_prediction-PRIOR_GP_UNC, ymax=default_prediction+PRIOR_GP_UNC), fill='red', alpha=0.3) +
  #		geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE),col='red') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'GP-prior.pdf')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)
#######################################################



optim_pars <- optim_res$par
pars_cov_mat <- myGPmodel$get_pars_cov_mat(optim_pars)
exp_cov_mat <- myGPmodel$get_exp_cov_mat(optim_pars) + Diagonal(x=expDt[,UNC^2])

Sigma_12 <- pars_cov_mat %*% t(model$jac)

# Going different routes with the mathematics I arrive at the same conditional
# covariance matrix
# Sigma_{beta|y} = (K_beta^{-1} + J^T Sigma_exp^{-1} J)^{-1} =  K_beta J^T (Sigma_exp + J K_beta J^T)^{-1} J K_beta
# but seemingly different conditional means
# (not sure, maybe I'm just missing something)
# originally I just use the conditional mean of the joint pdf, which gives
# mu_{beta|y} = K_beta J^T (Sigma_exp + J K_beta J^T)^{-1} y
# using the formulas presented in
# Lindholm et al., MACHINE LEARNING A First Course for Engineers and Scientists, Draft version: April 30, 2021
# theorem 9.4 and corollary 9.1 I get instead
# mu_{beta|y} = Sigma_{beta|y} J^T Sigma_exp^{-1} y = (K_beta^{-1} + J^T Sigma_exp^{-1} J)^{-1} J^T Sigma_exp^{-1} y 

# this mean vector is probably incorrect?!
pars_mean <- Sigma_12 %*% solve(exp_cov_mat,expDt[,RESIDUAL])
model$parsDt[,CONDITIONAL:=as.vector(pars_mean)]

#pars_cov <- pars_cov_mat - Sigma_12 %*% solve(exp_cov_mat,t(Sigma_12))
pars_cov <- as.matrix(pars_cov_mat) - as.matrix(Sigma_12 %*% solve(exp_cov_mat,t(Sigma_12))) # smaller in memory and faster to compute
# tmp_mat <- mult_xt_invCov_x(x=t(Sigma_12), D=Diagonal(x=expDt[,UNC^2]), S=model$jac, P=pars_cov_mat)
# pars_cov <- pars_cov_mat - tmp_mat

# this should be the correct one derived from Bayessian linear regression
#pars_mean_new <- pars_cov %*% t(model$jac) %*% solve(Diagonal(x=expDt[,UNC^2]), expDt[,RESIDUAL])
#model$parsDt[,CONDITIONAL_new:=as.vector(pars_mean_new)]

model$parsDt[,CONDITIONAL_UNC:=sqrt(diag(pars_cov))]

residual_predictions <- model$jac %*% pars_mean # this is the old, (maybe incorrect, but probably not) version
#residual_predictions_new <- model$jac %*% pars_mean_new # this is the new mean based on Bayessian linear regression
predictions_cov <- model$jac %*% pars_cov %*% t(model$jac)

# this is a way to compute the covariance matrix at the experiments without first computing the
# the very large model covariance matrix, this is much faster and saves memory (if it is not already computed)!
# prior_exp <- model$jac %*% pars_cov_mat %*% t(model$jac)
# predictions_cov <- prior_exp - mult_xt_invCov_x(prior_exp, D=Diagonal(x=expDt[,UNC^2]), S=model$jac, P=pars_cov_mat)

expDt[,RESIDUAL_GP:=as.vector(residual_predictions)]
#expDt[,RESIDUAL_GP_NEW:=as.vector(residual_predictions)]
expDt[,RESIDUAL_GP_UNC:=sqrt(diag(predictions_cov))]

ggp <- ggplot(data=expDt) +
  # geom_point(aes(x=L1,y=DATA, col=EXPID)) +
  geom_point(aes(x=L1,y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
  #geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP_NEW),col='cyan') +
  geom_line(aes(x=L1,y=default_prediction),col='green') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'exp-to-20MeV-new.png')
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

ggp <- ggplot(data=expDt[L1<5]) +
  geom_point(aes(x=L1,y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[L1<5],col='green') +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[L1<5],col='red') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[L1<5] , fill='red', alpha=0.3) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp


filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'mod-without-unc-matern-kernel-to-20MeV-new.pdf')
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