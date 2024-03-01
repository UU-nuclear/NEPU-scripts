
rm(list=ls())
gc()

library(mvnfast)

source("config/config-Fe56.R")
source("defect-model/defect-model.R")

# prepare the experimental data
expDt <- read_object(3, "expDt")
expDt[,UNC:=ORIG_UNC] # UNC was replaced using hetGP: undo

# reduce the energy range
setorder(expDt,REAC,L1)

minE <- 3
maxE <- 18.985
expDt <- expDt[L1<maxE & L1>minE]
expDt[,IDX:=seq_len(.N)]

# retrieve systematic uncertainties in the format of nucdataBaynet
sysDt <- read_object(4, "updSysDt")
sysDt[, ADJUSTABLE:=NULL]

sysDt <- sysDt[grepl("EXPID-",EXPID)] # remove the stuff related to the gp fit in step4
sysDt[,IDX:=seq_len(.N)]

# prepare the handlers to map systematic uncertainties of the experiments
normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

exp_sys_cov <- sysCompHandler$cov(sysDt, ret.mat = TRUE) # U matrix 
exp_sys_map <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE) # S matrix
exp_stat_unc <- Diagonal(x = expDt[,UNC^2]) # D matrix

expDt[,TOT_UNC:=sqrt(diag(exp_sys_map %*% exp_sys_cov %*% t(exp_sys_map) + exp_stat_unc))]
# create a list containing all the experimental uncertainties in the format
# of nucdataBaynet to pass to the GP model
exp_cov_mat_list <- list(D=exp_stat_unc, S=exp_sys_map, U=exp_sys_cov)


##########################################
# energy grid for the defect model
energies <- expDt[EXPID==22316003,sort(L1)]
energies <- c(energies,expDt[L1>max(energies),L1])
energies <- unique(round(energies*1e4))*1e-04

# take every second energy
last_energy <- energies[length(energies)] 
energies <- energies[rep(c(TRUE,FALSE),length.out=length(energies))]

if(energies[1] > expDt[,min(L1)]) {
  if(abs(energies[1] - expDt[,min(L1)]) < min(diff(energies))) {
    energies[1] <- expDt[,min(L1)]
  } else {
    energies <- c(expDt[,min(L1)],energies)
  }
}
if(energies[length(energies)] < expDt[,max(L1)]) energies <- c(energies, expDt[,max(L1)])

##########################################
# Initiate the defect model
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/talys-fit-fe56"
model <- defect_model(energies, talys_calc_dir, expDt)
pars <- model$parsDt[,V1]

# and map the experimental energies, to check that all looks ok
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
########################################

# energy grid for the sigma hyper-parameter
minE <- min(energies) - 1e-04
maxE <- max(energies) + 1e-04
hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out =30)

#############################################################################
# definition of the GP modelling object
hetGP_model <- function(hyper_par_energy_grid, hyper_pars, model, observed, observed_cov_mat, nugget=1e-04) {
  # observed_cov_mat should be a list holding named entries
  # D - a diagonal matrix with statistical uncertainties
  # U - a matrix holding the systematic components
  # S - a matrix that maps the systematic components in U to the data vector
  
  # define for convenience
  n_sigmas <- length(hyper_par_energy_grid)
  
  correlation_func <- function(x, length_scale) {
    d_over_l <- abs(outer(x,x,"-"))/length_scale
    
    idx <- which(d_over_l<1, arr.ind=TRUE)
    d_over_l <- d_over_l[idx]
    data <- (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)
    
    sparseMatrix(i=idx[,1], j=idx[,2], x=data)
  }
  
  make_cov_mat <- function(hyper_pars, model, nugget=1e-04) {
    sigmas <- hyper_pars[1:n_sigmas]
    len <- hyper_pars[n_sigmas + 1]
    
    model_energies <- model$parsDt[,unique(L1)]
    cov_mat_inner <-  correlation_func(model_energies, length_scale=len)
    
    sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
    sigma_scale <- outer(sigmas_tmp,sigmas_tmp)
    
    cov_mat_inner*sigma_scale + diag(x=nugget, nrow=nrow(cov_mat_inner))
  }
  
  update_cov_mat <- function(hyper_pars) {
    # if the hyper-parameters have changed from the last call
    # re-build covariance matrix
    if(!all(hyper_pars==cur_hyper_pars)) {
      cur_hyper_pars <<- hyper_pars
      cur_cov_mat <<- make_cov_mat(cur_hyper_pars, model)
    }
  }
  
  latent_cov_mat <- function(sigma, len) {
    sigma^2*correlation_func(hyper_par_energy_grid, length_scale=len)
  }
  
  log_likelihood <- function(hyper_pars) {
    
    # update the covariance matrix if hyper-parameters have changed
    update_cov_mat(hyper_pars)
    
    UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    SS <- cbind(observed_cov_mat$S, mapping_matrix)
    cholZ <- makeCholZ(observed_cov_mat$D, SS, UU) # from nucdataBaynet
    # Could this be done faster since only the second block of U is updated each time?
    
    chi_square <- chisquare(observed, observed_cov_mat$D, SS, UU, cholZ=cholZ) # chisquare from nucdataBaynet utilizing mult_xt_invCov_x()
    log_det_cov <- logDetCov(observed_cov_mat$D, SS, UU, cholZ=cholZ) # from nucdataBaynet
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
    chi_square <- (sigmas - sqrt(sigmas_mean)) %*% solve(sigmas_prior_cov_mat, (sigmas-sqrt(sigmas_mean)))
    
    P <- as.vector(-0.5*(length(sigmas)*log(2*pi) + log_det_sigmas_prior_cov_mat + chi_square))
    
    n_calls <<- n_calls + 1
    if((n_calls %% 25) == 0) print(paste("------- log_likelihood calls:",n_calls,":",L+P))
    
    L + P
  }
  
  get_pars_cov_mat <- function(hyper_pars) {
    # update the covariance matrix if hyper-parameters have changed
    update_cov_mat(hyper_pars)
    
    return(cur_cov_mat)
  }
  
  conditional_mean_exp <- function(hyper_pars) {
    # returns the conditional mean of parameters mapped to the experimental
    # energies. 
    
    update_cov_mat(hyper_pars)

    UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    SS <- cbind(observed_cov_mat$S, mapping_matrix)

    Sigma_yy_inv_y <- mult_invCov_x(observed, observed_cov_mat$D, SS, UU)
    
    return(mapping_matrix %*% cur_cov_mat %*% t(mapping_matrix) %*% Sigma_yy_inv_y)
  }
  
  conditional_covariance_exp <- function(hyper_pars) {
    # returns the conditional covariance matrix of parameters mapped to the
    # experimental energies. 
    update_cov_mat(hyper_pars)
    
    UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    SS <- cbind(observed_cov_mat$S, mapping_matrix)
    
    SxKxST <- mapping_matrix %*% cur_cov_mat %*% t(mapping_matrix)
    
    return(SxKxST - mult_xt_invCov_x(SxKxST,exp_cov_mat_list$D, SS, UU))
  }
  
  conditional_mean_pars <- function(hyper_pars) {
    # returns the conditional mean of parameters
    update_cov_mat(hyper_pars)
    
    UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    SS <- cbind(observed_cov_mat$S, mapping_matrix)
    
    return(mapping_matrix_pars  %*% cur_cov_mat %*% t(mapping_matrix) %*%  mult_invCov_x(observed, exp_cov_mat_list$D, SS, UU))
  }
  
  conditional_covariance_pars <- function(hyper_pars) {
    # returns the conditional covariance matrix of parameters
    update_cov_mat(hyper_pars)
    
    UU <- bdiag(observed_cov_mat$U, cur_cov_mat)
    SS <- cbind(observed_cov_mat$S, mapping_matrix)
    
    SxKxST <- mapping_matrix_pars %*% cur_cov_mat %*% t(mapping_matrix_pars)
    
    return(SxKxST - mult_xt_invCov_x(SxKxST,exp_cov_mat_list$D, SS, UU))
  }
  
  # member variables used by both log_likelihood() and grad_log_likelihood()
  cur_hyper_pars <- hyper_pars
  cur_cov_mat <- make_cov_mat(cur_hyper_pars, model)
  
  # create the mapping matrix
  # maps between the parameters of the model and the covariance matrix to match the energies
  # it also takes into account the scaling, to convert the value in the covariance matrix to
  # a relative variance. So multiplying this matrix (from the left) with a unit vector
  # will result in the default TALYS prediction for the exclusive channels at the energy grid
  # of the defect model (i.e. the column V1 in model$parsDt), multiplying this result with
  # model$jac will thereby result in the default TALYS prediction at the experimental points in
  # model$expDt. Sandwiching the matrix created by make_cov_mat() will result in the full
  # parameter covariance matrix K_beta
  model_energies <- model$parsDt[,unique(L1)]
  target_dt <- data.table(L1 = model$parsDt[,unique(L1)], ValueIndex = seq_along(model_energies))
  map_dt <- model$parsDt[target_dt, on = .(L1), nomatch = 0L, .(IDX, IDX2 = i.ValueIndex, V1)]
  mapping_matrix_pars <- sparseMatrix(map_dt[,IDX],map_dt[,IDX2],x=map_dt[,V1])
  
  # finally multiply together the mapping matrix and the model Jacobian so it can be used
  # to effidiently calculate the product r %*% (J %*% K_beta %*% J^T + Sigma_exp)^-1 %*% r^T
  # using the woodbury trick from nucdataBaynet
  mapping_matrix <- model$jac %*% mapping_matrix_pars
  
  # the mean of the GP on the sigma hyper-paramter is set to the estimate of the
  # standard deviation of the residual/default prediction, we use the more robust
  # estimate of the standard deviation obtained from 1.4826*MAD(x)
  #sigmas_mean <- 1.4826*mad(observed/(model$jac %*% model$parsDt[,V1]))
  #sigmas_mean <- 0.1
  sigmas_mean <- 0
  
  n_calls <- 0
  
  list(marginal_likelihood=log_likelihood,
       get_pars_cov_mat=get_pars_cov_mat,
       make_cov_mat = make_cov_mat,
       mapping_matrix = mapping_matrix,
       conditional_mean_exp = conditional_mean_exp,
       conditional_covariance_exp = conditional_covariance_exp,
       conditional_mean_pars=conditional_mean_pars,
       conditional_covariance_pars=conditional_covariance_pars,
       mapping_matrix_pars=mapping_matrix_pars)
}

# Test to use the GP model on data
guessSigma <- seq(from=0.2, to=0.1, length.out=length(hyper_par_energy_grid))
guessLen <- min(diff(energies))*3

guess_latent_sigma <- 0.4
guess_latent_len <- 4
hyper_pars_init <- c(guessSigma, guessLen, guess_latent_sigma, guess_latent_len)

# the observed data to optimize against: here the residual between the TALYS calculation
# used to create the defect model and the experimental cross-section data
y_data <-  expDt[,RESIDUAL]

myGPmodel <- hetGP_model(hyper_par_energy_grid,hyper_pars_init, model, y_data, exp_cov_mat_list, nugget=1e-04)

########################################################
# Create a graphical representation of the guess parameters

# prior covariance matrix mapped to the experiments
prior_cov_exp <- myGPmodel$mapping_matrix %*% myGPmodel$get_pars_cov_mat(hyper_pars_init) %*% t(myGPmodel$mapping_matrix)
expDt[,PRIOR_GP_UNC:=sqrt(diag(prior_cov_exp))]

# Draw a sample of the prior to visualise the length scale
# need to add an extra nugget for it to work (a bit vorying)
# If I don't cholesky decomposition fails
sample <- rmvn(1, rep(0,nrow(expDt)), prior_cov_exp + diag(x=1e-04,nrow=nrow(prior_cov_exp)), ncores = 1, isChol = FALSE)
# sample <- rmvn(1, rep(0,nrow(expDt)), prior_cov_exp, ncores = 1, isChol = FALSE)
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

############################################################
# now for the time consuming part: the hyper-parameter optimization

# limits of the hyper-parameters:
# It may be important to set proper box constraints on the length hyper-parameter
# This is because this affects the sparseness of the covariance matrix, and thereby
# the computation time (and memory use), if the optimizer starts to move out to 
# try very long length scales the calculation of the log-likelihood will slow down
# and the program may even crash if it runs out of memory. The proper length scale
# should reproduce the structure (few keV wide) so it is in the order of magnitude
# of the distance of the energy grid. Smaller than this will result in a diagonal
# cov matrix (so no smoothing), much larger than this will smooth out the structure.
# So the limits are set to l_min = min(diff(energies)) and l_max = 100 keV
min_hyper_pars <- c(rep(1e-03,length(guessSigma)),min(diff(energies)), 0.01, 0.1)
max_hyper_pars <- c(rep(1,length(guessSigma)),0.1, 10, 10)

# if you restart: hyper_pars_init <- optim_res$par
optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-6,length(hyper_par_energy_grid)),1e-07,1e-06,1e-06), maxit=100, factr=1e9)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)
# The convergence is quite slow requiring several hundred iterations

# converged with factr = 1e9
# X = 0.138694 0.122157 0.0963145 0.0737109 0.0594644 0.0519538 0.0458768 0.0390492 0.0333865 0.0296797 0.0265543 0.0235924 0.0243843 0.0338693 0.048549 0.05547 0.0482722 0.0347104 0.0267511 0.0295745 0.0381037 0.0441097 0.0481323 0.0565748 0.072904 0.0971988 0.126311 0.149135 0.151725 0.128329 0.00872765 0.0430673 4.18825 
# F = 21989.5
# final  value 21989.534292 
# converged

# with factr = 1e8
# LINE SEARCH 0 times; norm of step = 0.000323135
# X = 0.135261 0.122533 0.0975302 0.0740968 0.0592017 0.0514808 0.0457237 0.039037 0.0329877 0.0292485 0.0268311 0.0251733 0.0272793 0.0354216 0.0454496 0.0501119 0.0463224 0.0383574 0.0330552 0.0337998 0.0379449 0.0417224 0.0468316 0.0570177 0.0725147 0.0899733 0.10561 0.115193 0.114572 0.101618 0.00876602 0.0371622 4.3695 
# G = -5.85916 29.4242 -3.81785 13.5817 -20.5146 30.1048 3.0307 -24.6951 5.33694 20.122 -31.0319 40.3134 -52.3113 52.8382 -27.8476 -4.16474 -0.897169 -13.1489 20.5437 -17.1654 -18.3887 16.3177 1.86366 -3.30497 -14.3662 -2.83999 14.4279 -3.42281 -2.04791 4.22077 36.0481 4.99902 -17.1852 
# final  value 21983.950280 
# stopped after 101 iterations
###############################################################################
# create  graphical representation of the hyper-parameter optimized prior
hyper_pars_optim <- optim_res$par

prior_cov_exp <- myGPmodel$mapping_matrix %*% myGPmodel$get_pars_cov_mat(hyper_pars_optim) %*% t(myGPmodel$mapping_matrix)
expDt[,PRIOR_GP_UNC:=sqrt(diag(prior_cov_exp))]

# Draw a sample of the prior to visualise the length scale
# need to add an extra nugget for it to work (a bit vorying)
# If I don't cholesky decomposition fails
sample <- rmvn(1, rep(0,nrow(expDt)), prior_cov_exp + diag(x=1e-04,nrow=nrow(prior_cov_exp)), ncores = 1, isChol = FALSE)
# sample <- rmvn(1, rep(0,nrow(expDt)), prior_cov_exp, ncores = 1, isChol = FALSE)
expDt[,GP_SAMPLE:=sample[1,]]

ggp <- ggplot(data=expDt[L1<=max(energies)]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey', width = 0.2) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  #geom_point(aes(x=L1,y=DATA)) +
  geom_line(aes(x=L1,y=default_prediction),col='green') +
  geom_ribbon(aes(x=L1,ymin=default_prediction-PRIOR_GP_UNC, ymax=default_prediction+PRIOR_GP_UNC), fill='green', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction+GP_SAMPLE),col='red') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp
###############################################################################
# Now that we have the hyper-parameters chosen we do the conditioning on data

expDt[,RESIDUAL_GP:=as.vector(myGPmodel$conditional_mean_exp(hyper_pars_optim))]
expDt[,RESIDUAL_GP_UNC:=sqrt(diag(myGPmodel$conditional_covariance_exp(hyper_pars_optim)))]

expDt[,RESIDUAL_GP_UNC:=sqrt(diag(test))]

ggp <- ggplot(data=expDt[L1<=max(energies)]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey', width = 0.2) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),col='green') +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
  geom_ribbon(aes(x=L1, ymin=default_prediction+RESIDUAL_GP-RESIDUAL_GP_UNC, ymax=default_prediction+RESIDUAL_GP+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

###############################################################################
# make a fake expDt to plot the GP model at all energies of the model
# a bit of a trick:
# first create a fake expDt that has all the channels of the original expDt
# but at all energies of the model
new_energy_grid = energies
fake_expDt <- data.table(
  REAC=rep(expDt[,unique(REAC)], each=length(new_energy_grid)),
  L1=rep(new_energy_grid, times=length(expDt[,unique(REAC)]))
)

fake_expDt[,IDX:=seq_len(.N)]

# then create a new defect_model instance to get the Jacobian matrix that
# maps the model parameters to the new fake experimental energies
# it is vital that this defect_model instance has exactly the same energy grid
# as the original one
smooth_model <- defect_model(energies, talys_calc_dir, fake_expDt)

# We need the TALYS predicitions at this new energy grid too
# get the default model parameters (talys prediction)
pars <- smooth_model$parsDt[,V1]
# and map the fake experimental energies
predictions <- smooth_model$jac %*% pars
fake_expDt[,default_prediction:=as.vector(predictions)]

# finally we get the conditional mean and covariance matrix from the GP model
# and map them to the new energy grid using the Jacobian of the just created
# defect_model
cond_mean_pars <- myGPmodel$conditional_mean_pars(hyper_pars_optim)
fake_expDt[,RESIDUAL_GP:=as.vector(smooth_model$jac %*% cond_pars_mean)]

cond_cov_pars <- myGPmodel$conditional_covariance_pars(hyper_pars_optim)
cond_cov_smooth <- smooth_model$jac %*% cond_cov_pars %*% t(smooth_model$jac)
fake_expDt[,RESIDUAL_GP_UNC:=sqrt(diag(cond_cov_smooth))]


ggp <- ggplot(data=expDt) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='green') +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt,col='red') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='red', alpha=0.3) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp

filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'Fe56-with-full-covariance.pdf')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

ggp_LE <- ggplot(data=expDt[L1<5]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[L1<5],col='green') +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[L1<5],col='red') +
  #geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[L1<5] , fill='red', alpha=0.3) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp_LE

filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'Fe56-with-full-covariance-LE.pdf')
ggsave(filepath, ggp_LE, width = 29.7, height = 21.0, units = "cm", dpi = 300)
