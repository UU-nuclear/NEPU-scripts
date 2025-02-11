
rm(list=ls())
gc()

library(mvnfast)

source("config/config-Fe56.R")
source("defect-model/defect-model.R")
source("defect-model/structureGPmodel.R")

# prepare the experimental data
expDt <- read_object(3, "expDt")
expDt[,UNC:=ORIG_UNC] # UNC was replaced using hetGP: undo

# reduce the energy range
setorder(expDt,REAC,L1)

expDt_full <- copy(expDt)

minE <- 1
maxE <- 10
expDt <- expDt[L1<maxE & L1>minE]
expDt[,IDX:=seq_len(.N)]

##########################################
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
# retrieve talys parameter uncertainties

# the experimental data.table used in the TALYS parameter optimization
optExpDt <- read_object(6, "optExpDt")
# list of results from TALYS parameter optimization
talys_optRes <- read_object(10, "optRes")
# The jacobian of TALYS at the parameter optimum
talys.jac <- talys_optRes$jac
# The TALYS parameter prior cov. mat.
# just get the parameter covariance calculated during the fit
talysParCovmat <- talys_optRes$parCovLM

# check that it is positive definite!
if(!isSymmetric(talysParCovmat, tol=1e-8)) {
  print("the talys parameter covariance matrix does not appear to be symmetric")
}
# symmetrize the parameter covariance matrix
talysParCovmat <- (talysParCovmat + t(talysParCovmat)) / 2
if(!all(eigen(talysParCovmat)$values >=0)) {
  print("the talys parameter covariance matrix is not positive definite")
}
# Since we are doing the current fit on a reordered experimental data.table, which is 
# also a subset of the experimental data.table used for the TALYS fit we need to
# rearrange talys.jac accordingly, to map the data as in expDt
optExpDt <- optExpDt[L1<maxE & L1>minE]

# order the experimental data.table as the expDt data.table above
setorder(optExpDt,REAC,L1)
optExpDt[,NEW_IDX:=seq_len(.N)]

# sanity check: these should be the same tables:
stopifnot(all(optExpDt[,c('EXPID','REAC','L1','DATA','DATAREF')] == expDt[,c('EXPID','REAC','L1','DATA','DATAREF')]))

# reorder the talys jacobian so that it maps the data in expDt[,DATA]
reorder_idcs <- optExpDt[,IDX]
talys.jac <- talys.jac[reorder_idcs,]


# Now we add the TALYS parameter uncertainty as a systematic error component in
# residual, i.e. exp_data - TALYS_prediction

exp_cov_mat_list$S <- cbind(exp_cov_mat_list$S, talys.jac)
exp_cov_mat_list$U <- bdiag(exp_cov_mat_list$U, talysParCovmat)

# TEST: replace the systematic experimental uncertainty with the TALYS uncertainty
# exp_cov_mat_list$S <- talys.jac
# exp_cov_mat_list$U <- talysParCovmat


# TEST2: Take only the prior of the TALYS parameters
# optParamDt <- read_object(5, "optParamDt")
# talysParCovmat <- Diagonal(x = optParamDt[ADJUSTABLE==TRUE,PARUNC^2])
# 
# exp_cov_mat_list$S <- cbind(exp_cov_mat_list$S, talys.jac)
# exp_cov_mat_list$U <- bdiag(exp_cov_mat_list$U, talysParCovmat)
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

# first a model that will give me the default talys prediction for all channels
# so i can calculate the residual
model_full <- defect_model(energies, talys_calc_dir, expDt)
pars <- model_full$parsDt[,V1]

# and map the experimental energies, to check that all looks ok
predictions <- model_full$jac %*% pars # this will be the default TALYS prediction

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
# by providing a vector of reaction channel IDs in the argument include_channels
# only entries in the jocabian that maps to these experimental data will be non-zero
# so for example here I choose to include (n,inl) and (n,el) from the experimental
# data. Only (n,inl), (n,el) channels correspond to parameters in the model so by setting
# all parameters to a value 500 these two channels will be predicted to 500, while the
# (n,tot) channel (which is the sum of all channels) will be predicted at 1000 = 500 + 500.
# In this way when modeling the residual I can chose to ignore exclusive channels that I don't
# want to impose fluctuations on with the GP modeling.
include_channels <- c("CS/REAC/100000/TOT","CS/EL")
model <- defect_model(energies, talys_calc_dir, expDt, include_channels=include_channels)
default_pars <- model$parsDt[,V1]
pars <- rep(500, length.out=length(default_pars))

# and map the experimental energies, to check that all looks ok
predictions <- model$jac %*% pars # this will be the default TALYS prediction

expDt[,default_prediction_alt:=default_prediction + as.vector(predictions)]

ggp <- ggplot(data=expDt) +
  geom_point(aes(x=L1,y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction), col='red') +
  geom_line(aes(x=L1,y=default_prediction_alt), col='green') +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp
############################################

expDt[,RESIDUAL:=DATA-default_prediction]
# energy grid for the sigma hyper-parameter
minE <- min(energies) - 1e-04
maxE <- max(energies) + 1e-04
#hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out =30)
hyper_par_energy_grid <- c(seq(from=minE, to=5, length.out = 12), seq(from=5.5, to=maxE, length.out=9))

# Test to use the GP model on data
guessSigma <- seq(from=0.5, to=0.1, length.out=length(hyper_par_energy_grid))
guessLen <- min(diff(energies))*3

guess_latent_sigma <- 0.4
guess_latent_len <- 4
hyper_pars_init <- c(guessSigma, guessLen, guess_latent_sigma, guess_latent_len)

# the observed data to optimize against: here the residual between the TALYS calculation
# used to create the defect model and the experimental cross-section data
y_data <-  expDt[,RESIDUAL]

myGPmodel <- structureGPmodel(hyper_par_energy_grid,hyper_pars_init, model, y_data, exp_cov_mat_list, nugget=1e-04)

##################################################

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

# optimized parameters from an earlier run
hyper_pars_init <- c(0.414672, 0.356994, 0.299251, 0.247072, 0.202301, 0.165821, 0.138185, 0.1192, 0.106387, 0.0949331, 0.0818516, 0.0683738, 0.0554554, 0.049242, 0.0425604, 0.0338519, 0.028479, 0.0293926, 0.0324041, 0.0325962, 0.0278308, 0.00255047, 0.220547, 7.11532)
cat("log-likelihood at hyper-parameter guess:", myGPmodel$marginal_likelihood(hyper_pars_init), "\n")
# if you restart: hyper_pars_init <- optim_res$par
optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-5,length(hyper_par_energy_grid)),1e-06,1e-05,1e-05), maxit=300, factr=1e9)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)

# iterations 2
# function evaluations 7
# segments explored during Cauchy searches 25
# BFGS updates skipped 0
# active bounds at final generalized Cauchy point 0
# norm of the final projected gradient 2.88449
# final function value 65737.7
# 
# X = 0.414711 0.356971 0.299297 0.247056 0.202353 0.16581 0.138242 0.119258 0.10638 0.0949926 0.0818461 0.0683693 0.0555177 0.0492389 0.0426235 0.0338497 0.0285431 0.0293907 0.0324681 0.032594 0.0278949 0.00254438 0.220533 7.11551 
# F = 65737.7
# final  value 65737.731441 
# converged

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
smooth_model_full <- defect_model(energies, talys_calc_dir, fake_expDt)
smooth_model <- defect_model(energies, talys_calc_dir, fake_expDt, include_channels = include_channels)

# We need the TALYS predicitions at this new energy grid too
# get the default model parameters (talys prediction)
pars <- smooth_model_full$parsDt[,V1]
# and map the fake experimental energies
predictions <- smooth_model_full$jac %*% pars
fake_expDt[,default_prediction:=as.vector(predictions)]

# finally we get the conditional mean and covariance matrix from the GP model
# and map them to the new energy grid using the Jacobian of the just created
# defect_model
cond_mean_pars <- myGPmodel$conditional_mean_pars(hyper_pars_optim)
######### TEST SOMETHING ###########
# something seems strange, the structure seems to only be conditioned on the (n,tot)
# channel. However, a test to blow up the statistical uncertainty on all (n,tot) data
# shows that in this case the structure is taken only from the (n,inl) channel. So, 
# the effect that originally the structure seems to only be conditioned on the (n,tot)
# channel is due to the much smaller uncertainty in these data.
# exp_cov_mat_list_test <- copy(exp_cov_mat_list)
# expDt_tmp <- copy(expDt)
# expDt_tmp[REAC=="(26-FE-56(N,TOT),,SIG)",UNC:=10*UNC]
# exp_cov_mat_list_test$D <- Diagonal(x = expDt_tmp[,UNC^2])
# 
# myGPmodel_test <- structureGPmodel(hyper_par_energy_grid,hyper_pars_init, model, y_data, exp_cov_mat_list_test, nugget=1e-04)
# cond_mean_pars <- myGPmodel_test$conditional_mean_pars(hyper_pars_optim)
####################################
fake_expDt[,RESIDUAL_GP:=as.vector(smooth_model$jac %*% cond_mean_pars)]

cond_cov_pars <- myGPmodel$conditional_covariance_pars(hyper_pars_optim)
cond_cov_smooth <- smooth_model$jac %*% cond_cov_pars %*% t(smooth_model$jac)
fake_expDt[,RESIDUAL_GP_UNC:=sqrt(diag(cond_cov_smooth))]

####
ggp <- ggplot(data=expDt) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt,col='green') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='green', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='red') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp

min_plot_en <- 1.75
max_plot_en <- 2
# channels <- c("(26-FE-56(N,TOT),,SIG)","(26-FE-56(N,INL)26-FE-56,,SIG)","(26-FE-56(N,EL)26-FE-56,,SIG)")
channels <- c("(26-FE-56(N,TOT),,SIG)","(26-FE-56(N,INL)26-FE-56,,SIG)")
ggp <- ggplot(data=expDt[REAC %in% channels][L1>min_plot_en & L1<max_plot_en]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[REAC %in% channels][L1>min_plot_en & L1<max_plot_en],col='green') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[REAC %in% channels][L1>min_plot_en & L1<max_plot_en] , fill='green', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[REAC %in% channels][L1>min_plot_en & L1<max_plot_en],col='red') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp
# filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'Fe56-with-full-covariance.pdf')
# ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

ggp_LE <- ggplot(data=expDt[L1<2]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[L1<2],col='red') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[L1<2] , fill='red', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[L1<2],col='green') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp_LE


# filepath <- file.path(getwd(), 'defect-model', 'test_hetGP-figs', 'Fe56-with-full-covariance-LE.pdf')
# ggsave(filepath, ggp_LE, width = 29.7, height = 21.0, units = "cm", dpi = 300)

################################################################################
# Generate samples from the posterior, and write to output files
################################################################################
# To get samples of all possible channels I need to create a fake expDt which
# has entries for all channels that I'm interested in. Here I have only put the
# structure on the (n,el) and (n,inl) channels, which propagates to the (n,tot)
# channel. This mean I can make a fake expDt with only these three channels to
# generate samples. The energy grid of this fake expDt should be the one from
# the model I fitted, i.e. energies

channels_to_sample <- c("(26-FE-56(N,TOT),,SIG)","(26-FE-56(N,EL)26-FE-56,,SIG)","(26-FE-56(N,INL)26-FE-56,,SIG)")
sampling_expDt <- data.table(
  REAC=rep(channels_to_sample, each=length(energies)),
  L1=rep(new_energy_grid, times=length(channels_to_sample))
)
sampling_expDt[,IDX:=seq_len(.N)]

# sampling_model_full is not needed for the sampling,
# it is only used for the graphical representation
# to get the default talys prediction
sampling_model_full <- defect_model(energies, talys_calc_dir, sampling_expDt)
pars <- smooth_model_full$parsDt[,V1]
# and map the fake experimental energies
predictions <- sampling_model_full$jac %*% pars
sampling_expDt[,default_prediction:=as.vector(predictions)]

sampling_model <- defect_model(energies, talys_calc_dir, sampling_expDt, include_channels = include_channels)
# cond_cov_pars is the conditional covariance matrix of 'parameters' from above
cond_cov_sampling_model <- sampling_model$jac %*% cond_cov_pars %*% t(sampling_model$jac)
cond_mean_sampling_model <- as.vector(sampling_model$jac %*% cond_mean_pars)

# check that our covariance matrix is positive definite!
if(!isSymmetric(cond_cov_sampling_model, tol=1e-10)) {
  print("the talys parameter covariance matrix does not appear to be symmetric")
}
# symmetrize the parameter covariance matrix
cond_cov_sampling_model <- (cond_cov_sampling_model + t(cond_cov_sampling_model)) / 2
if(!all(eigen(cond_cov_sampling_model)$values >=0)) {
  print("the talys parameter covariance matrix is not positive definite")
}

samples <- rmvn(1000, cond_mean_sampling_model, cond_cov_sampling_model + diag(x=1e-04,nrow=nrow(cond_cov_sampling_model)), ncores = 1, isChol = FALSE)
residual_samples_dt <- data.table(t(samples))
random_residual_col_names <- paste("sample", seq_len(nrow(samples)),sep="")
setnames(residual_samples_dt, random_residual_col_names)

sampling_expDt <- cbind(sampling_expDt,residual_samples_dt)

cols <- rainbow(5)
ggp <- ggplot(data=sampling_expDt) +
  geom_point(aes(x=L1,y=DATA),data=expDt) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='black') +
  geom_line(aes(x=L1,y=default_prediction + sample1),col=cols[1]) +
  geom_line(aes(x=L1,y=default_prediction + sample2),col=cols[2]) +
  geom_line(aes(x=L1,y=default_prediction + sample3),col=cols[3]) +
  geom_line(aes(x=L1,y=default_prediction + sample4),col=cols[4]) +
  geom_line(aes(x=L1,y=default_prediction + sample5),col=cols[5]) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp

################################################################################
# write random files in TALYS format with the updated (n,el), (n,inl) and (n,tot)
# cross-sections. This seems like unnecessary work. Better to just write a big
# data table with samples that Erik can use to merge with the random files from
# talys. Writing files with the structure of talys files fro the residual does not
# seem to give anythig extra but more work for Erik to read many files.

# write a single file with random residuals
fwrite(sampling_expDt[,c("REAC","L1",random_residual_col_names), with=FALSE], 
       file = "defect-model/random-residuals.txt", 
       sep = ",", 
       dec = ".", # Decimal separator
       quote = FALSE, # Do not quote the output
       eol = "\n", # End of line character
       na = "NA", # How to denote missing values
       append = FALSE) # Do not append to the file, overwrite

################################################################################
#
################################################################################
# Some problem related to the reconstruction of the talys parameter covariance
# Error in .local(A, ...) : 
#   leading principal minor of order 2387 is not positive
# In addition: Warning message:
#   In .local(A, ...) :
#   CHOLMOD warning 'not positive definite' at file '../Cholesky/t_cholmod_rowfac.c', line 430

# Not a positive definite covariance matrix, can it be due to the fact that we in principle use the
# experimental covariance twice and the TALYS parameter uncertainty is then just a linear map of the same
# covariance matrix?
# probably not. If I ignore the experimental systematic uncertainty and only include TALYS parameter covariance,
# I still get non-positive definite covariance.

# If I take only the prior talys parameter cov. it works. This may mean that the problem is due to the double use
# of the experimental covariance, it may also mean that the reorder talys.jac is somehow incorrectly ordered. Since
# the ordering does not matter much if the parameter cov. mat. is diagonal

# An easy way to check is to calculate the predictions at the experiments using y = y0 + J(p - p0)

talys.pars <- talys_optRes$par
talys.func <- talys_optRes$fn[reorder_idcs]
expDt[,test_pred:=talys.func] # this appears to be correct!

talys.pars.mod <- copy(talys.pars)
# propagate one sigma unc according to prior uncertainty
talys.pars.mod <- talys.pars.mod + optParamDt[ADJUSTABLE==TRUE,PARUNC]
expDt[,test_pred:= talys.func + as.vector(talys.jac %*% (talys.pars.mod - talys.pars) )]

ggp <- ggplot(data=expDt) +
  geom_point(aes(x=L1,y=DATA)) +
  geom_line(aes(x=L1,y=default_prediction),col='green') +
  geom_line(aes(x=L1,y=test_pred),col='blue') +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

# it looks reasonable
# It seems I was doing something wrong when recalculating the covariance matrix
# of talys parameters. If I just retrieve it from the LM result it works fine