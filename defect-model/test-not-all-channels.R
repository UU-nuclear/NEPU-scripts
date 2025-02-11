
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
# hyper_pars_init <- c(0.135261, 0.122533, 0.0975302, 0.0740968, 0.0592017, 0.0514808, 0.0457237, 0.039037, 0.0329877, 0.0292485, 0.0268311, 0.0251733, 0.0272793, 0.0354216, 0.0454496, 0.0501119, 0.0463224, 0.0383574, 0.0330552, 0.0337998, 0.0379449, 0.0417224, 0.0468316, 0.0570177, 0.0725147, 0.0899733, 0.10561, 0.115193, 0.114572, 0.101618, 0.00876602, 0.0371622, 4.3695)

# if you restart: hyper_pars_init <- optim_res$par
optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-5,length(hyper_par_energy_grid)),1e-06,1e-05,1e-05), maxit=300, factr=1e9)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)

# LINE SEARCH 1 times; norm of step = 0.00799819
# X = 0.158417 0.115876 0.0930501 0.0776181 0.0614138 0.0517356 0.0497783 0.0406824 0.030559 0.0333295 0.031768 0.0200824 0.0157851 0.0238626 0.0420141 0.0545437 0.0527597 0.0424949 0.0475184 0.0725973 0.112394 0.14996 0.175622 0.196789 0.223271 0.265886 0.307425 0.311326 0.291844 0.283697 0.00911704 0.192767 5.12629 
# G = 2.97149 -48.7286 3.78161 52.9458 71.415 -16.5414 74.5622 188.695 -317.533 322.168 -20.6438 -173.948 146.638 -146.201 127.797 -287.591 262.765 -211.889 88.7664 -48.0199 -130.583 136.317 -100.372 132.606 -121.597 -54.4999 172.626 -39.4694 -105.285 61.0044 319.005 88.5357 -2.76949 
# final  value 22025.618286 
# stopped after 101 iterations

# LINE SEARCH 0 times; norm of step = 0.00372558
# X = 0.403902 0.363483 0.3072 0.245941 0.19086 0.151112 0.131027 0.123588 0.118088 0.107788 0.0908141 0.0734865 0.0583572 0.0486824 0.0383457 0.0322005 0.0293511 0.0233261 0.0243724 0.0315238 0.0286667 0.0025142 0.196584 7.11467 
# G = -426.758 963.879 -773.522 109.047 634.838 -1330.44 72.2141 785.248 -851.504 1398.63 -436.049 5.592 221.231 100.731 -230.182 -296.36 436.457 -233.062 -418.517 270.928 -40.8138 -8999.79 -13.5726 -0.198714 

# iterations 166
# function evaluations 82
# segments explored during Cauchy searches 89
# BFGS updates skipped 0
# active bounds at final generalized Cauchy point 0
# norm of the final projected gradient 2.88468
# final function value 65759.8
# 
# X = 0.414672 0.356994 0.299251 0.247072 0.202301 0.165821 0.138185 0.1192 0.106387 0.0949331 0.0818516 0.0683738 0.0554554 0.049242 0.0425604 0.0338519 0.028479 0.0293926 0.0324041 0.0325962 0.0278308 0.00255047 0.220547 7.11532 
# F = 65759.8
# final  value 65759.834569 
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

cur_cov_mat <- myGPmodel$get_pars_cov_mat(hyper_pars_optim)
mapping_matrix <- myGPmodel$mapping_matrix
mapping_matrix_pars <- myGPmodel$mapping_matrix_pars

UU <- bdiag(exp_cov_mat_list$U, cur_cov_mat)
SS <- cbind(exp_cov_mat_list$S, mapping_matrix)

# SxKxST <- mapping_matrix_pars %*% cur_cov_mat %*% t(mapping_matrix)
SxKxST <- mapping_matrix %*% cur_cov_mat %*% t(mapping_matrix_pars)

cond_cov_pars <- mapping_matrix_pars %*% cur_cov_mat %*% t(mapping_matrix_pars) - mult_xt_invCov_x(SxKxST,exp_cov_mat_list$D, SS, UU)

####


ggp <- ggplot(data=expDt) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt,col='red') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='red', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='green') +
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
# Try to set all experimental stat. unc equal to 3%

expDt_mod <- copy(expDt)

expDt_mod[, TOT_UNC:=sqrt(TOT_UNC^2 - UNC^2)]
expDt_mod[,UNC:=0.10*DATA]
expDt_mod[, TOT_UNC:=sqrt(TOT_UNC^2 + UNC^2)]

exp_cov_mat_list_mod <- copy(exp_cov_mat_list)
exp_cov_mat_list_mod$D <- Diagonal(x = expDt_mod[,UNC^2])

myGPmodel_test <- structureGPmodel(hyper_par_energy_grid,hyper_pars_init, model, y_data, exp_cov_mat_list_test, nugget=1e-04)
cond_mean_pars_test <- myGPmodel_test$conditional_mean_pars(hyper_pars_optim)

fake_expDt[,RESIDUAL_GP:=as.vector(smooth_model$jac %*% cond_mean_pars)]

ggp_LE <- ggplot(data=expDt_mod[L1<2]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[L1<2],col='red') +
  #geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[L1<2] , fill='red', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[L1<2],col='green') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp_LE
