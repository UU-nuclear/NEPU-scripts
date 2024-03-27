
rm(list=ls())
gc()

library(mvnfast)

source("config/config-Fe56.R")
source("defect-model/defect-model.R")
source("defect-model/structureGPmodel.R")

plotPath <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/figures"
dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
# prepare the experimental data
expDt <- read_object(3, "expDt")
expDt[,UNC:=ORIG_UNC] # UNC was replaced using hetGP: undo

# reduce the energy range
setorder(expDt,REAC,L1)

expDt_full <- copy(expDt)

minE <- 1
maxE <- 5
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

############################################

expDt[,RESIDUAL:=DATA-default_prediction]
# energy grid for the sigma hyper-parameter
minE <- min(energies) - 1e-04
maxE <- max(energies) + 1e-04
#hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out =30)
hyper_par_energy_grid <- seq(from=minE, to=maxE, length.out = 10)

# Test to use the GP model on data
guessSigma <- seq(from=0.5, to=0.1, length.out=length(hyper_par_energy_grid))
guessLen <- min(diff(energies))*3

guess_latent_sigma <- 0.4
guess_latent_len <- 4
hyper_pars_init <- c(guessSigma, guessLen, guess_latent_sigma, guess_latent_len)

# optimized values from a former run of the code
hyper_pars_init <- c(0.406154382, 0.369970063, 0.333784037, 0.297779673, 0.261661200, 0.225591640, 0.189486373, 0.153439282, 0.117371474, 0.081253550, 0.002479588, 0.326652088, 4.880907326)
# the observed data to optimize against: here the residual between the TALYS calculation
# used to create the defect model and the experimental cross-section data
y_data <-  expDt[,RESIDUAL]

myGPmodel <- structureGPmodel(hyper_par_energy_grid,hyper_pars_init, model, y_data, exp_cov_mat_list, nugget=1e-04)

##################################################

min_hyper_pars <- c(rep(1e-03,length(guessSigma)),min(diff(energies)), 0.01, 0.1)
max_hyper_pars <- c(rep(1,length(guessSigma)),0.05, 10, 10)

# optimized parameters from an earlier run
cat("log-likelihood at hyper-parameter guess:", myGPmodel$marginal_likelihood(hyper_pars_init), "\n")
# if you restart: hyper_pars_init <- optim_res$par
optim_control <- list(fnscale=-1, trace=6, REPORT=1, ndeps=c(rep(1e-5,length(hyper_par_energy_grid)),1e-06,1e-05,1e-05), maxit=100, factr=1e10)
optim_res <- optim(hyper_pars_init, fn=myGPmodel$marginal_likelihood, method = "L-BFGS-B", lower=min_hyper_pars, upper=max_hyper_pars, control=optim_control)

hyper_pars_optim <- optim_res$par
# iterations 11
# function evaluations 20
# segments explored during Cauchy searches 23
# BFGS updates skipped 0
# active bounds at final generalized Cauchy point 0
# norm of the final projected gradient 3.92602
# final function value 51877.1
# 
# X = 0.406154 0.36997 0.333784 0.29778 0.261661 0.225592 0.189486 0.153439 0.117371 0.0812536 0.00247959 0.326652 4.88091 
# F = 51877.1
# final  value 51877.128049 
# converged
##################################################
prior_cov_exp <- myGPmodel$mapping_matrix %*% myGPmodel$get_pars_cov_mat(hyper_pars_optim) %*% t(myGPmodel$mapping_matrix)
expDt[,PRIOR_GP_UNC:=sqrt(diag(prior_cov_exp))]

ggp <- ggplot(data=expDt) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),col='orange') +
  geom_ribbon(aes(x=L1,ymin=default_prediction-PRIOR_GP_UNC, ymax=default_prediction+PRIOR_GP_UNC), fill='orange', alpha=0.5) +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggp

filepath <- file.path(plotPath, 'Fe-56_optimized_prior.png')
ggsave(filepath, ggp)
filepath <- file.path(plotPath, 'Fe-56_optimized_prior.pdf')
ggsave(filepath, ggp)

################################################################

expDt[,RESIDUAL_GP:=as.vector(myGPmodel$conditional_mean_exp(hyper_pars_optim))]
conditional_cov_mat_exp <- myGPmodel$conditional_covariance_exp(hyper_pars_optim)
# conditional_corr_mat_exp <- cov2cor(conditional_cov_mat_exp)
expDt[,RESIDUAL_GP_UNC:=sqrt(diag(conditional_cov_mat_exp))]

# I can't calculate the conditional correlation matrix, because the is a zero along
# the diagonal in the conditional covariance matrix. This originates from the fact 
# that our prior covariance is a relative uncertainty. The zero in the covariance
# matrix appears for an experimental point of the (n,p) channel at 3.43 MeV, which
# is below the threshold for this reaction. At least the threshold determined from
# interpolation on the talys result. The theoretical threshold is 2.967 MeV.
# This is kind of a weakness of how the sum-rule model is constructed only from
# the TALYS calculation. It could set its first point to always be at the threshold
# for a reaction - calculated from tabulated masses.

ggp <- ggplot(data=expDt[L1<=max(energies)]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction),col='green') +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
  geom_ribbon(aes(x=L1, ymin=default_prediction+RESIDUAL_GP-RESIDUAL_GP_UNC, ymax=default_prediction+RESIDUAL_GP+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggp
################################################################
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
cond_cov_pars <- myGPmodel$conditional_covariance_pars(hyper_pars_optim)

fake_expDt[,RESIDUAL_GP:=as.vector(smooth_model$jac %*% cond_mean_pars)]
fake_expDt[,RESIDUAL_GP_UNC:=sqrt(diag(smooth_model$jac %*% cond_cov_pars %*% t(smooth_model$jac)))]

ggp <- ggplot(data=expDt) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt,col='green') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='green', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='red') +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggp

filepath <- file.path(plotPath, 'Fe-56_1.00-5.00MeV.pdf')
ggsave(filepath, ggp,  width=0.5*210, height=0.5*297, units='mm')
filepath <- file.path(plotPath, 'Fe-56_1.00-5.00MeV.png')
ggsave(filepath, ggp,  width=0.5*210, height=0.5*297, units='mm')

channels <- c("(26-FE-56(N,EL)26-FE-56,,SIG)", "(26-FE-56(N,INL)26-FE-56,,SIG)", "(26-FE-56(N,TOT),,SIG)" )
ggp <- ggplot(data=expDt[REAC %in% channels]) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[REAC %in% channels],col='green') +
  geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[REAC %in% channels] , fill='green', alpha=0.3) +
  geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[REAC %in% channels], col='red') +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw() +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
ggp

filepath <- file.path(plotPath, 'Fe-56_1.00-5.00MeV-not-all.pdf')
ggsave(filepath, ggp,  width=0.5*210, height=0.5*297, units='mm')
filepath <- file.path(plotPath, 'Fe-56_1.00-5.00MeV-not-all.png')
ggsave(filepath, ggp,  width=0.5*210, height=0.5*297, units='mm')

minE <- 1 # MeV
maxE <- minE + 0.25 # MeV
while(minE<expDt[,max(L1)]) {
  ggp <- ggplot(data=expDt[L1>minE & L1<maxE & REAC %in% channels]) +
    geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
    geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
    geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt[L1>minE & L1<maxE & REAC %in% channels],col='green') +
    geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt[L1>minE & L1<maxE & REAC %in% channels] , fill='green', alpha=0.3) +
    geom_line(aes(x=L1,y=default_prediction),data=fake_expDt[L1>minE & L1<maxE & REAC %in% channels],col='red') +
    facet_wrap(~ REAC, ncol=1, scales="free_y") +
    labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
    theme_bw() +
    theme(panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank())
  print(ggp)
  
  filepath <- file.path(plotPath, paste0('Fe-56_',format(minE, nsmall=3),'-',format(maxE, nsmall=3),'MeV.pdf'))
  ggsave(filepath, ggp,  width=0.5*210, height=0.5*297, units='mm')
  
  
  filepath <- file.path(plotPath, paste0('Fe-56_',format(minE, nsmall=3),'-',format(maxE, nsmall=3),'MeV.png'))
  ggsave(filepath, ggp,  width=0.5*210, height=0.5*297, units='mm')
  
  minE <- maxE
  maxE <- minE + 0.25
}

################################################################################
# plot covariance matrices
library(reshape2)
library(gplots)
library(RColorBrewer)

smooth_cov_mat <- smooth_model$jac %*% cond_cov_pars %*% t(smooth_model$jac)

# this covariance matrix contains zeros along the diagonal since we treat the
# uncertainty as relative only in the GP, therefore we should remove entries 
# corresponding to those where the TALYS prediction is zero

# test of the covariance matrix for (n,tot) channel
idcs <- fake_expDt[REAC=="(26-FE-56(N,TOT),,SIG)" & RESIDUAL_GP!=0, IDX]
sub_cov_mat <- smooth_cov_mat[idcs,idcs]

ncolors <- 256
color_palette <- colorRampPalette(c("red","white","blue"))(ncolors)

heatmap.2(as.matrix(cov2cor(sub_cov_mat)),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

image(as.matrix(cov2cor(sub_cov_mat)), useRaster=TRUE)

################################################################################

# test of the covariance matrix for (n,tot) channel
idcs <- expDt[REAC=="(26-FE-56(N,TOT),,SIG)" & RESIDUAL_GP!=0, IDX]
sub_cov_mat <- conditional_cov_mat_exp[idcs,idcs]

image(as.matrix(cov2cor(sub_cov_mat)), useRaster=TRUE)
