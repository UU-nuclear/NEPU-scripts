
rm(list=ls())
gc()

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
maxE <- 1.05
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

################################################################################

expDt_test <- copy(expDt)
sysDt_test <- copy(sysDt)

# test to only leave the inl channel
# expDt_test <- expDt_test[grepl("N,INL",REAC)]
# expDt_test[,IDX:=seq_len(.N)]
# expids <- paste0("EXPID-",expDt_test[,unique(EXPID)])
# sysDt_test <- sysDt_test[EXPID %in% expids]
# sysDt_test[,IDX:=seq_len(.N)]

exp_sys_cov_test <- sysCompHandler$cov(sysDt_test, ret.mat = TRUE) # U matrix 
exp_sys_map_test <- sysCompHandler$map(expDt_test, sysDt_test, ret.mat = TRUE) # S matrix
exp_stat_unc_test <- Diagonal(x = expDt_test[,UNC^2]) # D matrix

expDt_test[,TOT_UNC:=sqrt(diag(exp_sys_map_test %*% exp_sys_cov_test %*% t(exp_sys_map_test) + exp_stat_unc_test))]
# create a list containing all the experimental uncertainties in the format
# of nucdataBaynet to pass to the GP model

exp_cov_mat_list_test <- list(D=exp_stat_unc_test, S=exp_sys_map_test, U=exp_sys_cov_test)

model_test <- defect_model(energies, talys_calc_dir, expDt_test)

hyper_par_energy_grid <- c(min(energies)-0.01,0.5*(min(energies)+max(energies)),max(energies)+0.01)
hyper_pars_init <- c(rep(0.5,length(hyper_par_energy_grid)),0.005,0.5,1.0)
testGPmodel <- structureGPmodel(hyper_par_energy_grid,hyper_pars_init, model_test, expDt_test[,RESIDUAL], exp_cov_mat_list_test, nugget=1e-06)

ggp <- ggplot(data=expDt) +
  geom_point(aes(x=L1,y=DATA), data=expDt_test, col='green') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction), col='red') +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

cond_mean_pars_test <- testGPmodel$conditional_mean_pars(hyper_pars_init)
model_test$parsDt[, COND:=as.vector(cond_mean_pars_test)]

################################################################################
# THIS IS WRONG!!!! IF there is only data in n,inl no other exclusive channels
# should be affected by the conditioning!!!!!

# try with a simple test, no systematics included

K_beta0 <- testGPmodel$get_pars_cov_mat(hyper_pars_init)
K_beta <- testGPmodel$mapping_matrix_pars %*% K_beta0 %*% t(testGPmodel$mapping_matrix_pars)

Sigma_beta_delta <- K_beta %*% t(model_test$jac)

Sigma_delta_delta <- Diagonal(x=expDt_test[,UNC^2]) + model_test$jac %*% K_beta %*% t(model_test$jac)

conditional_mean <-  Sigma_beta_delta %*% solve(Sigma_delta_delta,expDt_test[,RESIDUAL])

model_test$parsDt[, TEST:=as.vector(conditional_mean)]

################################################################################
# Try with the "manual" construction of K_beta

correlation_func <- function(x, length_scale) {
  d_over_l <- abs(outer(x,x,"-"))/length_scale
  
  idx <- which(d_over_l<1, arr.ind=TRUE)
  d_over_l <- d_over_l[idx]
  data <- (2 + cos(2*pi*d_over_l))/3 * (1 - d_over_l) + sin(2*pi*d_over_l)/(2*pi)
  
  sparseMatrix(i=idx[,1], j=idx[,2], x=data)
}

make_cov_mat <- function(hyper_pars, model, n_sigmas, nugget=1e-04) {
  
  sigmas <- hyper_pars[1:n_sigmas]
  len <- hyper_pars[n_sigmas + 1]
  
  priors <- list()
  for(reac in model$parsDt[,unique(REAC)]) {
    
    ###########
    model_energies <- model$parsDt[REAC==reac,L1]
    model_default <- model$parsDt[REAC==reac,V1]
    cov_mat <- correlation_func(model_energies, length_scale=len)
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

K_beta_test <- make_cov_mat(hyper_pars_init, model_test, length(hyper_par_energy_grid))

Sigma_beta_delta <- K_beta_test %*% t(model_test$jac)

Sigma_delta_delta <- Diagonal(x=expDt_test[,UNC^2]) + model_test$jac %*% K_beta_test %*% t(model_test$jac)

conditional_mean_pars <-  Sigma_beta_delta %*% solve(Sigma_delta_delta,expDt_test[,RESIDUAL])

model_test$parsDt[, TEST2:=as.vector(conditional_mean_pars)]

# with this way of doing it, it works! So the mapping must be wrong somehow!

conditional_mean_exp <- model_test$jac %*% conditional_mean_pars

expDt_test[, RESIDUAL_GP:=as.vector(conditional_mean_exp)]
expDt[,RESIDUAL_GP:=as.vector(model_full$jac %*% conditional_mean_pars)]

ggp <- ggplot(data=expDt) +
  geom_point(aes(x=L1,y=DATA), data=expDt_test, col='green') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction), col='red') +
  geom_line(aes(x=L1,y=default_prediction + RESIDUAL_GP), col='blue') +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

################################################################################
# Let's try to plot the covariance matrices to find out what goes wrong

library(reshape2)
library(gplots)
library(RColorBrewer)

# This is the manually constructed matrix
ncolors <- 256
color_palette <- colorRampPalette(c("red","white","blue"))(ncolors)

heatmap.2(as.matrix(cov2cor(K_beta_test)),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# This is the matrix constructed using the mapping
heatmap.2(as.matrix(cov2cor(K_beta)),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# This is the mapping matrix
tmp_mat <- testGPmodel$mapping_matrix_pars
tmp_mat[tmp_mat!=0] <- 1
heatmap.2(as.matrix(tmp_mat),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# test
mapping <- testGPmodel$mapping_matrix_pars
mapping[mapping!=0] <- 1
tmp_mat <- mapping %*% cov2cor(K_beta0)
heatmap.2(as.matrix(tmp_mat),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

################################################################################

hyper_pars <- hyper_pars_init
n_sigmas <- length(hyper_par_energy_grid)
model <- model_test

sigmas <- hyper_pars[1:n_sigmas]
len <- hyper_pars[n_sigmas + 1]

model_energies <- model$parsDt[,unique(L1)]
corr_mat <- correlation_func(model_energies, length_scale=len)

heatmap.2(as.matrix(corr_mat),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

n_channels <- length(model$parsDt[,unique(REAC)])

sigmas_tmp <- approx(hyper_par_energy_grid,sigmas, xout=model_energies)$y
sigma_scale <- outer(sigmas_tmp,sigmas_tmp)
cov_mat <- corr_mat*sigma_scale

heatmap.2(as.matrix(cov_mat),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# mapping matrix: TALYS predictions for each channel at each energy, along diagonal
target_dt <- data.table(L1 = model_energies, ValueIndex = seq_along(model_energies))
map_dt <- model$parsDt[target_dt, on = .(L1), nomatch = 0L, .(IDX, IDX2 = IDX, V1)]
mapping_matrix_pars <- sparseMatrix(map_dt[,IDX],map_dt[,IDX2],x=map_dt[,V1])

heatmap.2(as.matrix(mapping_matrix_pars),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

cov_mat_full <- bdiag(replicate(n_channels, cov_mat))

heatmap.2(as.matrix(cov_mat_full),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")


heatmap.2(as.matrix(mapping_matrix_pars %*% cov_mat_full %*% t(mapping_matrix_pars)),
          Rowv=NA,Colv=NA,symm=TRUE,col=color_palette,scale="none",
          margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# calculate conditional mean
observed <- expDt_test[,RESIDUAL]

mapping_matrix <- model$jac %*% mapping_matrix_pars

UU <- bdiag(exp_cov_mat_list_test$U, cov_mat_full)
SS <- cbind(exp_cov_mat_list_test$S, mapping_matrix)

Sigma_yy_inv_y <- mult_invCov_x(observed, exp_cov_mat_list_test$D, SS, UU)

conditional_mean_exp <- mapping_matrix %*% cov_mat_full %*% t(mapping_matrix) %*% Sigma_yy_inv_y
expDt_test[,RESIDUAL_GP_NEW:=as.vector(conditional_mean_exp)]

conditional_mean_pars <- mapping_matrix_pars  %*% cov_mat_full %*% t(mapping_matrix) %*%  mult_invCov_x(observed, exp_cov_mat_list_test$D, SS, UU)

ggp <- ggplot(data=expDt) +
  geom_point(aes(x=L1,y=DATA), data=expDt_test, col='green') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=default_prediction), col='red') +
  geom_line(aes(x=L1,y=default_prediction + RESIDUAL_GP_NEW), data=expDt_test, col='blue') +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp
