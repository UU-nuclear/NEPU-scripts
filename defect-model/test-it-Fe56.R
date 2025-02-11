
rm(list=ls())
gc()

source("config/config-Fe56.R")
source("defect-model/defect-model.R")

expDt <- read_object(3, "expDt")
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc-Fe56"

# test on reduced energy range
setorder(expDt,REAC,L1)
expDt <- expDt[L1<3]
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
energies <- expDt[,sort(L1)]
# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
energies <- unique(round(energies*1e4))*1e-04

# seems there is an error when model energy grid falls exactly on experimental energy grid
# needs to be fixed!

# a simple not very dense energy grid to test
# energies <- seq(from=expDt[,min(L1)],to=expDt[,max(L1)],length.out=50)

# a simple dense grid
# energies <- seq(from=expDt[,min(L1)],to=expDt[,max(L1)], by=1e-03)
# if(max(energies)<expDt[,max(L1)]) energies <- c(energies,expDt[,max(L1)]+1e-03)


model <- defect_model(energies, talys_calc_dir, expDt)

# set all exclusive cross-sections to 1
# pars <- rep(1,nrow(model$parsDt))
pars <- model$parsDt[,V1]

# and map the experimental energies
predictions <- model$jac %*% pars

expDt[,test_prediction:=as.vector(predictions)]
expDt[,RESIDUAL:=DATA-test_prediction]


ggp <- ggplot(data=expDt[grepl("(N,INL)",REAC)]) +
		geom_line(aes(x=L1,y=test_prediction))
		#geom_line(aes(x=L1,y=DATA),col='red')
ggp

ggp <- ggplot(data=expDt[grepl("(N,TOT)",REAC)]) +
		geom_line(aes(x=L1,y=test_prediction)) #+
		#geom_line(aes(x=L1,y=DATA),col='red')
ggp

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1, y=DATA)) +
		geom_line(aes(x=L1,y=test_prediction),col='red') +
		facet_wrap(~ REAC, ncol=4)
ggp

# now try to create a GP prior covariance matrix on the parameters
# of the model

sigma2 <- (600)^2
l2 <- (3e-3)^2

cov_func <- function(x1,x2) {
	sigma2*exp(-0.5*(x1 - x2)^2/l2)
}

priors <- list()
for(reac in model$parsDt[,unique(REAC)]) {
	model_energies <- model$parsDt[REAC==reac,L1]
	cov_mat <- outer(model_energies,model_energies,cov_func)

	cov_mat <- Matrix(cov_mat, sparse = TRUE)

	priors <- append(priors,list(cov_mat))
}

full_prior_cov_mat <- bdiag(priors)
# each block in this matrix correspond to an exclusive reaction channel

library(gplots)

color_palette<-colorRampPalette(c("white","red","blue"),bias=1)
# heatmap.2(as.matrix(full_prior_cov_mat),Rowv=NA,Colv=NA,symm=TRUE,col=color_palette(256),scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# now test to propated this to the experiments

prior_exp <- model$jac %*% full_prior_cov_mat %*% t(model$jac)

# becomes to big to plot
# heatmap.2(as.matrix(prior_exp),Rowv=NA,Colv=NA,symm=TRUE,col=color_palette(256),scale="none",margins=c(8,8),trace = "none",dendrogram="none",density.info = "none")

# make inference with this GP
# cond_mean_of_pars = prior_mean_of_pars + Sigma_12 %*% Sigma_22 %*% (observation - prior_mean_of_exp)
# observation - experimental values
# prior_mean_of_pars = 0 - GP centered at 0
# prior_mean_of_exp = 0 - GP centered at 0
# Sigma_22 = prior_cov_mat at experiment + experiment_cov_mat
# Sigma_12 = 

Sigma_12 <- full_prior_cov_mat %*% t(model$jac)

#pars_mean <- Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC]),expDt[,DATA])
pars_mean <- Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC]),expDt[,RESIDUAL])

# note that it is the prior on the experimental energies that is acctually inverted so it should not matter much how large
# the energy grid of the model is?! <- It seems so, the calculation of the parameter covariance matrix is long though

pars_cov <- full_prior_cov_mat - Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC]),t(Sigma_12))

predictions <- model$jac %*% pars_mean
predictions_cov <- model$jac %*% pars_cov %*% t(model$jac)

expDt[,test_prediction:=as.vector(predictions)]
expDt[,test_prediction_unc:=sqrt(diag(predictions_cov))]


ggp <- ggplot(data=expDt) +
		geom_line(aes(x=L1,y=RESIDUAL)) +
		geom_line(aes(x=L1,y=test_prediction),col='red') +
		#geom_ribbon(aes(x=L1,ymin=test_prediction-test_prediction_unc, ymax=test_prediction+test_prediction_unc),fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		xlim(1,1.5)
ggp

model$parsDt[,POST:=as.vector(pars_mean)]
ggp <- ggplot(data=model$parsDt) +
		geom_line(aes(x=L1,y=POST),col='red') +
		#geom_ribbon(aes(x=L1,ymin=test_prediction-test_prediction_unc, ymax=test_prediction+test_prediction_unc),fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		xlim(1,1.5)
ggp

filepath <- file.path(getwd(), 'defect-model', 'gp-regression-Fe56.png')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

# create a smooth model that is not mapping to experimentla energies but another energy grid

new_energy_grid = seq(from=expDt[,min(L1)],to=expDt[,max(L1)], by = 1e-04)
fake_expDt <- data.table(
	REAC=rep(expDt[,unique(REAC)], each=length(new_energy_grid)),
	L1=rep(new_energy_grid, times=length(expDt[,unique(REAC)]))
	)

fake_expDt[,IDX:=seq_len(.N)]

smooth_model <- defect_model(energies, talys_calc_dir, fake_expDt)

smooth_predictions <- smooth_model$jac %*% pars_mean
smooth_predictions_cov <- smooth_model$jac %*% pars_cov %*% t(smooth_model$jac)

fake_expDt[,prediction:=as.vector(smooth_predictions)]
fake_expDt[,prediction_unc:=sqrt(diag(smooth_predictions_cov))]

L1_min <- 0
L1_max <- 20

ggp <- ggplot(data=expDt[L1<=L1_max & L1>=L1_min]) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_line(aes(x=L1,y=prediction), data=fake_expDt[L1<=L1_max & L1>=L1_min & prediction>=0], col='red') +
		geom_ribbon(aes(x=L1,ymin=prediction-prediction_unc, ymax=prediction+prediction_unc), data=fake_expDt[L1<=L1_max & L1>=L1_min & prediction>=0], fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=4, scales="free_y") 
ggp

filepath <- file.path(getwd(), 'defect-model', 'gp-regression-Fe56-smooth.png')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

# There appears to be a bug with open channels below threshold
# this seems to come from TALYS which for example predicts
# the (n,alpha) cross section up to about 2.6 MeV to be 1e-07 barn insted of 0