
rm(list=ls())
gc()

source("config/config-Cr52-hetGP.R")
source("defect-model/defect-model.R")

expDt <- read_object(3, "expDt")
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc"

# test on reduced energy range
expDt <- expDt[L1>1]
setorder(expDt,REAC,L1)
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
energies <- expDt[,sort(L1)]
# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
#energies <- unique(round(energies*1e4))*1e-04

# a simple not very dense energy grid to test
energies <- seq(from=expDt[,min(L1)],to=expDt[,max(L1)],length.out=50)

model <- defect_model(energies, talys_calc_dir, expDt)

# set all exclusive cross-sections to 1
pars <- rep(1,nrow(model$parsDt))

# and map the experimental energies
predictions <- model$jac %*% pars

expDt[,test_prediction:=as.vector(predictions)]


ggp <- ggplot(data=expDt[REAC=="(24-CR-52(N,INL)24-CR-52,,SIG)"]) +
		geom_line(aes(x=L1,y=test_prediction))
		#geom_line(aes(x=L1,y=DATA),col='red')
ggp

ggp <- ggplot(data=expDt[REAC=="(24-CR-52(N,TOT),,SIG)"]) +
		geom_line(aes(x=L1,y=test_prediction)) #+
		#geom_line(aes(x=L1,y=DATA),col='red')
ggp

ggp <- ggplot(data=expDt[REAC=="(24-CR-52(N,INL)24-CR-52,,SIG)"]) +
		geom_line(aes(x=L1,y=DATA))
ggp

ggp <- ggplot(data=expDt[REAC=="(24-CR-52(N,X)2-HE-4,,SIG)"]) +
		geom_line(aes(x=L1,y=test_prediction)) +
		geom_line(aes(x=L1,y=DATA),col='red')
ggp

ggp <- ggplot(data=expDt) +
		geom_line(aes(x=L1,y=test_prediction)) +
		#geom_line(aes(x=L1,y=DATA),col='red')
		facet_wrap(~ REAC, ncol=4)
ggp

# now try to create a GP prior covariance matrix on the parameters
# of the model

sigma2 <- (600)^2
l2 <- (3)^2

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

pars_mean <- Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC]),expDt[,DATA])

pars_cov <- full_prior_cov_mat - Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC]),t(Sigma_12))

predictions <- model$jac %*% pars_mean
predictions_cov <- model$jac %*% pars_cov %*% t(model$jac)

expDt[,test_prediction:=as.vector(predictions)]
expDt[,test_prediction_unc:=sqrt(diag(predictions_cov))]


ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_line(aes(x=L1,y=test_prediction),col='red') +
		geom_ribbon(aes(x=L1,ymin=test_prediction-test_prediction_unc, ymax=test_prediction+test_prediction_unc),fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=4, scales="free_y")
ggp

filepath <- file.path(getwd(), 'defect-model', 'gp-regression.png')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)