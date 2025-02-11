
rm(list=ls())
gc()

source("config/config-Fe56.R")
source("defect-model/defect-model.R")

expDt <- read_object(3, "expDt")
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc-Fe56"

# test on reduced energy range
setorder(expDt,REAC,L1)
expDt <- expDt[L1<40]
expDt[,IDX:=seq_len(.N)]

# a simple not very dense grid
energies <- seq(from=expDt[,min(L1)],to=expDt[,max(L1)], by=0.3)
if(max(energies)<expDt[,max(L1)]) energies <- c(energies,expDt[,max(L1)]+1e-03)

model <- defect_model(energies, talys_calc_dir, expDt)

# get the default model parameters (talys prediction)
pars <- model$parsDt[,V1]

# and map the experimental energies
predictions <- model$jac %*% pars

expDt[,default_prediction:=as.vector(predictions)]
expDt[,RESIDUAL:=DATA-as.vector(predictions)]

ggp <- ggplot(data=expDt) +
		geom_point(aes(x=L1, y=DATA)) +
		geom_line(aes(x=L1,y=default_prediction),col='red') +
		facet_wrap(~ REAC, scales="free_y", ncol=4)
ggp

# now try to create a GP prior covariance matrix on the parameters
# of the model

sigma2 <- (1./3.)^2
l2 <- (5)^2

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
# each block in this matrix correspond to an exclusive reaction channel

# now propagate this to the experiments

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

pars_mean <- Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC^2]),expDt[,RESIDUAL])
# pars_mean_alt <- mult_invCov_x(expDt[,RESIDUAL],D=Diagonal(x=expDt[,UNC^2]),S=model$jac,P=full_prior_cov_mat)

# note that it is the prior on the experimental energies that is acctually inverted so it should not matter much how large
# the energy grid of the model is?! <- It seems so, the calculation of the parameter covariance matrix is long though

pars_cov <- full_prior_cov_mat - Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt[,UNC^2]),t(Sigma_12))


residual_predictions <- model$jac %*% pars_mean
predictions_cov <- model$jac %*% pars_cov %*% t(model$jac)

expDt[,RESIDUAL_GP:=as.vector(residual_predictions)]
expDt[,RESIDUAL_GP_UNC:=sqrt(diag(predictions_cov))]


ggp <- ggplot(data=expDt) +
		geom_line(aes(x=L1,y=DATA)) +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=2, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

# model$parsDt[,POST:=as.vector(pars_mean)]
# ggp <- ggplot(data=model$parsDt) +
# 		geom_line(aes(x=L1,y=POST),col='red') +
# 		#geom_ribbon(aes(x=L1,ymin=test_prediction-test_prediction_unc, ymax=test_prediction+test_prediction_unc),fill='red', alpha=0.3) +
# 		facet_wrap(~ REAC, ncol=1, scales="free_y") +
# 		xlim(1,1.5)
# ggp

filepath <- file.path(getwd(), 'defect-model', 'gp-regression-Fe56.png')
ggsave(filepath, ggp, width = 29.7, height = 21.0, units = "cm", dpi = 300)

####################################################################################################
# create a smooth model that is not mapping to experimentla energies but another energy grid

new_energy_grid = seq(from=expDt[,min(L1)],to=expDt[,max(L1)], length.out = 380)
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
		theme(axis.text=element_text(size=12),
        		axis.title=element_text(size=18)) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1,ymin=DATA-UNC,ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP), data=fake_expDt, col='red') +
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), data=fake_expDt, fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=5, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

ggp <- ggplot(data=expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"]) +
		theme(axis.text=element_text(size=12),
        		axis.title=element_text(size=18)) +
		geom_point(aes(x=L1,y=DATA)) +
		geom_errorbar(aes(x=L1,ymin=DATA-UNC,ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP), data=fake_expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"], col='red') +
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), data=fake_expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"], fill='red', alpha=0.3) +
		geom_line(aes(x=L1,y=default_prediction), data=fake_expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"], col='green') +
		facet_wrap(~ REAC, ncol=5, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)", color="Legend Title") +
		scale_color_manual(name="Legend Title", values=c("red", "green", "black"),
                     labels=c("Line 1", "Line 2", "Line 3"))
ggp


ggp <- ggplot(data=expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"]) +
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=18),
        legend.text=element_text(size=18)) +
  geom_point(aes(x=L1, y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC, color="exp. data")) +
  geom_line(aes(x=L1, y=default_prediction, color='default TALYS'), data=fake_expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"]) +
  geom_line(aes(x=L1, y=default_prediction+RESIDUAL_GP, color='GP-regression'), data=fake_expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"], linetype='dashed') +
  geom_ribbon(aes(x=L1, ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), data=fake_expDt[REAC!="(26-FE-56(N,N+P)25-MN-55,,SIG)"], alpha=0.3, fill='red') +
  facet_wrap(~ REAC, ncol=5, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)", color="Legend Title") +
  scale_color_manual(name="", values=c('GP-regression'='red', 'default TALYS'='green', 'exp. data'='black'))

print(ggp)

filepath <- file.path(getwd(), 'defect-model', 'gp-regression-Fe56-smooth-4.pdf')
ggsave(filepath, ggp, width = 29.7, height = (0.6)*21.0, units = "cm", dpi = 300)

# There appears to be a bug with open channels below threshold
# this seems to come from TALYS which for example predicts
# the (n,alpha) cross section up to about 2.6 MeV to be 1e-07 barn insted of 0
# acctually the (n,alpha) channel is open - it has positive Q-value

##########################################################################################
# now do a GP in the low energy range wher onyl few channels are open


expDt_LE <- copy(expDt[L1<3.4])
setorder(expDt_LE,REAC,L1)
expDt_LE[,IDX:=seq_len(.N)]


# get energies from all experiments
energies <- expDt_LE[REAC=="(26-FE-56(N,TOT),,SIG)",sort(L1)]
# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
energies <- unique(round(energies*1e4))*1e-04

model_LE <- defect_model(energies, talys_calc_dir, expDt_LE)

# get the default model parameters (talys prediction)
pars <- model_LE$parsDt[,V1]

# and map the experimental energies
predictions <- model_LE$jac %*% pars

expDt_LE[,default_prediction:=as.vector(predictions)]
expDt_LE[,RESIDUAL:=DATA-as.vector(predictions)]


ggp <- ggplot(data=expDt_LE) +
		geom_point(aes(x=L1, y=DATA)) +
		geom_line(aes(x=L1,y=default_prediction),col='red') +
		facet_wrap(~ REAC, scales="free_y", ncol=4)
ggp

sigma2 <- (1.)^2
l2 <- (5e-04)^2

cov_func <- function(x1,x2) {
	sigma2*exp(-0.5*(x1 - x2)^2/l2)
}

priors <- list()
for(reac in model_LE$parsDt[,unique(REAC)]) {
	model_energies <- model_LE$parsDt[REAC==reac,L1]
	model_default <- model_LE$parsDt[REAC==reac,V1]
	cov_mat <- outer(model_energies,model_energies,cov_func)

	scale <- outer(model_default,model_default)

	cov_mat <- Matrix(cov_mat*scale, sparse = TRUE)

	priors <- append(priors,list(cov_mat))
}

full_prior_cov_mat <- bdiag(priors)

prior_exp <- model_LE$jac %*% full_prior_cov_mat %*% t(model_LE$jac)

Sigma_12 <- full_prior_cov_mat %*% t(model_LE$jac)

pars_mean <- Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt_LE[,UNC^2]),expDt_LE[,RESIDUAL])
model_LE$parsDt[,CONDITIONAL:=as.vector(pars_mean)]

# note that it is the prior on the experimental energies that is acctually inverted so it should not matter much how large
# the energy grid of the model is?! <- It seems so, the calculation of the parameter covariance matrix is long though

pars_cov <- full_prior_cov_mat - Sigma_12 %*% solve(prior_exp + Diagonal(x=expDt_LE[,UNC^2]),t(Sigma_12))
model_LE$parsDt[,CONDITIONAL_UNC:=sqrt(diag(pars_cov))]

residual_predictions <- model_LE$jac %*% pars_mean
predictions_cov <- model_LE$jac %*% pars_cov %*% t(model_LE$jac)

expDt_LE[,RESIDUAL_GP:=as.vector(residual_predictions)]
expDt_LE[,RESIDUAL_GP_UNC:=sqrt(diag(predictions_cov))]


ggp <- ggplot(data=expDt_LE[L1<2]) +
		geom_point(aes(x=L1,y=DATA)) +
		#geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),col='red') +
		geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC), fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

# I'm limiting the plot to cross sections that are larger than 
ggp <- ggplot(data=model_LE$parsDt[L1<2]) +
		geom_line(aes(x=L1,y=V1 + CONDITIONAL),col='red') +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp

new_energy_grid = energies
fake_expDt <- data.table(
	REAC=rep(expDt_LE[,unique(REAC)], each=length(new_energy_grid)),
	L1=rep(new_energy_grid, times=length(expDt_LE[,unique(REAC)]))
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

ggp <- ggplot(data=expDt_LE) +
		geom_point(aes(x=L1,y=DATA)) +
		#geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
		geom_line(aes(x=L1,y=default_prediction),data=fake_expDt,col='green') +
		geom_line(aes(x=L1,y=default_prediction+RESIDUAL_GP),data=fake_expDt,col='red') +
		#geom_ribbon(aes(x=L1,ymin=(default_prediction+RESIDUAL_GP)-RESIDUAL_GP_UNC, ymax=(default_prediction+RESIDUAL_GP)+RESIDUAL_GP_UNC),data=fake_expDt , fill='red', alpha=0.3) +
		facet_wrap(~ REAC, ncol=1, scales="free_y") +
		labs(x="neutron energy (MeV)", y="cross section (mbarn)")
ggp


filepath <- file.path(getwd(), 'defect-model', 'gp-regression-Fe56-low-energy.pdf')
ggsave(filepath, ggp, width = 0.9*29.7, height = (0.6)*21.0, units = "cm", dpi = 300)
