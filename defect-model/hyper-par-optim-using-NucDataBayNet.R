
rm(list=ls())
gc()

source("config/config-Fe56.R")
source("defect-model/defect-model.R")
source("res-struct-test/sysCompHandler_GP_sparsecov.R")

expDt <- read_object(3, "expDt")
talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc-Fe56"

# test on reduced energy range
setorder(expDt,REAC,L1)
# expDt <- expDt[L1<3]
expDt <- expDt[L1<1.3]
expDt[,IDX:=seq_len(.N)]

# get energies from all experiments
energies <- expDt[,sort(L1)]
# thin out the grid to have minimum 0.1 keV = 1e-04 MeV
energies <- unique(round(energies*1e4))*1e-04

model <- defect_model(energies, talys_calc_dir, expDt)

pars <- model$parsDt[,V1]

# and map the experimental energies
predictions <- model$jac %*% pars # this will be the default TALYS prediction

expDt[,default_prediction:=as.vector(predictions)]
expDt[,RESIDUAL:=DATA-default_prediction]

# now try to create a GP prior covariance matrix on the parameters
# of the model


gpHandler <- createSysCompGPHandler()