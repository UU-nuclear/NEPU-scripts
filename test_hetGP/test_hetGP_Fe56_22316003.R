
args = commandArgs(trailingOnly=TRUE)

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}

scriptnr <- 21L
overwrite <- FALSE

dir.create(file.path(outdataPath,scriptnr), recursive=TRUE, showWarnings=FALSE)

library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)

expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

experiment <- "22316003"
curExpDt <- expDt[EXPID==experiment]
#curExpDt <- curExpDt[L1<3] # just to test the script on my laptop

reac <- curExpDt[,unique(REAC)]
curModDt <- modDt[REAC==reac]


startTime <- Sys.time()
cat("Started calculations at", as.character(startTime), "\n")
#First the Matern5_2

#================= first talys ===========================
X <- matrix(curModDt$L1, ncol = 1)
Z <- curModDt$DATA
nvar <- 1
model_hetGP_talys_Mat5_2 <- mleHomGP(X = X, Z = Z, covtype = "Matern5_2")
mod_length_scale <- model_hetGP_talys_Mat5_2$theta

save_output_objects(scriptnr, "model_hetGP_talys_Mat5_2", overwrite)

Xgrid <- matrix(curModDt$L1, ncol=1)
pred_hetGP_talys_Mat5_2 <- as.data.table(predict(x = Xgrid, object = model_hetGP_talys_Mat5_2))
pred_hetGP_talys_Mat5_2[,L1:=Xgrid]

save_output_objects(scriptnr, "pred_hetGP_talys_Mat5_2", overwrite)

#==================== then experiment ======================================

X <- matrix(curExpDt$L1, ncol = 1)
Z <- curExpDt$DATA
nvar <- 1


settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training

# model the data with a heteroschedastic GP
model_hetGP_exp_Mat5_2 <- mleHetGP(X = X, Z = Z,
	covtype = "Matern5_2", settings = settings, known = list(theta=mod_length_scale))

save_output_objects(scriptnr, "model_hetGP_exp_Mat5_2", overwrite)

Xgrid <- matrix(curExpDt$L1, ncol=1)
pred_hetGP_exp_Mat5_2 <- as.data.table(predict(x = Xgrid, object = model_hetGP_exp_Mat5_2))
pred_hetGP_exp_Mat5_2[,L1:=Xgrid]

save_output_objects(scriptnr, "pred_hetGP_exp_Mat5_2", overwrite)

# Now the square exp kernel

#================= first talys ===========================

X <- matrix(curModDt$L1, ncol = 1)
Z <- curModDt$DATA
nvar <- 1
model_hetGP_talys_SqrExp <- mleHomGP(X = X, Z = Z, covtype = "Gaussian")
mod_length_scale <- model_hetGP_talys_SqrExp$theta

save_output_objects(scriptnr, "model_hetGP_talys_Mat5_2", overwrite)

Xgrid <- matrix(curModDt$L1, ncol=1)
pred_hetGP_talys_SqrExp <- as.data.table(predict(x = Xgrid, object = model_hetGP_talys_SqrExp))
pred_hetGP_talys_SqrExp[,L1:=Xgrid]

save_output_objects(scriptnr, "pred_hetGP_talys_Mat5_2", overwrite)

#==================== then experiment ======================================

X <- matrix(curExpDt$L1, ncol = 1)
Z <- curExpDt$DATA
nvar <- 1

settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training

# model the data with a heteroschedastic GP
model_hetGP_exp_SqrExp <- mleHetGP(X = X, Z = Z,
	covtype = "Gaussian", settings = settings, known = list(theta=mod_length_scale))

save_output_objects(scriptnr, "model_hetGP_exp_SqrExp", overwrite)

Xgrid <- matrix(curExpDt$L1, ncol=1)
pred_hetGP_exp_SqrExp <- as.data.table(predict(x = Xgrid, object = model_hetGP_exp_SqrExp))
pred_hetGP_exp_SqrExp[,L1:=Xgrid]

save_output_objects(scriptnr, "pred_hetGP_exp_SqrExp", overwrite)

stopTime <- Sys.time()
cat("Finished calculations at", as.character(stopTime), "\n")