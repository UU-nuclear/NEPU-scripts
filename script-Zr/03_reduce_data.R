#
# DESCRIPTION OF STEP
#
# In order to aviod overestimation of the cross sections channels
# affected by resonance structure, due to the fully correlated
# systematic uncertainties, the data is reduced in this step.
# This means that the original experimental data is averaged onto
# a new energy grid, closely related to the model calculation grid.
# not all channels are affected by this.
#
# The script replaces
# 03_extract_experimental_uncertainties.R.
# But executes it as a first step.
# Only after that step do we have access to the 
# statistcal uncertainties, which are needed in this step.
# 
# The current script modifies the subents object. This is not ideal
# and should be changed in the future. The subents hold all the
# relevant information from exfor and modifying this could lead
# to confusion since it will no longer hold the original experimental
# data. However, all mapping between experimental data and model
# calculations depend on the subents object. So, in the current version
# of the pipeline I see no way around this.
#

# execute 03_extract_experimental_uncertainties
#source("script/03_extract_experimental_uncertainties.R")
# and clear the memory
rm(list = ls())

library(plyr)
#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 3L
overwrite <- TRUE 

# overwrite is necessary since we need to change
# the expDt object, and also the subents

#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

# define the objects that should be returned
# we also store the original output of
# 03_extract_experimental_uncertainties.R,
# which is overwritten, under different names
outputObjectNames <- c("expDt", "modList", "modDt", "Smodexp",
                       "orig_expDt", "orig_modList", "orig_modDt", "orig_Smodexp")
#check_output_objects(scriptnr, outputObjectNames)

print("-----------------------------------------------------")
print("------------script 03_reduce_data.R------------------")
print("-----------------------------------------------------")

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

orig_subents <- read_object(1, "subents")
extNeedsDt <- read_object(2, "extNeedsDt")
orig_modList <- read_object(3, "modList")
orig_modDt <- read_object(3, "modDt") # default model prediction mapped to experimental energies
orig_expDt <- read_object(3, "expDt")
orig_Smodexp <- read_object(3, "Smodexp")

# Keep copies of the original versions of expDt and subents
expDt <- copy(orig_expDt)
subents <- copy(orig_subents)

# Define which reaction channels that needs this step
reaction_channels <- unique(expDt$REAC)
reaction_channels <- reaction_channels[
						grepl("(N,TOT)",reaction_channels) | 
						grepl("(N,INL)",reaction_channels)]
# a bit of a hack to only include the (n,tot) and (n,inl) channels

covMat_list <- list()
# loop over the relevant reaction channels
for(Reac in reaction_channels) {

	ExpIds <- orig_expDt[REAC == Reac, unique(EXPID)]
	# loop experiment by experiment
	for(curExpId in ExpIds) {
		curExpDt <- orig_expDt[EXPID==curExpId]
		curExpDt[, IDX := seq_len(.N)]
		if(nrow(curExpDt)<=1) {
			covMat_list[[curExpId]] = Diagonal(x=curExpDt$UNC^2)
			next # don't reduce data-sets with a single energy point
		}

		minE <- min(curExpDt$L1)
		maxE <- max(curExpDt$L1)

		grid_idcs <- which(energyGridrandomFiles>minE & energyGridrandomFiles<maxE)
		if(length(grid_idcs)>=nrow(curExpDt)) {
			covMat_list[[curExpId]] = Diagonal(x=curExpDt$UNC^2)
			next # don't reduce data-sets which is already more sparse than the computation grid
		}

		grid_idcs <- c(grid_idcs,max(grid_idcs)+1)

		suggested_energygrid <- energyGridrandomFiles[grid_idcs]

		#min_bin_width <- 0.5
		#grid_length <- max(2,floor((maxE-minE)/min_bin_width + 1))

		#suggested_energygrid <- seq(minE,maxE,length.out=grid_length)

		#points_per_bin <- count(cut(curExpDt$L1,breaks=suggested_energygrid,labels=FALSE))

		#E1 <- suggested_energygrid[1]
		E1 <- minE - 1e-06
		new_energyGrid <- c(E1)

		for(i in seq(2,length(suggested_energygrid))) {
			E2 <- suggested_energygrid[i]
			n_experiments <- nrow(curExpDt[L1>E1 & L1<=E2])
			if(n_experiments==1 & i>2) {
				Eexp <- curExpDt[L1>E1 & L1<=E2]$L1
				new_energyGrid <- c(new_energyGrid,Eexp)
				E1 <- Eexp
			} else if(n_experiments!=0) {
				new_energyGrid <- c(new_energyGrid,E2)
				E1 <- E2
			}
		}

		exp_energies <- curExpDt$L1
		n_exp_points <- length(exp_energies)
		n_interpolated <- length(new_energyGrid)
		SensitivityMatrix <- Matrix(nrow = n_exp_points, ncol = n_interpolated, data = 0, sparse = TRUE)

		for(j in seq(1,n_exp_points)) {
			E_pred <- exp_energies[j]
			for(i in seq(1,n_interpolated)) {
				if(i<length(new_energyGrid)) {
					E1 <- new_energyGrid[i]
					E2 <- new_energyGrid[i+1]
					if(E_pred>E1 & E_pred<=E2) {
						SensitivityMatrix[j,i] <- (E2-E_pred)/(E2-E1)
					}
				}
				if(i>1) {
					E1 <- new_energyGrid[i]
					E0 <- new_energyGrid[i-1]
					if(E_pred>E0 & E_pred<=E1) {
						SensitivityMatrix[j,i] <- (E_pred-E0)/(E1-E0)
					}
				}
			}
		}

		# calculate the reduced data, considering only the statistical uncertainties
		fact1 <- t(SensitivityMatrix) %*% solve(Diagonal(x=curExpDt$UNC^2)) %*% SensitivityMatrix
		fact2 <- t(SensitivityMatrix) %*% solve(Diagonal(x=curExpDt$UNC^2)) %*% getDt_DATA(curExpDt)

		cur_exp_reduced <- solve(fact1) %*% fact2
		cur_covMat <- solve(fact1)

		# then we add the variance of the mean
		#curFit_dt <- data.table(L1=new_energyGrid,DATA=as.vector(cur_exp_reduced))
		curExpDt[,FIT:=as.vector(SensitivityMatrix %*% as.vector(cur_exp_reduced))]
		curExpDt[,RESIDUAL:=DATA - FIT]

		# calculate the sample variance for each intervall on the energy grid
		curExpDt[,UNC_ORIG:=UNC]
		curExpDt[,BIN_NUM:=cut(L1,breaks=new_energyGrid,labels=FALSE)]

		#for(i in seq_len(100)) {
		for(i in seq_len(1)) {
			curExpDt[,NEWUNC:=sqrt(sum(RESIDUAL^2)/(.N)),by=BIN_NUM]
			curExpDt[,UNC:=pmax(UNC_ORIG,NEWUNC)]

			fact1 <- t(SensitivityMatrix) %*% solve(Diagonal(x=curExpDt$UNC^2)) %*% SensitivityMatrix
			fact2 <- t(SensitivityMatrix) %*% solve(Diagonal(x=curExpDt$UNC^2)) %*% getDt_DATA(curExpDt)

			new_exp_reduced <- solve(fact1) %*% fact2
			new_covMat <- solve(fact1)

			change <- new_exp_reduced - cur_exp_reduced
			chi2_of_change <- as.vector(t(change) %*% fact1 %*% change)
			#print(as.vector(chi2_of_change))
			if(chi2_of_change < 1.E-04) {
				break
			}
			if(i==100) {
				cat("too many iterations: chi2_of_change = ", chi2_of_change, "\n")
			}

			cur_exp_reduced <- new_exp_reduced
			cur_covMat <- new_covMat
			curExpDt[,FIT:=as.vector(SensitivityMatrix %*% as.vector(cur_exp_reduced))]
			curExpDt[,RESIDUAL:=DATA - FIT]

		}

		covMat_list[[curExpId]] = cur_covMat

		# remove the columns created while doing the fit
		curExpDt[,FIT:=NULL]
		curExpDt[,RESIDUAL:=NULL]
		curExpDt[,UNC_ORIG:=NULL]
		curExpDt[,BIN_NUM:=NULL]
		curExpDt[,NEWUNC:=NULL]

		# create a temporary data.table with the new xs at the reduced energy grid
		curExpDt <- curExpDt[1:length(cur_exp_reduced),]
		curExpDt[,L1:=new_energyGrid]
		curExpDt[,DATA:=as.vector(cur_exp_reduced)]
		curExpDt[,UNC:=sqrt(diag(cur_covMat))]
		curExpDt[,EXPID:=curExpId]
		curExpDt[,DIDX:=seq_len(.N)]

		# remove the original experimental data and replace with the reduced data set
		expDt <- expDt[EXPID!=curExpId]
		expDt <- rbindlist(list(expDt, curExpDt))

		# unfortunately it seems I need to modify the exfor subentries
		# in order for the rest of the pipeline to work
		subent_data_table <- data.table(
			"EN" = new_energyGrid,
			"DATA" = as.vector(cur_exp_reduced),
			"DATA-ERR" = sqrt(diag(cur_covMat))
			)
		subents[sapply(subents, function(y) curExpId %in% y)][[1]]$DATA$TABLE <- subent_data_table
	}
}

# Augment expDt with the model prediction
# this is easier said than done!
# in order to map the talys prediction to my modified experiments
# I need to have a subentry added in my exfor data base
# In order to not count the modified experiments twice I also need to remove the subents
# of the original data from my exfor database, this is done inside the loop above.

# IMPORTANT: the GLS fit to the energy grid results in a non-diagonal covariance matrix
# that sould be retained, it is no longer possible to sepparate "random" and "systematic"
# uncertainty in the same way as Georg was doing 

expDt[,IDX:= seq_len(.N)]
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
setkey(expDt, IDX)
setkey(extNeedsDt, IDX)
expDt[, DATAREF := as.vector(Sexp %*% extNeedsDt$V1)]
expDt[, REFDATA := 0]

# safeguard against very small and very large variations
# both may cause trouble in the GLS procedure
expDt[, UNC := pmin(1000, pmax(UNC, 1))]

# model predictions mapped to new expgrid

modList <- expDt[, list(SUBENT=list(createSubentStub(.BY[["REAC"]], energyGrid))), by="REAC"]
modDt <- exforHandler$extractData(modList$SUBENT, ret.values = FALSE)
modDt <- exforHandler$map(modDt, extNeedsDt, modList$SUBENT)
Smodexp <- exforHandler$getJac(modDt, extNeedsDt, modList$SUBENT)

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)

# also change the subents object from step01
save_output_objects(1L, c("subents","orig_subents"), TRUE)

ggp <- ggplot(orig_expDt[REAC=="(24-CR-52(N,INL)24-CR-52,,SIG)"]) +
	geom_point(aes(x=L1,y=DATA,col=EXPID)) + 
	geom_line(aes(x=L1,y=DATA,col=EXPID),data=expDt[REAC=="(24-CR-52(N,INL)24-CR-52,,SIG)"])

ggp <- ggplot(orig_expDt[REAC=="(24-CR-52(N,TOT),,SIG)"]) +
	geom_point(aes(x=L1,y=DATA,group=EXPID)) + 
	geom_line(aes(x=L1,y=DATA,group=EXPID),data=expDt[REAC=="(24-CR-52(N,TOT),,SIG)"],col='red') +
	geom_ribbon(aes(x=L1,ymin=DATA-UNC,ymax=DATA+UNC,group=EXPID),data=expDt[REAC=="(24-CR-52(N,TOT),,SIG)"],fill='red',alpha=0.3)