LM#################################################
#       SCRIPT Setup
##################################################
args = commandArgs(trailingOnly=TRUE)

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}


#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 4L
overwrite <- FALSE

library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)

outputObjectNames <- c("expDt_upd_unc")
check_output_objects(scriptnr, outputObjectNames)

# a script to test the idea of inflating the random uncertainty from resonance like data
# 
# We are trying to fit the average (smooth) cross section. The data in the low energy region
# is however not representative of this mean. They follow a different distribution, due to the resonance-
# like structure. Therefore, their weight in the fit should be representative of that the experimental
# points are samples from this distribution.
# If we do not do this the propagated uncertainty will be an underestimation of the true uncertainty.
# Additionally, the resonance range constrain the average cross-section too much, leading to a poor fit/

# The simple way of doing this would be to calculate the standard deviation
# from a model (GLS fit) in each energy interval and inflate the uncertainty of the data points to this,
# this is what Helgesson did in his thesis.
# Here is a different way of achieving this by using a heteroschedastic gaussian process
# 1) First we "train" a GP to fit the default talys calculation, giving us a representative length scale of
#    the mean cross section.
# 2) We then fit a het. GP to the data using this length scale. The het. GP estimates the variance of the 
#    data around a mean (with the length scale extracted from talys)
# 3) This can then be used to replace the uncertainty of the experimental data points, to be representative
#    of the distribution around the mean.
# So what we are doing is to estimate, using the GP, the distribution of data around a representative mean.
# We do not move the data points, but only inflate their "uncertainty" by considering them as samples from
# the determined distribution. This will then reflect in a larger uncertainty in the mean cross section once
# we fit the talys model


expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

setorder(expDt,REAC,EXPID,L1)

expDt_upd_unc <- copy(expDt)
expDt_upd_unc[,ORIG_UNC:=UNC]

reac_channels <- expDt[,unique(REAC)]
for(reac in reac_channels) {
	curModDt <- modDt[REAC==reac]

	# train a GP on the talys model, we expect the length scale of this
	# GP to be representative of how fast we expect the average cross-section
	# to vary with energy, the same length scale will be used for modeling
	# the data in the next step
	X <- matrix(curModDt$L1, ncol = 1)
	Z <- curModDt$DATA
	nvar <- 1

	model_hetGP <- mleHomGP(X = X, Z = Z, covtype = "Matern5_2")

	mod_length_scale <- model_hetGP$theta

	experiments <- expDt[reac==REAC,unique(EXPID)]
	for(experiment in experiments) {
		curExpDt <- expDt[EXPID==experiment]

		if(nrow(curExpDt)<=2) next

		# the detection of the abrupt change of binning is more difficult than I thought
		# I will leave it for now and treat all data the same

		# I select data sets to apply this correction to based on the minimum distance between
		# measured energies. The data set is considered "high resolution data" is the minimum
		# energy distance is smaller/equal Ec = 0.05 = 50 keV
		# setorder(curExpDt,L1)
		energies <- sort(curExpDt$L1)
		nn <- length(energies)
		dE <- energies[2:nn]-energies[1:nn-1]
		dE <- c(dE[1],dE) # the distance for the first point is to the next point
		if(min(dE) > 0.05) next

		X <- matrix(curExpDt$L1, ncol = 1)
		Z <- curExpDt$DATA
		nvar <- 1

		settings <- list(return.hom = TRUE) # To keep homoskedastic model used for training

		# model the data with a heteroschedastic GP
		model_hetGP_exp <- mleHetGP(X = X, Z = Z,
			covtype = "Matern5_2", settings = settings, known = list(theta=model_hetGP$theta))

		# create a prediction, this is just for the visualization
		Xgrid <- matrix(curExpDt$L1, ncol=1)
		pred_exp <- as.data.table(predict(x = Xgrid, object = model_hetGP_exp))
		pred_exp[,L1:=Xgrid]

		ggp <- ggplot(data=curExpDt)
		ggp <- ggp + theme_bw() + theme(legend.position="none")
		ggp <- ggp + theme(axis.text=element_text(size=4),
		                   axis.title=element_text(size=4),
		                   strip.text=element_text(size=3))
		ggp <- ggp + geom_line(aes(x=L1, y=DATA))
		ggp <- ggp + geom_line(aes(x=L1, y=DATA),data=curModDt[L1<=max(curExpDt$L1) & L1>=min(curExpDt$L1)],col='red')
		ggp <- ggp + geom_ribbon(aes(x=L1, ymin=qnorm(0.05,mean,sqrt(sd2)), ymax=qnorm(0.95,mean,sqrt(sd2))),data=pred_exp, col='green', alpha=0., linetype = "dashed")
		ggp <- ggp + geom_ribbon(aes(x=L1, ymin=mean-sqrt(nugs), ymax=mean+sqrt(nugs)),data=pred_exp, col='blue', alpha=0., linetype = "dashed")
		ggp <- ggp + geom_line(aes(x=L1, y=mean),data=pred_exp, col='blue')

		# store the fitted nuggets as updated stat. unc. in expDt
		expDt_upd_unc[EXPID==experiment,UNC:=pmax(sqrt(pred_exp$nugs),ORIG_UNC)]

		# save a figure
		dir.create(plotPath, recursive=TRUE, showWarnings=FALSE)
		filepath <- file.path(plotPath, paste0('hetGPfit_',reac,'_',experiment,'.png'))
		ggsave(filepath, ggp, width = 16*0.75, height = 9*0.75, units = "cm", dpi = 300)

	}
}

# replace the expDt data.table from step 03 with the corrected one
expDt <- expDt_upd_unc
save_output_objects(3L, "expDt", TRUE)

save_output_objects(scriptnr, outputObjectNames, overwrite)
