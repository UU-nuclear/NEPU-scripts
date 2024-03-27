#################################################
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

# a script to test the idea of inflating the random uncertainty from resonance like data
# 
# We are trying to fit the average (smooth) cross section. The data in the low energy region
# is however not representative of this mean. They follow a different distribution, due to the resonance-
# like structure. Therefore, their weight in the fit should be representative of that these are samples
# from this distribution. If we do not do this the propagated uncertainty will be an underestimation of
# the true uncertainty.

# The simple way of doing this would be to calculate the standard deviation
# from a model (GLS fit) in each energy interval and inflate the uncertainty of the data points to this,
# this is what Helgesson did in his thesis.
# Here is a different way of achieving this by using a heteroschedastic gaussian process
# 1) First we "train" a GP to fit the default talys calculation, giving us a representative length scale of
#    the mean cross section.
# 2) We then fit a het. GP to the data using this length scale. The het. GP estimates the variance of the 
#    data around a mean (with the length scale extracted from talys)
# 3) This can then be used to replace the uncertainty of the experimental data points, to be represtative
#    of the distribution around the mean.
# So what we are doing is to estimate, using the GP, the distribution of data around a representative mean.
# We do not move the data points, but only inflate their uncertainty by considering them as samples from
# the determined distribution. This will then reflect in a larger uncertainty in the mean cross section once
# we fit the talys model


expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

expDt_upd_unc <- copy(expDt)

reac_channels <- expDt[,unique(REAC)]
for(reac in reac_channels) {
	print(reac)

	experiments <- expDt[reac==REAC,unique(EXPID)]
	for(experiment in experiments) {
		curExpDt <- expDt[EXPID==experiment]

		if(nrow(curExpDt)<=2) next

		# the detection of the abrupt change of binning is more difficult than I thought
		# I will leave it for now and treat all data the same

		# I select data sets to apply this correction to based on the minimum distance between
		# measured energies. The data set is considered "high resolution data" is the minimum
		# energy distance is smaller/equal Ec = 0.05 = 50 keV
		setorder(curExpDt,L1)
		energies <- curExpDt$L1
		nn <- length(energies)
		dE <- energies[2:nn]-energies[1:nn-1]
		dE <- c(dE[1],dE) # the distance for the first point is to the next point
		if(min(dE) > 0.05) next

		print(experiment)
		
	}
}

