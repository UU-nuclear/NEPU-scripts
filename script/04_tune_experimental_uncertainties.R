#
# DESCRIPTION OF STEP
#
# If the experimental data is inconsistent
# (taking into account all uncertainties),
# an additional systematic uncertainty is
# introduced on this step using maximum
# likelihood optimization. 
# The underlying model of the 'true' cross
# section function is a Gaussian process 
# that penalizes large values of the
# second derivative.
#

#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)


if (length(args)==0) {
  source("./config/config.R")
  stop("No config file supplied, using default file config.R", call.=FALSE)
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}


#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 4L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
modList <- read_object(3, "modList")
modDt <- read_object(3, "modDt")
expDt <- read_object(3, "expDt")
sysUncDt <- read_object(3, "sysUncDt")

##################################################
#       START OF SCRIPT
##################################################

# define objects to be returned
outputObjectNames <- c("origSysDt", "updSysDt")
check_output_objects(scriptnr, outputObjectNames)

# instantiate handlers 
reacHandler <- createSysCompReacHandler(c(subents, modList$SUBENT))
normHandler <- createSysCompNormHandler("DATAREF")

# define normalization error of experiments retrieved from EXFOR
# and with some rule based correction/penalization
# (sysUncDt$UNC will be enlarged in this step if necessary using MLO)
normHandler$addSysUnc("EXPID", sysUncDt$EXPID, 
                      0, sysUncDt$UNC, rel = (sysUncDt$ERRTYPE == "sys-rel"))

# define the second derivative prior on the reaction cross sections
reacHandler$addMap("pw", pwMap)

# We define for every reaction channel a prior 
# that incorporates the knowledge that cross sections
# are smooth functions. We achieve this by defining an
# independent Gaussian prior on the second derivatives of the 
# cross section. Because thresholds are different for different 
# reactions, we locate the threshold and cut away the energy  
# grid points below that value.
# Also the prior on the 2nd derivative of the cross section
# in 'curUncs' is modified depending on the reaction channel.
# As a rule of thumb: Reaction channels with larger cross sections
# should permit bigger changes in the 2nd derivative of the cross
# section.
# Meaning of elements in curUncs:
# 1st element: prior uncertainty of value at threshold
# 2nd element: prior uncertainty of the slope of the xs at threshold
# Other elements: prior uncertainty of the 2nd derivative of the xs

curReac <- "(26-FE-56(N,INL)26-FE-56,,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 0.1)
curUncs <- c(1e4, 1e4, rep(2, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

curReac <- "(26-FE-56(N,P)25-MN-56,,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 0.1)
curUncs <- c(1e4, 1e4, rep(10, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

curReac <- "(26-FE-56(N,2N)26-FE-55,,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 0.1)
curUncs <- c(1e4, 1e4, rep(0.5, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

curReac <- "(26-FE-56(N,TOT),,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 0.1)
curUncs <- c(1e4, 1e4, rep(2, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))


curReac <- "(26-FE-56(N,EL)26-FE-56,,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 0.1)
curUncs <- c(1e4, 1e4, rep(2, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs,
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))


curReac <- "(26-FE-56(N,A)24-CR-53,,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 1)
curUncs <- c(1e4, 1e4, rep(50, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs, 
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))


curReac <- "(26-FE-56(N,D)25-MN-55,,SIG)"
curEnGrid <- seq(getThresEn(curReac, modDt, defaultThresEn), 32, by = 1)
curUncs <- c(1e4, 1e4, rep(5, length(curEnGrid)-2))
reacHandler$assignMapToReac("pw", curReac,
                            vals = rep(0, length(curEnGrid)),
                            uncs = curUncs, 
                            opts = list(ens = curEnGrid,
                                        order = 3, outsideZero = TRUE))

# Setup of all uncertainties done
# Initialize master handler and
# construct the associated datatable with systematic components

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(reacHandler)
sysCompHandler$addHandler(normHandler)

sysDt <- sysCompHandler$createSysDt()
sysDt[, REFDATA := 0]  
origSysDt <- copy(sysDt)

# define for MLO optimization that we want to optimize the
# experimental normalization uncertainties and not the uncertainties
# on the second derivatives of the cross section prior
# nor the statistical uncertainties of the experiments

sysDt[, ADJUSTABLE := FALSE]
expDt[, ADJUSTABLE := FALSE]
sysDt[grepl("^EXPID-", EXPID), ADJUSTABLE := TRUE]

# Fix some experimental normalization uncertainties 
# to help MLO find a good solution

# (26-FE-56(N,INL)26-FE-56,,SIG)"
curExpId <- "23134005"
sysDt[paste0("EXPID-",curExpId) == EXPID, UNC := 0.005]
expDt[curExpId == EXPID, UNC := pmax(DATAREF * 0.05, 1)]
sysDt[paste0("EXPID-",curExpId) == EXPID, ADJUSTABLE := FALSE]

curExpId <- "23171003"
sysDt[paste0("EXPID-",curExpId) == EXPID, UNC := 0.005]
expDt[curExpId == EXPID, UNC := pmax(DATAREF * 0.05, 1)]
sysDt[paste0("EXPID-",curExpId) == EXPID, ADJUSTABLE := FALSE]

# start the MLO optimization for the individual channels

optfuns <- createMLOptimFuns()
setkey(sysDt, IDX)
sysDt[, ORIGIDX := IDX]

# store modified uncertainties in updSysDt
updSysDt <- copy(sysDt)

# set the seed for random number generator
# to have reproducible results
set.seed(tuneExpUncSeed)

for (curReac in unique(expDt$REAC)) {

    cat("CURRENT REACTION: ", curReac, "\n")

    curExpIds <- expDt[REAC == curReac, unique(EXPID)]

    # create expDt only containing experiments of
    # current reaction channel
    curExpDt <- expDt[EXPID %in% curExpIds,]
    curExpDt[, IDX := seq_len(.N)]

    # create sysDt only containing systematic errors of
    # experiments in current reaction channel
    curSysDt <- sysDt[EXPID %in% paste0("EXPID-",curExpIds) | 
                      grepl("^REACEXP", EXPID),]   
    curSysDt[, IDX := seq_len(.N)]

    # set up the current optimization problem
    # upper limit is 50% for relative uncertainty components and
    # 200 mBarns for absolute uncertainty components
    optfuns$setDts(curExpDt, curSysDt, sysCompHandler = sysCompHandler)
    lowerLims <- curSysDt[ADJUSTABLE == TRUE, UNC]
    upperLims <- curSysDt[ADJUSTABLE == TRUE, ifelse(ERRTYPE=="sys-rel", 0.5, 200)]
    initUncs <- lowerLims + runif(length(lowerLims)) * (upperLims-lowerLims)

    # print logLike associated with reference specification
    # following from rule based approach
    cat("logLike before ML optimization: ", optfuns$logLike(lowerLims), "\n")
    # do the bare-bone optimization here
    optRes <- optim(initUncs, optfuns$logLike, optfuns$gradLogLike, 
                    method = "L-BFGS-B", lower = lowerLims, upper = upperLims,
                    control = list(fnscale = -1))  # fnscale -1 to maximize instead of minimize 

    # print logLike after correction of systematic uncertainties
    cat("logLike after ML optimization: ", optRes$value, "\n")

    # save the modified systematic uncertainties
    updSysDt[J(curSysDt[ADJUSTABLE==TRUE, ORIGIDX]),
             UNC := curSysDt[ADJUSTABLE==TRUE, optRes$par]]
}

# Quick before/after comparison
# updSysDt$UNC - sysDt$UNC


# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)



