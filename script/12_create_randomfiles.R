#
# DESCRIPTION OF STEP
#
# Use the 2nd order Taylor approximation of
# the posterior to create a sample of 
# TALYS parameter sets. Extend the paramters
# to a new energy grid (necessary for the 
# energy dependent parameters)
# Perform the respective
# calculations and save the results.
#

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

scriptnr <- 12L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

extNeedsDt <- read_object(2, "extNeedsDt")
optParamDt <- read_object(10, "optParamDt")
needsDt <- read_object(1, "needsDt")
optGpDt <- read_object(6, "optGpDt")
Sexp <- read_object(10, "Sexp")
mask <- read_object(10, "mask")
optSysDt_allpars <- read_object(10, "optSysDt_allpars")
finalPars <- read_object(11, "finalPars")
finalParCovmat <- read_object(11, "finalParCovmat")

##################################################
#       START OF SCRIPT
##################################################
print("-----------------------------------------------------")
print("----------------------script 12----------------------")
print("-----------------------------------------------------")

# define objects to be returned
outputObjectNames <- c("allParsets", "allResults")
check_output_objects(scriptnr, outputObjectNames)

# create the data table of final parameters
finalParamDt <- copy(optParamDt)
finalParamDt[PARNAME %in% optSysDt_allpars$PARNAME, ADJUSTABLE:=TRUE]
setkey(finalParamDt, IDX)
finalParamDt[ADJUSTABLE == TRUE, POSTVAL := paramTrafo$fun(finalPars)]
# IMPORTANT NOTE: POSTUNC is still with respect to transformed parameters
finalParamDt[ADJUSTABLE == TRUE, POSTUNC := sqrt(diag(finalParCovmat))]

# extend the energy range of energy dependent parameter grid
# we use the same increment as the last increment in the
# existing grid for energy dependent parameters
nn <- length(energyGridForParams)
dE <- energyGridForParams[nn] - energyGridForParams[nn-1] 
Elow = energyGridForParams[nn] + dE
Eup = energyGridrandomFiles[length(energyGridrandomFiles)]
energyGridForParams_extension <- seq(Elow,Eup+dE-Eup%%dE,by=dE)
# This is probably overkill:
# (1) If the length scale of the GP is short we only need additional energy grid close to 
#     the end of the optimized interval (say 3*l), so we could construct an energy grid
#     which is dense up to Elow + 3l and puts a final value at Eup
# (2) If the length scale is long, a dense energy grid is not needed in order to interpolate.
#     So we could set the energy grid density (by parameter in seq) to be a fraction of l
#     So far we have used that l>=3 and the energy grid density is 2, so an energy grid density
#     of 2/3*l should suffice
#  => This means we need several energy grids (one per parameter)
#     which will complicate the script a bit

# retrieve the names of the parameters to extend
energy_dep_par_names <- unique(str_remove(optParamDt[grepl("\\(.+\\)",PARNAME)]$PARNAME,"\\(.+\\)"))
energy_dep_par_names <- str_split(energy_dep_par_names," ",simplify=TRUE) # col 1 = names, col 2 = projectiles

projectiles <- unique(energy_dep_par_names[,2])
energy_dep_par_names <- unique(energy_dep_par_names[,1])

# add the energy dependent parameter specifications
endepParname_extension <- expand.grid(energy_dep_par_names,'(',energyGridForParams_extension,') ', projectiles)
endepParname_extension <- do.call(paste0, endepParname_extension)
endepInpList_extension <- as.list(rep(1, length(endepParname_extension)))
names(endepInpList_extension) <- endepParname_extension

# create a data table with the extended parameters
extParamDt <- data.table(PROJECTILE = unique(optParamDt$PROJECTILE),
                      ELEMENT = unique(optParamDt$ELEMENT),
                      MASS = unique(optParamDt$MASS),
                      PARNAME = names(endepInpList_extension),
                      PARVAL = endepInpList_extension,
                      ADJUSTABLE = TRUE,
                      PARUNC = 0.5)

# merge the existing and extended parameter data tables
allParamDt <- rbind(finalParamDt,extParamDt,fill=TRUE)
allParamDt[,IDX := seq_len(.N)]

# create a data table of systematics for the extended parameters
# THIS PART IS NOT CORRECT YET
par_names <- optParamDt[grepl("\\(.+\\)",PARNAME)]$PARNAME

tmpPar <- unique(str_remove(par_names,"\\(.+"))
tmpProj <- unique(str_remove(par_names,".+\\) "))
extParnames <- expand.grid(tmpPar,'(',energyGridForParams_extension,') ', tmpProj) 
extParnames <- do.call(paste0, extParnames)
expIDs <- paste0("TALYS-",str_remove(extParnames,"\\(.+\\)"))
energies <- unlist(regmatches(extParnames,gregexpr("\\(.+\\)",extParnames)))
energies <- as.numeric(substring(energies,2,nchar(energies)-1))

extSysDt <- data.table( IDX=NA,
                        EXPID=expIDs,
                        DIDX=NA,
                        ERRTYPE="talyspar_endep",
                        GPTYPE="sqrexp",
                        DATA=1,
                        UNC=0.5,
                        REFDATA=1,
                        ORIGIDX=NA,
                        EN=energies,
                        PARNAME=extParnames
  )
extSysDt[,IDX := seq_len(.N)+nrow(optSysDt_allpars)]
extSysDt[,DIDX := seq_len(.N)+nrow(optSysDt_allpars)]

# merge the existing and extended systematics data tables
allSysDt <- rbind(optSysDt_allpars,extSysDt,fill=TRUE)

# create data table with needs specified by the new energy grid
extNeedsDt <- needsDt[,{
  stopifnot(all(L2 == 0) & all(L3 == 0))
  list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
       L2 = 0, L3 = 0, V1 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]

extNeedsDt[, IDX := seq_len(.N)]


# now consider also the insensitive parameters
# available for variations
allParamDt[PARNAME %in% allSysDt$PARNAME, ADJUSTABLE:=TRUE]

# set energy grid for random-files
parval <- as.list(allParamDt[, PARVAL])
parval[[1]] <- as.vector(energyGridrandomFiles)
allParamDt[, PARVAL:= parval]

# see step 07_tune_talyspars.R for more explanation
# about setting up the talys handler
talysHnds <- createTalysHandlers()
talys <- talysHnds$talysOptHnd
talys$setPars(allParamDt)
talys$setParTrafo(paramTrafo$fun, paramTrafo$jac)
talys$setNeeds(extNeedsDt)
talys$setSexp(Sexp) # not sure if this is correct!!!
talys$setMask(mask)
if(!exists("talys_finite_diff")) talys_finite_diff <- 0.01
talys$setEps(talys_finite_diff)

# set the seed for the random number generator
# to have a reproducible creation of TALYS parameter sets
set.seed(talysFilesSeed)

#############################################################################
# extend the finalPars and finalParCovmat to the extended parameter set
#############################################################################

# start by recreating the prior covariance matrix for all parameters
# including the extended energy range

gpHandler <- createSysCompGPHandler()
sysCompHandler <- createSysCompHandler()
sysCompHandler$addGPHandler(gpHandler)

print("Building the covariance matrices...")
# Prior covariance matrix including extended parameters
P <- sysCompHandler$cov(allSysDt, optGpDt, ret.mat = TRUE)
expSel <- allSysDt[, !grepl("TALYS-", EXPID)]
talysSel <- !expSel
P0_all <- P[talysSel, talysSel]

# update the block in the covariance matrix related to the optimized parameters
# this is the posterior covariance matrix
optpars_indices <- optSysDt_allpars[, sort(IDX)]
#P0_all[optpars_indices,optpars_indices] <- finalParCovmat
P0_opt_opt <- P0_all[optpars_indices,optpars_indices]

# block of the prior covariance matrix
# correlating optimized and extended parameters
P0_opt_ext <- P0_all[-optpars_indices,optpars_indices]

# block of the prior covariance matrix correlating
# related to the extended parameters
P0_ext_ext <- P0_all[-optpars_indices,-optpars_indices]


# update extended energy range parameters

# prior mean of optimized parameters
p0_opt <- unlist(finalParamDt[ADJUSTABLE==TRUE, PARVAL])

# prior mean of extra parameters
p0_extra <- unlist(extParamDt[, PARVAL])

# conditional covariance matrix for extra parameters (doesn't depend on the value of p_opt)
P1_ext_ext <- P0_ext_ext - P0_opt_ext %*% solve(P0_opt_opt,t(P0_opt_ext))

# calculate the inverse of P0_opt_opt, which will be used many times during the sampling
P0_opt_opt_inv <- solve(P0_opt_opt)

# symmetrize P1_ext_ext
stopifnot(max(abs(P1_ext_ext - t(P1_ext_ext))) < 1.e-09)
P1_ext_ext <- 0.5*(P1_ext_ext + t(P1_ext_ext))

print("...done!")

print("sampling the optimized parameters...")
# create samples of parameter sets
# first sample the optimized parameters from the posterior mean & covariance matrix
opt_par_samples <- sample_mvn(numTalysFiles,finalPars,finalParCovmat) # this takes some time (few seconds)
print("...done!")


# now sample the extended parameters from the prior MVN conditioned on each
# sample from the posterior
sample_extended_pars <- function(optParset)
{
  # this thing is exceedingly slow!!!
  # It may be faster to use this package : https://cran.r-project.org/web/packages/mvnfast/mvnfast.pdf
  # ext_par_mean <- p0_extra + P0_opt_ext %*% solve(P0_opt_opt,optParset-p0_opt)
  ext_par_mean <- p0_extra + P0_opt_ext %*% P0_opt_opt_inv %*% (optParset-p0_opt) # this should be a bit faster
  sample_mvn(1,ext_par_mean,P1_ext_ext)
}

print("sampling the extended parameters...")
# this step is very slow even for just one sample
ext_par_samples <- apply(opt_par_samples,2,sample_extended_pars)
print("...done!")

# get the mean of the sampled extended parameters
ext_pars_mean <- matrix(rowMeans(ext_par_samples),ncol=1)

# bind the things together so we can do the talys calculations
variedParsets <- rbind(opt_par_samples,ext_par_samples)
optParset <- rbind(finalPars,ext_pars_mean)
allParsets <- cbind(optParset, variedParsets)

# perform calculations and save the result
#talysHnds$remHnd$ssh$execBash(paste0("mkdir -p '", pathTalys, "'; echo endofcommand"))
dir.create(savePathTalys, showWarnings=TRUE)
print(paste0("Storing talys results in: ", savePathTalys))

print("performing talys calculations...")
allResults <- talys$fun(allParsets, applySexp = FALSE, ret.dt=FALSE, saveDir = savePathTalys)
print("...done!")

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
