#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 6L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(2, "extNeedsDt")
modList <- read_object(3, "modList")
expDt <- read_object(3, "expDt")
updSysDt <- read_object(4, "updSysDt")
fullSensDt <- read_object(5, "fullSensDt") 
refInpList <- read_object(2, "refInpList")

##################################################
# define objects to be returned
outputObjectNames <- c("endep_SensDt","fullSensDt","optParamDt")
check_output_objects(scriptnr, outputObjectNames)

##################################################
# convert the sparse matrix given as data.table 
# into a spase matrix type as defined in package Matrix

Spar <- with(fullSensDt,
             sparseMatrix(i = IDX1, j = IDX2, x = X,
                          dims = c(nrow(extNeedsDt), nrow(refParamDt))))
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
Sglob <- Sexp %*% Spar 

# safeguard
stopifnot(all(dim(Sglob) == c(nrow(expDt), nrow(refParamDt))))

# convert the sparse matrix Sglob into a datatable
SglobDt <- as.data.table(summary(Sglob))
setnames(SglobDt, c("IDX1", "IDX2", "X"))

# we (linearly) propagate all parameter values equal one
# to the model predictions
imp1 <- as.vector(Spar %*% rep(1, nrow(refParamDt)))
# we propagate hypothetical experimental values equal one
# to the model prediction
imp2 <- as.vector(t(Sexp) %*% rep(1, nrow(Sexp)))
# then we select observables on the model grid that
# are affected by both the backpropagation from the
# experiment and the forward propagation of model parameters
impIdx <- which(imp1 * imp2 != 0)

optSparDt <- copy(fullSensDt)
setkey(optSparDt, IDX1)
optSparDt <- optSparDt[J(impIdx)]

paramImpactDt <- SglobDt[, list(IMP = max(abs(X))), by = "IDX2"]
paramImpactDt <- paramImpactDt[order(IMP, decreasing = TRUE)]
selParIdcs <- paramImpactDt[IMP >= 1, IDX2]

setkey(optSparDt, IDX2)
mask <- optSparDt[J(selParIdcs), list(DSTIDX = IDX1, SRCIDX = IDX2)]
setkey(mask, SRCIDX, DSTIDX)
adjParIdcs <- unique(mask$SRCIDX)

# make a copy of the reference parameter datatable 
# and define what and what not we want to optimize
optParamDt <- copy(refParamDt)
setkey(optParamDt, IDX)
optParamDt[, ADJUSTABLE := FALSE]
optParamDt[J(adjParIdcs), ADJUSTABLE := TRUE]

# safeguard
stopifnot(sum(optParamDt$ADJUSTABLE) == length(adjParIdcs))

##########################################################
# Now add the Gaussian Process for parameters

# augment parameter list with energy dependent parameters
# these are specified in the config script
# This will be done later after the initial sensitivity evaluation

# to do that, first remove the associated global specifications
endepParamDt <- optParamDt[ADJUSTABLE==TRUE & PARNAME %in% with(enParDt,paste(par,proj))]
optParamDt <- optParamDt[!(ADJUSTABLE==TRUE & PARNAME %in% with(enParDt,paste(par,proj)))]
cat("number of energy dependent parameters: ", nrow(endepParamDt) , "\n")
cat("out of: ", nrow(enParDt), " potential ones\n")

tmp <- str_split(endepParamDt[,PARNAME]," ",simplify=TRUE)
# add the energy dependent parameter specifications
endepParnames <- paste0(rep(tmp[,1],each=length(energyGridForParams)),"(",energyGridForParams,") ",rep(tmp[,2],each=length(energyGridForParams)))

endepParamDt <- endepParamDt[rep(seq_len(nrow(endepParamDt)), each=length(energyGridForParams))]
endepParamDt[,PARNAME:=endepParnames]
endepParamDt[,IDX:=NA_integer_]

optParamDt <- rbind(optParamDt,endepParamDt)

optParamDt[,OLD_IDX:=IDX]
optParamDt[,IDX := seq_len(.N)]

##########################################################
# Make a sensitivity calculation for the added parameters

# create the mask to specify for which parameters we want to
# calculate the Jacobian
# Use the functionality in TALYSeval to create a default mask
# that takes all parameters that have the keyword
# ADUSTABLE==TRUE, supply only the newly create parameters
calculation_mask <- createDefaultMask(optParamDt[is.na(OLD_IDX)],extNeedsDt)

optParamDt[,ADJUSTABLE:=TRUE]
optParamDt[1:2,ADJUSTABLE:=FALSE]
# create the parameter transformation object
paramTrafo <- parameterTransform(
                  x0 = unlist(optParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = optParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])

# We let the talys_wrapper take care of creating the jacobian
talysHnds <- createTalysHandlers()
talys <- talysHnds$talysOptHnd

talys$setPars(optParamDt)
# the optimization will be undertaken using transformed
# TALYS parameters to make sure that parameters remain
# within given limits.
talys$setParTrafo(paramTrafo$fun, paramTrafo$jac)
# define the required predictions from TALYS in order to
# map to the grid of the experimental data given in EXFOR
talys$setNeeds(extNeedsDt)
# the finite difference used to calculate numerically the
# Jacobian matrix
if(!exists("talys_finite_diff")) talys_finite_diff <- 0.01
talys$setEps(talys_finite_diff)
# Since we have already calculated most of the parameter
# sensitivities we exclude these from this calculation
talys$setMask(calculation_mask)

cat("added number of parameters ",nrow(endepParamDt),"\n")
cat("Started calculations at", as.character(Sys.time()), "\n")
endep_SensDt <- talys$jac(unlist(optParamDt[ADJUSTABLE==TRUE,PARVAL]),returnDt=TRUE)
cat("Finished calculations at", as.character(Sys.time()), "\n")

# save the needed files for reference
save_output_objects(scriptnr, "endep_SensDt", overwrite)
#endep_SensDt <- read_object(6,"endep_SensDt")

# only keep the parameters that have more than one
# data sensitive energy point 
endep_SensDt <- endep_SensDt[duplicated(endep_SensDt,by="IDX2")]

# copy fullSensDt with the old indexing 
fullSensDt_old <- copy(fullSensDt)

# make a mapping between the old and new index
# and map the old and the new index in fullSensDt
IDXmapping <- optParamDt[,IDX,OLD_IDX]
setkey(fullSensDt,IDX2)
setkey(IDXmapping,OLD_IDX)

tmp <- IDXmapping[fullSensDt]
tmp[,OLD_IDX:=NULL]
names(tmp) <- c("IDX2","IDX1","X")

# where is.na(idx)==TRUE the parameter has been
# replaced by energy dependence
fullSensDt <- tmp[!is.na(IDX2)] 
rm(tmp)

# bind together the energy dependent and independent
# sensitivity data.tables
fullSensDt <- rbind(endep_SensDt,fullSensDt)

#################################################################
# Now we test the sensitivty again
Spar <- with(fullSensDt,
             sparseMatrix(i = IDX1, j = IDX2, x = X,
                          dims = c(nrow(extNeedsDt), nrow(optParamDt))))
Sexp <- exforHandler$getJac(expDt, extNeedsDt, subents)
Sglob <- Sexp %*% Spar 

# safeguard
stopifnot(all(dim(Sglob) == c(nrow(expDt), nrow(optParamDt))))

# convert the sparse matrix Sglob into a datatable
SglobDt <- as.data.table(summary(Sglob))
setnames(SglobDt, c("IDX1", "IDX2", "X"))

# we (linearly) propagate all parameter values equal one
# to the model predictions
imp1 <- as.vector(Spar %*% rep(1, nrow(optParamDt)))
# we propagate hypothetical experimental values equal one
# to the model prediction
imp2 <- as.vector(t(Sexp) %*% rep(1, nrow(Sexp)))
# then we select observables on the model grid that
# are affected by both the backpropagation from the
# experiment and the forward propagation of model parameters
impIdx <- which(imp1 * imp2 != 0)

optSparDt <- copy(fullSensDt)
setkey(optSparDt, IDX1)
optSparDt <- optSparDt[J(impIdx)]

paramImpactDt <- SglobDt[, list(IMP = max(abs(X))), by = "IDX2"]
paramImpactDt <- paramImpactDt[order(IMP, decreasing = TRUE)]
selParIdcs <- paramImpactDt[IMP >= 1, IDX2]

setkey(optSparDt, IDX2)
mask <- optSparDt[J(selParIdcs), list(DSTIDX = IDX1, SRCIDX = IDX2)]
setkey(mask, SRCIDX, DSTIDX)
adjParIdcs <- unique(mask$SRCIDX)

# make a copy of the reference parameter datatable 
# and define what and what not we want to optimize
optParamDt_old <- copy(optParamDt)
setkey(optParamDt, IDX)
optParamDt[, ADJUSTABLE := FALSE]
optParamDt[J(adjParIdcs), ADJUSTABLE := TRUE]

# safeguard
stopifnot(sum(optParamDt$ADJUSTABLE) == length(adjParIdcs))

# Voila! Now we have optParamDt as a data table with energy dependent
# parameters, and wether or not they should be adjustable

# save_output_objects(scriptnr, outputObjectNames, overwrite)