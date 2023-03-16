#
# DESCRIPTION OF STEP
#
# Select parameters to include in the LM optimization
# this is the same step as the creation of the reference
# Jacobian, but is put in it's own script because the
# Jacobian computation is so resource intense. We do
# want to execute it more than once if something goes
# wrong in the parameter selection.
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

scriptnr <- 5
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

#subents <- read_object(1, "subents")
subents <- read_object(4, "fake_subents")
refParamDt <- read_object(2, "refParamDt")
extNeedsDt <- read_object(4, "extNeedsDt")
fullSensDt <- read_object(5, "fullSensDt") 
#expDt <- read_object(3, "expDt")
expDt <- read_object(4, "fake_expDt")

##################################################
#      define objects to be returned 
##################################################

outputObjectNames <- c("optParamDt","Sexp","mask")
check_output_objects(scriptnr, outputObjectNames)

##################################################
#      
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

# find all energy dependent parameters which have at least one point that is adjustbale
adjustable_endep_par_names <- optParamDt[ADJUSTABLE==TRUE]$PARNAME[grepl("\\(.+\\)",optParamDt[ADJUSTABLE==TRUE]$PARNAME)]
adjustable_endep_par_names <- unique(str_remove(adjustable_endep_par_names,"\\(.+\\)"))
optParamDt$tmp = str_remove(optParamDt$PARNAME,"\\(.+\\)")
optParamDt[tmp %in% adjustable_endep_par_names]$ADJUSTABLE=TRUE
optParamDt[,tmp:=NULL] # remove the temporary column from the data table

cat("number of adjustable parameters: ",sum(optParamDt[,ADJUSTABLE])," / ",sum(refParamDt[,ADJUSTABLE]),"\n")

save_output_objects(scriptnr, outputObjectNames, overwrite)
