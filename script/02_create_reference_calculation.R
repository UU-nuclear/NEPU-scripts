#################################################
#  IMPORTANT NOTE
#  The cluster must be set up properly 
#  before the execution of this script.
#################################################

#
# DESCRIPTION OF STEP
#
# 1) create a parameter specification 'refInpList'
# 2) create an extended output specification 'extNeedsDt'
#    also containing results of reference calculation
#
# The same information as in 'refInpList' is also
# stored in 'refParamDt' as a datatable.
# Also the object 'rawRes' is stored which contains
# in addition to the results more details about the
# calculations, such as an excerpt of the main
# TALYS output file.
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

scriptnr <- 2L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEP
##################################################

subents <- read_object(1, "subents") 
expDt <- read_object(1, "expDt")
needsDt <- read_object(1, "needsDt")

##################################################
#       START OF SCRIPT
##################################################

# define the objects that should be returned
outputObjectNames <- c("refInpList", "refParamDt", "extNeedsDt", "rawRes") 
check_output_objects(scriptnr, outputObjectNames) 

################################################################
# define the data structures containing the reference parameters
################################################################

paramTemplate <- readLines(param_template_path)
paramTemplate <- removeLinesFromTalysTemplate(
    c("rescuefile", "best", "ompenergyfile"), paramTemplate
)

refInpHeader <- list(projectile = "n",
                  element = "Fe",
                  mass = 56L,
                  energy = energyGrid,
# FIXME: temporary to speed up calculation
                  endf = "y")
#                 template = paramTemplate)

adjParList <- extractAdjustedTalysParameters(paramTemplate)
if (isTRUE(fewParameterTest)) adjParList <- adjParList[1:3]

# augment parameter list with energy dependent parameters
tmpPar <- paste0(c('v1','d1','w1','vso1','wso1','rc'),'adjust')
tmpProj <- c('n','p','d','t','h','a')

if (isTRUE(fewParameterTest)) {
    tmpPar <- tmpPar[1:3]
    tmpProj <- tmpProj[1]
}

# to do that, first remove the associated global specifications
enParDt <- expand.grid(par = tmpPar, proj = tmpProj)
enParReg <- paste0(with(enParDt, paste0('^ *',par,' +',proj,' *$')), collapse='|')
idcsToRemove <- grep(enParReg, names(adjParList))
if (length(idcsToRemove) > 0)
    adjParList <- adjParList[-idcsToRemove]

# add the energy dependent parameter specifications
endepParnames <- expand.grid(tmpPar,'(',energyGridForParams,') ', tmpProj) 
endepParnames <- do.call(paste0, endepParnames)
endepInpList <- as.list(rep(1, length(endepParnames)))
names(endepInpList) <- endepParnames

# voila! Here is the complete reference parameter list
refInpList <- c(refInpHeader, adjParList, endepInpList)

# FIXME: functions such as createInputsForJacobian are picky
#        PROJECTILE, ELEMENT fields must contain uppercase strings
#        and MASS must be integer. This should be fixed in the future
refParamDt <- data.table(PROJECTILE = toupper(refInpList$projectile),
                      ELEMENT = toupper(refInpList$element),
                      MASS = as.integer(refInpList$mass),
                      PARNAME = names(refInpList),
                      PARVAL = refInpList)

# remove the redundant information about the reaction system
refParamDt <- refParamDt[! PARNAME %in% c("projectile", "element", "mass")]

# only parameters with "adjust" in the parameter name 
# are considered as adjustable
refParamDt[, ADJUSTABLE := grepl("adjust", PARNAME)]
refParamDt[ADJUSTABLE == TRUE, PARUNC := unlist(PARVAL)*0.1]

# index is needed for the sensitivity calculations later
refParamDt[, IDX := seq_len(.N)]

# Note: function convertToInput(paramDt, needsDt) produces a 
#                list of refInpList from a paramDt and needsDt 
#                specification which can even cover different
#                reaction systems (i.e. require several TALYS
#                calculations)

################################################################
# extend the energy points in needsDt to a regular grid
# for the calculation. We always can go back to the 
# experimental energies using linear interpolation
################################################################

extNeedsDt <- needsDt[,{
    stopifnot(all(L2 == 0) & all(L3 == 0))
    list(L1 = defineEnergyGrid(L1, energyGrid, enPolicy="compgrid"),
         L2 = 0, L3 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]

extNeedsDt[, IDX := seq_len(.N)]


################################################################
# perform the reference calculation and save the output
################################################################

talysHnds <- createTalysHandlers()
talysHnd <- talysHnds$talysHnd

runObj <- talysHnd$run(list(refInpList), extNeedsDt)
# save the information about the job for 
# later recovery if something goes wrong
save_output_objects(scriptnr, "runObj", overwrite)

cat("Started calculations at", as.character(Sys.time()), "\n")
cat("Waiting for termination...\n")
while (talysHnd$isRunning(runObj)) { 
    Sys.sleep(pollTime)
}
rawRes <- talysHnd$result(runObj)
talysHnds$clustHnd$closeCon()

extNeedsDt <- rawRes[[1]]$result 
extNeedsDt[, IDX := seq_len(.N)]

# save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
