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

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}

talysHnds <- createTalysHandlers()

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

print("-----------------------------------------------------")
print("----------------------script 02----------------------")
print("-----------------------------------------------------")

# define the objects that should be returned
outputObjectNames <- c("refInpList", "refParamDt", "extNeedsDt", "rawRes") 
check_output_objects(scriptnr, outputObjectNames) 

################################################################
# define the data structures containing the reference parameters
################################################################

tendl_mass_str <- sprintf("%03d",tendl_mass)
param_template_path <- paste0(workdir, "/indata/", tendl_projectile, "_", tendl_element, "_", tendl_mass_str, "_", tendl_year, ".inp")
if(!file.exists(param_template_path)) {
    # the file is not on disk: download from tendl
    # note: at present only works for tendl_projectile = n
    tendl_url <- paste0("https://tendl.web.psi.ch/tendl_", tendl_year, "/neutron_file/", tendl_element, "/", tendl_element, tendl_mass_str, "/input/talys.inp.0000")
    if(download.file(tendl_url, destfile=param_template_path, method="curl")) {
        stop(paste("not able to find parameter template on disk, nor at:",tendl_url))
    }
}
paramTemplate <- readLines(param_template_path)
paramTemplate <- removeLinesFromTalysTemplate(
    c("rescuefile", "best", "ompenergyfile"), paramTemplate
)

# create the energy grid
ExpEn_min <- min(expDt[,L1])
ExpEn_max <- max(expDt[,L1])

#energyGrid <- energyGridrandomFiles[energyGridrandomFiles>=minExpEn & energyGridrandomFiles<=maxExpEn]
tmp_idx <- which(energyGridrandomFiles>=ExpEn_min & energyGridrandomFiles<=ExpEn_max)
tmp_idx <- c(max(tmp_idx[1]-2,1),max(tmp_idx[1]-1,1),tmp_idx,tmp_idx[length(tmp_idx)]+1)
tmp_idx <- unique(tmp_idx)
energyGrid <- energyGridrandomFiles[tmp_idx]

refInpHeader <- list(projectile = tendl_projectile,
                  element = tendl_element,
                  mass = tendl_mass,
                  energy = energyGrid,
# FIXME: temporary to speed up calculation
                  #endf = "y",
                  outexcitation = "y",
                  channels = "y",
                  filechannels = "y",
                  filetotal = "y",
                  fileresidual = "y",
                  bins = 60)
#                 template = paramTemplate)

adjParList <- extractAdjustedTalysParameters(paramTemplate)
if (isTRUE(fewParameterTest)) {
    if(exists("max_number_of_free_pars")) {
        print(paste("limiting the number of free parameters to",max_number_of_free_pars))
        adjParList <- adjParList[1:max_number_of_free_pars]
    } else {
        print("limiting the number of free parameters to 3")
        adjParList <- adjParList[1:3]
    }
}

# to do that, first remove the associated global specifications
# enParDt <- expand.grid(par = tmpPar, proj = tmpProj)
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

# set the prior parameter uncertainties according to Nuclear Data Sheets 113 (2012) 2841–2934
# multiplied with a factor (here 2) and maximum prior uncertainty 0.5
library(stringi)
unc_multiplier <- 2.0
par_unc <- read.csv("indata/par_unc.txt",comment.char='#')
for(i in 1:nrow(par_unc)) {
    name <- as.character(par_unc$NAME[i])
    particle <- as.character(par_unc$PARTICLE[i])
    #unc <- par_unc$UNC[i]*0.01*unc_multiplier

    # set a minimum uncertainty on the prior parameter to 10%
    unc <- max(par_unc$UNC[i]*0.01*unc_multiplier,0.1)
    # the uncertainy should be smaller than the parameter transformation
    # unc <- min(unc,0.8) 
    if(is.na(particle)) {
        refParamDt[grepl(name,refParamDt$PARNAME),PARUNC := min(unc,0.5)]
    } else {
        refParamDt[grepl(name,refParamDt$PARNAME) & stri_sub(refParamDt$PARNAME,-1)==particle ,PARUNC := min(unc,0.5)]
    }
}

# index is needed for the sensitivity calculations later
refParamDt[, IDX := seq_len(.N)]

adjParNames <- refParamDt[ADJUSTABLE==TRUE,PARNAME]
for( i in 1:nrow(parRanges) ) {
  refParamDt[grepl(parRanges[i]$keyword,PARNAME),PARMIN:=parRanges[i]$min]
  refParamDt[grepl(parRanges[i]$keyword,PARNAME),PARMAX:=parRanges[i]$max]
}

# convert the parameter specification (specifiaclly the uncertainty) to the internal parameter space
paramTrafo <- parameterTransform(
                  x0 = unlist(refParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = refParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])

par_unc_int <- refParamDt[ADJUSTABLE==TRUE,PARUNC]/diag(paramTrafo$jac(refParamDt[ADJUSTABLE==TRUE,unlist(PARVAL)]))
par_val_int <- paramTrafo$invfun(refParamDt[ADJUSTABLE==TRUE,unlist(PARVAL)])

refParamDt[ADJUSTABLE==TRUE,PARUNC_EXT:=PARUNC]
#refParamDt[ADJUSTABLE==TRUE,PARUNC:=par_unc_int]
refParamDt[ADJUSTABLE==TRUE,PARUNC:=pmin(par_unc_int,unlist(PARVAL)-PARMIN)]
refParamDt[ADJUSTABLE==TRUE,PARVAL_EXT:=PARVAL]
refParamDt[ADJUSTABLE==TRUE,PARVAL:=as.list(par_val_int)]

print(refParamDt)

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

talysHnd <- talysHnds$talysHnd

startTime <- Sys.time()
cat("Started calculations at", as.character(startTime), "\n")
runObj <- talysHnd$run(list(refInpList), extNeedsDt)
stopTime <- Sys.time()
cat("Finished calculations at", as.character(stopTime), "\n")

exec_time <- as.double(stopTime-startTime,units="mins")

cat("talys execution time: ",exec_time,"\n")

nbr_free_pars <- nrow(refParamDt[ADJUSTABLE==TRUE])

cat("estimated coretime for full Jacobian: ", exec_time*nbr_free_pars/60, " corehours\n")
cat("number of free parameters: ", nbr_free_pars,"\n")
# save the information about the job for 
# later recovery if something goes wrong
save_output_objects(scriptnr, "runObj", overwrite)

rawRes <- talysHnd$result(runObj)
#talysHnds$clustHnd$closeCon()#

extNeedsDt <- rawRes[[1]]$result

# in order to use := extNeedsDt must be a data.table, seems that in the newer version of R
# it does not recognize that extNeedsDt is a data.table, forcing it with the following line
# solves the problem:
# Check that is.data.table(DT) == TRUE. Otherwise, := and `:=`(...) are defined for use in j, once only and in particular ways. See help(":=").
extNeedsDt <- data.table(extNeedsDt) ## Just add this line compared to the original one.
extNeedsDt[, IDX := seq_len(.N)]#

print(extNeedsDt)

## save the needed files for reference
save_output_objects(scriptnr, outputObjectNames, overwrite)
