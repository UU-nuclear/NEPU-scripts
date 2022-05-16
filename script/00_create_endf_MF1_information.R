#
# DESCRIPTION OF STEP
#
# Print pertinent information to a txt-file 
# explaining the pipeline and what data was 
# used. Also include paramter priors and 
# posterior values, with uncertainties.
#
library(stargazer)
#################################################
#       SCRIPT Setup
##################################################

args = commandArgs(trailingOnly=TRUE)
defaultConfig <- "config.R"


if (length(args)==0) {
  source(defaultConfig)
  print(paste0("No config file supplied, using default file: ", defaultConfig))
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}


#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 8L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

optExpDt <- read_object(6, "optExpDt")
optSysDt_optpars <- read_object(7, "optSysDt_optpars")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")
finalPars <- read_object(8, "finalPars")
finalParCovmat <- read_object(8, "finalParCovmat")
##################################################
#       START OF SCRIPT
##################################################

endfTxtPath <- file.path(getwd(), "endftxt") 
pilelineInfoPath <-file.path(endfTxtPath, "pipeline_information")
genTablesPath <-file.path(endfTxtPath, "generated_tables")
talystxtPath <-file.path(endfTxtPath, "talystxt")

paths <- list( endfTxtPath, endfTxtPath, genTablesPath, talystxtPath)

for (path in paths) {
  if(!dir.exists(path)){
    dir.create(file.path(path), showWarnings = FALSE, recursive = TRUE)
  }
}

##################################################
# Reactions table 
##################################################
reacNumTab <- file.path(genTablesPath, "reacNumTab.txt")

reacNumDt <- optExpDt[, .N, by="REAC"][order(-N)]
reacNumDt <- setNames(reacNumDt, c("EXFOR REACTION STRING", "NUM PTS"))
stargazer(reacNumDt, type = "text", digit.separator="",  summary=FALSE, rownames=FALSE, out = reacNumTab)

##################################################
#Access numbers
##################################################

exforAccesNumsTab <- file.path(genTablesPath, "exforAccesNumsTab.txt")

exforAccesNums <- optExpDt[, as.numeric(unique(EXPID))]
exforAccesNumsMat <- matrix(data=exforAccesNums, nrow=24, ncol = 4)


stargazer(exforAccesNumsMat, type = "text", digit.separator="",  summary=FALSE, rownames=FALSE, out = exforAccesNumsTab)


##################################################
# Parameter table
##################################################
paramtxt <- file.path(pilelineInfoPath, "paramtxt.txt")
paramTab <- file.path(genTablesPath, "paramTab.txt")

finalParamDt <- optSysDt_allpars
finalParamDt[,ADJUSTABLE:=FALSE]
optpars_indices <- optSysDt_optpars[, sort(IDX)]
finalParamDt[IDX %in% optpars_indices, ADJUSTABLE:=TRUE]
finalParamDt[, DATA := finalPars]
finalParamDt[, POSTUNC := sqrt(diag(finalParCovmat))]

paramDt <- finalParamDt[IDX > 2, 
                        data.table("Name"=unlist(PARNAME), 
                                   "Prior"=unlist(REFDATA),
                                   "PriorUnc"=unlist(UNC), 
                                   "Updt"=unlist(DATA), 
                                   "UpdtUnc"=unlist(POSTUNC), 
                                   "ADJUSTABLE"=unlist(ADJUSTABLE))
                        ]
cols <- c("Prior","PriorUnc","Updt","UpdtUnc")
paramDt[,(cols) := signif(.SD,4) , .SDcols=cols]
paramAdjDt <- paramDt[ADJUSTABLE==TRUE]
paramAdjDt <- paramDt[, "Adjust." := ifelse(ADJUSTABLE==TRUE, "LM", "POSTUPDT")]

stargazer(paramDt[ , !"ADJUSTABLE"], type = "text", summary=FALSE, rownames=FALSE, out = paramTab)



##################################################
# Concatenate all files 
##################################################

talystxt <- file.path(talystxtPath, "endf_n.txt") 
endftxt <- file.path(endfTxtPath, "endftxt.txt")
endftxtfold <- file.path(endfTxtPath, "endftxtfold.txt")
totaltxt <- file.path(pilelineInfoPath, "totaltxt.txt")

system(paste0("cp ", totaltxt," ", endftxt))

cmdStrInsert <- paste0("sed -e '/TABLEREAC/ {' -e 'r ", reacNumTab,"' -e 'd' -e '}' -i ", endftxt)
system(cmdStrInsert)

cmdStrInsert <- paste0("sed -e '/TABLEEXPID/ {' -e 'r ", exforAccesNumsTab,"' -e 'd' -e '}' -i ", endftxt)
system(cmdStrInsert)

cmdStrInsert <- paste0("sed -e '/TABLEPAR/ {' -e 'r ", paramTab,"' -e 'd' -e '}' -i ", endftxt)
system(cmdStrInsert)

cmdStrInsert <- paste0("sed -e '/TENDLINFO/ {' -e 'r ", talystxt,"' -e 'd' -e '}' -i ", endftxt)
system(cmdStrInsert)

cmdStrFold <- paste0("fold -w 67 -s ", endftxt, " > ", endftxtfold)
system(cmdStrFold)
