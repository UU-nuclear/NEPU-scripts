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

if(length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
} else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
}


#################################################
#       SCRIPT PARAMETERS
##################################################

scriptnr <- 0L
overwrite <- FALSE

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

optExpDt <- read_object(6, "optExpDt")
optSysDt_optpars <- read_object(7, "optSysDt_optpars")
optSysDt_allpars <- read_object(7, "optSysDt_allpars")
optParamDt <- read_object(10, "optParamDt")
finalPars <- read_object(11, "finalPars")
finalParCovmat <- read_object(11, "finalParCovmat")
finalParamDt <- read_object(11, "finalParamDt")
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

exforAccesNums <- optExpDt[, unique(EXPID)]
# pad exforAccesNums with " " so that its length is divisible by 4
exforAccesNums <- c(exforAccesNums,rep(" ", 4*ceiling(length(exforAccesNums)/4) - length(exforAccesNums)))
exforAccesNumsMat <- matrix(data=exforAccesNums, ncol = 4, byrow=TRUE)

stargazer(exforAccesNumsMat, type = "text", digit.separator="",  summary=FALSE, rownames=FALSE, out = exforAccesNumsTab)


##################################################
# Parameter table
##################################################
paramtxt <- file.path(pilelineInfoPath, "paramtxt.txt")
paramTab <- file.path(genTablesPath, "paramTab.txt")

UpdtUnc <- paste0("[+",signif(finalParamDt[ADJUSTABLE==TRUE,POSTUNC_UP],2),
                    ",-",signif(finalParamDt[ADJUSTABLE==TRUE,POSTUNC_Low],2),"]")

finalParamDt[,ADJUSTABLE:=optParamDt[,ADJUSTABLE]]

paramDt <- finalParamDt[IDX > 2, 
                        data.table("Name"=unlist(PARNAME), 
                                   "Prior"=unlist(PARVAL),
                                   "PriorUnc"=unlist(PARUNC), 
                                   "Updt"=unlist(POSTVAL), 
                                   "UpdtUnc"=UpdtUnc, 
                                   "ADJUSTABLE"=unlist(ADJUSTABLE))
                        ]
cols <- c("Prior","PriorUnc","Updt")
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
