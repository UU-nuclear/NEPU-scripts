#
# DESCRIPTION OF STEP
#
# 1) extract a list of subentries 'subents' to be 
#    used in the evaluation
# 2) extract a datatable 'expDt' summarizing and
#    indexing information in subents falling
#    in the energy range between 'minExpEn'
#    and 'maxExpEn'
# 3) create a datatable 'needsDt' containing the
#    information about required TALYS calculations
#    and the relevant output quantities
#

#################################################
#       SCRIPT Setup

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

scriptnr <- 1L
overwrite <- FALSE

#################################################
#       START OF SCRIPT
##################################################

print("-----------------------------------------------------")
print("----------------------script 01----------------------")
print("-----------------------------------------------------")

outputObjectNames <- c("subents", "expDt", "needsDt","expDtFull")
check_output_objects(scriptnr, outputObjectNames)

db <- connectExfor(mongo_colname, mongo_dbname, "mongodb://localhost")

# target reaction strings matching this regular expression

# reacPat <- "\\(26-FE-56\\(N,[^)]+\\)[^,]*,,SIG\\)"
# moved to config file

# construct the query string for mongodb

queryStr <- makeQueryStr(and(
    paste0('BIB.REACTION: { $regex: "', reacPat, '", $options: ""}'),
    paste0('DATA.TABLE.DATA: { $exists: true }'),
    paste0('DATA.TABLE.EN: { $exists: true }')
))

# create a datatable with an overview of the data
# filter out reaction expressions where:
#   * subexpressions are not cross sections (not SIG) 
#   * the reaction is not a character vector of length 1

reacOverviewDt <- db$find(queryStr, {
    if (is.list(BIB$REACTION) || length(BIB$REACTION) != 1)
       NULL
    else
    {
        reacExpr <- parseReacExpr(BIB$REACTION)
        unlistedReacExpr <- unlist(reacExpr, recursive = TRUE)
        if (any(unlistedReacExpr[grepl("quantspec", names(unlistedReacExpr))] != ",SIG"))
            NULL
        else
            list(ID = ID,
                REAC = reacStrucToStr(parseReacExpr(BIB$REACTION)),
                COLS = paste0(DATA$DESCR, collapse=","))
   } 
})

# retrieve the exfor subentries corresponding to
# the IDs in the summary datatable

idsStr <- paste0(paste0('"', reacOverviewDt$ID, '"'), collapse=",") 

queryStr <- makeQueryStr(
    paste0('ID: { $in: [', idsStr, ']}')                              
)

it <- exforIterator(queryStr)
subentList <- list()
while (!is.null((curSub <- it$getNext()))) {
    #if (curSub$ID %in% c("23313002", "23313003"))
    #    next # because measured at angles (wrong EXFOR classification)
    if (curSub$ID %in% c("13840005","12830008")) {
        # subentry 13840005 is obviously wrong, it comes from the same experiment as subentry 13840002,
        # but is, according to the entry, from a different run. It shows a completely different XS.
        # subentry 12830008 agrees with 13840005
        next
    }
    subentList <- c(subentList, list(curSub))
}

# just keep the subentries that can be mapped to
# TALYS predictions by package talysExforMapping
# and have valid uncertainty specifications

isMapable <- exforHandler$canMap(subentList, quiet=TRUE)
hasUncertainties <- hasValidUncertainties(subentList)

subents <- subentList[isMapable & hasUncertainties]
expDt <- exforHandler$extractData(subents, ret.values = TRUE)
expDtFull <- copy(expDt)
expDt <- expDt[L1 > minExpEn & L1 <= maxExpEn] do not introduce the cut-off energy
expDt <- expDt[!(duplicated(expDt$EXPID) &  duplicated(expDt$L1)),] # removing duplicates same energy from same experiment
expDt[, IDX := seq_len(.N)]

needsDt <- exforHandler$needs(expDt, subents)

# save the relevant objects for further processing
save_output_objects(scriptnr, outputObjectNames, overwrite)
