#!/usr/bin/env Rscript
##################################################
#
#   GENERATE NEPU CONFIG FROM input.inp
#
#   usage:
#     Rscript make_config.R [input.inp] [output.R]
#
#   - reads KEY = VALUE pairs from input.inp
#     (full lines starting with '#' are comments;
#      '#' inside a value is kept verbatim, needed for
#      R expressions with trailing comments)
#   - the config template is embedded in this script
#     (no separate template file required)
#   - writes config-<ELEMENT><MASS>.R, or [output.R] if
#     supplied as the second argument
#   - REFUSES to overwrite: if the output file already
#     exists the script stops with an error
#
#   requires R >= 4.0  (for raw string literals)
#
##################################################

if (getRversion() < "4.0.0")
    stop("R >= 4.0.0 required (raw string literals used to embed template). ",
         "Current version: ", as.character(getRversion()))

args <- commandArgs(trailingOnly = TRUE)
inpFile <- if (length(args) >= 1) args[1] else "input.inp"
outArg  <- if (length(args) >= 2) args[2] else NULL

if (!file.exists(inpFile)) stop("input file not found: ", inpFile)

# ---- parse the input file -------------------------------------------------

inpLines <- readLines(inpFile, warn = FALSE)
vals <- list()
for (i in seq_along(inpLines)) {
    ln <- inpLines[i]
    s <- sub("^[ \t]+", "", ln)
    if (nchar(s) == 0L) next                  # empty line
    if (substr(s, 1, 1) == "#") next          # full-line comment
    eq <- regexpr("=", ln, fixed = TRUE)
    if (eq < 0L)
        stop("cannot parse line ", i, " of ", inpFile, ": ", ln)
    key <- gsub("^[ \t]+|[ \t]+$", "", substr(ln, 1, eq - 1L))
    val <- substr(ln, eq + 1L, nchar(ln))
    val <- gsub("^[ \t]+|[ \t]+$", "", val)   # trim, keep inner '#'
    if (nchar(key) == 0L)
        stop("empty key on line ", i, " of ", inpFile)
    vals[[key]] <- val
}

req <- function(key) {
    if (is.null(vals[[key]]))
        stop("missing required key in ", inpFile, ": ", key)
    vals[[key]]
}

# ---- derived quantities ---------------------------------------------------

element <- req("ELEMENT")
mass    <- req("MASS")
z       <- req("Z")
proj    <- req("PROJECTILE")

reacTriple <- paste0(z, "-", toupper(element), "-", mass)  # e.g. 40-ZR-90
projUC <- toupper(proj)

# comma list -> R vector literal with single quotes, e.g. c('v1','d1')
quoteVec <- function(csv) {
    parts <- gsub("^[ \t]+|[ \t]+$", "", strsplit(csv, ",", fixed = TRUE)[[1]])
    parts <- parts[nchar(parts) > 0L]
    if (length(parts) == 0L) stop("empty list value")
    paste0("c(", paste0("'", parts, "'", collapse = ","), ")")
}

# exclude_exfor_entries: empty -> commented example line, otherwise active line
excl <- if (is.null(vals[["EXCLUDE_EXFOR_ENTRIES"]])) "" else vals[["EXCLUDE_EXFOR_ENTRIES"]]
if (nchar(gsub("[ \t]", "", excl)) == 0L) {
    exclLine <- '#exclude_exfor_entries <- c("23313002", "23313003", "22433007") # (n,inel) measured at angles (wrong EXFOR classification)'
} else {
    ids <- gsub("^[ \t]+|[ \t]+$", "", strsplit(excl, ",", fixed = TRUE)[[1]])
    ids <- ids[nchar(ids) > 0L]
    exclLine <- paste0('exclude_exfor_entries <- c(',
                       paste0('"', ids, '"', collapse = ", "), ')')
}

# ---- placeholder table ----------------------------------------------------

repl <- c(
    "@WORKDIR@"                  = req("WORKDIR"),
    "@REAC_TRIPLE@"              = reacTriple,
    "@PROJ_UC@"                  = projUC,
    "@FEW_PARAMETER_TEST@"       = req("FEW_PARAMETER_TEST"),
    "@POLL_TIME@"                = req("POLL_TIME"),
    "@MONGO_DBNAME@"             = req("MONGO_DBNAME"),
    "@MONGO_COLNAME@"            = req("MONGO_COLNAME"),
    "@MIN_ENERGY@"               = req("MIN_ENERGY"),
    "@MAX_ENERGY@"               = req("MAX_ENERGY"),
    "@EXCLUDE_EXFOR_LINE@"       = exclLine,
    "@ENERGY_GRID_RANDOM_FILES@" = req("ENERGY_GRID_RANDOM_FILES"),
    "@DEFAULT_THRES_EN@"         = req("DEFAULT_THRES_EN"),
    "@PARAM_GRID_STRIDE@"        = req("PARAM_GRID_STRIDE"),
    "@ENDEP_PARS_C@"             = quoteVec(req("ENDEP_PARS")),
    "@ENDEP_PROJ_C@"             = quoteVec(req("ENDEP_PROJ")),
    "@TENDL_YEAR@"               = req("TENDL_YEAR"),
    "@TENDL_ELEMENT@"            = req("TENDL_ELEMENT"),
    "@TENDL_MASS@"               = req("TENDL_MASS"),
    "@TENDL_PROJECTILE@"         = req("TENDL_PROJECTILE"),
    "@TALYS_FINITE_DIFF@"        = req("TALYS_FINITE_DIFF"),
    "@TUNE_EXP_UNC_SEED@"        = req("TUNE_EXP_UNC_SEED"),
    "@MAXIT_LM@"                 = req("MAXIT_LM"),
    "@RELTOL_LM@"                = req("RELTOL_LM"),
    "@OUTDATA_DIR@"              = req("OUTDATA_DIR"),
    "@TALYS_FILES_SEED@"         = req("TALYS_FILES_SEED"),
    "@NUM_TALYS_FILES@"          = req("NUM_TALYS_FILES"),
    "@PATH_TALYS@"               = req("PATH_TALYS"),
    "@AUTO_SYS_UNC_FLOOR@"       = req("AUTO_SYS_UNC_FLOOR"),
    "@AUTO_STAT_UNC_FLOOR@"      = req("AUTO_STAT_UNC_FLOOR"),
    "@PPP_COMP@"                 = req("PPP_COMP"),
    "@ABS_ERR_MIN@"              = req("ABS_ERR_MIN"),
    "@ABS_ERR_MAX@"              = req("ABS_ERR_MAX"),
    "@REL_ERR_MIN@"              = req("REL_ERR_MIN"),
    "@REL_ERR_MAX@"              = req("REL_ERR_MAX")
)

# ---- output filename, refuse to overwrite ---------------------------------

outFile <- if (!is.null(outArg)) outArg else paste0("config-", element, mass, ".R")

if (file.exists(outFile))
    stop("output file already exists, refusing to overwrite: ", outFile,
         "\n  remove it first, or pass a different output name as the second argument.")

# ---- embedded template ----------------------------------------------------

template <- r"---[##################################################
#
#       CONFIGURATION OF PIPELINE
#
##################################################

# add a user library where we can install additional packages
# userLib <- "/TMC/alf/pipeline/eval-Cr-isotopes/NEPU-scripts/R-libs-user"
#userLib <- "/home/jinba562/NEPU-apptainer/NEPU-scripts/R-libs-user"
# .libPaths( c( .libPaths(), userLib) )


# working directory
workdir <- "@WORKDIR@"
setwd(workdir)

source("config/required_packages.R")
source("config/required_sourcefiles.R")

tmp_dir <- file.path("/dev/shm",Sys.getenv("SLURM_JOB_ID"))

createTalysHandlers <- function() {

    # Initialize the talysR mpi interface

    # Important note: 1) Scripts that will run talys under mpi should call this function at the very
    #                    begining of the script. The execution of this function, blocks any futher
    #                    execution of the R-scripts on the slave ranks (rank-id>0). The R code in the
    #                    scripts should be executed only on the main thread, so call this as early
    #                    as possible!
    #                 2) TMPDIR = "/dev/shm" is an important specification because /dev/shm usually
    #                    resides in main memory. TALYS produces many thousand files per run
    #                    and normal disks and shared file systems cannot deal with this load
    #                    so it is a good idea to store them in main memory.
    #                 3) maxNumCPU set the number of requested talys workers. The number is an upper
    #                    limit on the number of workers. If maxNumCPU=0 the number of workers will be
    #                    the number of availible workers as given by the MPI interface.
    runOpts <- list(TMPDIR = tmp_dir)
    talysHnd <- initTALYSmpi(runOpts = runOpts, maxNumCPU=0)

    # initialize an alternative TALYS handler
    talysOptHnd <- createTalysFun(talysHnd, TMPDIR=tmp_dir)

    # Difference between talysHnd and talysOptHnd:
    #   talysHnd is a lower-level interface that provides
    #            the functions run, isRunning, and result.
    #            The input specification is passed as a list
    #            with input keywords and values and the output
    #            specification as a datatable enumerating the
    #            observables of interest

    #   talysOptHnd provides the functions fun and jac which 
    #               take a vector x as input and return either
    #               a vector of observables (fun) or the Jacobian 
    #               matrix (jac). Default parameter values and
    #               which values are present in x is specified
    #               via additional setter functions. Functions
    #               provided by talysOptHnd rely on those 
    #               provided by talysHnd.

    list(talysHnd = talysHnd,
         talysOptHnd = talysOptHnd)
}

# specify the reaction(s) to extract data fromthe EXFOR data base
# target reaction strings matching this regular expression
reacPat <- "\\(@REAC_TRIPLE@\\(@PROJ_UC@,[^)]+\\)[^,]*,,SIG\\)"
#reacPat <- "\\(@REAC_TRIPLE@\\(@PROJ_UC@,(TOT|INL|P|2N|EL)\\)[^,]*,,SIG\\)"
# should pipeline be executed with a very 
# small number of adjustable model parameters
# for testing purposes
#Commenting and changing here JB
#fewParameterTest <- FALSE
fewParameterTest <- @FEW_PARAMETER_TEST@
# time interval to check for completed
# TALYS calculation in seconds
pollTime <- @POLL_TIME@

# settings to retrieve EXFOR entries from
# the MongoDb database
mongo_dbname <- "@MONGO_DBNAME@"
mongo_colname <- "@MONGO_COLNAME@"

# only use experimental data in that energy range
minExpEn <- @MIN_ENERGY@
maxExpEn <- @MAX_ENERGY@

# exclude exfor entries if needed
@EXCLUDE_EXFOR_LINE@

# Specify energy grid for the final random files created in step 9.
# The grid used during the fit is based on this one, but limited to the range of
# the experimental data. Therefore the following must hold
# energyGridrandomFiles[1] < minExpEn
# energyGridrandomFiles[length(energyGridrandomFiles)] > maxExpEn
# Note that for proper error propagation the energy dependent paramters should cover the same range.

energyGridrandomFiles <- @ENERGY_GRID_RANDOM_FILES@
# TESTING ONLY - much coarser grid
#energyGridrandomFiles <- c(seq(1.0, 20, by=2.0), seq(22, 50, by=5.0))
# moved the creation of the energy grid used in the fit to script 02, in order to limit it to where there is data


# default threshold energy for reaction channels
# if automatic determination fails
# (because of vanishing reaction cross section at all energies)
defaultThresEn <- @DEFAULT_THRES_EN@ 

# energy grid for energy-dependent TALYS parameters
# copy energyGridrandomFiles 
energyGridForParams <- energyGridrandomFiles
# take every 3rd point for the parameters
energyGridForParams <- energyGridForParams[seq(@PARAM_GRID_STRIDE@,length(energyGridForParams),by=@PARAM_GRID_STRIDE@)]
# limit at the maximum energy, including the energy grid-point above maxExpEn
energyGridForParams <- energyGridForParams[1:(max(which(energyGridForParams < maxExpEn))+1)]
# add the energy at "zero" (must be a positive energy, otherwise talys ignores it)
energyGridForParams <- c(1.e-06,energyGridForParams)

# specify which paramters to make energy dependent
# a data.frame object named enParDt containing the columns par and proj
# is used for this, for example
# > enParDt
#           par proj
# 1    v1adjust    n
# 2    d1adjust    n
tmpPar <- paste0(@ENDEP_PARS_C@,'adjust')
tmpProj <- @ENDEP_PROJ_C@
enParDt <- data.table(expand.grid(par = tmpPar, proj = tmpProj))
#enParDt <- enParDt[!(par=='rcadjust' & proj=='n')] # remove Coloumb radius for the neutron

# specification of the TALYS input file used as template
# param_template_path <- file.path(workdir,"indata/n_Fe_056.inp")
# I think that it should not really matter that this file is for Fe-56, only the parameters are extracted from the file
# the target and projectile are specified sepparately
# the input will be searched for in the indata directory, if not found there, it will be downloaded from
# https://tendl.imperial.ac.uk/
# use the following keywords to specify which nuclide and projectile
tendl_year <- @TENDL_YEAR@
tendl_element <- "@TENDL_ELEMENT@"
tendl_mass <- @TENDL_MASS@
tendl_projectile <- "@TENDL_PROJECTILE@"

#by jb
# FORCE use of local TALYS input file - must come AFTER tendl_ vars
#param_template_path <- file.path(workdir, "indata/n_Zr_090_2019.inp")

# Verify file exists and stop if it doesn't
#if (!file.exists(param_template_path)) {
#  stop(paste("TALYS input file not found at:", param_template_path))
#}
#cat("Using TALYS input file:", param_template_path, "\n")
##############################

# optional argument to set the finite difference used by talys to calculate the Jacobian
# default value is talys_finite_diff <- 0.01
talys_finite_diff <- @TALYS_FINITE_DIFF@

# random generator seed for optimization of experimental uncertainties
# impacts the initial extra uncertainties in the optimization setup
tuneExpUncSeed <- @TUNE_EXP_UNC_SEED@

# set up the handlers to map TALYS results to EXFOR entries
subentHandler <- createSubentHandler(createDefaultSubentHandlerList())
exforHandler <- createExforHandler(subentHandler)
# abuAgent <- createAbuAgent("talys/structure/abundance/")
# subentHandler$getHandlerByName("handler_ntot_nat")$configure(list(abuAgent = abuAgent))

# maximum number of iterations for Levenberg-Marquardt algorithm
maxitLM <- @MAXIT_LM@

# if the relative difference between subsequent iterations of the
# Levenberg-Marquardt algorithm falls below this value,
# the optimization procedure terminates
reltolLM <- @RELTOL_LM@

# where to save output data
outdataPath <- file.path(workdir, "@OUTDATA_DIR@")
dir.create(outdataPath, recursive=TRUE, showWarnings=FALSE)

# specify the directory were status information and plots during the 
# optimization using the Levenberg-Marquardt algorithm should be stored
savePathLM <- file.path(outdataPath, "/LMalgo")

# random seed to create TALYS randomfiles
talysFilesSeed <- @TALYS_FILES_SEED@

# number of TALYS randomfiles to be created
numTalysFiles <- @NUM_TALYS_FILES@

# where to store the TALYS results
# content of TALYS result directories is stored as tar archives
# needed for the creation of ENdf randomfiles using modified TASMAN
#pathTalys <- file.path(outdataPath, "random-files")

pathTalys <- "@PATH_TALYS@"
savePathTalys <- pathTalys

# where to save plots produced by the scripts in eval-fe56/script/visualization
plotPath <- file.path(outdataPath, '/plots')

autoSysUncFloor  <- @AUTO_SYS_UNC_FLOOR@
autoStatUncFloor <- @AUTO_STAT_UNC_FLOOR@
pppcomp <- "@PPP_COMP@"
absErrMin <- @ABS_ERR_MIN@
absErrMax <- @ABS_ERR_MAX@
relErrMin <- @REL_ERR_MIN@
relErrMax <- @REL_ERR_MAX@
]---"

# ---- substitute and write -------------------------------------------------

for (ph in names(repl)) {
    if (!grepl(ph, template, fixed = TRUE))
        stop("placeholder not found in embedded template: ", ph)
    template <- gsub(ph, repl[[ph]], template, fixed = TRUE)
}

leftover <- regmatches(template, gregexpr("@[A-Z_]+@", template))[[1]]
if (length(leftover) > 0L)
    stop("unreplaced placeholders in template: ",
         paste(unique(leftover), collapse = ", "))

writeLines(strsplit(template, "\n", fixed = TRUE)[[1]], outFile)
cat("wrote", outFile, "\n")
