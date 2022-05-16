##################################################
#       CONVENIENCE FUNCTIONS NEEDED IN ALL STEPS
##################################################

check_output_objects <- function(scriptnr, objnames, overwrite=FALSE) {
    outpath <- file.path(outdataPath, sprintf("%02d", scriptnr))
    outfiles <- file.path(outpath, paste0(objnames, ".rda"))
    print(outfiles)
    dir.create(outpath, showWarnings=FALSE, recursive=TRUE)
    stopifnot(dir.exists(outpath))
    if (any(file.exists(outfiles))  && !overwrite)
    stop(paste0("AT STEP ", scriptnr, ": some of the result files ", paste0(outfiles, collapse=", "), " already exists\n"))
}

save_output_objects <- function(scriptnr, objnames, overwrite=FALSE) {
    outpath <- file.path(outdataPath, sprintf("%02d", scriptnr))
    for (objname in objnames) {
        curObj <- get(objname, envir=parent.frame()) 
        outfile <- file.path(outpath, paste0(objname, ".rda"))
        if (!file.exists(outfile) || overwrite)
            saveRDS(curObj, outfile)
    } 
}

save_objects <- function(scriptnr, objnames) {
  outpath <- file.path(outdataPath, sprintf("%02d", scriptnr))
  
  for (objname in objnames) {
    curObj <- get(objname, envir=parent.frame()) 
    outfile <- file.path(outpath, paste0(objname, ".rda"))
    
    #if (!file.exists(outfile) || overwrite)
    saveRDS(curObj, outfile)
  } 
}

read_object <- function(scriptnr, objname, outdata_path=outdataPath) {
  
  outpath <- file.path(outdata_path, sprintf("%02d", scriptnr))
  outfile <- file.path(outpath, paste0(objname, ".rda"))
  print(outfile)
  readRDS(outfile)
}



