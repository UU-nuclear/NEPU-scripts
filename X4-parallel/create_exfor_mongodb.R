#  The MIT License
#  
#  Copyright (c) 2019 Georg Schnabel 
#  
#  Permission is hereby granted, free of charge, 
#  to any person obtaining a copy of this software and 
#  associated documentation files (the "Software"), to 
#  deal in the Software without restriction, including 
#  without limitation the rights to use, copy, modify, 
#  merge, publish, distribute, sublicense, and/or sell 
#  copies of the Software, and to permit persons to whom 
#  the Software is furnished to do so, 
#  subject to the following conditions:
#  
#  The above copyright notice and this permission notice 
#  shall be included in all copies or substantial portions of the Software.
#  
#  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, 
#  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES 
#  OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. 
#  IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR 
#  ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, 
#  TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE 
#  SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#  


# load necessary packages and script modules

args = commandArgs(trailingOnly=TRUE)

cat("args = ", args, "\n")

if (length(args)==0) {
  stop("you need to provide the path to where the exfor data files are.", call.=FALSE)
}  else if (length(args) > 1) {
  stop("Script only accepts one argument.", call.=FALSE)
} else {
  print(paste0("looking for files in: ", args[1]))
}

library(mongolite)
library(data.table)
library(exforParser)

exforFiles <- list.files(args[1],"*\\.x4",recursive=TRUE,full.names=TRUE)

m <- mongo("entries", db=basename(dirname(args[1])), url="mongodb://localhost")

errorCounter <- 0
errorFiles <- character(0)

# loop over entries
for (idx in seq_along(exforFiles)) {

    cat("read file ", idx, " of ", length(exforFiles), "\n")
    curFile <- exforFiles[idx]
    curText <- try(readChar(curFile, file.info(curFile)$size), silent=TRUE)
    
    # remove carriage return '\r', present in *.x4 files downloaded from https://www-nds.iaea.org/exfor-master/x4toc4/
    curText <- gsub("\\\r","",curText) 

    if ("try-error" %in% class(curText)) 
    {
        cat("Problems reading ", curFile, "\n")
        errorCounter <- errorCounter + 1
        errorFiles <- c(errorFiles, curFile)
        next
    }
    curEntry <- parseEntry(curText)
    firstSub <- NULL

    # loop over subentries
    for (idx2 in seq_along(curEntry$SUBENT)) {
        curSub <- curEntry$SUBENT[[idx2]]
        if (idx2==1) 
        {
            firstSub <- curSub
        }
        else
        {
            curSub <- transformSubent(firstSub,curSub)
        }
        jsonObj <- convToJSON(curSub)
        m$insert(jsonObj)
    }
}
