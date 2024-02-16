# I will take as input a default talys calculation , done with the keywords
# channels = "y", filechannels="y"
# in order to produce all exclusive cross sections in separate files

#' A simple model of cross-section sum-rules 
#' 
#' @param energies a vector of numbers representing the energies at which the
#' cross-section should be provided.
#' @param talys_calc_dir a string representing the output directory of a TALYS
#' calculation of the cross-sections.
#' @param expDt a data.table object of experimental data as produced by the
#' package talysExforMapping
#' @returns A list consisting of a data.table object of model parameters and 
#' a sparse Jacobain matrix to map between the model parameters and the
#' experimental data in expDt.
#' 
defect_model <- function(energies, talys_calc_dir, expDt) {

	expDt <- copy(expDt) # create a copy to not change the input one by reference
	# We read results of a talys calculation to get the open exclusive reaction channels
	# and there thresholds

	linear_interpolation <- function(curDefectDt, curExpDt) {
		# x is the energy grid - defectDt[,L1]
		# xi is the inetrpolated energies - expDt[,L1]

		x <- curDefectDt[,L1]
		xi <- curExpDt[,L1]

		dx <- diff(x)
		if (any(dx == 0)) {
			stop("creation of interpolation_matrix: x values must be distinct")
		}

		# Find the indices for xi in x
		ix <- findInterval(xi, x, rightmost.closed=TRUE)

		# Calculate weights for interpolation
		wx <- (xi - x[ix]) / dx[ix]
		#wx <- pmax(pmin(wx, 1), 0)  # Ensure weights are between 0 and 1

		# each row has two entries
		# rows <- rep(seq_len(length(xi)), each = 2)
		rows <- rep(curExpDt[,IDX], each = 2)
		# the coloumns correspond to the indicies of x which are interpolated between
		cols <- as.vector(rbind(curDefectDt[ix,IDX], curDefectDt[ix+1,IDX]))
		# the values is the partial derivative of the interpolated function value
		# with respect to the y-value at the energy grid x
		values <- as.vector(rbind(1-wx, wx))

		# remove zeros
		rows <- rows[values!=0]
		cols <- cols[values!=0]
		values <- values[values!=0]

		data.table(
			rows = rows,
			cols = cols,
			values = values
		)
	}

	read_excl_xs <- function(xs_file) {
		dt <- as.data.table(read.table(file.path(talys_calc_dir, xs_file), col.names=c("E","xs","gamma xs","xs/res.prod xs"))[1:2])
		reac_id <- gsub("\\D", "", xs_file)
		reac_id <- paste0("CS/REAC/",reac_id,"/TOT")
		dt[,ID := reac_id]
	}

	talys_files <- list.files(talys_calc_dir)
	exclusive_xs <- talys_files[grepl("xs\\d{6}.tot",talys_files)]

	xs_tables <- lapply(exclusive_xs,FUN=read_excl_xs)
	exl_xs_dt <- rbindlist(xs_tables)

	# the sum of all of these at each energy is the total non-elastic xs
	# we need to add also the elastic xs to these
	elastic_xs_dt <- as.data.table(read.table(file.path(talys_calc_dir, "elastic.tot"), col.names=c("E","xs")))
	elastic_xs_dt[,ID:="CS/EL"]
	xs_dt <- rbind(exl_xs_dt,elastic_xs_dt)

	# create data.table for the model defect
	defectDt <- data.table(
		REAC = rep(xs_dt[,unique(ID)], each=length(energies)),
		L1 = rep(energies, times=length(xs_dt[,unique(ID)])),
		V1 = 0
	)

	# set numeric identifiers for the number of n, p, d, t, h, a in the outgoing channel
	defectDt[grepl("/REAC/",REAC),N:=as.numeric(substring(REAC,9,9))]
	defectDt[grepl("/REAC/",REAC),P:=as.numeric(substring(REAC,10,10))]
	defectDt[grepl("/REAC/",REAC),D:=as.numeric(substring(REAC,11,11))]
	defectDt[grepl("/REAC/",REAC),T:=as.numeric(substring(REAC,12,12))]
	defectDt[grepl("/REAC/",REAC),H:=as.numeric(substring(REAC,13,13))]
	defectDt[grepl("/REAC/",REAC),A:=as.numeric(substring(REAC,14,14))]

	# set N, P, D, T, H, A to zero for NA values
	cols_to_replace <- c("N", "P", "D", "T", "H", "A")
	defectDt[, (cols_to_replace) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = cols_to_replace]

	# remove parameters below the threshold of the reaction channel
	for(reac in defectDt[,unique(REAC)]) {
		xs_threshold <- xs_dt[ID==reac & xs==0, suppressWarnings(max(E))]
		threshold_idx <- findInterval(xs_threshold,defectDt[REAC==reac,L1])
		if(threshold_idx<defectDt[REAC==reac,.N]) {
			# if the larget energy in the model energygrid below the threshold
			# is the very last energy the channel is below the threshold for all
			# energies in the energy grid and should not be included in the model
			# at all.
			# without this if statement I allways get one free parameter for each channel
			# at the highest energygrid point
			defectDt[REAC==reac][seq(threshold_idx,defectDt[REAC==reac,.N])]$V1 <- 1
		}
	}
	defectDt <- defectDt[V1!=0,]

	# create an index column to use for the creation of the Jacobian
	defectDt[, IDX:=seq_len(.N)]

	# finally set the parameters from the default TALYS calculation
	for(reac in defectDt[,unique(REAC)]) {
		defectDt[REAC==reac,V1:=approx(x=xs_dt[ID==reac,E], y=xs_dt[ID==reac,xs], xout=L1)$y]
	}

	# now I can get out derived cross section in the following way
	# total cross-section
	# plot(defectDt[,sum(V1),by=L1])
	#
	# (n,x alpha) cross-section
	# select by at least one alpha in the outgoing channel
	# weight by the number of alphas in the outgoing channel
	# plot(defectDt[,sum(A*V1,na.rm=TRUE),by=L1])
	#
	# (n,x proton) cross-section
	# weight by the number of proton in the outgoing channel
	# plot(defectDt[,sum(P*V1,na.rm=TRUE),by=L1])
	#
	# (n,x neutron) cross-section
	# weight by the number of proton in the outgoing channel
	# plot(defectDt[,sum(N*V1,na.rm=TRUE),by=L1])

	# add the talys readction identifier to the expDt data.table
	generateTalysParticleStr <- function(processStr) {

	  if (length(processStr) != 1)
	    stop("processStr must be of length 1")

	  processStr <- toupper(processStr)
	  processStr <- gsub(" ","", processStr)
	  numXpart <- strsplit(processStr, "+", fixed = TRUE)[[1]]
	  pat <- "^([0-9]?)([NPDTA]|HE3)$"
	  if (all(grepl(pat, numXpart))) {

	    regRes <- regexec(pat, numXpart)
	    regStrs <- regmatches(numXpart, regRes)
	    particleStrs <- sapply(regStrs, function(x) x[3])
	    numStrs <- sapply(regStrs, function(x) x[2])
	    numAry <- rep(0, 6)
	    names(numAry) <- c("N","P","D","T","HE3","A")
	    numAry[particleStrs] <- ifelse (numStrs=="", 1, as.numeric(numStrs))
	    paste0(numAry, collapse="")
	  } else if(processStr=="X") {
	    "XXXXXX"
	  } else NULL
	}

	for(reac in expDt[,unique(REAC)]) {

		reacStruc <- parseReacExpr(reac)
		curProj <- reacStruc$projectile
		curElem <- reacStruc$target$sym
		curMass <- reacStruc$target$A
		curProc <- reacStruc$process

		particleStr <- generateTalysParticleStr(sub("INL", curProj, curProc))
		talysReacStr <- if (is.null(particleStr)) {
		  paste0("CS/", curProc)
		} else if(particleStr=="XXXXXX") {
		  particleStr <- paste0(sprintf("%03d", reacStruc$residual$Z),sprintf("%03d", reacStruc$residual$A))
		  paste0("CS/PROD/", particleStr, "/TOT")
		} else {
		  paste0("CS/REAC/", particleStr, "/TOT")
		}

		expDt[REAC==reac,REACID:=talysReacStr]
	}


	excl_reacs_exp <- expDt[REACID %in% defectDt[,unique(REAC)],unique(REACID)]
	incl_reacs_exp <- expDt[!(REACID %in% excl_reacs_exp),unique(REACID)]
	incl_reacs_exp <- incl_reacs_exp[incl_reacs_exp!="CS/TOT"]

	# loop over all talys reaction identifiers for exclusive channels in expDt
	JacobiansDt <- data.table()
	for(reacid in excl_reacs_exp) {
		curDefectDt <- defectDt[REAC==reacid]

		# strip off energies below threshold, they will be predicted as zero by the model
		minE <- curDefectDt[,min(L1)]
		maxE <- curDefectDt[,max(L1)]
		curExpDt <- expDt[REACID==reacid & L1>=minE & L1<=maxE]

		# create the data.table for the linear interpolation
		# Jacobian of the current channel

		JacobiansDt <- rbind(JacobiansDt,linear_interpolation(curDefectDt,curExpDt))
	}

	# production cross-sections
	for(reac_exp in incl_reacs_exp) {

		npdtha <- ""
		if(reac_exp=="CS/PROD/000001/TOT") { # neutron production
			this_defectDt <- defectDt[N>0]
			npdtha <- "N"
		} else if(reac_exp=="CS/PROD/001001/TOT") { # proton production
			this_defectDt <- defectDt[P>0]
			npdtha <- "P"
		} else if (reac_exp=="CS/PROD/001002/TOT") { # deuteron production
			this_defectDt <- defectDt[D>0]
			npdtha <- "D"
		} else if (reac_exp=="CS/PROD/001003/TOT") { # triton production
			this_defectDt <- defectDt[T>0]
			npdtha <- "T"
		} else if (reac_exp=="CS/PROD/002003/TOT") { # He-3 production
			this_defectDt <- defectDt[H>0]
			npdtha <- "H"
		} else if (reac_exp=="CS/PROD/002004/TOT") { # alpha production
			this_defectDt <- defectDt[A>0]
			npdtha <- "A"
		}

		if(npdtha=="") {
			cat("defect_model does not know how to map ",reac_exp,"\n")
			stop()
		}

		for(reacid in this_defectDt[,unique(REAC)]) { # loop to calculate block by block
			curDefectDt <- this_defectDt[REAC==reacid]
			multiplicity <- curDefectDt[,unique(get(npdtha))]
			# cat(reacid," : ", npdtha, ":", multiplicity, "\n")
			stopifnot(length(multiplicity)==1)

			if(multiplicity==0) {
				next
			}

			# strip off energies below threshold, they will be predicted as zero by the model
			minE <- curDefectDt[,min(L1)]
			maxE <- curDefectDt[,max(L1)]
			curExpDt <- expDt[REACID==reac_exp & L1>=minE & L1<=maxE]

			JacobiansDt <- rbind(JacobiansDt,linear_interpolation(curDefectDt,curExpDt))
			
		}
	}

	# the total cross-section
	if("CS/TOT" %in% expDt[,unique(REACID)]) {
		for(reacid in defectDt[,unique(REAC)]) { # loop to calculate block by block
			curDefectDt <- defectDt[REAC==reacid]
			# strip off energies below threshold, they will be predicted as zero by the model
			minE <- curDefectDt[,min(L1)]
			maxE <- curDefectDt[,max(L1)]
			curExpDt <- expDt[REACID=="CS/TOT" & L1>=minE & L1<=maxE]

			JacobiansDt <- rbind(JacobiansDt,linear_interpolation(curDefectDt,curExpDt))
		}
	}

	matrix_rows <- nrow(expDt)
	matrix_cols <- nrow(defectDt)

	jacobian_matrix <- sparseMatrix(
			i = JacobiansDt$rows,
			j = JacobiansDt$cols,
			x = JacobiansDt$values,
			dims = c(matrix_rows, matrix_cols)
		)

	return(list(parsDt=defectDt,jac=jacobian_matrix))
}