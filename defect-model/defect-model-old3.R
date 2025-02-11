# I will take as input a default talys calculation , done with the keywords
# channels = "y", filechannels="y"
# in order to produce all exclusive cross sections in separate files



defect_model <- function(energies, talys_calc_dir, expDt) {

	create_interpolation_matrix <- function(x, xi) {

		# first strip of the values in xi that are out of the range of x
		xi_inside <- xi[xi>=min(x) & xi<=max(x)]

		# count the number of striped values to the left of min(x)
		n_below <- length(xi[xi<min(x)])

		dx <- diff(x)
		if (any(dx == 0)) {
			stop("creation of interpolation_matrix: x values must be distinct")
		}

		# Find the indices for xi in x
		ix <- findInterval(xi_inside, x, rightmost.closed=TRUE)

		# Calculate weights for interpolation
		wx <- (xi_inside - x[ix]) / dx[ix]
		#wx <- pmax(pmin(wx, 1), 0)  # Ensure weights are between 0 and 1

		# Create the sparse interpolation matrix
		matrix_rows <- length(xi)
		matrix_cols <- length(x)

		# each row has two entries, shift the rows to the left by the number of points out of range
		rows <- rep(seq_len(length(xi_inside)), each = 2) + n_below
		# the coloumns correspond to the indicies of x which are interpolated between
		cols <- as.vector(rbind(ix, ix+1))
		# the values is the partial derivative of the interpolated function value
		# with respect to the y-value at the energy grid x
		values <- as.vector(rbind(1-wx, wx))

		interpolation_matrix <- sparseMatrix(
			i = rows,
			j = cols,
			x = values,
			dims = c(matrix_rows, matrix_cols)
		)

		return(interpolation_matrix)

	}

	# We read results of a talys calculatino to get the open exclusive reaction channels
	# and there thresholds

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
	# defectDt[grepl("/REAC/",REAC),N:=as.numeric(substring(REAC,9,9))]
	# defectDt[grepl("/REAC/",REAC),P:=as.numeric(substring(REAC,10,10))]
	# defectDt[grepl("/REAC/",REAC),D:=as.numeric(substring(REAC,11,11))]
	# defectDt[grepl("/REAC/",REAC),T:=as.numeric(substring(REAC,12,12))]
	# defectDt[grepl("/REAC/",REAC),H:=as.numeric(substring(REAC,13,13))]
	# defectDt[grepl("/REAC/",REAC),A:=as.numeric(substring(REAC,14,14))]

	# interpolate on the talys calculation to get out the thresholds of the reaction channels
	# I use constant interpolation with f=1, to get left-continuous approximation, i.e.
	# the highest energy on the talys grid with 0 cross-section will get a positive number
	for(reac in defectDt[,unique(REAC)]) {
		xs <- approx(xs_dt[ID==reac,E],xs_dt[ID==reac,xs],defectDt[REAC==reac,L1],method="constant",f=1)
		defectDt[REAC==reac,V1:=xs$y]
	}

	# then I can remove everything that has VAL==0 to be approximately under the threshold
	defectDt <- defectDt[V1!=0,]

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


	# start by creating linear interpolation matrices for all the exclusive reactions
	excl_channel_Jacs <- list()
	exp_energies <- expDt[,L1] # each model must be able to map to all experimental energies
	for(reac in defectDt[,unique(REAC)]) {
		model_energy_grid <- defectDt[REAC==reac,L1]

		Jac <- create_interpolation_matrix(model_energy_grid,exp_energies)
		excl_channel_Jacs <- append(excl_channel_Jacs, Jac)
	}
	names(excl_channel_Jacs) <- defectDt[,unique(REAC)]

	JacobiansDt <- data.table(
		REAC = defectDt[,unique(REAC)],
		JAC = excl_channel_Jacs
		)

	# set numeric identifiers for the number of n, p, d, t, h, a in the outgoing channel
	JacobiansDt[grepl("/REAC/",REAC),n:=as.numeric(substring(REAC,9,9))]
	JacobiansDt[grepl("/REAC/",REAC),p:=as.numeric(substring(REAC,10,10))]
	JacobiansDt[grepl("/REAC/",REAC),d:=as.numeric(substring(REAC,11,11))]
	JacobiansDt[grepl("/REAC/",REAC),t:=as.numeric(substring(REAC,12,12))]
	JacobiansDt[grepl("/REAC/",REAC),h:=as.numeric(substring(REAC,13,13))]
	JacobiansDt[grepl("/REAC/",REAC),a:=as.numeric(substring(REAC,14,14))]

	# set N, P, D, T, H, A to zero for NA values
	cols_to_replace <- c("n", "p", "d", "t", "h", "a")
	JacobiansDt[, (cols_to_replace) := lapply(.SD, function(x) ifelse(is.na(x), 0, x)), .SDcols = cols_to_replace]

	# to acctually get the matrix I need to do
	# JacobiansDt[channel=="CS/REAC/100000/TOT",Jacobian[[1]]]
	# note the [[1]]

	# now create the Jacobians for the inclusive reactions
	# based on the sum of the other ones
	# for the total cross-section the entry in the Jacobian
	# is the sum of all exclusive channels
	# This just becomes stacking all the individual Jacobains
	# horizontally, using cbind(J1,J2,...)
	# The same happens for inclusive channels like (n,xa) but for
	# exclusive channels that do not produce an outgoing alpha
	# its respective Jacobian should be replaced by a same sized
	# null-matrix

	incl_channel_Jacs <- list()
	exp_reacs <- expDt[,unique(REACID)]
	exp_reacs_incl <- exp_reacs[!(exp_reacs %in% defectDt[,unique(REAC)])]
	for(reac in exp_reacs_incl) {
		if(reac=="CS/TOT") {
			incl_channel_Jacs <- append(incl_channel_Jacs,do.call("cbind",JacobiansDt[,JAC]))
		} else if(grepl("CS/PROD/\\d{6}/TOT",reac)) { # production cross/section
			ZZZ <- as.numeric(substring(reac,9,11))
			AAA <- as.numeric(substring(reac,12,14))

			particle <- "unknown"

			if(AAA==1 & ZZZ==0) {
				particle <- "n"
			} else if(AAA==1 & ZZZ==1) {
				particle <- "p"
			} else if(AAA==2 & ZZZ==1) {
				particle <- "d"
			} else if(AAA==3 & ZZZ==1) {
				particle <- "t"
			} else if(AAA==3 & ZZZ==2) {
				particle <- "h"
			} else if(AAA==4 & ZZZ==2) {
				particle <- "a"
			}

			if(particle=="unknown") {
				cat("defect_model doesn't know how to map, ", reac,
					"!!!!\nThis may lead to undefined behaiour. Proceed with caution, or remove ",
					reac, " from expDt.\n")
			}

			# create a list of jacobians for each exclusive channel by
			# multiplying by the number of n,p,d,t,h or a in the outgoing channel
			jac_lst <- mapply('*',JacobiansDt[,JAC],as.vector(JacobiansDt[, particle, with = FALSE][[particle]]))
			# then create a Full Jacobian by stacking them
			incl_channel_Jacs <- append(incl_channel_Jacs,do.call("cbind",jac_lst))
		}
	}

	# Finally, build the full Jacobian matrix mapping between the parameters of the defect_model
	# and the experimental data points

	# first build a block matrix for the exclusive channels
	FullJacobian <- do.call("cbind", JacobiansDt)

	# then attach the inclusive channels

#	fun <- function(x) {
#		# function to calculate all cross-sections from sum rules, given cross-sections
#		# of all exclusive channels at each energy of the energy-grid
#		# x is a vector with all the exclusive + elastic cross section defects
#
#		defectDt[,V1:=x]
#
#		# the function must return the 'predicted' cross section in the correct order
#		# as given by the experimental data, this should be done by an Sexp sensitivity
#		# matrix as is done in the talys-wrapper
#
#
#	}



#	jac <- function(x) {
#		this$jac
#	}

	return(list(parsDt=defectDt,jac=jac))
}