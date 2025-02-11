# I will take as input a default talys calculation , done with the keywords
# channels = "y", filechannels="y"
# in order to produce all exclusive cross sections in separate files



defect_model <- function(energies, talys_calc_dir, expDt) {

	# needsDt should be the output of the exfor_Handler, as in step01.
	# i.e. it should have the actual experimental energies in the L1 column

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

#	projectile <- needsDt[,unique(PROJECTILE)]
#	element <- needsDt[,unique(ELEMENT)]
#	mass <- needsDt[,unique(MASS)]
	# create data.table for the model defect
	defectDt <- data.table(
		#PROJECTILE = projectile,
		#ELEMENT = element,
		#MASS = mass,
		REAC = rep(xs_dt[,unique(ID)], each=length(energies)),
		L1 = rep(energies, times=length(xs_dt[,unique(ID)])),
		#L2 = 0,
		#L3 = 0,
		V1 = 0
	)

	# set numeric identifiers for the number of n, p, d, t, h, a in the outgoing channel
	defectDt[grepl("/REAC/",REAC),N:=as.numeric(substring(REAC,9,9))]
	defectDt[grepl("/REAC/",REAC),P:=as.numeric(substring(REAC,10,10))]
	defectDt[grepl("/REAC/",REAC),D:=as.numeric(substring(REAC,11,11))]
	defectDt[grepl("/REAC/",REAC),T:=as.numeric(substring(REAC,12,12))]
	defectDt[grepl("/REAC/",REAC),H:=as.numeric(substring(REAC,13,13))]
	defectDt[grepl("/REAC/",REAC),A:=as.numeric(substring(REAC,14,14))]

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

	prediction <- function(j) {
		# compute the cross sections at the experiments by varying parameter j
		# returns at data table with the IDX of the paramter, the experimental data points and
		# the predicted cross section

		Emin <- ifelse(defectDt[j-1,REAC]==defectDt[j,REAC],defectDt[j-1,L1],defectDt[j,L1])
		Emax <- ifelse(defectDt[j+1,REAC]==defectDt[j,REAC],defectDt[j+1,L1],defectDt[j,L1])
		if(j==1) {
			Emin <- defectDt[j,L1]
		}
		if(j==nrow(defectDt)) {
			Emax <- defectDt[j,L1]
		}

		defectDt[j,V1:=1]

		reac <- defectDt[j,REAC]
		expDt_sub <- expDt[L1<=Emax & L1>=Emin]
		expDt_sub[REACID==reac, PRED:=approx(defectDt[REAC==reac & L1<=Emax & L1>=Emin,L1],defectDt[REAC==reac & L1<=Emax & L1>=Emin,V1],L1)$y]

#		reacs <- expDt_sub[,unique(REACID)]
#		incl_reacs <- reacs[!grepl("CS/REAC/\\d{6}/TOT|CS/EL",reacs)]
#		for(incl_reac in incl_reacs) {
#			# this could probably be done in a better way without looping and if statements
#			if(incl_reac=="CS/TOT") {
#				xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(V1),by=L1]
#				expDt_sub[REACID==incl_reac,PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]
#			} else if(incl_reac=="CS/PROD/001001/TOT") { # (nx,p)
#				xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(P*V1,na.rm=TRUE),by=L1]
#				expDt_sub[REACID==incl_reac,PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]
#			} else if(incl_reac=="CS/PROD/002004/TOT") { # (nx,alpha)
#				xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(A*V1,na.rm=TRUE),by=L1]
#				expDt_sub[REACID==incl_reac,PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]
#			}
#		}

		# total xs
		xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(V1),by=L1]
		expDt_sub[REACID=="CS/TOT",PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]

		# (nx,p) xs
		xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(P*V1,na.rm=TRUE),by=L1]
		expDt_sub[REACID=="CS/PROD/001001/TOT",PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]

		# (nx,alpha) xs
		xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(A*V1,na.rm=TRUE),by=L1]
		expDt_sub[REACID=="CS/PROD/002004/TOT",PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]

		defectDt[j,V1:=0]

		resultDt <- expDt_sub[PRED!=0,c("IDX","PRED")]
		resultDt[,IDX2:=defectDt[j,IDX]]

		return(resultDt)
	}

	# compute a data table to be used to compute the sparse Jacobian matrix
	setorder(defectDt,REAC,L1)
	defectDt[,IDX:=seq_len(.N)]
	defectDt[,V1:=0]
	expDt[,PRED:=0]
	pars <- rep(0,nrow(defectDt))
	
	jacDt <- Reduce(function(...) merge(..., all = TRUE), lapply(seq_along(pars),prediction))
	names(jacDt) <- c("IDX1","X","IDX2")

	jac <- sparseMatrix(i=jacDt[,IDX1],j=jacDt[,IDX2],x=jacDt[,X])
	rm(jacDt)
	# IDX1: index in expDt
	# IDX2: index in defectDt
	# X: Sensitivity, change in predicted cross-section at IDX1 due to a unit change in defectDt at IDX2

	# The computation of the Jacobian could probably be done using linearmap_piecewise.R in nucdataBaynet
	# for each exclusive cross-section there will be a block that is just a mapping matrix for a piecewise
	# linear.
	# for the total and production cross sections 

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