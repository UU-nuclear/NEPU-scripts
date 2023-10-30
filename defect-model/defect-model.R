defect_model <- function(energies, talys_pred, compound_A, compound_Z) {

	if(length(talys_pred)!=1) {
		cat("error in create_defect_model, length(talys_pred) != 1\n")
		return(list())
	}

	talys_pred <- talys_pred[[1]]

	#reaction_channels <- talys_pred$result[,unique(REAC)]
	#reaction_channels <- c(reaction_channels[reaction_channels=="CS/EL"],reaction_channels[grepl("REAC",reaction_channels)])
	#reaction_channels[reaction_channels=="CS/EL"] <- "CS/REAC/000000/TOT"

	# create a data.table with the cross section model-defect for the exlusive channels
	defectDt <- data.table(
			REAC = rep(talys_pred$result[,unique(REAC)], each=length(energies)),
			L1 = rep(energies, times=length(talys_pred$result[,unique(REAC)])),
			VAL = 0
		)

	pars <- defectDt[grepl("REAC|CS/EL",REAC),VAL]


	#paramDt <- defectDt[grepl("REAC|CS/EL",REAC)]

	this <- list(
		jac <- ,
		...
	)

	fun <- function(x) {
		# function to calculate all cross-sections from sum rules, given cross-sections
		# of all exclusive channels at each energy of the energy-grid

		#paramDt[,VAL:=x]

		# exclusive cross sections - the regular expression "REAC|CS/EL" gives all exclusive plus elastic cross-sections
		defectDt[grepl("REAC|CS/EL",REAC),VAL:=x]

		# total cross section
		defectDt[REAC=="CS/TOT",VAL:=defectDt[grepl("REAC|CS/EL",REAC),sum(VAL),by=L1]$V1]

		excl_reacs <- defectDt[grepl("REAC",REAC), unique(REAC)]
		particleIDs <- gsub("\\D","",str_extract(excl_reacs,"REAC/.+/TOT"))
		strsplit(particleIDs, "")

		particleA <- c(1,1,2,3,4,5) # the mass numbers for NPDTHA
		particleZ <- c(0,1,1,1,2,2) # the stomic numbers for NPDTHA

		# loop over the exclusive reactions and calculate all residuals that each produce,
		# store this in some good structure
		for(particleID in particleIDs) {
			strsplit(particleIDs, "")[[1]]

		}

		# This should be done in the reverse order:
		# looping over the exclusive cross-sections and producing all possible production cross-sections
		# this will be much more straight forward.
		# However, we need first to extract all exclusive cross sections from TALYS not only the ones
		# where there is experimental data

		for(reac in defectDt[grepl("PROD",REAC),unique(REAC)]) {
			particles <- gsub("\\D","",str_extract(reac,"PROD/.+/TOT")) # extract the number identifer of the particle identification string

			particles <- strsplit(particles, "")[[1]] # split the string
			for(position in seq_along(particles)) {
				print(particles[position])

				# now find exclusive cross-sections for each position and number
				# the exlusive reaction string is NPDTHA, where
				# N is number of neutron
				# P is number of protons
				# D is number of deuterons
				# T is number of tritons
				# H is number of He-3
				# A is number of alphas
				# the production reaction string is ZZZAAA (i think, check it)
				# for example PROD/002004 is alpha production, i.e. alpha + any other particle
				# so I should find all exclusive reactions that have an alpha in the outgoing channel
				# REAC/000001/TOT - exclusive alpha
				# REAC/000002/TOT - exclusive 2 alpha
				# REAC/100001/TOT - exclusive alpha + neutron
				# so all REACTION strings that do not have 0 in the A position
				# another example is PROD/001001 - proton production
				# so extract all cross-section with a number different from 0 in the proton position (2)
				defectDt[grepl("REAC",REAC),unique(REAC)]
			}

		}
		# production cross sections: todo
		# for particle production (proton, alpha, etc) it should be easy
		# e.g. proton production is just the sum over all exclusive channels
		# weighted by the number of protons in the outgoing channel.
		# for residual nucleus it is more cumbersome, I need also the target
		# nucleus and add up all channels that leads to a specific residual
		#for(channel in defectDt[grepl("/PROD/",REAC),unique(REAC)]) {
		#}

		defectDt[,VAL]
	}

	# calculate the Jacobian for this model
	Jacobian <- Matrix(rep(0,length.out = nrow(defectDt[grepl("REAC|CS/EL",REAC)])*nrow(defectDt)),nrow=nrow(defectDt))
	for(column in seq_len(defectDt[grepl("REAC|CS/EL",REAC),.N])) {
		#print(IDX1)
		pars <- rep(0,length.out=nrow(defectDt[grepl("REAC|CS/EL", REAC)]))
		pars[column] <- 1
		Jacobian[,column] <- fun(pars)
	}

	jac <- function(x) {
		this$jac
	}

	return(list())
}