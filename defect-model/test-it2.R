
energy_step = 0.1
energies <- seq(from=extNeedsDt[,min(L1)], to=extNeedsDt[,max(L1)], by=energy_step)

talys_calc_dir <- "/home/alf/projects/NucDat/NAISS/UPPMAX/2023-22-58/eval-fe56-scripts/defect-model/test_calc"

talys_files <- list.files(talys_calc_dir)

exclusive_xs <- talys_files[grepl("xs\\d{6}.tot",talys_files)]

# read_excl_xs <- function(xs_file) {
# 	as.data.table(read.table(file.path(talys_calc_dir, xs_file), col.names=c("E",xs_file,"gamma xs","xs/res.prod xs"))[1:2])
# }
# 
# xs_tables <- lapply(exclusive_xs,FUN=read_excl_xs)
# 
# exclusive_xs_ids <- gsub("\\D", "", exclusive_xs)


########################
# merge(xs_tables[[1]],xs_tables[[2]],by="E")

# merge the exclusive xs data into a single data.table
# exl_xs_dt <- Reduce(function(...) merge(..., by="E"),xs_tables)

# ggp <- ggplot(data=exl_xs_dt) + geom_line(aes(x=E, y=DATA))

######################

read_excl_xs <- function(xs_file) {
	dt <- as.data.table(read.table(file.path(talys_calc_dir, xs_file), col.names=c("E","xs","gamma xs","xs/res.prod xs"))[1:2])
	reac_id <- gsub("\\D", "", xs_file)
	reac_id <- paste0("CS/REAC/",reac_id,"/TOT")
	dt[,ID := reac_id]
}

xs_tables <- lapply(exclusive_xs,FUN=read_excl_xs)
exl_xs_dt <- rbindlist(xs_tables)

# the sum of all of these at each energy is the total non-elastic xs
# we need to add also the elastic xs to these
elastic_xs_dt <- as.data.table(read.table(file.path(talys_calc_dir, "elastic.tot"), col.names=c("E","xs")))
elastic_xs_dt[,ID:="CS/EL"]
xs_dt <- rbind(exl_xs_dt,elastic_xs_dt)

ggplot(data=xs_dt[xs>0]) + geom_line(aes(x=E,y=xs,col=ID))

ggplot(data=xs_dt[A>0 & xs>0]) + geom_line(aes(x=E,y=xs,col=ID))

xs_dt[grepl("/REAC/",ID),N:=substring(ID,9,9)]
xs_dt[grepl("/REAC/",ID),P:=substring(ID,10,10)]
xs_dt[grepl("/REAC/",ID),D:=substring(ID,11,11)]
xs_dt[grepl("/REAC/",ID),T:=substring(ID,12,12)]
xs_dt[grepl("/REAC/",ID),H:=substring(ID,13,13)]
xs_dt[grepl("/REAC/",ID),A:=substring(ID,14,14)]

# test to see that it works
xs_tot_dt <- xs_dt[,sum(xs),by=E]
xs_tot_file <- as.data.table(read.table(file.path(talys_calc_dir, "totalxs.tot"), col.names=c("E","xs")))

ggplot() + geom_line(aes(x=E,y=V1),data=xs_tot_dt) + geom_line(aes(x=E,y=xs),data=xs_tot_file,color='red')

#########

excl_reacs <- exl_xs_dt[,unique(ID)]

# loop over all exclusive reaction channels and create production channels
for(reac in excl_reacs) {
	npdtha <- strsplit(gsub("\\D","",str_extract(reac,"REAC/.+/TOT")),"")[[1]]
	print(npdtha)
}


########################
extNeedsDt <- read_object(2, "extNeedsDt")
talys_pred <- read_object(2, "rawRes")
expDt <- read_object(3, "expDt")
yexp <- read_object(7, "yexp")

energy_step = 0.02
energies <- seq(from=extNeedsDt[,min(L1)], to=extNeedsDt[,max(L1)], by=energy_step)

defectDt <- data.table(
		REAC = rep(xs_dt[,unique(ID)], each=length(energies)),
		L1 = rep(energies, times=length(xs_dt[,unique(ID)])),
		VAL = 0
	)

defectDt[grepl("/REAC/",REAC),N:=as.numeric(substring(REAC,9,9))]
defectDt[grepl("/REAC/",REAC),P:=as.numeric(substring(REAC,10,10))]
defectDt[grepl("/REAC/",REAC),D:=as.numeric(substring(REAC,11,11))]
defectDt[grepl("/REAC/",REAC),T:=as.numeric(substring(REAC,12,12))]
defectDt[grepl("/REAC/",REAC),H:=as.numeric(substring(REAC,13,13))]
defectDt[grepl("/REAC/",REAC),A:=as.numeric(substring(REAC,14,14))]

defectDt[,unique(REAC)]

# interpolate on the talys calculation to get out the thresholds of the reaction channels
# I use constant interpolation with f=1, to get left-continuous approximation, i.e.
# the highest energy on the talys grid with 0 cross-section will get a positive number
# then I can remove everything that has the number 0 to be approximately under the threshold
for(reac in defectDt[,unique(REAC)]) {
	xs <- approx(xs_dt[ID==reac,E],xs_dt[ID==reac,xs],defectDt[REAC==reac,L1],method="constant",f=1)
	defectDt[REAC==reac,VAL:=xs$y]
}

defectDt <- defectDt[VAL!=0,]

ggplot(data=defectDt[A>0]) + geom_line(aes(x=L1,y=VAL,col=REAC))# + geom_line(aes(x=E,y=xs,col=ID),data=xs_dt)

# the total cross section can be calculated as
defectDt[,sum(VAL),by=L1]

# total cross-section
plot(defectDt[,sum(VAL),by=L1])

# (n,x alpha) cross-section
# select by at least one alpha in the outgoing channel
# weight by the number of alphas in the outgoing channel
plot(defectDt[,sum(A*VAL,na.rm=TRUE),by=L1])

# (n,x proton) cross-section
# weight by the number of proton in the outgoing channel
plot(defectDt[,sum(P*VAL,na.rm=TRUE),by=L1])

# (n,x neutron) cross-section
# weight by the number of proton in the outgoing channel
plot(defectDt[,sum(N*VAL,na.rm=TRUE),by=L1])

# (n,gamma) cross section
#plot(defectDt[REAC=="CS/REAC/000000/TOT",VAL,by=L1])

needsDt01 <- read_object(1, "needsDt")


##################

setkey(expDt, EXPID)
subentDt <- data.table(EXPID = sapply(subents, function(x) x$ID),
                       LISTPOS = seq_along(subents), key="EXPID")
subentDt <- unique(subentDt, by="EXPID")
joinedDt <- merge(expDt, subentDt)

spec <- joinedDt[, {
      subMat <- subentHandler$getJac(subents[[LISTPOS[1]]], needsDt, DIDX)
      df <- summary(subMat)
      df$i <- IDX[df$i]
      df[df$x > .Machine$double.eps*100, ]
    } ,by="EXPID"]

sparseMatrix(i=spec$i, j=spec$j, x=spec$x,
                               dims=c(nrow(expDt), nrow(needsDt)))


######################

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

exp_energies <- expDt[REAC==reac,L1]
mod_energies <- defectDt[REAC==talysReacStr,L1]

idx1 <- findInterval(exp_energies, mod_energies)
#stopifnot(idx1 > 0, idx1 <= length(modRes$L1))
isOnEdge <- idx1==length(mod_energies)
idx1[isOnEdge] <- idx1[isOnEdge] - 1  # special case if point on right edge

idx2 <- idx1 + 1
en1 <- mod_energies[idx1]
en2 <- mod_energies[idx2]
len <- en2 - en1

rowi <- rep(seq_along(exp_energies), 2)
colj <- c(idx1, idx2)  # colj is only with respect to the selected rows
# has to be converted to the global index in extout later
valx <- c((en2 - exp_energies) / len,
        (exp_energies - en1) / len)
rowi <- rowi[valx > .Machine$double.eps * 100]
colj <- colj[valx > .Machine$double.eps * 100]
valx <- valx[valx > .Machine$double.eps * 100]

sparseMatrix(i=rowi, j=defectDt[REAC==talysReacStr]$IDX[colj], x=valx,
                   dims=c(length(rowidcs), nrow(extout)))


pars <- rep(0,nrow(defectDt))
#pars <- rep(0,50)
IDX1 <- c()
IDX2 <- c()
X <- c()
start <- Sys.time()

expDt[,PRED:=0]
for(j in seq_len(50)) { # this for loop is extremely slow, could probably be done in a smarter way
	cur_pars <- pars
	cur_pars[j] <- 1
	prediction(cur_pars)
	predictions <- expDt[,PRED]

	IDXexp <- expDt[PRED!=0,IDX]
	VAL <- expDt[PRED!=0,PRED]
	IDXpar <- rep(j,length(IDXexp))

	IDX1 <- c(IDX1,IDXexp)
	IDX2 <- c(IDX2,IDXpar)
	X <- c(X,VAL)

	expDt[,PRED:=0]
}

stop <- Sys.time()

print(stop-start)

########################

pars <- rep(0,nrow(defectDt))
defectDt[,V1:=pars]
#pars <- rep(0,50)
IDX1 <- c()
IDX2 <- c()
X <- c()
start <- Sys.time()

for(j in seq_len(300)) { # this for loop is extremely slow, could probably be done in a smarter way
	prediction(j)
	predictions <- expDt[,PRED]

	IDXexp <- expDt[PRED!=0,IDX]
	VAL <- expDt[PRED!=0,PRED]
	IDXpar <- rep(j,length(IDXexp))

	IDX1 <- c(IDX1,IDXexp)
	IDX2 <- c(IDX2,IDXpar)
	X <- c(X,VAL)

	expDt[,PRED:=0]
}

stop <- Sys.time()

print(stop-start)

####################################

pars <- rep(0,nrow(defectDt))
defectDt[,V1:=pars]

result_dt <- NULL

start <- Sys.time()
for(j in seq_len(300)) { # this for loop is extremely slow, could probably be done in a smarter way
	merge(result_dt,prediction_alt(j))
}
stop <- Sys.time()
print(stop-start)


######################################

pars <- rep(0,nrow(defectDt))
defectDt[,V1:=pars]

res <- lapply(seq_len(300),prediction_alt)

Reduce(merge, res)

########33
setorder(defectDt,REAC,L1)
defectDt[,IDX:=seq_len(.N)]
expDt[,PRED:=0]
pars <- rep(0,nrow(defectDt))
defectDt[,V1:=pars]
start <- Sys.time()
resDt <- Reduce(function(...) merge(..., all = TRUE), lapply(seq_len(300),prediction_alt))
stop <- Sys.time()
print(stop-start)

names(resDt) <- c("IDX1","X","IDX2")

jac <- sparseMatrix(i=resDt[,IDX1],j=resDt[,IDX2],x=resDt[,X])
####

prediction <- function(j) {
		# compute the cross sections at the experiments by varying parameter j
		defectDt[j,V1:=1]

		Emin <- ifelse(defectDt[j-1,REAC]==defectDt[j,REAC],defectDt[j-1,L1],defectDt[j,L1])
		Emax <- ifelse(defectDt[j+1,REAC]==defectDt[j,REAC],defectDt[j+1,L1],defectDt[j,L1])
		if(j==1) {
			Emin <- defectDt[j,L1]
		}
		if(j==nrow(defectDt)) {
			Emax <- defectDt[j,L1]
		}

		reac <- defectDt[j,REAC]
		expDt[REACID==reac & L1<=Emax & L1>=Emin, PRED:=approx(defectDt[REAC==reac & L1<=Emax & L1>=Emin,L1],defectDt[REAC==reac & L1<=Emax & L1>=Emin,V1],L1)$y]

		reacs <- expDt[,unique(REACID)]
		incl_reacs <- reacs[!grepl("CS/REAC/\\d{6}/TOT|CS/EL",reacs)]
		for(incl_reac in incl_reacs) {
			# this could probably be done in a better way without looping and if statements
			if(incl_reac=="CS/TOT") {
				xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(V1),by=L1]
				expDt[REACID==incl_reac & L1<=Emax & L1>=Emin,PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]
			} else if(incl_reac=="CS/PROD/001001/TOT") { # (nx,p)
				xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(P*V1,na.rm=TRUE),by=L1]
				expDt[REACID==incl_reac & L1<=Emax & L1>=Emin,PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]
			} else if(incl_reac=="CS/PROD/002004/TOT") { # (nx,alpha)
				xs_table <- defectDt[L1<=Emax & L1>=Emin,sum(A*V1,na.rm=TRUE),by=L1]
				expDt[REACID==incl_reac & L1<=Emax & L1>=Emin,PRED:=approx(xs_table$L1,xs_table$V1,L1)$y]
			}
		}

		defectDt[j,V1:=0]
	}