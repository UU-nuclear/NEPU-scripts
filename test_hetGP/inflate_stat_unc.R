library(hetGP)
library(MASS)
library(ggplot2)
library(data.table)


expDt <- read_object(3,"expDt")
modDt <- read_object(3, "modDt")

for(expID in expDt[,unique(EXPID)]) {
	curExpDt <- expDt[EXPID==expID]

	# maybe this should not be done
	# it is difficult to estimate which region has resonances and which don't
	#  limit to the region of resonance structure
	# Ec <- 0.100 # 100 keV
	# # data with energy spacing below Ec are considered as resonance data
	# # data with energy spacing larger than Ec are considered as average cross-section
	# setorder(curExpDt,L1)
	# energies <- curExpDt$L1
	# nn <- length(energies)
	# dE <- energies[2:nn]-energies[1:nn-1]
	# dE <- c(dE[1],dE) # the distance for the first point is to the next point
	# curExpDt[,dE := dE]
	# curExpDt <- curExpDt[dE<Ec]

	reac_channel <- curExpDt[,unique(REAC)]
	curModDt <- modDt[REAC==reac_channel]
	setorder(curModDt,L1)

	rows_to_keep <- which(curModDt[,L1]>min(curExpDt[,L1]) & curModDt[,L1]<max(curExpDt[,L1]))
	rows_to_keep <- c(min(rows_to_keep)-1,rows_to_keep,max(rows_to_keep)+1)
	curModDt <- curModDt[rows_to_keep]

	# fit a homoschedastic GP to the talys model to get the length scale
	talysGP <- mleHomGP(X = matrix(curModDt$L1, ncol = 1), Z = curModDt$DATA, covtype = "Matern5_2")

	# fit a heteroschedastic GP to the experimental data, using the length scale from the talys model
	# expGP <- mleHetGP(X = matrix(curExpDt$L1, ncol = 1),
	# 					Z = curExpDt$DATA,
	# 					covtype = "Matern5_2",
	# 					known = list(theta=talysGP$theta),
	# 					settings = list(linkThetas='none'))

	expGP <- mleHetGP(X = matrix(curExpDt$L1, ncol = 1),
						Z = curExpDt$DATA,
						covtype = "Matern5_2",
						known = list(theta=talysGP$theta))
}

ggp <- ggplot(data=curExpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=4),
                   axis.title=element_text(size=4),
                   strip.text=element_text(size=3))
ggp <- ggp + geom_line(aes(x=L1, y=DATA))
ggp

Xgrid <- matrix(curExpDt$L1, ncol=1)
pred_exp <- as.data.table(predict(x = Xgrid, object = expGP))
pred_exp[,L1:=Xgrid]
ggp <- ggp + geom_ribbon(aes(x=L1, ymin=mean-sqrt(nugs), ymax=mean+sqrt(nugs)),data=pred_exp, col='blue', alpha=0., linetype = "dashed")
ggp <- ggp + geom_line(aes(x=L1, y=mean),data=pred_exp, col='blue')

ggp + geom_errorbar(aes(x=L1,ymin=DATA-UNC,ymax=DATA+UNC)) 