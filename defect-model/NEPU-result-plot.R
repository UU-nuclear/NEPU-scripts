
allResults <- read_object(12, "allResults")
expDt <- read_object(3, "expDt")
expDt[,UNC:=ORIG_UNC] # UNC was replaced using hetGP: undo

# reduce the energy range
setorder(expDt,REAC,L1)

expDt_full <- copy(expDt)

minE <- 1
maxE <- 10
expDt <- expDt[L1<maxE & L1>minE]
expDt[,IDX:=seq_len(.N)]
##########################################
# retrieve systematic uncertainties in the format of nucdataBaynet
sysDt <- read_object(4, "updSysDt")
sysDt[, ADJUSTABLE:=NULL]

sysDt <- sysDt[grepl("EXPID-",EXPID)] # remove the stuff related to the gp fit in step4
sysDt[,IDX:=seq_len(.N)]

# prepare the handlers to map systematic uncertainties of the experiments
normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

# create global handler and register the individual handlers
sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

exp_sys_cov <- sysCompHandler$cov(sysDt, ret.mat = TRUE) # U matrix 
exp_sys_map <- sysCompHandler$map(expDt, sysDt, ret.mat = TRUE) # S matrix
exp_stat_unc <- Diagonal(x = expDt[,UNC^2]) # D matrix

expDt[,TOT_UNC:=sqrt(diag(exp_sys_map %*% exp_sys_cov %*% t(exp_sys_map) + exp_stat_unc))]

#####################################################
reactions <- expDt[,unique(REAC)]

# create a data.table with the model energy-grid
# and the reaction channel specification

needsDt <- read_object(1, "needsDt")
extNeedsDt <- needsDt[,{
  stopifnot(all(L2 == 0) & all(L3 == 0))
  list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
       L2 = 0, L3 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]
extNeedsDt[, IDX := seq_len(.N)]

# create the grid which to interpolate to

# create the interpolation matrix using the exforHandler
modSubents <- lapply(reactions, createSubentStub, en=energyGridrandomFiles) # fake exfor sub-entries
modDt <- exforHandler$extractData(modSubents, ret.values=FALSE) # fake experimental data
Smod <- exforHandler$getJac(modDt, extNeedsDt, modSubents)
setkey(modDt, IDX, DIDX)

modDt[,L2:=NULL]
modDt[,L3:=NULL]

reorder <- function(x)
{
  as.vector(Smod %*% x)
}
allResults <- apply(allResults,2,reorder)

talys_cov <- cov(t(allResults[,2:ncol(allResults)]))
talys_mode <- allResults[,1]

modDt[,mode:=talys_mode]
modDt[,UNC:=sqrt(diag(talys_cov))]

modDt_red <- modDt[L1<=10 & L1>=1.0]

ggp <- ggplot(data=expDt) +
  #geom_point(aes(x=L1,y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=mode),data=modDt_red, col='red') +
  geom_ribbon(aes(x=L1,ymin=mode-UNC, ymax=mode+UNC),data=modDt_red, fill='red', alpha=0.2) +
  facet_wrap(~ REAC, ncol=2, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp

min_plot_en <- 1
max_plot_en <- 2
channels <- c("(26-FE-56(N,TOT),,SIG)","(26-FE-56(N,INL)26-FE-56,,SIG)")
ggp <- ggplot(data=expDt[REAC %in% channels][L1>min_plot_en & L1<max_plot_en]) +
  #geom_point(aes(x=L1,y=DATA)) +
  geom_errorbar(aes(x=L1, ymin=DATA-TOT_UNC, ymax=DATA+TOT_UNC), col='grey') +
  geom_errorbar(aes(x=L1, ymin=DATA-UNC, ymax=DATA+UNC)) +
  geom_line(aes(x=L1,y=mode),data=modDt_red[REAC %in% channels][L1<=2.0], col='red') +
  geom_ribbon(aes(x=L1,ymin=mode-UNC, ymax=mode+UNC),data=modDt_red[REAC %in% channels][L1<=2.0], fill='red', alpha=0.2) +
  facet_wrap(~ REAC, ncol=1, scales="free_y") +
  labs(x="neutron energy (MeV)", y="cross section (mbarn)") +
  theme_bw()
ggp
