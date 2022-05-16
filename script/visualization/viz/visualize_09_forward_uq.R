

source("config.R")
library(ggplot2)
#outdataPathRun <- paste0(outdataPath, "_2-200_MeV_100_files_small_lenght_good")
outdataPathRun <- paste0(outdataPath)
plotPath <- paste0(outdataPathRun, "/plots")
crossSectionSamples<- read_object(9, 'allResults', outdata_path = outdataPathRun)
optExpDt <- read_object(6, "optExpDt", outdata_path = outdataPathRun)
updSysDt <- read_object(4, "updSysDt", outdata_path = outdataPathRun)
needsDt <- read_object(1, "needsDt", outdata_path = outdataPathRun)
optRes <- read_object(7, "optRes", outdata_path = outdataPathRun)
subents <- read_object(1, "subents", outdata_path = outdataPathRun)
optGpDt <- read_object(6, "optGpDt", outdata_path = outdataPathRun)
mapAssignment <- read_object(6, "mapAssignment", outdata_path = outdataPathRun)

#plot_path <- file.path(rootpath, outdataDir, "plots")

normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

S <- sysCompHandler$map(optExpDt, updSysDt, ret.mat = TRUE)
updX <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(optExpDt)
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))
optExpDt[, UPDUNC := updUnc]






#energyGridrandomFiles <- energyGrid
extNeedsDt <- needsDt[,{
  stopifnot(all(L2 == 0) & all(L3 == 0))
  list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy="compgrid"),
       L2 = 0, L3 = 0, V1 = 0)
}, by=c("PROJECTILE", "ELEMENT", "MASS", "REAC")]

extNeedsDt[, IDX := seq_len(.N)]


Sexp <- exforHandler$getJac(optExpDt, extNeedsDt, subents)


reacMap <- optExpDt[, {
  
  # Construct mask for columns in extNeedsDt that correspond to the current reaction in optExpDt
  # drop=FALSE ensures that Sexp dimensions are retained if length(IDX) == 1
  mask <- apply(Sexp[IDX, , drop=FALSE], 2, function(x) !all(x == 0) ) 
  
  # Apply the mask to extNeedsDt and retrive the corresponding unique talys reaction codes 
  talysReac <- extNeedsDt[mask, unique(REAC)]
  
  # Stop if the mapping exforReac -> talysReac is not unique
  stopifnot(length(talysReac) == 1)
  
  # Return a data table containing the unique mapping
  data.table("REACTALYS" = talysReac)
}, by=REAC]

reacMap
#meanCov <- cov.wt(t(crossSectionSamples),rep(1, dim(crossSectionSamples)[2]))
meanCov <- cov.wt(t(crossSectionSamples))
centerVector <- meanCov$center
modelUnc <- sqrt(diag(meanCov$cov))

extNeedsDt[,V1:=crossSectionSamples[,1]]
reacExpidDt <- setnames(optExpDt[,as.vector(unique(EXPID))[1],by=REAC], "V1", "EXPID") 

talysExpDt <- data.table(IDX=seq(1, nrow(reacExpidDt)*length(energyGridrandomFiles)),
                         EXPID =rep(NA, each=length(energyGridrandomFiles)),
                         DIDX=rep(seq(1,length(energyGridrandomFiles)), times = nrow(reacExpidDt)),
                         REACTALYS=extNeedsDt$REAC,
                         REAC= "NA",
                         L1=extNeedsDt$L1,
                         L2=NA,
                         L3=NA,
                         DATALM=crossSectionSamples[,1],
                         DATA=centerVector,
                         DATAREF=NA,
                         REFDATA=NA,
                         UNC=modelUnc
)
talysExpDt[,REAC:=reacMap[REACTALYS==.BY, REAC], by=REACTALYS]

reactions <- talysExpDt[,unique(REAC)]
#reactions <- c("(26-FE-56(N,P)25-MN-56,,SIG)")

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

for (curReac in reactions) {
  
  curExpDt <- talysExpDt[REAC == curReac]
  
  
  curGpDt <-  optGpDt[EXPID %in% mapAssignment[REAC==curReac, EXPID]]
  sig <- curGpDt[PARNAME=="sigma", PARVAL]
  len <- curGpDt[PARNAME=="len", PARVAL]
  nug <- curGpDt[PARNAME=="nugget", PARVAL]
  lableString <- sprintf("sigma == %3.2f~~lambda == %3.2f ", sig, len)
  #lableString <- sprintf("sigma == %3.2f~~lambda == %3.2f~~tau ==  %3.2f  ", sig, len, nug)
  
  ggp <- ggplot(curExpDt) + theme_bw(base_family = "Times") + guides(col = FALSE)
  ggp <- ggp + theme(axis.text=element_text(size=16),
                     axis.title=element_text(size=18))
  
  ggp <- ggp + xlab(expression("Incident Energy"~italic(E)~"(MeV)")) + 
    ylab(expression("Angle Integrated Cross Section"~sigma~"(mb)"))
  
  ggp <- ggp + ggtitle(curReac)
  
  
  #print(reacMap[REACTALYS==curReac, REAC])
  ggp <- ggp + geom_point(data = optExpDt[REAC==curReac], 
                          aes(x = L1, y = DATA), 
                          col = "black", alpha=0.5, size = 0.4)
  ggp <- ggp + geom_errorbar(data = optExpDt[REAC==curReac], 
                             aes(x = L1, ymin = DATA - UPDUNC, 
                                 ymax= DATA + UPDUNC), 
                             color="gray", size = 0.4, alpha=0.5)  
  
  ggp <- ggp + geom_vline(xintercept=minExpEn, linetype="dashed", 
                          color="gray70", size=0.5)
  ggp <- ggp + geom_vline(xintercept=maxExpEn, linetype="dashed", 
                          color="gray70", size=0.5)
  
  ggp <- ggp + geom_ribbon(data = curExpDt, 
                           aes(x = L1, ymin = DATA -2*UNC, ymax = DATA +2*UNC), 
                           size = 0.2, 
                           fill = "lightskyblue2", 
                           alpha=0.5)   
  
  ggp <- ggp + geom_ribbon(data = curExpDt, 
                           aes(x = L1, ymin = DATA - UNC, ymax = DATA + UNC), 
                           size = 0.2, 
                           fill = "lightskyblue2", 
                           alpha=0.8) 

#  ggp <- ggp + geom_ribbon(data = curExpDt, 
#                           aes(x = L1, ymin = qnorm(0.05, DATALM, UNC), ymax = qnorm(0.95, DATALM, UNC)), 
#                           size = 0.2, 
#                           fill = "lightskyblue2", 
#                           alpha=0.5)   

#  ggp <- ggp + geom_ribbon(data = curExpDt, 
#                           aes(x = L1, ymin = qnorm(0.32, DATALM, UNC), ymax = qnorm(0.68, DATALM, UNC)), 
#                           size = 0.2, 
#                           fill = "skyblue2", 
#                           alpha=0.5) 
  
#  ggp <- ggp + geom_ribbon(data = curExpDt, 
#                           aes(x = L1, ymin = DATALM-UNC, ymax =DATALM+UNC), 
#                           size = 0.2, 
#                           fill = "lightskyblue2", 
#                           alpha=0.5)    
  ggp <- ggp + geom_line(data = curExpDt, aes(x = L1, y = DATA), size = 0.4, color="black")  
  ggp <- ggp + geom_line(data = curExpDt, aes(x = L1, y = DATALM), size = 0.4, color="red")  
  #ggp <- ggp + scale_y_continuous(labels=scientific_10)
  #ggp <- ggp + scale_x_log10() + scale_y_log10()
  #ggp <- ggp + scale_x_log10()
  ggp <- ggp + annotate("text", x = Inf, y = Inf, label = lableString, hjust =1.1, vjust = 3.1, size=5, parse = TRUE)   
  ggp <- ggp + xlim(2, 200)
  print(ggp)
  
  filepath <- file.path(plotPath,"09", "forward_prop_gp_prior_obs") 
  filename <- paste0(curReac, "_forward_prop_gp_prior_obs.pdf")
  if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  
  path_to_file <- file.path(filepath, filename)
  print(path_to_file)
  ggsave(path_to_file, ggp, width = 36, height = 24, units = "cm", dpi = 300)
  
  #filepath <- file.path(plot_path, "default_forward_prop_pdf") 
  #filename <- paste0(curReac, "_observable_mean_stderr.pdf")
  #if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  #ggsave(file.path(filepath, filename), ggp, width = 36, height = 24, units = "cm", dpi = 300)  
  #filepath <- file.path(plot_path, "default_forward_prop_png") 
  #filename <- paste0(curReac, "_observable_mean_stderr.png")
  #if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
  #ggsave(file.path(filepath, filename), ggp, width = 36, height = 24, units = "cm", dpi = 300)  
}
