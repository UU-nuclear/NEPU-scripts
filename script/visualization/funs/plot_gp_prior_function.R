plot_gp_prior <- function(expDt, optGpDt){
  expDt[,{  
    
    curReac <- .BY
    
    curGpDt <-  optGpDt[EXPID %in% mapAssignment[REAC==.BY, EXPID]]
    sig <- curGpDt[PARNAME=="sigma", PARVAL]
    len <- curGpDt[PARNAME=="len", PARVAL]
    nug <- curGpDt[PARNAME=="nugget", PARVAL]
    
    lableString <- sprintf("sigma == %3.2f~~lambda == %3.2f~~tau ==  %3.2f  ", sig, len, nug)
    
    ggr <- ggplot(.SD) + theme_bw() + ggtitle(.BY)
    
    r <- DATA-LMFIT 
    #r <- ORIGDATA
    
    rpred <- DATA
    #rpred <- GPMEAN
    
    
    
    alpha <- 1 - .N /nrow(optExpDt)
    
    
    
    #ggr <- ggr + geom_line(aes(x = L1, y = DATAREF), color="black", size = 0.5)
    ggr <- ggr + geom_ribbon(aes(x = L1, ymin = 0-2*GPPRIOR, ymax = 0+2*GPPRIOR), fill="olivedrab", alpha=0.4) 
    ggr <- ggr + geom_ribbon(aes(x = L1, ymin = 0-GPPRIOR, ymax = 0+GPPRIOR), fill="olivedrab", alpha=0.4)  
    ggr <- ggr + geom_line(aes(x = L1, y = 0), color="black", linetype="dashed", size = 0.5) 
    
    
    ggr <- ggr + geom_errorbar(aes(x = L1, ymin = r - UPDUNC, ymax= r + UPDUNC), color="gray", size = 0.8, alpha=alpha)  
    ggr <- ggr + geom_point(aes(x = L1, y = r), size = 0.5)      
    ggr <- ggr + annotate("text", x = Inf, y = Inf, label = lableString, hjust =1.1, vjust = 3.1, size=5, parse = TRUE)   
    
    print(ggr)
    filepath <- file.path(plotPath,"06", "OPT_OBS_GP_PRIOR") 
    filename <- paste0(curReac, "_GP_PRIOR.pdf")
    if(!dir.exists(filepath)) dir.create(file.path(filepath), showWarnings = FALSE, recursive = TRUE)
    
    path_to_file <- file.path(filepath, filename)
    print(path_to_file)
    ggsave(path_to_file, ggr, width = 36, height = 24, units = "cm", dpi = 300)
    
  }
  , by=REAC]
}