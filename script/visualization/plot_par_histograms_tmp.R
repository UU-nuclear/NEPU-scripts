
library(ggridges)
library(ggplot2)

allParamDt <- read_object(12, "allParamDt")
allParsets <- read_object(12, "allParsets")
optParamDt <- read_object(10,"optParamDt")

# transform the sampled parameters to the external parameter space
paramTrafo <- parameterTransform(
                  x0 = unlist(allParamDt[ADJUSTABLE==TRUE,PARVAL]),
                  delta = allParamDt[ADJUSTABLE==TRUE,unlist(PARVAL) - PARMIN])

allParsets_ext <- paramTrafo$fun(allParsets)


parnames <- allParamDt[ADJUSTABLE==TRUE,PARNAME]
selection <- grepl("\\(.+\\)",parnames)
enindep_par_names <- parnames[!selection]
endep_par_names <- parnames[selection]

adj_parnames <- optParamDt[ADJUSTABLE==TRUE,PARNAME]
selection <- grepl("\\(.+\\)",adj_parnames)
adj_enindep_par_names <- adj_parnames[!selection]
adj_endep_par_names <- adj_parnames[selection]

adj_enindep_par_names <- adj_enindep_par_names[c(7,8,9,10,14,18)]

plotDt <- data.table(PARNAME=parnames,parval=as.vector(allParsets_ext[,2:ncol(allParsets)]/allParsets_ext[,1]))

# first the parameters that are not energy dependent and that have been adjusted in LM
histplot1 <- ggplot(data=plotDt[PARNAME %in% adj_enindep_par_names],
					aes(y=PARNAME, x=parval,  fill=PARNAME)) +
			geom_density_ridges(alpha=0.6, stat="binline", bins=100) +
		    theme_ridges() +
		    theme(
		      legend.position="none",
		      panel.spacing = unit(0.1, "lines"),
		      strip.text.x = element_text(size = 8),
		      panel.background = element_rect(fill='white'),
		      plot.background = element_rect(fill='white')
		    ) +
		    xlab(bquote(p[ext]/p[0])) +
		    ylab("")
print(histplot1)

cur_plotPath <- file.path(plotPath, 'parameters/')
dir.create(cur_plotPath, recursive=TRUE, showWarnings=FALSE)
filename <- file.path(cur_plotPath, paste0("parameters_ext.png"))
ggsave(filename, histplot1, width = 29.7/2, height = 21.0/2, units = "cm", dpi = 300)

# the same energy independent parameters but in the internal parameter space
plotDt_int <- data.table(PARNAME=parnames,parval=as.vector(allParsets[,2:ncol(allParsets)]/allParsets_ext[,1]))

# first the parameters that are not energy dependent and that have been adjusted in LM
histplot3 <- ggplot(data=plotDt_int[PARNAME %in% adj_enindep_par_names],
					aes(y=PARNAME, x=parval,  fill=PARNAME)) +
			geom_density_ridges(alpha=0.6, stat="binline", bins=100) +
		    theme_ridges() +
		    theme(
		      legend.position="none",
		      panel.spacing = unit(0.1, "lines"),
		      strip.text.x = element_text(size = 8),
		      panel.background = element_rect(fill='white'),
		      plot.background = element_rect(fill='white')
		    ) +
		    xlab(bquote(p[int]/p[0])) +
		    ylab("")
print(histplot3)

cur_plotPath <- file.path(plotPath, 'parameters/')
dir.create(cur_plotPath, recursive=TRUE, showWarnings=FALSE)
filename <- file.path(cur_plotPath, paste0("parameters_int.png"))
ggsave(filename, histplot3, width = 29.7/2, height = 21.0/2, units = "cm", dpi = 300)

# now the parameters that are not energy dependent and that were not adjusted in LM
names <- enindep_par_names[!(enindep_par_names %in% adj_enindep_par_names)]
histplot2 <- ggplot(data=plotDt[PARNAME %in% names],
					aes(y=PARNAME, x=parval,  fill=PARNAME)) +
			geom_density_ridges(alpha=0.6, stat="binline", binwidth=0.01) +
		    theme_ridges() +
		    theme(
		      legend.position="none",
		      panel.spacing = unit(0.1, "lines"),
		      strip.text.x = element_text(size = 8),
		      panel.background = element_rect(fill='white'),
		      plot.background = element_rect(fill='white')
		    ) +
		    xlab("") +
		    ylab("Probability (%)")
print(histplot2)
