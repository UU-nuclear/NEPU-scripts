
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

plotDt <- data.table(PARNAME=parnames,parval=as.vector(allParsets_ext[,2:ncol(allParsets)]))

# first the parameters that are not energy dependent and that have been adjusted in LM
histplot1 <- ggplot(data=plotDt[PARNAME %in% adj_enindep_par_names],
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
		    xlab("parameter value") +
		    ylab("")
print(histplot1)

cur_plotPath <- file.path(plotPath, 'parameters/')
dir.create(cur_plotPath, recursive=TRUE, showWarnings=FALSE)
filename <- file.path(cur_plotPath, paste0("parameters.png"))
ggsave(filename, histplot1, width = 0.5*29.7, height = 21.0, units = "cm", dpi = 300)

# the same energy independent parameters but in the internal parameter space
plotDt_int <- data.table(PARNAME=parnames,parval=as.vector(allParsets[,2:ncol(allParsets)]))

# first the parameters that are not energy dependent and that have been adjusted in LM
histplot3 <- ggplot(data=plotDt_int[PARNAME %in% adj_enindep_par_names],
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
		    xlab("parameter value") +
		    ylab("")
print(histplot3)

cur_plotPath <- file.path(plotPath, 'parameters/')
dir.create(cur_plotPath, recursive=TRUE, showWarnings=FALSE)
filename <- file.path(cur_plotPath, paste0("parameters_internal.png"))
ggsave(filename, histplot3, width = 0.5*29.7, height = 21.0, units = "cm", dpi = 300)

# now the parameters that are not energy dependent and that were not adjusted in LM
names <- enindep_par_names[!enindep_par_names %in% adj_enindep_par_names]
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

cur_plotPath <- file.path(plotPath, 'parameters/')
dir.create(cur_plotPath, recursive=TRUE, showWarnings=FALSE)
filename <- file.path(cur_plotPath, paste0("parameters_not_adjusted.png"))
ggsave(filename, histplot1, width = 0.5*29.7, height = 2*21.0, units = "cm", dpi = 300)

energies <- str_extract(adj_endep_par_names,"\\(.+\\)")
energies <- str_sub(energies,2,-2)

plotDt[,EN:=NULL]
plotDt[PARNAME %in% adj_endep_par_names,EN:=rep(energies,times=.N/length(energies))]
plotDt_endep <- plotDt[PARNAME %in% adj_endep_par_names]

plotDt_endep[,PARAMETER:=str_replace(PARNAME,"\\(.+\\)"," ")]

cur_plotPath <- file.path(plotPath, 'parameters/endep/')
dir.create(cur_plotPath, recursive=TRUE, showWarnings=FALSE)

for(par in unique(plotDt_endep[,PARAMETER])) {
	print(par)

	histplot3 <- ggplot(data=plotDt_endep[PARAMETER == par][order(EN)],
					aes(y=as.numeric(EN), x=parval, fill=EN, group=EN)) +
			geom_density_ridges(alpha=0.6, stat="binline", binwidth=0.01) +
		    theme_ridges() +
		    theme(
		      legend.position="none",
		      panel.spacing = unit(0.1, "lines"),
		      strip.text.x = element_text(size = 8),
		      panel.background = element_rect(fill='white'),
		      plot.background = element_rect(fill='white')
		    ) +
		    xlab(par) +
		    ylab("Energy (MeV)")

	filename <- file.path(cur_plotPath, paste0(str_replace(par," ","_"),".png"))
	ggsave(filename, histplot3, width = 29.7, height = 21.0, units = "cm", dpi = 300)
	
}