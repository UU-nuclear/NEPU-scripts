#source("config/config-Fe56.R")
source("config/config-Cr52-hetGP.R")


expDt <- read_object(3, "expDt")
#updSysDt <- read_object(4, "updSysDt")

plotDir <- file.path(plotPath, 'before_MLO_uncertainties')
dir.create(plotDir, recursive=TRUE, showWarnings=FALSE)
for(channel in unique(expDt[,REAC])) {
	cur_expDt <- expDt[REAC==channel]
	#cur_SysDt <- updSysDt[REAC==channel]

	ggp1 <- ggplot(cur_expDt) + theme_bw() +
	theme(text=element_text(size=2),
			axis.text=element_text(size=6),
		    axis.title=element_text(size=12),
		    plot.title=element_text(size=3),
		    legend.text=element_text(size=2),
		    legend.title=element_text(size=2)) +
	xlab("energy (MeV)") + ylab("cross section (mbarn)") +
	guides(col = "none") +
	#geom_errorbar(aes(x = L1, ymin = DATA - UPDUNC, ymax = DATA + UPDUNC), col = "black",
	#                               linewidth = 0.1, width = 0.3) +
	geom_errorbar(aes(x = L1, ymin = DATA - ORIG_UNC, ymax = DATA + ORIG_UNC, col = EXPID),
	                           linewidth = 0.1, width = 0.2) +
	geom_point(aes(x = L1, y = DATA, col = EXPID),size=0.2) +
	xlim(2,23) +
	ylim(0,160)

	#new_scale_colour() +
	#geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN, col='orig. sys. unc.'), linewidth = 0.1) +
	#geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_STAT, col='only stat. unc.'), linewidth = 0.1) +
	#geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_MLO, col='extra sys. unc.'), linewidth = 0.1) +
	#geom_line(data = cur_SysDt[ERRTYPE=="pw"],aes(x = EN_POLY_CHAIN, y = XS_POLY_CHAIN_PRIOR, col='prior'), linewidth = 0.1) +
	#scale_color_manual(name='Regression Model',
	#                     breaks=c('orig. sys. unc.', 'only stat. unc.', 'extra sys. unc.','prior'),
	#                     values=c('extra sys. unc.'='green', 'only stat. unc.'='red', 'orig. sys. unc.'='black','prior'='blue')) +
	#ylim(0,NA)

	filename <- file.path(plotDir, paste0(channel,'.png'))
	ggsave(filename, ggp1, width = 0.2*53, height = 0.2*29.8125, units = "cm", dpi = 300)
}