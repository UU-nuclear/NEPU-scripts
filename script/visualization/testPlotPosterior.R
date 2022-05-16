
library(matrixStats)
library(ggplot2)

allResults <- read_object(9, 'allResults')

# total xs is allResults[1:length(energyGrid),n_calc]
# all total xs is allResults[1:length(energyGrid),n_calc]
# plot(rep(length(energyGrid)),allResults[1:length(energyGrid),]) plots all calculations total xs

all_total_xs <- allResults[1:length(energyGrid),] # matrix[length(energyGrid),n_calc]
mean_total_xs <- rowMeans(all_total_xs)
sds_total_xs <- rowSds(all_total_xs)

# I belive that the best paramter estimate is calc 001
best_total_xs <- allResults[1:length(energyGrid),1]

total_xs_df <- data.frame(energyGrid,mean_total_xs,sds_total_xs,best_total_xs)

ggp <- ggplot(total_xs_df,aes(x=energyGrid,y=mean_total_xs))
ggp <- ggp + geom_line()
ggp <- ggp + geom_ribbon(aes(ymin=mean_total_xs-sds_total_xs,ymax=mean_total_xs+sds_total_xs,alpha=0.3))
ggp <- ggp + geom_line(aes(y=best_total_xs),col='red')