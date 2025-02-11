source("config/config-Cr52-hetGP.R")


extNeedsDt <- read_object(2, "extNeedsDt")
talys_pred <- read_object(2, "rawRes")

energy_step = 0.1
energies <- seq(from=extNeedsDt[,min(L1)], to=extNeedsDt[,max(L1)], by=energy_step)

projectile_A <- 1
projectile_Z <- 0
target_A <- 52
target_Z <- 24

compound A <- 53
compound_Z <- 24

##

dt0 <- data.table()

dt0 <- rbind(dt0,data.table(
	a = c(1,2,3),
	b = c(3,4,5)
))

dt0 <- rbind(dt0,dt2 <- data.table(
	a = c(6,7,8),
	b = c(9,10,11)
))


rbind(dt1,dt2)