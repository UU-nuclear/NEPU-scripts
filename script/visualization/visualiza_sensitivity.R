optRes <- read_object(10, "optRes")
optParamDt <- read_object(10,"optParamDt")

optExpDt <- read_object(6, "optExpDt")


optParamDt <- optParamDt[ADJUSTABLE==TRUE]
J <- optRes$jac

tmpDt <- copy(optExpDt)

# try with par#1 == "v2adjust n"

par_num <- 57

tmpDt[, JAC := J[,par_num]]

ggp <- ggplot(data=tmpDt)
ggp <- ggp + theme_bw() + theme(legend.position="none")
ggp <- ggp + theme(axis.text=element_text(size=9),
                   axis.title=element_text(size=10),
                   strip.text=element_text(size=8))
ggp <- ggp + xlab('enegy [MeV]') + ylab('dsigma/dp [mbarn]')

ggp <- ggp + geom_point(aes(x=L1, y=JAC), size=0.2, alpha=0.6)
ggp <- ggp + geom_line(aes(x=L1, y=JAC), size=0.2, alpha=0.6)
ggp <- ggp + facet_wrap(~REAC, scales='free_y')

print(ggp)

model <- lm(JAC ~ L1, data = tmpDt[REAC=="(26-FE-56(N,TOT),,SIG)"])

x <- energyGrid
y <- approx(tmpDt[REAC=="(26-FE-56(N,TOT),,SIG)"]$L1,tmpDt[REAC=="(26-FE-56(N,TOT),,SIG)"]$JAC,xout=x)
plot(y[[1]],y[[2]])

plot(tmpDt[REAC=="(26-FE-56(N,TOT),,SIG)"]$L1,tmpDt[REAC=="(26-FE-56(N,TOT),,SIG)"]$JAC)
points(y[[1]],y[[2]], col='red', pch=19)

JacDt <- data.table(energyGrid)
reactions <- unique(tmpDt$REAC)
for(reaction in reactions) {
	print(reaction)
	model <- lm(JAC ~ L1, data = tmpDt[REAC==reaction])
	JacDt[, new := approx(tmpDt[REAC==reaction]$L1,tmpDt[REAC==reaction]$JAC,xout=energyGrid)[[2]]]
}

