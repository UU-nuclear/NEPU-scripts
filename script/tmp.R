###############################
# perform the check
##############################

# convert the sparse matrix given as data.table 
# into a spase matrix type as defined in package Matrix
Spar <- with(fullSensDt,
             sparseMatrix(i = IDX1, j = IDX2, x = X,
                          dims = c(nrow(extNeedsDt), nrow(refParamDt))))
Sexp <- exforHandler$getJac(optExpDt, extNeedsDt, subents)
Sglob <- Sexp %*% Spar 

# safeguard
stopifnot(all(dim(Sglob) == c(nrow(optExpDt), nrow(refParamDt))))

# convert the sparse matrix Sglob into a datatable
SglobDt <- as.data.table(summary(Sglob))
setnames(SglobDt, c("IDX1", "IDX2", "X"))

# we (linearly) propagate all parameter values equal one
# to the model predictions
imp1 <- as.vector(Spar %*% rep(1, nrow(refParamDt)))
# we propagate hypothetical experimental values equal one
# to the model prediction
imp2 <- as.vector(t(Sexp) %*% rep(1, nrow(Sexp)))
# then we select observables on the model grid that
# are affected by both the backpropagation from the
# experiment and the forward propagation of model parameters
impIdx <- which(imp1 * imp2 != 0)

optSparDt <- copy(fullSensDt)
setkey(optSparDt, IDX1)
optSparDt <- optSparDt[J(impIdx)]

paramImpactDt <- SglobDt[, list(IMP = max(abs(X))), by = "IDX2"]
paramImpactDt <- paramImpactDt[order(IMP, decreasing = TRUE)]
selParIdcs <- paramImpactDt[IMP >= 1, IDX2]

setkey(optSparDt, IDX2)
mask <- optSparDt[J(selParIdcs), list(DSTIDX = IDX1, SRCIDX = IDX2)]
setkey(mask, SRCIDX, DSTIDX)
adjParIdcs <- unique(mask$SRCIDX)

# make a copy of the reference parameter datatable 
# and define what and what not we want to optimize
optParamDt_post <- copy(refParamDt)
setkey(optParamDt_post, IDX)
optParamDt_post[, ADJUSTABLE := FALSE]
optParamDt_post[J(adjParIdcs), ADJUSTABLE := TRUE]

# safeguard
stopifnot(sum(optParamDt_post$ADJUSTABLE) == length(adjParIdcs))

adjustable_par_names_post <- optParamDt_post[ADJUSTABLE==TRUE]$PARNAME
adjustable_par_names_prior <- optParamDt[ADJUSTABLE==TRUE]$PARNAME

if(!all(adjustable_par_names_post %in% adjustable_par_names_prior)) {
  print("There are parameters sensitive to the experimental data at the posterior paramter set that were not sensitive at the initial paramter set")
}