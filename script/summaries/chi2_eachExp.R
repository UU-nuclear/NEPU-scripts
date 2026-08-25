#################################################
#       Print Chi-square summary table
#       for all reaction channels
#################################################
args = commandArgs(trailingOnly=TRUE)

if (length(args)==1) {
  print(paste0("Setting as config file: ", args[1]))
  source(args[1])
}

library(data.table)

##################################################
#       OUTPUT FROM PREVIOUS STEPS
##################################################

subents <- read_object(1, "subents")
needsDt <- read_object(1, "needsDt")
extNeedsDt <- read_object(2, "extNeedsDt")
origSysDt <- read_object(4, "origSysDt")
updSysDt <- read_object(4, "updSysDt")
modDt <- read_object(3, "modDt")
optExpDt <- read_object(7, "optExpDt")
allResults <- read_object(12, 'allResults')

optParamDt <- read_object(10, "optParamDt")

# sort the experimental data table by reaction, energy
optExpDt <- optExpDt[order(REAC, L1)]
optExpDt[, IDX := seq_len(.N)]

reactions <- optExpDt[, unique(REAC)]

normHandler <- createSysCompNormHandler("DATAREF")
normHandler$addSysUnc("EXPID", "", 0, 0, TRUE)

sysCompHandler <- createSysCompHandler()
sysCompHandler$addHandler(normHandler)

S <- sysCompHandler$map(optExpDt, origSysDt, ret.mat = TRUE)
origX <- sysCompHandler$cov(origSysDt, ret.mat = TRUE)
updX <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
statUnc <- getDt_UNC(optExpDt)

origUnc <- sqrt(statUnc^2 + diag(S %*% origX %*% t(S)))
updUnc <- sqrt(statUnc^2 + diag(S %*% updX %*% t(S)))

setkey(optExpDt, IDX)
optExpDt[, ORIGUNC := origUnc]
optExpDt[, UPDUNC := updUnc]


RandomFilesNeedsDt <- needsDt[, {
    stopifnot(all(L2 == 0) & all(L3 == 0))
    list(L1 = defineEnergyGrid(L1, energyGridrandomFiles, enPolicy = "compgrid"),
         L2 = 0, L3 = 0)
}, by = c("PROJECTILE", "ELEMENT", "MASS", "REAC")]
RandomFilesNeedsDt[, IDX := seq_len(.N)]

allResults[is.na(allResults)] <- 0
sampledResults <- allResults[, 2:ncol(allResults)]

Sexp <- exforHandler$getJac(optExpDt, RandomFilesNeedsDt, subents)
transform_function <- function(x) { as.vector(Sexp %*% x) }
allResults_transform <- apply(allResults, 2, FUN = transform_function)
sampledResults_transform <- allResults_transform[, 2:ncol(allResults_transform)]

uncinfo_transform <- cov.wt(t(sampledResults_transform))
model_covariance <- uncinfo_transform$cov
model_predictions <- uncinfo_transform$center

optExpDt <- optExpDt[order(REAC, L1)]
optExpDt[, FIT_MEAN := model_predictions]
optExpDt[, FIT_MODE := allResults_transform[, 1]]

##################################################
#   Get TALYS default predictions at exp. points
##################################################

# modDt contains the TALYS reference (default/prior) calculation
# Map it to experimental points via interpolation
optExpDt[, TALYS_DEFAULT := NA_real_]
for (curReac in reactions) {
    curModDt <- modDt[REAC == curReac]
    curExpDt <- optExpDt[REAC == curReac]
    if (nrow(curModDt) > 1) {
        interp <- approx(curModDt$L1, curModDt$DATA, xout = curExpDt$L1, rule = 2)
        optExpDt[REAC == curReac, TALYS_DEFAULT := interp$y]
    }
}

##################################################
#   Chi-square functions
##################################################

getChi2 <- function(reaction_string, prediction_col) {
    reaction.expDt <- copy(optExpDt[REAC == reaction_string])
    reaction.expDt[, IDX := seq_len(.N)]
    reaction.stat_unc <- getDt_UNC(reaction.expDt)
    reaction.D <- Diagonal(x = reaction.stat_unc^2)
    reaction.S <- sysCompHandler$map(reaction.expDt, origSysDt, ret.mat = TRUE)
    reaction.P <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
    reaction.d <- reaction.expDt[, DATA - get(prediction_col)]
    chisquare(reaction.d, reaction.D, reaction.S, reaction.P)
}

getChi2_with_model_unc <- function(reaction_string, prediction_col) {
    reaction.expDt <- copy(optExpDt[REAC == reaction_string])
    reaction.model_cov <- model_covariance[reaction.expDt$IDX, reaction.expDt$IDX]
    reaction.expDt[, IDX := seq_len(.N)]
    reaction.stat_unc <- getDt_UNC(reaction.expDt)
    reaction.D <- Diagonal(x = reaction.stat_unc^2)
    reaction.S <- sysCompHandler$map(reaction.expDt, origSysDt, ret.mat = TRUE)
    reaction.P <- sysCompHandler$cov(updSysDt, ret.mat = TRUE)
    reaction.exp_cov <- reaction.S %*% reaction.P %*% t(reaction.S) + reaction.D
    reaction.tot_cov <- reaction.model_cov + reaction.exp_cov
    reaction.d <- reaction.expDt[, DATA - get(prediction_col)]
    tryCatch(
        as.vector(t(reaction.d) %*% solve(reaction.tot_cov, reaction.d)),
        error = function(e) NA_real_
    )
}

getNpoints <- function(reaction_string) {
    nrow(optExpDt[REAC == reaction_string])
}

##################################################
#   Compute all chi-squares
##################################################

cat("\n")
cat("================================================================\n")
cat("  Chi-square values\n")
cat("================================================================\n\n")

results <- data.table(
    Channel = character(),
    N_points = integer(),
    Chi2_default_exp = numeric(),
    Chi2_posterior_exp = numeric(),
    Chi2_posterior_exp_fit = numeric()
)

for (curReac in reactions) {
    n <- getNpoints(curReac)

    # Chi2 with TALYS default (prior) prediction, experimental covariance only
    chi2_default <- tryCatch(
        getChi2(curReac, "TALYS_DEFAULT"),
        error = function(e) NA_real_
    )

    # Chi2 with posterior mean, experimental covariance only
    chi2_post_exp <- getChi2(curReac, "FIT_MEAN")

    # Chi2 with posterior mean, experimental + model covariance
    chi2_post_expfit <- getChi2_with_model_unc(curReac, "FIT_MEAN")

    results <- rbind(results, data.table(
        Channel = curReac,
        N_points = n,
        Chi2_default_exp = round(as.numeric(chi2_default), 1),
        Chi2_posterior_exp = round(as.numeric(chi2_post_exp), 1),
        Chi2_posterior_exp_fit = round(as.numeric(chi2_post_expfit), 1)
    ))
}

# Add per-point values
results[, Chi2_default_per_n := round(Chi2_default_exp / N_points, 2)]
results[, Chi2_posterior_per_n := round(Chi2_posterior_exp / N_points, 2)]
results[, Chi2_postfit_per_n := round(Chi2_posterior_exp_fit / N_points, 2)]

# Print table
cat(sprintf("%-45s %6s %12s %8s %12s %8s %12s %8s\n",
    "Channel", "N", "Chi2_def", "/n", "Chi2_post", "/n", "Chi2_p+f", "/n"))
cat(paste(rep("-", 120), collapse = ""), "\n")

for (i in seq_len(nrow(results))) {
    r <- results[i]
    cat(sprintf("%-45s %6d %12.1f %8.2f %12.1f %8.2f %12.1f %8.2f\n",
        r$Channel, r$N_points,
        r$Chi2_default_exp, r$Chi2_default_per_n,
        r$Chi2_posterior_exp, r$Chi2_posterior_per_n,
        r$Chi2_posterior_exp_fit, r$Chi2_postfit_per_n))
}

cat(paste(rep("-", 120), collapse = ""), "\n")

cat("\n")
cat("Chi2_def     = chi-square using TALYS default (TENDL) vs data, experimental covariance\n")
cat("Chi2_post    = chi-square using posterior mean vs data, experimental covariance\n")
cat("Chi2_p+f     = chi-square using posterior mean vs data, experimental + model covariance\n")
cat("/n           = chi-square per data point (good fit ~ 1)\n")
cat("\n")

# Save to file
dir.create(plotPath, recursive = TRUE, showWarnings = FALSE)
write.csv(results, file.path(plotPath, "chi2_summary_new_eachexp_new.csv"), row.names = FALSE)
cat("Table saved to:", file.path(plotPath, "chi2_summary_new_eachexp_new.csv"), "\n")


##################################################
#   Per-experiment chi2 breakdown
#
#   HOW IT WORKS:
#   We use the SAME getChi2 and getChi2_with_model_unc
#   functions already defined above. The trick is:
#   temporarily replace optExpDt with only the rows
#   for one experiment, then call the same function.
#   This guarantees identical Cd and identical logic
#   as the channel-level table above.
#
#   The full channel Cd has off-diagonal blocks between
#   experiments (shared normalisation structure via S0XS0T).
#   When we isolate one experiment those cross-terms
#   disappear. So per-exp chi2 tells us each experiment's
#   individual misfit — not its exact share of channel chi2.
##################################################

cat("\n")
cat(paste(rep("=", 120), collapse=""), "\n")
cat("  Per-experiment chi-square breakdown\n")
cat("  Logic: same getChi2 functions, optExpDt temporarily\n")
cat("  restricted to one experiment at a time\n")
cat(paste(rep("=", 120), collapse=""), "\n")

full_optExpDt <- optExpDt   # save full table
exp_rows <- list()          # collect rows for CSV

for (curReac in reactions) {

    expids <- full_optExpDt[REAC == curReac, unique(EXPID)]
    n_chan  <- full_optExpDt[REAC == curReac, .N]

    # Print channel header
    cat(sprintf("\n  Channel: %s  (N=%d)\n", curReac, n_chan))
    cat(sprintf("  %-22s %4s %12s %7s %12s %7s %12s %7s\n",
        "EXPID", "N",
        "Chi2_def", "/n",
        "Chi2_post", "/n",
        "Chi2_p+f", "/n"))
    cat("  ", paste(rep("-", 100), collapse=""), "\n")

    for (curExpid in expids) {

        rows  <- full_optExpDt[REAC == curReac & EXPID == curExpid]
        n_exp <- nrow(rows)
        if (n_exp == 0) next

        # ── KEY STEP ──────────────────────────────────────────
        # Temporarily make optExpDt contain only this experiment.
        # getChi2 and getChi2_with_model_unc both read optExpDt
        # internally, so they will now compute chi2 for this
        # experiment only, using the exact same Cd structure:
        #   D  = diag(stat_unc^2)   for these n_exp rows
        #   SP = from dummy handler  ≈ 0
        # So chi2 = sum_i (d_i^2 / sigma_stat_i^2)
        # where d_i = DATA_i - prediction_i
        # ──────────────────────────────────────────────────────
        optExpDt <- rows
        optExpDt[, IDX := seq_len(.N)]

        # Chi2 vs TALYS default
        chi2_def <- tryCatch(
            as.numeric(getChi2(curReac, "TALYS_DEFAULT")),
            error = function(e) NA_real_)

        # Chi2 vs posterior mean
        chi2_post <- tryCatch(
            as.numeric(getChi2(curReac, "FIT_MEAN")),
            error = function(e) NA_real_)

        # Chi2 vs posterior mean + model covariance
        # NOTE: model_covariance[rows$IDX, rows$IDX] is the
        # sub-block of the full MC sample covariance matrix
        # for these specific data point indices.
        # This captures how much parameter uncertainty
        # contributes to the cross section uncertainty
        # at these experimental energies.
        chi2_pf <- tryCatch(
            as.numeric(getChi2_with_model_unc(curReac, "FIT_MEAN")),
            error = function(e) NA_real_)

        # Restore full table before next iteration
        optExpDt <- full_optExpDt

        cat(sprintf("  %-22s %4d %12.1f %7.2f %12.1f %7.2f %12.1f %7.2f\n",
            curExpid, n_exp,
            ifelse(is.na(chi2_def),  0, chi2_def),
            ifelse(is.na(chi2_def),  0, chi2_def  / n_exp),
            ifelse(is.na(chi2_post), 0, chi2_post),
            ifelse(is.na(chi2_post), 0, chi2_post / n_exp),
            ifelse(is.na(chi2_pf),   0, chi2_pf),
            ifelse(is.na(chi2_pf),   0, chi2_pf   / n_exp)))

        exp_rows[[length(exp_rows) + 1]] <- data.table(
            Channel          = curReac,
            EXPID            = curExpid,
            N_points         = n_exp,
            Chi2_def         = round(ifelse(is.na(chi2_def),  NA, chi2_def),  1),
            Chi2_def_per_n   = round(ifelse(is.na(chi2_def),  NA, chi2_def  / n_exp), 2),
            Chi2_post        = round(ifelse(is.na(chi2_post), NA, chi2_post), 1),
            Chi2_post_per_n  = round(ifelse(is.na(chi2_post), NA, chi2_post / n_exp), 2),
            Chi2_pf          = round(ifelse(is.na(chi2_pf),   NA, chi2_pf),   1),
            Chi2_pf_per_n    = round(ifelse(is.na(chi2_pf),   NA, chi2_pf   / n_exp), 2))
    }

    cat("  ", paste(rep("-", 100), collapse=""), "\n")
}

# Restore optExpDt
optExpDt <- full_optExpDt

# Combine all rows into one table
exp_chi2_dt <- rbindlist(exp_rows)

##################################################
#   Save per-experiment table to CSV
##################################################
write.csv(exp_chi2_dt,
          file.path(plotPath, "chi2_per_experiment_new.csv"),
          row.names = FALSE)
cat("\nPer-experiment table saved to:",
    file.path(plotPath, "chi2_per_experiment_new.csv"), "\n")

cat("\n")
cat("Column definitions:\n")
cat("  Chi2_def      = sum_i (DATA_i - TALYS_DEFAULT_i)^2 / sigma_stat_i^2\n")
cat("  Chi2_post     = sum_i (DATA_i - FIT_MEAN_i)^2 / sigma_stat_i^2\n")
cat("  Chi2_p+f      = d^T (C_model + D)^-1 d  where d = DATA - FIT_MEAN\n")
cat("  /n            = chi2 divided by number of data points in that experiment\n")
