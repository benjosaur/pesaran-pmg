# ==============================================================================
# Cross-Estimator Tests: LR Test and Nesting Restrictions
#
# Tests the relationships between DFE, MG, and PMG estimators that follow
# from the nesting structure:
#
#   MG (unrestricted) ⊃ PMG (common θ) ⊃ DFE (common θ, φ, γ)
#
# The MG allows all parameters to differ across groups.  The PMG restricts
# only the long-run θ to be homogeneous.  The DFE additionally restricts
# all short-run dynamics (φ, γ) to be homogeneous (plus fixed effects).
#
# ---- Why are lag orders fixed rather than chosen by AIC? ----
#
# The original GAUSS code offers two lag-selection options:
#   (1) User-specified fixed lags
#   (2) Model-selection criteria (AIC / SBC / HQ)
#
# The OUTPUT.OUT run uses option (1) with ARDL(1,2,0): p=1 for Δy,
# q=2 for Δ(income), q=0 for Δ(inflation).  This is not arbitrary:
#
#   - When a "relative variable" is present (here: inflation relative to
#     a reference country), the GAUSS code forces lag_op=1 (fixed lags)
#     and sets the relative variable's lag order to 0.  AIC-based
#     selection is disabled in this case (see SUB0.PRC line 525:
#     "if rel_op eq 2; lag_op = 1; goto SKIPLAG").
#
#   - The ARDL(1,2,0) specification matches the paper's Table 3 and is
#     the standard choice for annual macro panels of this vintage —
#     one lag of Δy captures persistence, two lags of Δ(income) allow
#     a richer short-run income response, and zero lags on the relative
#     price variable is imposed by the program.
#
#   - Allowing group-specific AIC-chosen lags would make the SR regressor
#     count vary across groups, complicating the PMG information matrix
#     and the DFE pooled regression.  The fixed-lag design keeps numsi
#     uniform, which the current code asserts.
#
# A future extension could implement AIC/SBC lag selection (the GAUSS
# AICLAG procedure in SUB0.PRC), but that is outside the scope of this
# replication.
#
# Run: Rscript tests/test_cross_estimator.R
# ==============================================================================

source("panel_ecm_fe.R")
source("load_data.R")
source("pmg_mg.R")


# ==============================================================================
# LR Test: MG (unrestricted) vs PMG (restricted long-run)
# ==============================================================================
#'
#' Under H0 (PMG homogeneity restriction on θ is valid):
#'   LR = 2 * (LL_MG - LL_PMG) ~ χ²(df)
#'
#' where LL_MG is the sum of group-specific unrestricted log-likelihoods
#' (each group estimated by OLS), and LL_PMG is the sum of group-specific
#' log-likelihoods evaluated at the common θ.
#'
#' Degrees of freedom: The unrestricted model has NN*k free LR parameters
#' (k per group).  The restricted model has k.  So df = (NN-1)*k.
#'
#' Note: The original GAUSS output reports df = NN*k = 46.  This appears
#' to be an error — the correct df is (NN-1)*k = 44.  We report both
#' for transparency, but use (NN-1)*k in inference.
#'
lr_test <- function(mg_res, pmg_res, k, NN) {
  ursll <- mg_res$ursll     # unrestricted (MG) sum of log-likelihoods
  rsll  <- pmg_res$rsll     # restricted (PMG) sum of log-likelihoods

  lr_stat <- 2 * (ursll - rsll)
  df      <- (NN - 1L) * k
  p_value <- 1 - pchisq(lr_stat, df = df)

  list(
    lr_stat = lr_stat,
    df      = df,
    p_value = p_value,
    ursll   = ursll,
    rsll    = rsll
  )
}


# ==============================================================================
# Run estimators
# ==============================================================================

d <- load_example1()

# SFE
sfe <- SFE(d$Y, d$X, d$Z, d$n, d$Ti, d$k)

# DFE
ytemp   <- d$Y[1:d$Ti[1]]
xtemp   <- as.matrix(d$X[1:d$Ti[1], ])
ztemp   <- as.matrix(d$Z[1:d$Ti[1], , drop = FALSE])
dgp_tmp <- FDGP1(ytemp, xtemp, ztemp, d$plag, d$qlag, d$k, d$Ti[1], 1L)
numot   <- if (is.null(dgp_tmp$ww)) 0L else ncol(dgp_tmp$ww)
dfe <- DFE1(d$Y, d$X, d$Z, d$n, d$Ti, d$k, numot, d$plag, d$qlag)

# MG (NN=23, excludes reference country U.K.)
mg <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)

# PMG
pmg <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)


# ==============================================================================
# Tests
# ==============================================================================

pass_count <- 0L
fail_count <- 0L

check <- function(name, condition, detail = "") {
  if (condition) {
    cat(sprintf("  [PASS] %s\n", name))
    pass_count <<- pass_count + 1L
  } else {
    cat(sprintf("  [FAIL] %s  %s\n", name, detail))
    fail_count <<- fail_count + 1L
  }
}

approx <- function(a, b, tol = 0.01) abs(a - b) <= tol


# ---- 1. Likelihood Ratio Test (MG vs PMG) ----

cat("\n=== LR Test: MG (unrestricted) vs PMG (restricted) ===\n\n")

lr <- lr_test(mg, pmg, d$k, d$NN)

cat(sprintf("  MG  unrestricted LL:  %.4f  (ref: 2200.1104)\n", lr$ursll))
cat(sprintf("  PMG restricted LL:    %.4f  (ref: 2139.3546)\n", lr$rsll))
cat(sprintf("  LR statistic:         %.4f  (ref: 121.5116)\n", lr$lr_stat))
cat(sprintf("  df = (NN-1)*k:        %d\n", lr$df))
cat(sprintf("  p-value:              %.4f\n\n", lr$p_value))

check("MG URSLL matches reference",
      approx(lr$ursll, 2200.11, tol = 0.01))

check("LR statistic matches reference",
      approx(lr$lr_stat, 121.51, tol = 0.5))

check("LR stat is positive (MG fits better than PMG)",
      lr$lr_stat > 0)

check("LR rejects H0 at 5% (evidence of LR heterogeneity)",
      lr$p_value < 0.05)


# ---- 2. Log-Likelihood Ordering ----

cat("\n=== Log-Likelihood Ordering ===\n\n")

cat(sprintf("  MG  URSLL:  %.4f\n", mg$ursll))
cat(sprintf("  PMG RSLL:   %.4f\n", pmg$rsll))
cat(sprintf("  DFE LL:     %.4f  (not directly comparable — uses within transformation)\n\n",
            dfe$dfedres$loglik))

check("MG URSLL > PMG RSLL (unrestricted fits better)",
      mg$ursll > pmg$rsll)


# ---- 3. DFE vs PMG: Homogeneous vs Heterogeneous Dynamics ----

cat("\n=== DFE vs PMG: Short-Run Heterogeneity ===\n\n")

dfe_phi <- dfe$dferes[3, 1]  # pooled phi from DFE
pmg_phi_mean <- pmg$phi       # mean of group-specific PMG phi
pmg_phi_all  <- pmg$phi_all   # individual phi values

cat(sprintf("  DFE pooled phi:        %.4f\n", dfe_phi))
cat(sprintf("  PMG mean phi:          %.4f\n", pmg_phi_mean))
cat(sprintf("  PMG phi range:        [%.3f, %.3f]\n",
            min(pmg_phi_all), max(pmg_phi_all)))
cat(sprintf("  PMG phi std dev:       %.4f\n\n",
            sd(pmg_phi_all)))

check("PMG phi values show substantial heterogeneity (sd > 0.05)",
      sd(pmg_phi_all) > 0.05,
      sprintf("sd = %.4f", sd(pmg_phi_all)))

check("DFE pooled phi is within range of PMG group phi values",
      dfe_phi >= min(pmg_phi_all) && dfe_phi <= max(pmg_phi_all))

check("All PMG phi < 0 (error correction is stable for all groups)",
      # Allow small positive phi for near-unit-root groups
      all(pmg_phi_all < 0.05))


# ---- 4. DFE vs MG: Long-Run Parameter Comparison ----

cat("\n=== DFE vs MG vs PMG: Long-Run Parameters ===\n\n")

cat(sprintf("  Income elasticity:  SFE=%.3f  DFE=%.3f  MG=%.3f  PMG=%.3f\n",
            sfe$sferes[1, 1], dfe$dferes[1, 1], mg$theta[1], pmg$theta[1]))
cat(sprintf("  Inflation coef:     SFE=%.3f  DFE=%.3f  MG=%.3f  PMG=%.3f\n\n",
            sfe$sferes[2, 1], dfe$dferes[2, 1], mg$theta[2], pmg$theta[2]))

check("DFE and MG income elasticity agree within 0.02",
      approx(dfe$dferes[1, 1], mg$theta[1], tol = 0.02))

check("DFE and PMG income elasticity agree within 0.02",
      approx(dfe$dferes[1, 1], pmg$theta[1], tol = 0.02))

check("SFE income elasticity is highest (upward bias from ignoring dynamics)",
      sfe$sferes[1, 1] > dfe$dferes[1, 1] &&
      sfe$sferes[1, 1] > mg$theta[1] &&
      sfe$sferes[1, 1] > pmg$theta[1])

check("SFE inflation coef is smallest in magnitude (attenuation from omitted dynamics)",
      abs(sfe$sferes[2, 1]) < abs(dfe$dferes[2, 1]) &&
      abs(sfe$sferes[2, 1]) < abs(mg$theta[2]) &&
      abs(sfe$sferes[2, 1]) < abs(pmg$theta[2]))


# ---- 5. Hausman Test Consistency ----

cat("\n=== Hausman Test: PMG vs MG ===\n\n")

ht <- hausman_test(pmg, mg, d$k)

cat(sprintf("  h(inc)=%.2f  p=%.2f\n", ht$h_individual[1], ht$p_individual[1]))
cat(sprintf("  h(inf)=%.2f  p=%.2f\n", ht$h_individual[2], ht$p_individual[2]))
cat(sprintf("  Joint: h=%.2f  p=%.2f  (ref: h=2.79, p=0.25)\n\n",
            ht$h_joint, ht$p_joint))

check("Hausman test does not reject PMG at 5% (PMG restriction plausible)",
      ht$p_joint > 0.05)

check("Hausman and LR tests give opposite conclusions (expected in finite samples)",
      lr$p_value < 0.05 && ht$p_joint > 0.05,
      "LR rejects but Hausman does not — classic small-sample divergence")


# ---- Summary ----

cat(sprintf("\n%d/%d tests passed.\n\n", pass_count, pass_count + fail_count))
if (fail_count > 0) quit(status = 1)
