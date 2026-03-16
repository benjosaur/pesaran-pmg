# ==============================================================================
# Replication Script for Pesaran, Shin & Smith (1999)
# "Pooled Mean Group Estimation of Dynamic Heterogeneous Panels"
# JASA, 94(446), 621-634.
#
# Replicates Table 3 (Example 1: OECD Consumption Function)
# ==============================================================================

# Source all R files
source("panel_ecm_fe.R")
source("load_data.R")
source("pmg_mg.R")


# ==============================================================================
# Example 1: OECD Consumption Function (24 countries, 1960-1993)
# ==============================================================================

cat("\n")
cat("==============================================================\n")
cat("  Pesaran, Shin & Smith (1999) - Replication\n")
cat("  Example 1: OECD Consumption Function\n")
cat("==============================================================\n\n")

# Load data
d <- load_example1()
cat(sprintf("Loaded: N=%d countries, k=%d regressors\n", d$n, d$k))
cat(sprintf("Time periods by group: %s\n",
            paste(d$Ti, collapse = " ")))
cat(sprintf("Total observations (raw): %d\n\n", sum(d$Ti)))


# --- 1. Static Fixed Effects (SFE) ---
cat("--- Static Fixed Effects (SFE) ---\n")
sfe <- SFE(d$Y, d$X, d$Z, d$n, d$Ti, d$k)
cat(sprintf("  Sample size: T(1)+...+T(N) = %d\n", sfe$qobs))
cat(sprintf("  inc:  %.4f  (se=%.4f, t=%.2f)\n",
            sfe$sferes[1, 1], sfe$sferes[1, 2], sfe$sferes[1, 3]))
cat(sprintf("  inf:  %.4f  (se=%.4f, t=%.2f)\n",
            sfe$sferes[2, 1], sfe$sferes[2, 2], sfe$sferes[2, 3]))
cat("\n")


# --- 2. Dynamic Fixed Effects (DFE) ---
cat("--- Dynamic Fixed Effects (DFE) ---\n")

# Compute numot (number of SR dynamics excl. intercept) from FDGP1
ytemp   <- d$Y[1:d$Ti[1]]
xtemp   <- as.matrix(d$X[1:d$Ti[1], ])
ztemp   <- as.matrix(d$Z[1:d$Ti[1], , drop = FALSE])
dgp_tmp <- FDGP1(ytemp, xtemp, ztemp, d$plag, d$qlag, d$k, d$Ti[1], 1L)
numot   <- if (is.null(dgp_tmp$ww)) 0L else ncol(dgp_tmp$ww)

dfe <- DFE1(d$Y, d$X, d$Z, d$n, d$Ti, d$k, numot, d$plag, d$qlag)
cat(sprintf("  Sample size: T(1)+...+T(N) = %d\n", dfe$qobs))
cat("  Long-Run Coefficients:\n")
cat(sprintf("    inc:  %.4f  (se=%.4f, t=%.2f)\n",
            dfe$dferes[1, 1], dfe$dferes[1, 2], dfe$dferes[1, 3]))
cat(sprintf("    inf:  %.4f  (se=%.4f, t=%.2f)\n",
            dfe$dferes[2, 1], dfe$dferes[2, 2], dfe$dferes[2, 3]))
cat("  Error Correction:\n")
cat(sprintf("    Phi:  %.4f  (se=%.4f, t=%.2f)\n",
            dfe$dferes[3, 1], dfe$dferes[3, 2], dfe$dferes[3, 3]))
cat("\n")


# --- 3. Mean Group (MG) ---
cat("--- Mean Group (MG) ---\n")
mg <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)
cat("  Long-Run Coefficients:\n")
cat(sprintf("    inc:  %.3f  (se=%.3f, t=%.3f)\n",
            mg$theta[1], mg$theta_se[1], mg$theta_t[1]))
cat(sprintf("    inf:  %.3f  (se=%.3f, t=%.3f)\n",
            mg$theta[2], mg$theta_se[2], mg$theta_t[2]))
cat("  Error Correction:\n")
cat(sprintf("    Phi:  %.3f  (se=%.3f, t=%.3f)\n",
            mg$phi, mg$phi_se, mg$phi_t))
cat("  Short-Run Coefficients:\n")
cat(sprintf("    inc:  %.3f  (se=%.3f, t=%.3f)\n",
            mg$beta[1], mg$beta_se[1], mg$beta_t[1]))
cat(sprintf("    inf:  %.3f  (se=%.3f, t=%.3f)\n",
            mg$beta[2], mg$beta_se[2], mg$beta_t[2]))
sr_names <- c("dinc", "dinc(-1)", "Inpt")
for (j in seq_along(mg$other)) {
  nm <- if (j <= length(sr_names)) sr_names[j] else sprintf("sr%d", j)
  cat(sprintf("    %s:  %.3f  (se=%.3f, t=%.3f)\n",
              nm, mg$other[j], mg$other_se[j], mg$other_t[j]))
}
cat("\n")


# --- 4. Pooled Mean Group (PMG) via Newton-Raphson ---
cat("--- Pooled Mean Group (PMG, Newton-Raphson) ---\n")
pmg <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN,
              verbose = FALSE)
cat(sprintf("  Converged: %s after %d iterations\n",
            pmg$converged, pmg$iterations))
cat(sprintf("  Restricted log-likelihood: %.4f\n", pmg$loglik))
cat("  Long-Run Coefficients:\n")
cat(sprintf("    inc:  %.4f  (se=%.4f, t=%.3f)\n",
            pmg$theta[1], pmg$theta_se[1], pmg$theta_t[1]))
cat(sprintf("    inf:  %.4f  (se=%.4f, t=%.3f)\n",
            pmg$theta[2], pmg$theta_se[2], pmg$theta_t[2]))
cat("  Error Correction:\n")
cat(sprintf("    Phi:  %.3f  (se=%.3f, t=%.3f)\n",
            pmg$phi, pmg$phi_se, pmg$phi_t))
cat("  Short-Run Coefficients:\n")
cat(sprintf("    inc:  %.3f  (se=%.3f, t=%.3f)\n",
            pmg$beta[1], pmg$beta_se[1], pmg$beta_t[1]))
cat(sprintf("    inf:  %.3f  (se=%.3f, t=%.3f)\n",
            pmg$beta[2], pmg$beta_se[2], pmg$beta_t[2]))
for (j in seq_along(pmg$other)) {
  nm <- if (j <= length(sr_names)) sr_names[j] else sprintf("sr%d", j)
  cat(sprintf("    %s:  %.3f  (se=%.3f, t=%.3f)\n",
              nm, pmg$other[j], pmg$other_se[j], pmg$other_t[j]))
}
cat("\n")


# --- 5. Hausman Test (PMG vs MG) ---
cat("--- Hausman Test (PMG vs MG) ---\n")
ht <- hausman_test(pmg, mg, d$k)
cat(sprintf("  inc:  h=%.2f, p=%.2f\n", ht$h_individual[1], ht$p_individual[1]))
cat(sprintf("  inf:  h=%.2f, p=%.2f\n", ht$h_individual[2], ht$p_individual[2]))
cat(sprintf("  Joint: h=%.2f, p=%.2f  (df=%d)\n",
            ht$h_joint, ht$p_joint, d$k))
cat("\n")


# --- Comparison Table ---
cat("==============================================================\n")
cat("  Comparison of Long-Run Coefficients\n")
cat("  (cf. Table 3 / OUTPUT.OUT from the paper)\n")
cat("==============================================================\n\n")
cat("              PMG              MG             SFE             DFE\n")
cat("        Coef.   (SE)     Coef.   (SE)    Coef.   (SE)    Coef.   (SE)\n")
cat(sprintf("  inc  %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)\n",
            pmg$theta[1], pmg$theta_se[1],
            mg$theta[1], mg$theta_se[1],
            sfe$sferes[1, 1], sfe$sferes[1, 2],
            dfe$dferes[1, 1], dfe$dferes[1, 2]))
cat(sprintf("  inf  %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)\n",
            pmg$theta[2], pmg$theta_se[2],
            mg$theta[2], mg$theta_se[2],
            sfe$sferes[2, 1], sfe$sferes[2, 2],
            dfe$dferes[2, 1], dfe$dferes[2, 2]))
cat("\nReference (from OUTPUT.OUT):\n")
cat("  inc   0.915 (0.008)   0.916 (0.028)   0.953 (0.007)   0.925 (0.019)\n")
cat("  inf  -0.501 (0.052)  -0.305 (0.128)  -0.095 (0.027)  -0.307 (0.076)\n")
cat("\n")
