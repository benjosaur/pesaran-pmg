# ==============================================================================
# Mock Panel ECM: Simulate and Estimate
#
# DGP: Error Correction Model with known parameters
#   Δy_it = phi_i * (y_{i,t-1} - theta * x_it) + gamma_i * Δx_it + mu_i + e_it
#
# True parameters:
#   theta = 0.8  (long-run elasticity)
#   phi_i ~ Uniform(-0.5, -0.1)  (heterogeneous adjustment speed)
#   gamma_i ~ Normal(0.3, 0.1)   (short-run response)
# ==============================================================================

set.seed(42)
source("panel_helpers.R")

# --- Simulation parameters ---
N  <- 3          # units: A, B, C
TT <- 40         # time points each
units <- LETTERS[1:N]

theta_true <- 0.8
phi_true   <- c(-0.35, -0.20, -0.45)  # heterogeneous phi
gamma_true <- c(0.30, 0.25, 0.40)     # heterogeneous SR response
mu_true    <- c(0.5, 1.0, 0.2)        # fixed effects
sigma_e    <- 0.02                     # error std dev

# --- Generate data ---
dates <- seq(as.Date("2020-01-01"), by = "month", length.out = TT)
date_str <- format(dates, "%Y-%m")

df_list <- list()

for (i in seq_len(N)) {
  # Exogenous regressor: random walk with drift
  x <- cumsum(rnorm(TT, mean = 0.01, sd = 0.05))

  # Generate y from ECM DGP
  y <- numeric(TT)
  y[1] <- mu_true[i] / (1 - 1 - phi_true[i]) + theta_true * x[1] + rnorm(1, 0, 0.1)
  # More stable initialization
  y[1] <- theta_true * x[1] + mu_true[i] + rnorm(1, 0, 0.05)

  for (t in 2:TT) {
    ec_term <- y[t-1] - theta_true * x[t]
    dy <- phi_true[i] * ec_term + gamma_true[i] * (x[t] - x[t-1]) + mu_true[i] * (1 + phi_true[i]) + rnorm(1, 0, sigma_e)
    y[t] <- y[t-1] + dy
  }

  df_list[[i]] <- data.frame(
    unit = units[i],
    date = date_str,
    y    = y,
    x    = x,
    stringsAsFactors = FALSE
  )
}

df <- do.call(rbind, df_list)

cat("=== Mock Panel Data ===\n\n")
cat(sprintf("Units: %s\n", paste(units, collapse = ", ")))
cat(sprintf("Periods: %s to %s (%d months)\n", date_str[1], date_str[TT], TT))
cat(sprintf("Total rows: %d\n\n", nrow(df)))
cat("Head:\n")
print(head(df, 6))
cat("\nTail:\n")
print(tail(df, 6))

cat("\n\nTrue DGP parameters:\n")
cat(sprintf("  theta (LR):  %.2f\n", theta_true))
cat(sprintf("  phi (EC):    %s\n", paste(sprintf("%.2f", phi_true), collapse = ", ")))
cat(sprintf("  gamma (SR):  %s\n", paste(sprintf("%.2f", gamma_true), collapse = ", ")))


# ==============================================================================
# Estimate via panel_helpers
# ==============================================================================

cat("\n\n=== DFE: ARDL(1,0,0) — plag=1, qlag=0 ===\n\n")
dfe_00 <- dfe_panel(df, y = "y", x = "x", unit = "unit", time = "date",
                     plag = 1, qlag = 0)
print(dfe_00)
cat("\nFull dferes matrix:\n")
print(dfe_00$dferes)

cat("\n\n=== DFE: ARDL(1,1,0) — plag=1, qlag=1 ===\n\n")
dfe_10 <- dfe_panel(df, y = "y", x = "x", unit = "unit", time = "date",
                     plag = 1, qlag = 1)
print(dfe_10)
cat("\nFull dferes matrix:\n")
print(dfe_10$dferes)

cat("\n\n=== DFE: ARDL(1,5,0) — plag=1, qlag=5 (your spec) ===\n\n")
dfe_50 <- dfe_panel(df, y = "y", x = "x", unit = "unit", time = "date",
                     plag = 1, qlag = 5)
print(dfe_50)
cat("\nFull dferes matrix:\n")
print(dfe_50$dferes)

cat("\n\n=== DFE: plag=0, qlag=5 (static y, 5 x lags — your reported spec) ===\n\n")
dfe_05 <- dfe_panel(df, y = "y", x = "x", unit = "unit", time = "date",
                     plag = 0, qlag = 5)
print(dfe_05)
cat("\nFull dferes matrix:\n")
print(dfe_05$dferes)
cat("\nNote: plag=0 forces phi=-1 (static model). theta = beta.\n")
cat("Row 1 (LR theta) and row 2 (beta on x) will be identical.\n")


# ==============================================================================
# MG and PMG for comparison
# ==============================================================================

cat("\n\n=== MG: ARDL(1,1,0) ===\n\n")
mg <- mg_panel(df, y = "y", x = "x", unit = "unit", time = "date",
               plag = 1, qlag = 1)
print(mg)
cat("\nGroup-specific phi:\n")
print(round(mg$group_phi, 3))

cat("\n\n=== PMG: ARDL(1,1,0) ===\n\n")
pmg <- pmg_panel(df, y = "y", x = "x", unit = "unit", time = "date",
                  plag = 1, qlag = 1)
print(pmg)
cat("\nGroup-specific phi:\n")
print(round(pmg$group_phi, 3))

cat("\n\n=== Summary ===\n\n")
cat(sprintf("%-20s %10s %10s %10s\n", "", "theta(LR)", "phi(EC)", "obs"))
cat(sprintf("%-20s %10.3f %10s %10s\n", "True", theta_true, "varies", ""))
cat(sprintf("%-20s %10.3f %10.3f %10d\n", "DFE ARDL(1,0,0)", dfe_00$lr$coef, dfe_00$dferes[2,1], dfe_00$qobs))
cat(sprintf("%-20s %10.3f %10.3f %10d\n", "DFE ARDL(1,1,0)", dfe_10$lr$coef, dfe_10$dferes[2,1], dfe_10$qobs))
cat(sprintf("%-20s %10.3f %10.3f %10d\n", "DFE ARDL(1,5,0)", dfe_50$lr$coef, dfe_50$dferes[2,1], dfe_50$qobs))
cat(sprintf("%-20s %10.3f %10s %10d\n", "DFE plag=0,qlag=5", dfe_05$lr$coef, "  -1 (forced)", dfe_05$qobs))
cat(sprintf("%-20s %10.3f %10.3f %10s\n", "MG ARDL(1,1,0)", mg$lr$coef, mg$phi, ""))
cat(sprintf("%-20s %10.3f %10.3f %10s\n", "PMG ARDL(1,1,0)", pmg$lr$coef, pmg$phi, ""))
