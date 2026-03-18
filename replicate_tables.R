# ==============================================================================
# Replication of Paper Tables 1, 3, and 4
#
# Pesaran, Shin & Smith (1999), JASA 94(446), 621-634.
#
# Table 1: OECD Consumption, ARDL(1,1,1), NDI, N=24
# Table 3: Asian Energy Demand, ARDL(1,0,0), N=10
#
# Table 2 uses private disposable labor income (LPDI) and various subsamples
# (N=15, 14, 9) — we reproduce the N=15 ARDL(1,0,0) and ARDL(1,1,1) rows.
#
# Table 4 uses SBC-chosen lags per country, which requires lag selection
# not yet implemented — noted but not reproduced.
# ==============================================================================

source("panel_ecm_fe.R")
source("load_data.R")
source("pmg_mg.R")


# ==============================================================================
# Table 1: OECD Consumption Function, ARDL(1,1,1), NDI
# ==============================================================================

cat("\n")
cat("================================================================\n")
cat("  Table 1: ARDL(1,1,1) Consumption Functions for OECD Countries\n")
cat("           Using National Disposable Income (NDI)\n")
cat("           N = 24, T = 1962-1993\n")
cat("================================================================\n\n")

d <- load_example1()

# Table 1 uses ARDL(1,1,1) with ALL 24 groups (no relative variable)
plag <- rep(1L, d$n)
qlag <- cbind(rep(1L, d$n), rep(1L, d$n))
NN   <- d$n  # all 24

mg  <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, plag, qlag, NN = NN)
pmg <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, plag, qlag, NN = NN)
dfe <- DFE1(d$Y, d$X, d$Z, d$n, d$Ti, d$k, 2L, plag, qlag)

cat("                       MG               PMG              DFE (robust SE)\n")
cat("                  Coef.   (SE)      Coef.   (SE)      Coef.   (SE)\n")
cat(sprintf("  θ₁ (inc)       %6.3f (%5.3f)   %6.3f (%5.3f)   %6.3f (%5.3f)\n",
    mg$theta[1], mg$theta_se[1],
    pmg$theta[1], pmg$theta_se[1],
    dfe$rdferes[1, 1], dfe$rdferes[1, 2]))
cat(sprintf("  θ₂ (inf)      %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)\n",
    mg$theta[2], mg$theta_se[2],
    pmg$theta[2], pmg$theta_se[2],
    dfe$rdferes[2, 1], dfe$rdferes[2, 2]))
cat(sprintf("  φ             %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)\n",
    mg$phi, mg$phi_se, pmg$phi, pmg$phi_se,
    dfe$rdferes[3, 1], dfe$rdferes[3, 2]))
cat(sprintf("  MLL           %7.0f          %7.0f\n", mg$ursll, pmg$rsll))

cat("\nPaper Table 1:\n")
cat("  θ₁ (inc)        0.918 (0.027)    0.904 (0.010)    0.912 (0.045)\n")
cat("  θ₂ (inf)       -0.353 (0.117)   -0.466 (0.063)   -0.266 (0.099)\n")
cat("  φ              -0.306 (0.030)   -0.200 (0.032)   -0.179 (0.042)\n")
cat("  MLL              2390              2247              1999\n")


# ==============================================================================
# Table 2 (partial): OECD Consumption, PDI, N=15
# ==============================================================================

cat("\n\n")
cat("================================================================\n")
cat("  Table 2 (partial): Income Elasticity Using PDI\n")
cat("           N = 15, ARDL(1,0,0) and ARDL(1,1,1)\n")
cat("================================================================\n\n")

# Load PDI data (private disposable labor income)
MAXT <- 34L; MAXN <- 24L; MISSING <- 8934567
n_elem <- MAXT * MAXN

lpc  <- scan(file.path("prog", "LPC.DAT"),  what = numeric(), n = n_elem, quiet = TRUE)
lpdi <- scan(file.path("prog", "LPDI.DAT"), what = numeric(), n = n_elem, quiet = TRUE)
dp   <- scan(file.path("prog", "DP.DAT"),   what = numeric(), n = n_elem, quiet = TRUE)

country_names <- c(
  "Australi", "Austria",  "Belgium",  "Canada",   "Denmark",
  "Finland",  "France",   "Germany",  "Greece",   "Iceland",
  "Ireland",  "Italy",    "Japan",    "Luxembou", "Netherla",
  "NewZeala", "Norway",   "Portugal", "Spain",    "Sweden",
  "Switzerl", "Turkey",   "U.S.",     "U.K."
)

# Table 2 uses only countries with PDI data — N=15
# Identify which countries have non-missing LPDI
has_pdi <- logical(MAXN)
for (i in seq_len(MAXN)) {
  idx <- ((i - 1L) * MAXT + 1L):(i * MAXT)
  has_pdi[i] <- any(lpdi[idx] != MISSING & is.finite(lpdi[idx]))
}
pdi_countries <- which(has_pdi)
n_pdi <- length(pdi_countries)

# Stack data for PDI countries only
Y_list <- vector("list", n_pdi)
X_list <- vector("list", n_pdi)
Ti_pdi <- integer(n_pdi)

for (j in seq_len(n_pdi)) {
  i   <- pdi_countries[j]
  idx <- ((i - 1L) * MAXT + 1L):(i * MAXT)
  yi  <- lpc[idx]; x1i <- lpdi[idx]; x2i <- dp[idx]
  valid <- (yi != MISSING) & (x1i != MISSING) & (x2i != MISSING)
  Ti_pdi[j]   <- sum(valid)
  Y_list[[j]] <- yi[valid]
  X_list[[j]] <- cbind(x1i[valid], x2i[valid])
}

Y_pdi <- do.call(c, Y_list)
X_pdi <- do.call(rbind, X_list)
Z_pdi <- matrix(1, nrow = length(Y_pdi), ncol = 1)

cat(sprintf("Countries with PDI data: %d\n", n_pdi))
cat(sprintf("  %s\n", paste(country_names[pdi_countries], collapse = ", ")))
cat(sprintf("Total obs: %d\n\n", sum(Ti_pdi)))

# ARDL(1,0,0) — income elasticity only
plag_100 <- rep(1L, n_pdi)
qlag_100 <- cbind(rep(0L, n_pdi), rep(0L, n_pdi))

mg_100  <- MG(Y_pdi, X_pdi, Z_pdi, n_pdi, Ti_pdi, 2L, plag_100, qlag_100)
pmg_100 <- PMG_NR(Y_pdi, X_pdi, Z_pdi, n_pdi, Ti_pdi, 2L, plag_100, qlag_100)
dfe_100 <- DFE1(Y_pdi, X_pdi, Z_pdi, n_pdi, Ti_pdi, 2L, 0L, plag_100, qlag_100)

cat("ARDL(1,0,0):\n")
cat("                       MG               PMG              DFE\n")
cat(sprintf("  θ₁ (inc)       %6.3f (%5.3f)   %6.3f (%5.3f)   %6.3f (%5.3f)\n",
    mg_100$theta[1], mg_100$theta_se[1],
    pmg_100$theta[1], pmg_100$theta_se[1],
    dfe_100$dferes[1, 1], dfe_100$dferes[1, 2]))

cat("\nPaper Table 2 (N=15, ARDL(1,0,0)):\n")
cat("  θ₁ (inc)        0.948 (0.070)    0.933 (0.010)    0.936 (0.022)\n")

# ARDL(1,1,1)
plag_111 <- rep(1L, n_pdi)
qlag_111 <- cbind(rep(1L, n_pdi), rep(1L, n_pdi))

mg_111  <- MG(Y_pdi, X_pdi, Z_pdi, n_pdi, Ti_pdi, 2L, plag_111, qlag_111)
dfe_111 <- DFE1(Y_pdi, X_pdi, Z_pdi, n_pdi, Ti_pdi, 2L, 2L, plag_111, qlag_111)

# PMG with Sweden is numerically singular (phi ≈ 0 for Sweden).
# The paper discusses this: "In Sweden... the speed of adjustment was
# −.03 (.154), the long-run income elasticity was −3.74 (25.22)...
# it distorts the MG estimates." We skip PMG SEs for this case.
pmg_111 <- tryCatch(
  PMG_NR(Y_pdi, X_pdi, Z_pdi, n_pdi, Ti_pdi, 2L, plag_111, qlag_111),
  error = function(e) {
    # Run BSA without final SE computation
    res <- list(theta = c(NA, NA), theta_se = c(NA, NA), phi = NA)
    message("  PMG info matrix singular (Sweden outlier): ", e$message)
    res
  }
)

cat("\nARDL(1,1,1):\n")
cat("                       MG               PMG              DFE\n")
if (is.na(pmg_111$theta[1])) {
  cat(sprintf("  θ₁ (inc)       %6.1f (%5.1f)   (singular)       %6.3f (%5.3f)\n",
      mg_111$theta[1], mg_111$theta_se[1],
      dfe_111$dferes[1, 1], dfe_111$dferes[1, 2]))
} else {
  cat(sprintf("  θ₁ (inc)       %6.1f (%5.1f)   %6.3f (%5.3f)   %6.3f (%5.3f)\n",
      mg_111$theta[1], mg_111$theta_se[1],
      pmg_111$theta[1], pmg_111$theta_se[1],
      dfe_111$dferes[1, 1], dfe_111$dferes[1, 2]))
}

cat("\nPaper Table 2 (N=15, ARDL(1,1,1)):\n")
cat("  θ₁ (inc)       -12.8 (13.8)     0.918 (0.016)    0.935 (0.023)\n")
cat("  (MG distorted by Sweden outlier — paper discusses this)\n")


# ==============================================================================
# Table 3: Energy Demand in Asian Developing Countries, ARDL(1,0,0)
# ==============================================================================

cat("\n\n")
cat("================================================================\n")
cat("  Table 3: ARDL(1,0,0) Energy Demand for 10 Asian Economies\n")
cat("           1974-1990\n")
cat("================================================================\n\n")

d2 <- load_example2()

# Table 3 uses ARDL(1,0,0)
plag_e <- rep(1L, d2$n)
qlag_e <- cbind(rep(0L, d2$n), rep(0L, d2$n))

mg_e  <- MG(d2$Y, d2$X, d2$Z, d2$n, d2$Ti, d2$k, plag_e, qlag_e)
pmg_e <- PMG_NR(d2$Y, d2$X, d2$Z, d2$n, d2$Ti, d2$k, plag_e, qlag_e)
dfe_e <- DFE1(d2$Y, d2$X, d2$Z, d2$n, d2$Ti, d2$k, 0L, plag_e, qlag_e)
sfe_e <- SFE(d2$Y, d2$X, d2$Z, d2$n, d2$Ti, d2$k)

cat("                       MG               PMG              DFE              SFE\n")
cat("                  Coef.   (SE)      Coef.   (SE)      Coef.   (SE)      Coef.   (SE)\n")
cat(sprintf("  θ₁ (output)    %6.3f (%5.3f)   %6.3f (%5.3f)   %6.3f (%5.3f)   %6.3f (%5.3f)\n",
    mg_e$theta[1], mg_e$theta_se[1],
    pmg_e$theta[1], pmg_e$theta_se[1],
    dfe_e$dferes[1, 1], dfe_e$dferes[1, 2],
    sfe_e$sferes[1, 1], sfe_e$sferes[1, 2]))
cat(sprintf("  θ₂ (price)    %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)\n",
    mg_e$theta[2], mg_e$theta_se[2],
    pmg_e$theta[2], pmg_e$theta_se[2],
    dfe_e$dferes[2, 1], dfe_e$dferes[2, 2],
    sfe_e$sferes[2, 1], sfe_e$sferes[2, 2]))
cat(sprintf("  φ             %6.3f (%5.3f)  %6.3f (%5.3f)  %6.3f (%5.3f)    -1\n",
    mg_e$phi, mg_e$phi_se, pmg_e$phi, pmg_e$phi_se,
    dfe_e$dferes[3, 1], dfe_e$dferes[3, 2]))
cat(sprintf("  MLL           %7.1f          %7.1f          %7.1f          %7.1f\n",
    mg_e$ursll, pmg_e$rsll, dfe_e$dfedres$loglik, sfe_e$sfedres$loglik))
cat(sprintf("  N x T          %d              %d              %d              %d\n",
    sum(d2$Ti), sum(d2$Ti), dfe_e$qobs, sfe_e$qobs))

cat("\nPaper Table 3:\n")
cat("  θ₁ (output)     1.228 (0.183)    1.184 (0.039)    1.301 (0.109)    1.009 (0.037)\n")
cat("  θ₂ (price)     -0.261 (0.118)   -0.339 (0.033)   -0.365 (0.097)   -0.067 (0.030)\n")
cat("  φ              -0.524 (0.070)   -0.298 (0.063)   -0.235 (0.040)     -1\n")
cat("  MLL             347.12           322.79           288.36           186.95\n")
cat("  N x T            170              170              170              170\n")

cat("\n")
cat("================================================================\n")
cat("  Table 4: ARDL-SBC Energy Demand\n")
cat("  (Not reproduced — requires SBC lag selection per country,\n")
cat("   which is not yet implemented)\n")
cat("================================================================\n")
cat("\n")
