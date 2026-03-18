# ==============================================================================
# Replication via data.frame API
#
# Converts the raw .DAT files into a standard data.frame, then runs
# all estimators through the panel_helpers.R convenience wrappers.
# Verifies that the data.frame path produces identical results to
# the raw stacked-vector path (replicate.R / replicate_tables.R).
# ==============================================================================

source("panel_helpers.R")


# ==============================================================================
# Build data.frame from raw .DAT files (Example 1)
# ==============================================================================

cat("=== Building data.frame from raw .DAT files ===\n\n")

MAXT <- 34L; MAXN <- 24L; MISSING <- 8934567

lpc  <- scan("prog/LPC.DAT",  quiet = TRUE, n = MAXT * MAXN)
lndi <- scan("prog/LNDI.DAT", quiet = TRUE, n = MAXT * MAXN)
dp   <- scan("prog/DP.DAT",   quiet = TRUE, n = MAXT * MAXN)

countries <- c(
  "Australi", "Austria",  "Belgium",  "Canada",   "Denmark",
  "Finland",  "France",   "Germany",  "Greece",   "Iceland",
  "Ireland",  "Italy",    "Japan",    "Luxembou", "Netherla",
  "NewZeala", "Norway",   "Portugal", "Spain",    "Sweden",
  "Switzerl", "Turkey",   "U.S.",     "U.K."
)

df <- data.frame(
  country = rep(countries, each = MAXT),
  year    = rep(1960:(1960 + MAXT - 1), MAXN),
  lpc     = lpc,
  lndi    = lndi,
  dp      = dp
)

# Replace sentinel with NA
df$lpc[df$lpc   == MISSING] <- NA
df$lndi[df$lndi == MISSING] <- NA
df$dp[df$dp     == MISSING] <- NA

cat(sprintf("data.frame: %d rows x %d cols\n", nrow(df), ncol(df)))
cat(sprintf("Countries: %d\n", length(unique(df$country)))  )
cat(sprintf("Years: %d-%d\n", min(df$year), max(df$year)))
cat(sprintf("NAs: lpc=%d, lndi=%d, dp=%d\n",
    sum(is.na(df$lpc)), sum(is.na(df$lndi)), sum(is.na(df$dp))))
cat("\nHead:\n")
print(head(df, 5))
cat("\nTail:\n")
print(tail(df, 5))


# ==============================================================================
# Test 1: Table 1 spec — ARDL(1,1,1), all 24 groups
# ==============================================================================

cat("\n\n=== Table 1: ARDL(1,1,1), N=24, via data.frame API ===\n\n")

pmg_t1 <- pmg_panel(df, y = "lpc", x = c("lndi", "dp"),
                     unit = "country", time = "year",
                     plag = 1, qlag = c(1, 1))
mg_t1  <- mg_panel(df, y = "lpc", x = c("lndi", "dp"),
                    unit = "country", time = "year",
                    plag = 1, qlag = c(1, 1))
dfe_t1 <- dfe_panel(df, y = "lpc", x = c("lndi", "dp"),
                     unit = "country", time = "year",
                     plag = 1, qlag = c(1, 1))
sfe_t1 <- sfe_panel(df, y = "lpc", x = c("lndi", "dp"),
                     unit = "country", time = "year")

cat("PMG:\n"); print(pmg_t1)
cat("\nMG:\n");  print(mg_t1)
cat("\nDFE:\n"); print(dfe_t1)
cat("\nSFE:\n"); print(sfe_t1)


# ==============================================================================
# Test 2: OUTPUT.OUT spec — ARDL(1,2,0), relative variable, NN=23
# ==============================================================================

cat("\n\n=== OUTPUT.OUT: ARDL(1,2,0), ref_group='U.K.', via data.frame API ===\n\n")

pmg_out <- pmg_panel(df, y = "lpc", x = c("lndi", "dp"),
                      unit = "country", time = "year",
                      plag = 1, qlag = c(2, 0), ref_group = "U.K.")
mg_out  <- mg_panel(df, y = "lpc", x = c("lndi", "dp"),
                     unit = "country", time = "year",
                     plag = 1, qlag = c(2, 0), ref_group = "U.K.")

cat("PMG:\n"); print(pmg_out)
cat("\nMG:\n");  print(mg_out)


# ==============================================================================
# Test 3: Verify data.frame results match raw API results
# ==============================================================================

cat("\n\n=== Verification: data.frame API vs raw API ===\n\n")

source("load_data.R")
d <- load_example1()

# Table 1 via raw API
plag_111 <- rep(1L, d$n)
qlag_111 <- cbind(rep(1L, d$n), rep(1L, d$n))
mg_raw   <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, plag_111, qlag_111, NN = d$n)
pmg_raw  <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, plag_111, qlag_111, NN = d$n)
sfe_raw  <- SFE(d$Y, d$X, d$Z, d$n, d$Ti, d$k)

pass <- 0L; fail <- 0L
check <- function(name, a, b, tol = 1e-6) {
  ok <- all(abs(a - b) <= tol)
  if (ok) { pass <<- pass + 1L } else { fail <<- fail + 1L }
  cat(sprintf("  [%s] %s: max diff = %.2e\n", if (ok) "PASS" else "FAIL",
              name, max(abs(a - b))))
}

cat("Table 1 (ARDL(1,1,1), N=24):\n")
check("MG theta",       mg_t1$theta,    mg_raw$theta)
check("MG theta_se",    mg_t1$theta_se, mg_raw$theta_se)
check("MG phi",         mg_t1$phi,      mg_raw$phi)
check("PMG theta",      pmg_t1$theta,   pmg_raw$theta)
check("PMG theta_se",   pmg_t1$theta_se, pmg_raw$theta_se)
check("PMG phi",        pmg_t1$phi,     pmg_raw$phi)
check("SFE coefficients", sfe_t1$lr$coef, as.numeric(sfe_raw$sferes[1:2, 1]))
check("SFE obs",        sfe_t1$qobs,    sfe_raw$qobs)

# OUTPUT.OUT via raw API
mg_raw2  <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)
pmg_raw2 <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)

cat("\nOUTPUT.OUT (ARDL(1,2,0), NN=23):\n")
check("MG theta",     mg_out$theta,    mg_raw2$theta)
check("MG theta_se",  mg_out$theta_se, mg_raw2$theta_se)
check("PMG theta",    pmg_out$theta,   pmg_raw2$theta)
check("PMG theta_se", pmg_out$theta_se, pmg_raw2$theta_se)
check("PMG phi",      pmg_out$phi,     pmg_raw2$phi)
check("PMG loglik",   pmg_out$loglik,  pmg_raw2$loglik)

cat(sprintf("\n%d/%d checks passed.\n", pass, pass + fail))
if (fail > 0) quit(status = 1)


# ==============================================================================
# Example 2: Energy Demand (show it works with different data)
# ==============================================================================

cat("\n\n=== Example 2: Asian Energy Demand via data.frame ===\n\n")

MAXT2 <- 18L; MAXN2 <- 10L
eda <- matrix(scan("prog/EDA.DAT", quiet = TRUE, n = MAXT2 * MAXN2 * 4),
              ncol = 4, byrow = TRUE)

# Try reading country names
j2names <- tryCatch(
  { nn <- trimws(readLines("prog/J2name.dat", warn = FALSE))
    nn[nchar(nn) > 0][1:MAXN2] },
  error = function(e) paste0("Country", 1:MAXN2)
)

df2 <- data.frame(
  country = rep(j2names, each = MAXT2),
  year    = rep(1974:(1974 + MAXT2 - 1), MAXN2),
  lec     = log(eda[, 3] / eda[, 2]),   # ln(TOTC/POPUL)
  lgdp    = log(eda[, 4]),              # ln(RGDPL)
  lprice  = log(eda[, 1])               # ln(AERP)
)

cat("Energy demand data.frame:\n")
print(str(df2))
cat("\n")

pmg_e <- pmg_panel(df2, y = "lec", x = c("lgdp", "lprice"),
                    unit = "country", time = "year",
                    plag = 1, qlag = c(0, 0))
print(pmg_e)

cat("\nPaper Table 3 reference:\n")
cat("  θ₁ (output)  1.184 (0.039)\n")
cat("  θ₂ (price)  -0.339 (0.033)\n")
cat("  φ           -0.298 (0.063)\n")
cat("\n")
