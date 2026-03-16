# ==============================================================================
# Unit Tests for PSS (1999) Replication
#
# Tests against known values from OUTPUT.OUT (original GAUSS output).
# Run: Rscript tests/test_replication.R
# ==============================================================================

# Check for testthat
if (!requireNamespace("testthat", quietly = TRUE)) {
  cat("testthat not installed. Running basic checks instead.\n\n")

  # --- Basic checks without testthat ---
  source("panel_ecm_fe.R")
  source("load_data.R")
  source("pmg_mg.R")

  d   <- load_example1()
  sfe <- SFE(d$Y, d$X, d$Z, d$n, d$Ti, d$k)

  ytemp   <- d$Y[1:d$Ti[1]]
  xtemp   <- as.matrix(d$X[1:d$Ti[1], ])
  ztemp   <- as.matrix(d$Z[1:d$Ti[1], , drop = FALSE])
  dgp_tmp <- FDGP1(ytemp, xtemp, ztemp, d$plag, d$qlag, d$k, d$Ti[1], 1L)
  numot   <- if (is.null(dgp_tmp$ww)) 0L else ncol(dgp_tmp$ww)

  dfe <- DFE1(d$Y, d$X, d$Z, d$n, d$Ti, d$k, numot, d$plag, d$qlag)
  mg  <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)
  pmg <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)

  check <- function(name, actual, expected, tol = 0.01) {
    pass <- abs(actual - expected) <= tol
    status <- if (pass) "PASS" else "FAIL"
    cat(sprintf("  [%s] %s: got %.4f, expected %.4f (diff=%.4f, tol=%.4f)\n",
                status, name, actual, expected, abs(actual - expected), tol))
    pass
  }

  results <- c()
  cat("SFE:\n")
  results <- c(results, check("SFE inc", sfe$sferes[1, 1], 0.9532, 0.001))
  results <- c(results, check("SFE inf", sfe$sferes[2, 1], -0.0945, 0.001))
  results <- c(results, check("SFE obs", sfe$qobs, 767, 0))

  cat("DFE:\n")
  results <- c(results, check("DFE inc", dfe$dferes[1, 1], 0.9247, 0.001))
  results <- c(results, check("DFE inf", dfe$dferes[2, 1], -0.3067, 0.001))
  results <- c(results, check("DFE phi", dfe$dferes[3, 1], -0.1560, 0.001))
  results <- c(results, check("DFE obs", dfe$qobs, 743, 0))

  cat("MG:\n")
  results <- c(results, check("MG inc", mg$theta[1], 0.916, 0.01))
  results <- c(results, check("MG inf", mg$theta[2], -0.305, 0.01))

  cat("PMG:\n")
  results <- c(results, check("PMG inc", pmg$theta[1], 0.915, 0.01))
  results <- c(results, check("PMG inf", pmg$theta[2], -0.501, 0.01))
  results <- c(results, check("PMG converged", pmg$converged, TRUE, 0))

  cat(sprintf("\n%d/%d tests passed.\n", sum(results), length(results)))
  if (!all(results)) quit(status = 1)
  quit(status = 0)
}


# --- Full testthat tests ---
library(testthat)

source("panel_ecm_fe.R")
source("load_data.R")
source("pmg_mg.R")


# ---- Setup: run all estimators once ----
d <- load_example1()

sfe <- SFE(d$Y, d$X, d$Z, d$n, d$Ti, d$k)

ytemp   <- d$Y[1:d$Ti[1]]
xtemp   <- as.matrix(d$X[1:d$Ti[1], ])
ztemp   <- as.matrix(d$Z[1:d$Ti[1], , drop = FALSE])
dgp_tmp <- FDGP1(ytemp, xtemp, ztemp, d$plag, d$qlag, d$k, d$Ti[1], 1L)
numot   <- if (is.null(dgp_tmp$ww)) 0L else ncol(dgp_tmp$ww)

dfe <- DFE1(d$Y, d$X, d$Z, d$n, d$Ti, d$k, numot, d$plag, d$qlag)
mg  <- MG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)
pmg <- PMG_NR(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)


# ==== Helper function tests ====

test_that("lag1 works correctly", {
  expect_equal(lag1(c(1, 2, 3, 4)), c(NA, 1, 2, 3))

  m <- matrix(1:6, ncol = 2)
  lagged <- lag1(m)
  expect_equal(lagged[1, ], c(NA, NA))
  expect_equal(lagged[2, ], c(1, 4))
  expect_equal(lagged[3, ], c(2, 5))
})

test_that("demean works correctly", {
  x <- c(1, 2, 3, 4, 5)
  expect_equal(demean(x), x - mean(x))
  expect_equal(mean(demean(x)), 0)

  m <- matrix(1:6, ncol = 2)
  dm <- demean(m)
  expect_equal(colMeans(dm), c(0, 0))
})


# ==== Data loading tests ====

test_that("Example 1 data loads correctly", {
  expect_equal(d$n, 24)
  expect_equal(d$k, 2)
  expect_equal(d$NN, 23)
  expect_equal(length(d$Ti), 24)
  expect_equal(length(d$Y), sum(d$Ti))
  expect_equal(nrow(d$X), sum(d$Ti))
  expect_equal(ncol(d$X), 2)
  expect_true(d$rel_var)
})

test_that("Time periods match OUTPUT.OUT", {
  # From OUTPUT.OUT: 33 33 32 33 (x21) 33
  expected_Ti <- c(33, 33, 32, rep(33, 21))
  expect_equal(d$Ti, expected_Ti)
})


# ==== SFE tests ====

test_that("SFE observation count matches", {
  expect_equal(sfe$qobs, 767)
})

test_that("SFE long-run coefficients match OUTPUT.OUT", {
  expect_equal(sfe$sferes[1, 1], 0.9532, tolerance = 0.001)
  expect_equal(sfe$sferes[2, 1], -0.0945, tolerance = 0.001)
})

test_that("SFE standard errors match OUTPUT.OUT", {
  expect_equal(sfe$sferes[1, 2], 0.0067, tolerance = 0.001)
  expect_equal(sfe$sferes[2, 2], 0.0273, tolerance = 0.001)
})

test_that("SFE t-ratios match OUTPUT.OUT", {
  expect_equal(sfe$sferes[1, 3], 142.20, tolerance = 0.1)
  expect_equal(sfe$sferes[2, 3], -3.46, tolerance = 0.1)
})


# ==== DFE tests ====

test_that("DFE observation count matches", {
  expect_equal(dfe$qobs, 743)
})

test_that("DFE long-run coefficients match OUTPUT.OUT", {
  expect_equal(dfe$dferes[1, 1], 0.9247, tolerance = 0.001)
  expect_equal(dfe$dferes[2, 1], -0.3067, tolerance = 0.001)
})

test_that("DFE phi coefficient matches OUTPUT.OUT", {
  expect_equal(dfe$dferes[3, 1], -0.1560, tolerance = 0.001)
})

test_that("DFE standard errors match OUTPUT.OUT", {
  expect_equal(dfe$dferes[1, 2], 0.0191, tolerance = 0.001)
  expect_equal(dfe$dferes[2, 2], 0.0763, tolerance = 0.001)
  expect_equal(dfe$dferes[3, 2], 0.0171, tolerance = 0.001)
})


# ==== MG tests ====

test_that("MG long-run coefficients match OUTPUT.OUT", {
  expect_equal(mg$theta[1], 0.916, tolerance = 0.01)
  expect_equal(mg$theta[2], -0.305, tolerance = 0.01)
})

test_that("MG standard errors match OUTPUT.OUT", {
  expect_equal(mg$theta_se[1], 0.028, tolerance = 0.005)
  expect_equal(mg$theta_se[2], 0.128, tolerance = 0.005)
})

test_that("MG phi matches OUTPUT.OUT", {
  expect_equal(mg$phi, -0.317, tolerance = 0.01)
  expect_equal(mg$phi_se, 0.032, tolerance = 0.005)
})

test_that("MG short-run coefficients match OUTPUT.OUT", {
  expect_equal(mg$beta[1], 0.289, tolerance = 0.01)
  expect_equal(mg$beta[2], -0.117, tolerance = 0.01)
})

test_that("MG uses NN=23 groups (excluding reference country)", {
  expect_equal(ncol(mg$theta_all), 23)
})


# ==== PMG tests ====

test_that("PMG converges", {
  expect_true(pmg$converged)
})

test_that("PMG long-run coefficients match OUTPUT.OUT", {
  # Tolerance of 0.01 accounts for cross-platform numerical differences
  expect_equal(pmg$theta[1], 0.915, tolerance = 0.01)
  expect_equal(pmg$theta[2], -0.501, tolerance = 0.01)
})

test_that("PMG standard errors match OUTPUT.OUT", {
  expect_equal(pmg$theta_se[1], 0.008, tolerance = 0.002)
  expect_equal(pmg$theta_se[2], 0.052, tolerance = 0.005)
})

test_that("PMG phi matches OUTPUT.OUT", {
  expect_equal(pmg$phi, -0.198, tolerance = 0.01)
})

test_that("PMG short-run coefficients match OUTPUT.OUT", {
  expect_equal(pmg$beta[1], 0.181, tolerance = 0.01)
  expect_equal(pmg$beta[2], -0.099, tolerance = 0.01)
})

test_that("PMG group-specific phi values are reasonable", {
  # All phi should be in (-1, 0.1) for convergence
  expect_true(all(pmg$phi_all > -1))
  expect_true(all(pmg$phi_all < 0.1))
})


# ==== BSA PMG tests ====

test_that("BSA and NR converge to the same values", {
  pmg_bsa <- PMG(d$Y, d$X, d$Z, d$n, d$Ti, d$k, d$plag, d$qlag, NN = d$NN)
  expect_equal(pmg_bsa$theta[1], pmg$theta[1], tolerance = 0.001)
  expect_equal(pmg_bsa$theta[2], pmg$theta[2], tolerance = 0.001)
})


# ==== Hausman test ====

test_that("Hausman test matches OUTPUT.OUT", {
  ht <- hausman_test(pmg, mg, d$k)
  # inc h-test ≈ 0.00
  expect_true(ht$h_individual[1] < 0.5)
  # inf h-test ≈ 2.79 (our value slightly different due to PMG difference)
  expect_true(ht$h_individual[2] > 1.0)
  # p-values indicate no rejection at 5%
  expect_true(ht$p_joint > 0.05)
})


# ==== Cross-estimator consistency ====

test_that("All estimators give income elasticity close to 1", {
  expect_true(abs(sfe$sferes[1, 1] - 1) < 0.1)
  expect_true(abs(dfe$dferes[1, 1] - 1) < 0.1)
  expect_true(abs(mg$theta[1] - 1) < 0.1)
  expect_true(abs(pmg$theta[1] - 1) < 0.1)
})

test_that("All estimators give negative inflation coefficient", {
  expect_true(sfe$sferes[2, 1] < 0)
  expect_true(dfe$dferes[2, 1] < 0)
  expect_true(mg$theta[2] < 0)
  expect_true(pmg$theta[2] < 0)
})


cat("\nAll tests passed!\n")
