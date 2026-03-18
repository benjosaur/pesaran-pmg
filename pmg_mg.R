# ==============================================================================
# Mean Group (MG) and Pooled Mean Group (PMG) Estimators
#
# Based on: Pesaran, Shin & Smith (1999)
# "Pooled Mean Group Estimation of Dynamic Heterogeneous Panels"
# Journal of the American Statistical Association, 94(446), 621-634.
#
# Translated from the original GAUSS implementation (SUB1.PRC) to R.
#
# Provides:
#   MG()    - Mean Group estimator
#   PMG()   - Pooled Mean Group estimator (Back-Substitution Algorithm)
#
# Requires: lag1() from panel_ecm_fe.R
# ==============================================================================


# --- Data Generation Helper --------------------------------------------------

#' DGP1: Construct estimation variables for group i (MG/PMG version)
#'
#' Unlike FDGP1 in panel_ecm_fe.R, this function includes ALL of Z
#' (including the intercept) in the short-run regressors ww. This is
#' needed for MG/PMG which do NOT use the within-transformation.
#'
DGP1 <- function(ytemp, xtemp, ztemp, plag, qlag, k, mxi, tis) {

  pp <- plag[mxi]

  y1temp <- lag1(ytemp)
  dytemp <- ytemp - y1temp
  x1temp <- lag1(xtemp)
  dxtemp <- xtemp - x1temp

  ww    <- NULL
  wflag <- 0L

  # Lagged Delta-y terms (for ARDL dynamics)
  if (pp >= 2L) {
    dy_lags <- matrix(NA_real_, tis, pp - 1L)
    tmp <- dytemp
    for (j in seq_len(pp - 1L)) {
      tmp <- lag1(tmp)
      dy_lags[, j] <- tmp
    }
    ww    <- dy_lags
    wflag <- 1L
  }

  # Current and lagged Delta-x terms
  dx_all <- NULL
  for (j in seq_len(k)) {
    q_ij <- qlag[mxi, j]
    if (q_ij >= 1L) {
      dx_j     <- dxtemp[, j]
      dx_block <- matrix(dx_j, ncol = 1)
      if (q_ij >= 2L) {
        tmp <- dx_j
        for (l in seq_len(q_ij - 1L)) {
          tmp      <- lag1(tmp)
          dx_block <- cbind(dx_block, tmp)
        }
      }
      dx_all <- cbind(dx_all, dx_block)
    }
  }
  if (!is.null(dx_all)) {
    ww    <- if (wflag == 0L) dx_all else cbind(ww, dx_all)
    wflag <- 1L
  }

  # ALL fixed regressors including intercept
  # (This is the key difference from FDGP1 which excludes the intercept)
  if (wflag == 0L) {
    ww <- ztemp
  } else {
    ww <- cbind(ww, ztemp)
  }
  wflag <- 1L

  # Sample adjustment: lag1() always loses row 1; additional lags lose more
  ltemp <- c(pp, qlag[mxi, ])
  mxlag <- max(ltemp)
  sini  <- max(2L, mxlag + 1L)

  y0 <- ytemp[sini:tis]
  y1 <- y1temp[sini:tis]
  dy <- dytemp[sini:tis]
  x0 <- xtemp[sini:tis, , drop = FALSE]
  ww <- ww[sini:tis, , drop = FALSE]

  list(y0 = y0, dy = dy, y1 = y1, x0 = x0, ww = ww, mxlag = mxlag)
}


# ==============================================================================
# MG: Mean Group Estimator
# ==============================================================================
#'
#' For each group i, runs unrestricted OLS on the ECM:
#'   Δy = φ_i y_{t-1} + β_i' x_t + γ_i' w_t + e_t
#'
#' Computes long-run parameters θ_i = -β_i / φ_i for each group,
#' then averages: θ̄ = (1/NN) Σ θ_i
#'
#' SE uses cross-sectional variance:
#'   SE(θ̄_j) = sqrt( Σ(θ_ij - θ̄_j)² / (NN(NN-1)) )
#'
#' @param Y     Stacked dependent variable
#' @param X     Stacked regressors, sum(Ti) x k
#' @param Z     Stacked fixed regressors (intercept last column)
#' @param n     Total number of groups
#' @param Ti    Time periods per group (length n)
#' @param k     Number of regressors in X
#' @param plag  Lag order for y (length n)
#' @param qlag  Lag orders for X (n x k matrix)
#' @param NN    Number of groups to use for MG (default n; set to n-1
#'              when a relative variable is present)
#'
#' @return List with theta, theta_se, theta_t, theta_all, phi, phi_se,
#'         phi_all, beta, beta_se, other, other_se, numsi, ursll
#'
MG <- function(Y, X, Z, n, Ti, k, plag, qlag, NN = n) {

  theta_all  <- matrix(0, k, NN)
  phi_all    <- numeric(NN)
  sigma2_all <- numeric(NN)
  numsi      <- integer(NN)

  # Store all short-run params for cross-sectional averaging
  beta_store  <- matrix(0, k, NN)
  other_store <- list()
  sumlike     <- 0

  ini <- 1L
  for (i in seq_len(NN)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    ytemp <- Y[ini:fini]
    xtemp <- as.matrix(X[ini:fini, , drop = FALSE])
    ztemp <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp <- DGP1(ytemp, xtemp, ztemp, plag, qlag, k, i, Tis)
    y0 <- dgp$y0;  dy <- dgp$dy;  y1 <- dgp$y1
    x0 <- dgp$x0;  ww <- dgp$ww

    pp <- plag[i]
    numsi[i] <- ncol(ww)
    Ti_eff <- length(dy)

    # Build regression: dep ~ reg
    if (pp == 0L) {
      reg <- cbind(x0, ww)
      dep <- y0
    } else {
      reg <- cbind(y1, x0, ww)
      dep <- dy
    }

    # OLS
    SRb  <- solve(crossprod(reg), crossprod(reg, dep))
    res  <- dep - reg %*% SRb
    sig2 <- as.numeric(crossprod(res) / Ti_eff)  # ML sigma^2 (divide by T)

    # Log-likelihood
    like <- -(Ti_eff / 2) * log(2 * pi * sig2) - Ti_eff / 2
    sumlike <- sumlike + like

    # Extract long-run parameters
    if (pp == 0L) {
      beta_i  <- as.numeric(SRb[1:k])
      theta_i <- beta_i
      phi_i   <- -1
      other_i <- as.numeric(SRb[(k + 1):length(SRb)])
    } else {
      phi_i   <- as.numeric(SRb[1])
      beta_i  <- as.numeric(SRb[2:(k + 1)])
      theta_i <- -beta_i / phi_i
      other_i <- as.numeric(SRb[(k + 2):length(SRb)])
    }

    theta_all[, i]  <- theta_i
    phi_all[i]      <- phi_i
    beta_store[, i] <- beta_i
    sigma2_all[i]   <- sig2
    other_store[[i]] <- other_i

    ini <- fini + 1L
  }

  # --- MG averages and cross-sectional SEs ---

  # Long-run theta
  theta_mg  <- rowMeans(theta_all)
  theta_dev <- theta_all - theta_mg
  theta_cov <- tcrossprod(theta_dev) / (NN * (NN - 1))
  theta_se  <- sqrt(diag(theta_cov))
  theta_t   <- theta_mg / theta_se

  # Error correction phi
  phi_mg  <- mean(phi_all)
  phi_dev <- phi_all - phi_mg
  phi_se  <- sqrt(sum(phi_dev^2) / (NN * (NN - 1)))
  phi_t   <- phi_mg / phi_se

  # Short-run beta (= coefficients on x in ECM)
  beta_mg  <- rowMeans(beta_store)
  beta_dev <- beta_store - beta_mg
  beta_se  <- sqrt(rowSums(beta_dev^2) / (NN * (NN - 1)))
  beta_t   <- beta_mg / beta_se

  # Other short-run params (lagged Δy, current/lagged Δx, intercept)
  n_other   <- length(other_store[[1]])
  other_mat <- matrix(0, n_other, NN)
  for (i in seq_len(NN)) other_mat[, i] <- other_store[[i]]
  other_mg  <- rowMeans(other_mat)
  other_dev <- other_mat - other_mg
  other_se  <- sqrt(rowSums(other_dev^2) / (NN * (NN - 1)))
  other_t   <- other_mg / other_se

  list(
    theta = theta_mg, theta_se = theta_se, theta_t = theta_t,
    theta_cov = theta_cov, theta_all = theta_all,
    phi = phi_mg, phi_se = phi_se, phi_t = phi_t, phi_all = phi_all,
    beta = beta_mg, beta_se = beta_se, beta_t = beta_t,
    other = other_mg, other_se = other_se, other_t = other_t,
    sigma2_all = sigma2_all, numsi = numsi, ursll = sumlike
  )
}


# ==============================================================================
# PMG: Pooled Mean Group Estimator (Back-Substitution Algorithm)
# ==============================================================================
#'
#' Constrains the long-run parameters θ to be the same across all groups,
#' while allowing φ_i, γ_i to differ.
#'
#' Iterative algorithm:
#'   Given θ^(j-1), for each group:
#'     1. Compute EC term: ξ = y_{t-1} - X θ^(j-1)
#'     2. Regress Δy on [ξ, W] → get φ̂_i, ŝ_i, σ̂²_i
#'     3. Update θ via concentrated ML normal equations
#'   Repeat until convergence.
#'
#' After convergence, SEs for θ come from the full information matrix.
#' SEs for φ, β, and other SR params use cross-sectional variance.
#'
#' @param Y         Stacked dependent variable
#' @param X         Stacked regressors, sum(Ti) x k
#' @param Z         Stacked fixed regressors (intercept last column)
#' @param n         Total number of groups
#' @param Ti        Time periods per group (length n)
#' @param k         Number of regressors in X
#' @param plag      Lag order for y (length n)
#' @param qlag      Lag orders for X (n x k matrix)
#' @param NN        Number of groups for PMG (default n)
#' @param ini_theta Initial values for θ (length k). If NULL, uses MG.
#' @param maxiter   Maximum iterations (default 1000)
#' @param tol       Convergence tolerance (default 1e-4)
#' @param verbose   Print iteration info (default FALSE)
#'
#' @return List with theta, theta_se, theta_t, phi, phi_se, phi_all,
#'         beta, beta_se, other, other_se, converged, iterations, loglik
#'
PMG <- function(Y, X, Z, n, Ti, k, plag, qlag, NN = n,
                ini_theta = NULL, maxiter = 1000L, tol = 1e-4,
                verbose = FALSE) {

  # --- Initial values for theta ---
  if (is.null(ini_theta)) {
    mg_res    <- MG(Y, X, Z, n, Ti, k, plag, qlag, NN)
    oldtheta  <- mg_res$theta
  } else {
    oldtheta <- ini_theta
  }

  # --- BSA iterations ---
  cflag   <- FALSE
  oldRSLL <- NULL
  niter   <- 0L

  for (iter in seq_len(maxiter)) {
    res   <- pmg_bsa_step(Y, X, Z, Ti, k, plag, qlag, NN, oldtheta)
    theta <- res$theta
    RSLL  <- res$sumlike

    if (iter == 1L) {
      oldRSLL  <- RSLL
      oldtheta <- theta
      next
    }

    convlike <- RSLL - oldRSLL
    convt    <- abs(theta - oldtheta)
    convtr   <- sum(convt)
    tol2     <- tol * k

    if (verbose) {
      cat(sprintf("Iter %d: RSLL=%.6f, dRSLL=%.6f, sum|dtheta|=%.8f\n",
                  iter, RSLL, convlike, convtr))
    }

    if (abs(convlike) < tol && convtr < tol2) {
      cflag <- TRUE
      niter <- iter
      break
    }

    oldRSLL  <- RSLL
    oldtheta <- theta
    niter    <- iter
  }

  if (!cflag) {
    warning("PMG did not converge after ", maxiter, " iterations")
    niter <- maxiter
  }

  # --- Final results with SEs ---
  final <- pmg_final(Y, X, Z, Ti, k, plag, qlag, NN, theta)
  final$converged  <- cflag
  final$iterations <- niter
  final$loglik     <- RSLL

  final
}


# --- BSA iteration step (internal) -------------------------------------------

pmg_bsa_step <- function(Y, X, Z, Ti, k, plag, qlag, NN, oldtheta) {

  sumden  <- matrix(0, k, k)
  sumnum  <- matrix(0, k, 1)
  sumlike <- 0

  ini <- 1L
  for (i in seq_len(NN)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    ytemp <- Y[ini:fini]
    xtemp <- as.matrix(X[ini:fini, , drop = FALSE])
    ztemp <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp <- DGP1(ytemp, xtemp, ztemp, plag, qlag, k, i, Tis)
    y0 <- dgp$y0;  dy <- dgp$dy;  y1 <- dgp$y1
    x0 <- dgp$x0;  ww <- dgp$ww

    pp     <- plag[i]
    Ti_eff <- length(dy)

    # --- Estimate phi and SR params given old theta ---
    if (pp == 0L) {
      faihat <- -1
      tmp0   <- y0 - x0 %*% oldtheta
      SRb    <- solve(crossprod(ww), crossprod(ww, tmp0))
      tmp2   <- tmp0 - ww %*% SRb
    } else {
      xi     <- y1 - x0 %*% oldtheta       # EC term
      reg    <- cbind(xi, ww)
      SRb    <- solve(crossprod(reg), crossprod(reg, dy))
      faihat <- as.numeric(SRb[1])
      tmp2   <- dy - reg %*% SRb
    }

    sig2 <- as.numeric(crossprod(tmp2) / Ti_eff)

    # Log-likelihood
    like    <- -(Ti_eff / 2) * log(2 * pi * sig2) - Ti_eff / 2
    sumlike <- sumlike + like

    # --- Concentrated ML normal equations for theta ---
    # H = I - ww(ww'ww)^{-1}ww'  (projection orthogonal to ww)
    ww_inv <- solve(crossprod(ww))
    tmp3   <- crossprod(x0) - crossprod(x0, ww) %*% ww_inv %*% crossprod(ww, x0)
    T1den  <- (faihat^2 / sig2) * tmp3

    tmp4 <- dy - faihat * y1
    tmp5 <- crossprod(x0, tmp4) -
            crossprod(x0, ww) %*% ww_inv %*% crossprod(ww, tmp4)
    T1num <- -(faihat / sig2) * tmp5

    sumden <- sumden + T1den
    sumnum <- sumnum + T1num

    ini <- fini + 1L
  }

  theta <- as.numeric(solve(sumden, sumnum))
  list(theta = theta, sumlike = sumlike)
}


# --- Final results after convergence (internal) ------------------------------
#'
#' Computes SEs from the full information matrix (matches PMLESR1 in SUB1.PRC).
#' The information matrix has blocks for:
#'   θ (k x k), φ_1..φ_NN (diagonal), s_1..s_NN (block-diagonal)
#'
pmg_final <- function(Y, X, Z, Ti, k, plag, qlag, NN, theta) {

  phi_all    <- numeric(NN)
  phi_fixed  <- logical(NN)   # TRUE when plag==0 forces phi=-1
  beta_all   <- matrix(0, k, NN)
  sigma2_all <- numeric(NN)
  numsi      <- integer(NN)
  sr_all     <- list()
  nfai0      <- 0L
  sumlike    <- 0

  # Information matrix blocks
  G11 <- matrix(0, k, k)       # θ block
  G13 <- matrix(0, k, NN)      # θ-φ cross
  G33 <- matrix(0, NN, NN)     # φ block (diagonal)

  G14_list <- vector("list", NN)
  G34_list <- vector("list", NN)
  G44_list <- vector("list", NN)

  ini <- 1L
  for (i in seq_len(NN)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    ytemp <- Y[ini:fini]
    xtemp <- as.matrix(X[ini:fini, , drop = FALSE])
    ztemp <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp <- DGP1(ytemp, xtemp, ztemp, plag, qlag, k, i, Tis)
    y0 <- dgp$y0;  dy <- dgp$dy;  y1 <- dgp$y1
    x0 <- dgp$x0;  ww <- dgp$ww

    pp     <- plag[i]
    Ti_eff <- length(dy)
    numsi[i] <- ncol(ww)

    # --- Estimate phi and SR params given converged theta ---
    if (pp == 0L) {
      faihat <- -1
      nfai0  <- nfai0 + 1L
      phi_fixed[i] <- TRUE
      tmp0   <- y0 - x0 %*% theta
      si     <- solve(crossprod(ww), crossprod(ww, tmp0))
      tmp2   <- tmp0 - ww %*% si
    } else {
      xi     <- y1 - x0 %*% theta
      reg    <- cbind(xi, ww)
      SRb    <- solve(crossprod(reg), crossprod(reg, dy))
      faihat <- as.numeric(SRb[1])
      si     <- SRb[-1, , drop = FALSE]
      tmp2   <- dy - reg %*% SRb
    }

    sig2 <- as.numeric(crossprod(tmp2) / Ti_eff)

    like    <- -(Ti_eff / 2) * log(2 * pi * sig2) - Ti_eff / 2
    sumlike <- sumlike + like

    phi_all[i]      <- faihat
    beta_all[, i]   <- -faihat * theta    # β_i = -φ_i θ
    sigma2_all[i]   <- sig2
    sr_all[[i]]     <- as.numeric(si)

    # EC term for information matrix
    ecm <- y1 - x0 %*% theta

    # θ block
    G11 <- G11 + (faihat^2 / sig2) * crossprod(x0)
    G13[, i] <- as.numeric(-(faihat / sig2) * crossprod(x0, ecm))

    # φ block
    if (pp == 0L) {
      G33[i, i] <- 0   # No variance when φ = -1
    } else {
      G33[i, i] <- as.numeric((1 / sig2) * crossprod(ecm))
    }

    # s_i blocks
    G14_list[[i]] <- -(faihat / sig2) * crossprod(x0, ww)
    G34_list[[i]] <- (1 / sig2) * crossprod(ecm, ww)
    G44_list[[i]] <- (1 / sig2) * crossprod(ww)

    ini <- fini + 1L
  }

  # --- Assemble full information matrix ---
  sumsi <- sum(numsi)

  G14 <- do.call(cbind, G14_list)   # k x sumsi

  G34 <- matrix(0, NN, sumsi)
  G44 <- matrix(0, sumsi, sumsi)
  col_s <- 1L
  for (i in seq_len(NN)) {
    col_e <- col_s + numsi[i] - 1L
    G34[i, col_s:col_e] <- G34_list[[i]]
    G44[col_s:col_e, col_s:col_e] <- G44_list[[i]]
    col_s <- col_e + 1L
  }

  g1 <- cbind(G11, G13, G14)
  g3 <- cbind(t(G13), G33, G34)
  g4 <- cbind(t(G14), t(G34), G44)
  GG <- rbind(g1, g3, g4)

  # Remove rows/columns for groups where φ = -1 (plag==0, no variance)
  if (nfai0 > 0L) {
    remove_idx <- k + which(phi_fixed)
    GG <- GG[-remove_idx, -remove_idx]
  }

  # Invert information matrix
  GGinv   <- solve(GG)
  PMLECOV <- GGinv[1:k, 1:k, drop = FALSE]

  # θ SE from information matrix
  theta_se <- sqrt(diag(PMLECOV))
  theta_t  <- theta / theta_se

  # --- Cross-sectional SEs for SR parameters (MG-style averaging) ---

  # φ
  phi_mg  <- mean(phi_all)
  phi_dev <- phi_all - phi_mg
  phi_se  <- sqrt(sum(phi_dev^2) / (NN * (NN - 1)))
  phi_t   <- phi_mg / phi_se

  # β_i = -φ_i θ (short-run coefficients on X)
  beta_mg  <- rowMeans(beta_all)
  beta_dev <- beta_all - beta_mg
  beta_se  <- sqrt(rowSums(beta_dev^2) / (NN * (NN - 1)))
  beta_t   <- beta_mg / beta_se

  # Other SR params (lagged Δy, current/lagged Δx, intercept)
  # Requires uniform lag specification across groups
  stopifnot(length(unique(numsi)) == 1L)
  n_other   <- numsi[1]
  other_mat <- matrix(0, n_other, NN)
  for (i in seq_len(NN)) other_mat[, i] <- sr_all[[i]]
  other_mg  <- rowMeans(other_mat)
  other_dev <- other_mat - other_mg
  other_se  <- sqrt(rowSums(other_dev^2) / (NN * (NN - 1)))
  other_t   <- other_mg / other_se

  list(
    theta = theta, theta_se = theta_se, theta_t = theta_t,
    theta_cov = PMLECOV,
    phi = phi_mg, phi_se = phi_se, phi_t = phi_t, phi_all = phi_all,
    beta = beta_mg, beta_se = beta_se, beta_t = beta_t,
    other = other_mg, other_se = other_se, other_t = other_t,
    sigma2_all = sigma2_all, numsi = numsi,
    rsll = sumlike
  )
}


# ==============================================================================
# Hausman Test: PMG vs MG
# ==============================================================================
#'
#' Tests H0: PMG and MG estimates are both consistent (PMG is efficient).
#' Under H0, the test statistic is Chi-squared with k degrees of freedom.
#'
#' h = (θ_MG - θ_PMG)' (V_MG - V_PMG)^{-1} (θ_MG - θ_PMG)
#'
hausman_test <- function(pmg_res, mg_res, k) {
  diff  <- mg_res$theta - pmg_res$theta
  V_diff <- mg_res$theta_cov - pmg_res$theta_cov

  # Individual tests
  h_individual <- diff^2 / diag(V_diff)
  p_individual <- 1 - pchisq(h_individual, df = 1)

  # Joint test
  h_joint <- as.numeric(t(diff) %*% solve(V_diff) %*% diff)
  p_joint <- 1 - pchisq(h_joint, df = k)

  list(
    h_individual = h_individual,
    p_individual = p_individual,
    h_joint = h_joint,
    p_joint = p_joint
  )
}


# ==============================================================================
# PMG_NR: Pooled Mean Group Estimator (Newton-Raphson Algorithm)
# ==============================================================================
#'
#' Same model as PMG() but uses Newton-Raphson optimization, which jointly
#' updates θ and φ_1,...,φ_NN using first and second derivatives of the
#' log-likelihood. This matches the GAUSS PMLENR1 procedure.
#'
#' @param Y,X,Z,n,Ti,k,plag,qlag,NN  Same as PMG()
#' @param ini_theta Initial θ (if NULL, uses MG)
#' @param maxiter,tol,verbose  Convergence parameters
#'
PMG_NR <- function(Y, X, Z, n, Ti, k, plag, qlag, NN = n,
                   ini_theta = NULL, maxiter = 1000L, tol = 1e-4,
                   verbose = FALSE) {

  # --- Initial theta ---
  if (is.null(ini_theta)) {
    mg_res   <- MG(Y, X, Z, n, Ti, k, plag, qlag, NN)
    oldtheta <- mg_res$theta
  } else {
    oldtheta <- ini_theta
  }

  # --- Initial phi given theta (matches INIPHI1 in SUB1.PRC) ---
  oldphi <- pmg_initial_phi(Y, X, Z, Ti, k, plag, qlag, NN, oldtheta)

  # --- NR iterations ---
  cflag   <- FALSE
  oldRSLL <- NULL
  niter   <- 0L

  for (iter in seq_len(maxiter)) {
    res <- pmg_nr_step(Y, X, Z, Ti, k, plag, qlag, NN, oldtheta, oldphi)
    theta <- res$theta
    phi   <- res$phi
    RSLL  <- res$sumlike

    if (iter == 1L) {
      oldRSLL  <- RSLL
      oldtheta <- theta
      oldphi   <- phi
      next
    }

    convlike <- RSLL - oldRSLL
    convt    <- abs(theta - oldtheta)
    convtr   <- sum(convt)
    tol2     <- tol * k

    if (verbose) {
      cat(sprintf("Iter %d: RSLL=%.6f, dRSLL=%.6f, sum|dtheta|=%.8f\n",
                  iter, RSLL, convlike, convtr))
    }

    if (abs(convlike) < tol && convtr < tol2) {
      cflag <- TRUE
      niter <- iter
      break
    }

    oldRSLL  <- RSLL
    oldtheta <- theta
    oldphi   <- phi
    niter    <- iter
  }

  if (!cflag) {
    warning("PMG_NR did not converge after ", maxiter, " iterations")
    niter <- maxiter
  }

  # --- Final results with SEs ---
  final <- pmg_final(Y, X, Z, Ti, k, plag, qlag, NN, theta)
  final$converged  <- cflag
  final$iterations <- niter
  final$loglik     <- RSLL

  final
}


# --- Initial phi estimates for NR (matches INIPHI1 in SUB1.PRC) ---

pmg_initial_phi <- function(Y, X, Z, Ti, k, plag, qlag, NN, theta) {

  phiini <- numeric(NN)

  ini <- 1L
  for (i in seq_len(NN)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    ytemp <- Y[ini:fini]
    xtemp <- as.matrix(X[ini:fini, , drop = FALSE])
    ztemp <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp <- DGP1(ytemp, xtemp, ztemp, plag, qlag, k, i, Tis)
    dy <- dgp$dy;  y1 <- dgp$y1;  x0 <- dgp$x0;  ww <- dgp$ww

    pp <- plag[i]

    # EC term given theta
    xi <- y1 - x0 %*% theta

    if (pp == 0L) {
      phiini[i] <- -1
    } else {
      # Concentrated OLS for phi: phi = (xi'H*dy) / (xi'H*xi)
      ww_inv <- solve(crossprod(ww))
      den <- as.numeric(crossprod(xi) -
               crossprod(xi, ww) %*% ww_inv %*% crossprod(ww, xi))
      num <- as.numeric(crossprod(xi, dy) -
               crossprod(xi, ww) %*% ww_inv %*% crossprod(ww, dy))
      phiini[i] <- num / den
    }

    ini <- fini + 1L
  }

  phiini
}


# --- NR iteration step (matches PMLENR1 in SUB1.PRC) ---

pmg_nr_step <- function(Y, X, Z, Ti, k, plag, qlag, NN, oldtheta, oldphi) {

  # Gradient and Hessian blocks
  sumden  <- matrix(0, k, k)     # Hessian: θ block
  sumnum  <- matrix(0, k, 1)     # Gradient: θ

  phifd   <- numeric(NN)         # Gradient: each φ_i
  phisd   <- matrix(0, NN, NN)   # Hessian: φ block (diagonal)

  thphisd <- matrix(0, k, NN)    # Hessian: θ-φ cross

  sumlike <- 0

  ini <- 1L
  for (i in seq_len(NN)) {
    faihat <- oldphi[i]
    Tis    <- Ti[i]
    fini   <- ini + Tis - 1L

    ytemp <- Y[ini:fini]
    xtemp <- as.matrix(X[ini:fini, , drop = FALSE])
    ztemp <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp <- DGP1(ytemp, xtemp, ztemp, plag, qlag, k, i, Tis)
    dy <- dgp$dy;  y1 <- dgp$y1;  x0 <- dgp$x0;  ww <- dgp$ww

    Ti_eff <- length(dy)

    # EC term and residual given old theta and old phi
    xi   <- y1 - x0 %*% oldtheta                        # EC term
    ndep <- dy - xi * faihat                             # dy - phi*ecm
    ww_inv <- solve(crossprod(ww))
    res  <- ndep - ww %*% ww_inv %*% crossprod(ww, ndep) # residual
    sig2 <- as.numeric(crossprod(res) / Ti_eff)

    like    <- -(Ti_eff / 2) * log(2 * pi * sig2) - Ti_eff / 2
    sumlike <- sumlike + like

    # --- Gradient ---
    # For theta:
    jnk1 <- crossprod(x0, res) -
             crossprod(x0, ww) %*% ww_inv %*% crossprod(ww, res)
    sumnum <- sumnum - (faihat / sig2) * jnk1

    # For phi_i:
    jnk2 <- as.numeric(crossprod(xi, res) -
              crossprod(xi, ww) %*% ww_inv %*% crossprod(ww, res))
    phifd[i] <- (1 / sig2) * jnk2

    # --- Hessian ---
    # theta-theta:
    jnk3 <- crossprod(x0) -
             crossprod(x0, ww) %*% ww_inv %*% crossprod(ww, x0)
    sumden <- sumden + (faihat^2 / sig2) * jnk3

    # phi_i-phi_i:
    jnk4 <- as.numeric(crossprod(xi) -
              crossprod(xi, ww) %*% ww_inv %*% crossprod(ww, xi))
    phisd[i, i] <- (1 / sig2) * jnk4

    # theta-phi_i:
    jnk5 <- crossprod(x0, xi) -
             crossprod(x0, ww) %*% ww_inv %*% crossprod(ww, xi)
    thphisd[, i] <- as.numeric(-(faihat / sig2) * jnk5)

    ini <- fini + 1L
  }

  # Assemble full gradient and Hessian
  TNUM  <- c(as.numeric(sumnum), phifd)              # (k+NN) x 1
  TDEN1 <- cbind(sumden, thphisd)                     # k x (k+NN)
  TDEN2 <- cbind(t(thphisd), phisd)                   # NN x (k+NN)
  TDEN  <- rbind(TDEN1, TDEN2)                        # (k+NN) x (k+NN)

  # NR update
  oldpara <- c(oldtheta, oldphi)
  para    <- oldpara + as.numeric(solve(TDEN, TNUM))
  theta   <- para[1:k]
  phi     <- para[(k + 1):length(para)]

  list(theta = theta, phi = phi, sumlike = sumlike)
}
