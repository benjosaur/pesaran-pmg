# ==============================================================================
# Panel Error Correction Model: Static and Dynamic Fixed Effects Estimators
#
# Based on: Pesaran, Shin & Smith (1999)
# "Pooled Mean Group Estimation of Dynamic Heterogeneous Panels"
# Journal of the American Statistical Association, 94(446), 621-634.
#
# Translated from the original GAUSS implementation to R.
#
# Provides:
#   SFE()   - Static Fixed Effects estimator of long-run parameters
#   DFE1()  - Dynamic Fixed Effects estimator (single regressor set)
#   DFE2()  - Dynamic Fixed Effects estimator (two regressor sets)
#   FDGP1() - Data generation helper for DFE1
#   FDGP2() - Data generation helper for DFE2
# ==============================================================================


# --- Helpers ------------------------------------------------------------------

#' Lag a vector (or each column of a matrix) by one period within a group
lag1 <- function(x) {
  if (is.matrix(x)) {
    rbind(NA, x[-nrow(x), , drop = FALSE])
  } else {
    c(NA, x[-length(x)])
  }
}

#' Within-group (FE) demeaning: z_tilde = z - bar{z}
#' This is the "within transformation" Q = I - iota*(iota'iota)^{-1}*iota'
#' applied to z, which projects out the group mean (= the fixed effect).
demean <- function(z) {
  if (is.matrix(z)) {
    sweep(z, 2, colMeans(z))
  } else {
    z - mean(z)
  }
}

#' Basic regression diagnostics (R^2, log-likelihood, information criteria)
#' Stands in for the GAUSS DIAGNOS() procedure (not provided in source).
diagnos <- function(dep, reg, res, yhat, nreg) {
  n_obs  <- length(dep)
  k      <- nreg
  ss_res <- sum(res^2)
  ss_tot <- sum((dep - mean(dep))^2)

  rsq     <- 1 - ss_res / ss_tot
  adj_rsq <- 1 - (1 - rsq) * (n_obs - 1) / (n_obs - k)
  sig2    <- ss_res / (n_obs - k)
  ll      <- -n_obs / 2 * (log(2 * pi) + log(sig2) + 1)
  aic     <- -2 * ll + 2 * k
  bic     <- -2 * ll + k * log(n_obs)
  hq      <- -2 * ll + 2 * k * log(log(n_obs))

  list(rsq = rsq, adj_rsq = adj_rsq, sigma2 = sig2,
       loglik = ll, aic = aic, bic = bic, hq = hq)
}


# ==============================================================================
# SFE: Static Fixed Effects Estimation of the Long-Run Parameters
# ==============================================================================
#'
#' Estimates the cointegrating (long-run) relationship directly from levels:
#'
#'   y_{it} = theta' x_{it} + alpha_i + u_{it}             [levels]
#'
#' Equivalently, this imposes phi = -1 in the ECM:
#'
#'   dy_{it} = -(y_{i,t-1} - theta' x_{it}) + alpha_i + e_{it}
#'
#' The within (FE) transformation removes alpha_i by demeaning per group.
#' The pooled within estimator is:
#'
#'   theta_hat = [sum_i X_i' Q_i X_i]^{-1} [sum_i X_i' Q_i y_i]
#'
#' @param Y   Stacked dependent variable, length sum(Ti)
#' @param X   Stacked regressors in the cointegrating relation, sum(Ti) x k
#' @param Z   Stacked fixed regressors; intercept must be the LAST column
#' @param n   Number of cross-sectional units (groups)
#' @param Ti  Integer vector of length n, time periods per group
#' @param k   Number of regressors in X entering the cointegrating vector
#'
#' @return List with sferes, sfeind, rsferes, sfedres, qobs, ftheta
#'
SFE <- function(Y, X, Z, n, Ti, k) {

  ntmp   <- k + ncol(Z) - 1L         # x-regressors + fixed regressors excl. intercept
  sumden <- matrix(0, ntmp, ntmp)
  sumnum <- matrix(0, ntmp, 1)
  meanY  <- numeric(n)
  meanX  <- matrix(0, n, ntmp)

  # Collectors for pooled demeaned data
  QY_list  <- vector("list", n)
  QX_list  <- vector("list", n)
  QY1_list <- vector("list", n)
  QDY_list <- vector("list", n)

  ini <- 1L

  for (i in seq_len(n)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    # --- Extract group-i data ---
    y0 <- Y[ini:fini]
    y1 <- lag1(y0)                             # y_{t-1}
    dy <- y0 - y1                              # Delta y
    x0 <- as.matrix(X[ini:fini, , drop = FALSE])
    z0 <- as.matrix(Z[ini:fini, , drop = FALSE])

    # Append any fixed regressors (excluding the intercept in last column)
    if (ncol(z0) > 1L) {
      x0 <- cbind(x0, z0[, 1:(ncol(z0) - 1L), drop = FALSE])
    }

    # Drop the first obs (lost to the lag)
    idx <- 2:Tis
    y0 <- y0[idx];  y1 <- y1[idx];  dy <- dy[idx]
    x0 <- x0[idx, , drop = FALSE]

    # --- Within-group demeaning (FE transformation) ---
    QY_list[[i]]  <- demean(y0)
    QX_list[[i]]  <- demean(x0)
    QY1_list[[i]] <- demean(y1)
    QDY_list[[i]] <- demean(dy)

    # Group means for recovering alpha_i later
    meanY[i]   <- mean(y0)
    meanX[i, ] <- colMeans(x0)

    # Accumulate normal equations: sum_i X_i'Q_i X_i and sum_i X_i'Q_i y_i
    Ti_eff <- length(y0)
    iota   <- matrix(1, Ti_eff, 1)
    PX     <- iota %*% solve(crossprod(iota), t(iota))   # projection onto constant
    Fden   <- crossprod(x0) - crossprod(x0, PX) %*% x0
    Fnum   <- crossprod(x0, y0) - crossprod(x0, PX) %*% y0

    sumden <- sumden + Fden
    sumnum <- sumnum + Fnum
    ini    <- fini + 1L
  }

  # Pool demeaned data
  QY  <- do.call(c, QY_list)
  QX  <- do.call(rbind, QX_list)
  QY1 <- do.call(c, QY1_list)
  QDY <- do.call(c, QDY_list)

  # ---- FE estimator of theta ----
  Ftheta <- solve(sumden, sumnum)

  Qres  <- QY - QX %*% Ftheta
  nreg  <- nrow(Ftheta) + n                    # regressors + n fixed effects
  Qobs  <- length(QY)
  Qsig2 <- as.numeric(crossprod(Qres) / (Qobs - nreg))

  # Conventional variance / SE / t
  Qtemp <- solve(crossprod(QX))
  Qvar  <- Qsig2 * Qtemp
  Qse   <- sqrt(diag(Qvar))
  Qt    <- as.numeric(Ftheta) / Qse
  sferes <- cbind(estimate = as.numeric(Ftheta), se = Qse, t = Qt)

  # ---- ECM-form diagnostics ----
  epara   <- c(-1, as.numeric(Ftheta))
  reg     <- cbind(QY1, QX)
  yhat    <- reg %*% epara
  sfedres <- diagnos(QDY, reg, as.numeric(QDY - yhat), yhat, nreg)

  # ---- Cluster-robust (White) standard errors ----
  meat <- matrix(0, ntmp, ntmp)
  tini <- 1L
  for (i in seq_len(n)) {
    diffT <- Ti[i] - 1L
    tfini <- tini + diffT - 1L
    Qu_i  <- QY[tini:tfini] - QX[tini:tfini, , drop = FALSE] %*% Ftheta
    Xi    <- QX[tini:tfini, , drop = FALSE]
    meat  <- meat + crossprod(Xi, Qu_i) %*% crossprod(Qu_i, Xi)
    tini  <- tfini + 1L
  }
  RQvar   <- Qtemp %*% meat %*% Qtemp
  RQse    <- sqrt(diag(RQvar))
  RQt     <- as.numeric(Ftheta) / RQse
  rsferes <- cbind(estimate = as.numeric(Ftheta), robust_se = RQse, robust_t = RQt)

  # ---- Individual fixed effects: alpha_i = mean(y_i) - mean(x_i)' theta ----
  Qind   <- numeric(n)
  Qindse <- numeric(n)
  for (i in seq_len(n)) {
    Qind[i]   <- meanY[i] - sum(meanX[i, ] * as.numeric(Ftheta))
    Qindse[i] <- sqrt(Qsig2 / Ti[i])
  }
  Qindt  <- Qind / Qindse
  sfeind <- cbind(alpha = Qind, se = Qindse, t = Qindt)

  ftheta_lr <- as.numeric(Ftheta[1:k])

  list(sferes  = sferes,
       sfeind  = sfeind,
       rsferes = rsferes,
       sfedres = sfedres,
       qobs    = Qobs,
       ftheta  = ftheta_lr)
}


# ==============================================================================
# FDGP1: Data Generation for DFE1
# ==============================================================================
#'
#' Constructs the variables needed for the dynamic FE regression for group i:
#'   y0 (current y), dy (Delta y), y1 (lagged y), x0 (current x),
#'   ww (short-run dynamics: lagged Delta y, current/lagged Delta x, fixed regs)
#'
FDGP1 <- function(ytemp, xtemp, ztemp, plag, qlag, k, tis, mxi) {

  pp <- plag[1]

  y1temp <- lag1(ytemp)
  dytemp <- ytemp - y1temp
  x1temp <- lag1(xtemp)
  dxtemp <- xtemp - x1temp

  ww    <- NULL
  wflag <- 0L

  # --- Lagged Delta-y terms (for ARDL dynamics) ---
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

  # --- Current and lagged Delta-x terms ---
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

  # --- Fixed regressors (excluding intercept in last column) ---
  if (ncol(ztemp) > 1L) {
    z_extra <- ztemp[, 1:(ncol(ztemp) - 1L), drop = FALSE]
    ww      <- if (wflag == 0L) z_extra else cbind(ww, z_extra)
    wflag   <- 1L
  }

  # --- Trim to valid sample (drop rows lost to lagging) ---
  ltemp <- c(pp, qlag[mxi, ])
  mxlag <- max(ltemp)
  sini  <- if (mxlag == 0L) 2L else mxlag + 1L

  y0 <- ytemp[sini:tis]
  y1 <- y1temp[sini:tis]
  dy <- dytemp[sini:tis]
  x0 <- xtemp[sini:tis, , drop = FALSE]
  if (wflag == 0L) ww <- NULL else ww <- ww[sini:tis, , drop = FALSE]

  list(y0 = y0, dy = dy, y1 = y1, x0 = x0,
       ww = ww, wflag = wflag, mxlag = mxlag)
}


# ==============================================================================
# DFE1: Dynamic Fixed Effects Estimation (single regressor set)
# ==============================================================================
#'
#' Estimates the panel ARDL / ECM:
#'
#'   dy_{it} = phi * y_{i,t-1} + beta' x_{it} + alpha_i
#'             + [lagged dy, current/lagged dx] * gamma + e_{it}
#'
#' The within transformation removes alpha_i.  After estimating the
#' short-run parameters (phi, beta, gamma) by pooled within-OLS,
#' the long-run (cointegrating) parameters are recovered:
#'
#'   theta = -beta / phi
#'
#' with standard errors via the delta method.
#'
#' @param Y     Stacked dependent variable
#' @param X     Stacked regressors in the LR relation, sum(Ti) x k
#' @param Z     Stacked fixed regressors (intercept last column)
#' @param n     Number of groups
#' @param Ti    Time periods per group (length n)
#' @param k     Number of x-regressors in the cointegrating relation
#' @param numot Number of "other" short-run terms (from W)
#' @param plag  Lag order for Delta y (scalar or vector)
#' @param qlag  Lag orders for Delta x (n x k matrix)
#'
#' @return List with dferes, dfeind, rdferes, dfedres, qobs, qtheta
#'
DFE1 <- function(Y, X, Z, n, Ti, k, numot, plag, qlag) {

  pp <- plag[1]

  dimind <- if (pp == 0L) k + numot else k + 1L + numot

  sumden  <- matrix(0, dimind, dimind)
  sumnum  <- matrix(0, dimind, 1)
  meandep <- numeric(n)
  meanind <- matrix(0, n, dimind)

  QDY_list  <- vector("list", n)
  QY0_list  <- vector("list", n)
  QY1_list  <- vector("list", n)
  QX0_list  <- vector("list", n)
  QWW_list  <- vector("list", n)
  QREG_list <- vector("list", n)

  global_wflag <- 0L
  global_mxlag <- 0L

  ini <- 1L

  for (i in seq_len(n)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    ytemp <- Y[ini:fini]
    xtemp <- as.matrix(X[ini:fini, , drop = FALSE])
    ztemp <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp   <- FDGP1(ytemp, xtemp, ztemp, plag, qlag, k, Tis, i)
    y0 <- dgp$y0;  dy <- dgp$dy;  y1 <- dgp$y1
    x0 <- dgp$x0;  ww <- dgp$ww
    wflag <- dgp$wflag;  mxlag <- dgp$mxlag

    global_wflag <- max(global_wflag, wflag)
    global_mxlag <- max(global_mxlag, mxlag)

    # Demeaning
    QY0_list[[i]]  <- demean(y0)
    QY1_list[[i]]  <- demean(y1)
    QDY_list[[i]]  <- demean(dy)
    QX0_list[[i]]  <- demean(x0)
    if (wflag == 1L) QWW_list[[i]] <- demean(ww)

    # Build full regressor matrix
    if (pp == 0L) {
      Freg <- if (wflag == 0L) x0 else cbind(x0, ww)
    } else {
      Freg <- if (wflag == 0L) cbind(y1, x0) else cbind(y1, x0, ww)
    }

    QREG_list[[i]] <- demean(Freg)

    # Accumulate normal equations
    Ti_eff  <- length(dy)
    iota    <- matrix(1, Ti_eff, 1)
    PX      <- iota %*% solve(crossprod(iota), t(iota))
    Fden    <- crossprod(Freg) - crossprod(Freg, PX) %*% Freg
    dep_var <- if (pp == 0L) y0 else dy
    Fnum    <- crossprod(Freg, dep_var) - crossprod(Freg, PX) %*% dep_var

    sumden <- sumden + Fden
    sumnum <- sumnum + Fnum

    meandep[i]  <- mean(dy)
    meanind[i,] <- colMeans(Freg)

    ini <- fini + 1L
  }

  # Pool demeaned data
  QDY  <- do.call(c, QDY_list)
  QY0  <- do.call(c, QY0_list)
  QY1  <- do.call(c, QY1_list)
  QX0  <- do.call(rbind, QX0_list)
  Qreg <- do.call(rbind, QREG_list)
  QWW  <- if (global_wflag == 1L) do.call(rbind, QWW_list) else NULL

  # ---- FE estimator of all short-run parameters ----
  Fpara <- solve(sumden, sumnum)

  dep_pooled <- if (pp == 0L) QY0 else QDY
  Qres  <- dep_pooled - Qreg %*% Fpara
  Qobs  <- length(QDY)
  nreg  <- nrow(Fpara) + n
  Qsig2 <- as.numeric(crossprod(Qres) / (Qobs - nreg))

  Qtemp <- solve(crossprod(Qreg))
  Qvar  <- Qsig2 * Qtemp
  Qse   <- sqrt(diag(Qvar))
  Qt    <- as.numeric(Fpara) / Qse
  dferes <- cbind(estimate = as.numeric(Fpara), se = Qse, t = Qt)

  # ---- Cluster-robust (White) SE ----
  meat <- matrix(0, dimind, dimind)
  tini <- 1L
  for (i in seq_len(n)) {
    diffT <- if (global_mxlag == 0L) Ti[i] - 1L else Ti[i] - global_mxlag
    tfini <- tini + diffT - 1L
    deptmp <- if (global_mxlag == 0L) QY0[tini:tfini] else QDY[tini:tfini]
    indtmp <- Qreg[tini:tfini, , drop = FALSE]
    Qu     <- deptmp - indtmp %*% Fpara
    meat   <- meat + crossprod(indtmp, Qu) %*% crossprod(Qu, indtmp)
    tini   <- tfini + 1L
  }
  RQvar   <- Qtemp %*% meat %*% Qtemp
  RQse    <- sqrt(diag(RQvar))
  RQt     <- as.numeric(Fpara) / RQse
  rdferes <- cbind(estimate = as.numeric(Fpara), robust_se = RQse, robust_t = RQt)

  # ---- Individual fixed effects ----
  Qind   <- numeric(n)
  Qindse <- numeric(n)
  for (i in seq_len(n)) {
    Qind[i]   <- meandep[i] - sum(meanind[i, ] * as.numeric(Fpara))
    Qindse[i] <- sqrt(Qsig2 / Ti[i])
  }
  Qindt  <- Qind / Qindse
  dfeind <- cbind(alpha = Qind, se = Qindse, t = Qindt)

  # ---- Long-run parameters ----
  if (pp == 0L) {
    # Static case: beta IS the LR parameter
    Qtheta   <- as.numeric(Fpara[1:k])
    setheta  <- Qse[1:k]
    ttheta   <- Qtheta / setheta
    thres    <- cbind(theta = Qtheta, se = setheta, t = ttheta)

    Rsetheta <- RQse[1:k]
    Rttheta  <- Qtheta / Rsetheta
    rthres   <- cbind(theta = Qtheta, robust_se = Rsetheta, robust_t = Rttheta)

  } else {
    # Dynamic case: theta = -beta / phi  (delta method for SE)
    phi    <- as.numeric(Fpara[1])
    beta   <- as.numeric(Fpara[2:(k + 1)])
    Qtheta <- -beta / phi

    # Jacobian of theta w.r.t. (phi, beta):
    #   d(theta_j)/d(phi)    = beta_j / phi^2
    #   d(theta_j)/d(beta_m) = -1/phi  if j==m, 0 otherwise
    dtheta1 <- beta / phi^2
    dtheta2 <- diag(-1 / phi, nrow = k)
    dtheta  <- cbind(dtheta1, dtheta2)              # k x (k+1)

    vtheta   <- dtheta %*% Qvar[1:(k+1), 1:(k+1)] %*% t(dtheta)
    setheta  <- sqrt(diag(vtheta))
    ttheta   <- Qtheta / setheta
    thres    <- cbind(theta = Qtheta, se = setheta, t = ttheta)

    Rvtheta  <- dtheta %*% RQvar[1:(k+1), 1:(k+1)] %*% t(dtheta)
    Rsetheta <- sqrt(diag(Rvtheta))
    Rttheta  <- Qtheta / Rsetheta
    rthres   <- cbind(theta = Qtheta, robust_se = Rsetheta, robust_t = Rttheta)
  }

  # ---- ECM-form diagnostics ----
  reg_diag   <- if (is.null(QWW)) cbind(QY1, QX0) else cbind(QY1, QX0, QWW)
  dep_diag   <- QDY
  Qpara_diag <- if (pp == 0L) c(-1, as.numeric(Fpara)) else as.numeric(Fpara)
  yhat       <- reg_diag %*% Qpara_diag
  dfedres    <- diagnos(dep_diag, reg_diag, as.numeric(dep_diag - yhat), yhat, nreg)

  # Stack LR results on top, then all SR parameters
  dferes_full  <- rbind(thres, dferes)
  rdferes_full <- rbind(rthres, rdferes)

  list(dferes  = dferes_full,
       dfeind  = dfeind,
       rdferes = rdferes_full,
       dfedres = dfedres,
       qobs    = Qobs,
       qtheta  = Qtheta)
}


# ==============================================================================
# FDGP2: Data Generation for DFE2 (two regressor sets)
# ==============================================================================
FDGP2 <- function(ytemp, x1temp, x2temp, ztemp, plag, q1lag, q2lag,
                  k1, k2, mxi, tis) {

  y1temp  <- lag1(ytemp)
  dytemp  <- ytemp - y1temp
  x11temp <- lag1(x1temp)
  dx1temp <- x1temp - x11temp
  x21temp <- lag1(x2temp)
  dx2temp <- x2temp - x21temp

  ww    <- NULL
  wflag <- 0L

  # Lagged Delta-y terms
  pp <- plag[mxi]
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

  # Current and lagged Delta-x1 terms
  dx1_all <- NULL
  for (j in seq_len(k1)) {
    q_ij <- q1lag[mxi, j]
    if (q_ij >= 1L) {
      dx_j     <- dx1temp[, j]
      dx_block <- matrix(dx_j, ncol = 1)
      if (q_ij >= 2L) {
        tmp <- dx_j
        for (l in seq_len(q_ij - 1L)) {
          tmp      <- lag1(tmp)
          dx_block <- cbind(dx_block, tmp)
        }
      }
      dx1_all <- cbind(dx1_all, dx_block)
    }
  }
  if (!is.null(dx1_all)) {
    ww    <- if (wflag == 0L) dx1_all else cbind(ww, dx1_all)
    wflag <- 1L
  }

  # Current and lagged Delta-x2 terms
  dx2_all <- NULL
  for (j in seq_len(k2)) {
    q_ij <- q2lag[mxi, j]
    if (q_ij >= 1L) {
      dx_j     <- dx2temp[, j]
      dx_block <- matrix(dx_j, ncol = 1)
      if (q_ij >= 2L) {
        tmp <- dx_j
        for (l in seq_len(q_ij - 1L)) {
          tmp      <- lag1(tmp)
          dx_block <- cbind(dx_block, tmp)
        }
      }
      dx2_all <- cbind(dx2_all, dx_block)
    }
  }
  if (!is.null(dx2_all)) {
    ww    <- if (wflag == 0L) dx2_all else cbind(ww, dx2_all)
    wflag <- 1L
  }

  # Fixed regressors (excluding intercept)
  if (ncol(ztemp) > 1L) {
    z_extra <- ztemp[, 1:(ncol(ztemp) - 1L), drop = FALSE]
    ww      <- if (wflag == 0L) z_extra else cbind(ww, z_extra)
    wflag   <- 1L
  }

  # Sample adjustment
  ltemp <- c(plag[mxi], q1lag[mxi, ], q2lag[mxi, ])
  mxlag <- max(ltemp)
  sini  <- if (mxlag == 0L) 2L else mxlag + 1L

  y0  <- ytemp[sini:tis]
  y1  <- y1temp[sini:tis]
  dy  <- dytemp[sini:tis]
  x10 <- x1temp[sini:tis, , drop = FALSE]
  x20 <- x2temp[sini:tis, , drop = FALSE]
  if (wflag == 0L) ww <- NULL else ww <- ww[sini:tis, , drop = FALSE]

  list(y0 = y0, dy = dy, y1 = y1, x10 = x10, x20 = x20,
       ww = ww, wflag = wflag, mxlag = mxlag)
}


# ==============================================================================
# DFE2: Dynamic Fixed Effects Estimation (two regressor sets)
# ==============================================================================
#'
#' Same as DFE1 but with two separate sets of cointegrating regressors x1, x2:
#'
#'   dy_{it} = phi * y_{i,t-1} + beta1' x1_{it} + beta2' x2_{it}
#'             + alpha_i + W * gamma + e_{it}
#'
#' Long-run: theta1 = -beta1/phi,  theta2 = -beta2/phi
#'
DFE2 <- function(Y, X1, X2, Z, n, Ti, k1, k2, numot, plag, q1lag, q2lag) {

  pp <- plag[1]

  dimind <- if (pp == 0L) k1 + k2 + numot else k1 + k2 + 1L + numot

  sumden  <- matrix(0, dimind, dimind)
  sumnum  <- matrix(0, dimind, 1)
  meandep <- numeric(n)
  meanind <- matrix(0, n, dimind)

  QDY_list  <- vector("list", n)
  QY0_list  <- vector("list", n)
  QY1_list  <- vector("list", n)
  QX10_list <- vector("list", n)
  QX20_list <- vector("list", n)
  QWW_list  <- vector("list", n)
  QREG_list <- vector("list", n)

  global_wflag <- 0L
  global_mxlag <- 0L

  ini <- 1L

  for (i in seq_len(n)) {
    Tis  <- Ti[i]
    fini <- ini + Tis - 1L

    ytemp  <- Y[ini:fini]
    x1temp <- as.matrix(X1[ini:fini, , drop = FALSE])
    x2temp <- as.matrix(X2[ini:fini, , drop = FALSE])
    ztemp  <- as.matrix(Z[ini:fini, , drop = FALSE])

    dgp <- FDGP2(ytemp, x1temp, x2temp, ztemp, plag, q1lag, q2lag,
                 k1, k2, i, Tis)
    y0 <- dgp$y0;  dy <- dgp$dy;  y1 <- dgp$y1
    x10 <- dgp$x10;  x20 <- dgp$x20
    ww <- dgp$ww;  wflag <- dgp$wflag;  mxlag <- dgp$mxlag

    global_wflag <- max(global_wflag, wflag)
    global_mxlag <- max(global_mxlag, mxlag)

    QY0_list[[i]]  <- demean(y0)
    QY1_list[[i]]  <- demean(y1)
    QDY_list[[i]]  <- demean(dy)
    QX10_list[[i]] <- demean(x10)
    QX20_list[[i]] <- demean(x20)
    if (wflag == 1L) QWW_list[[i]] <- demean(ww)

    if (pp == 0L) {
      Freg <- if (wflag == 0L) cbind(x10, x20) else cbind(x10, x20, ww)
    } else {
      Freg <- if (wflag == 0L) cbind(y1, x10, x20) else cbind(y1, x10, x20, ww)
    }

    QREG_list[[i]] <- demean(Freg)

    Ti_eff  <- length(dy)
    iota    <- matrix(1, Ti_eff, 1)
    PX      <- iota %*% solve(crossprod(iota), t(iota))
    Fden    <- crossprod(Freg) - crossprod(Freg, PX) %*% Freg
    dep_var <- if (pp == 0L) y0 else dy
    Fnum    <- crossprod(Freg, dep_var) - crossprod(Freg, PX) %*% dep_var

    sumden <- sumden + Fden
    sumnum <- sumnum + Fnum

    meandep[i]  <- mean(dy)
    meanind[i,] <- colMeans(Freg)

    ini <- fini + 1L
  }

  QDY  <- do.call(c, QDY_list)
  QY0  <- do.call(c, QY0_list)
  QY1  <- do.call(c, QY1_list)
  QX10 <- do.call(rbind, QX10_list)
  QX20 <- do.call(rbind, QX20_list)
  Qreg <- do.call(rbind, QREG_list)
  QWW  <- if (global_wflag == 1L) do.call(rbind, QWW_list) else NULL

  Fpara <- solve(sumden, sumnum)

  dep_pooled <- if (pp == 0L) QY0 else QDY
  Qres  <- dep_pooled - Qreg %*% Fpara
  Qobs  <- length(QDY)
  nreg  <- nrow(Fpara) + n
  Qsig2 <- as.numeric(crossprod(Qres) / (Qobs - nreg))

  Qtemp <- solve(crossprod(Qreg))
  Qvar  <- Qsig2 * Qtemp
  Qse   <- sqrt(diag(Qvar))
  Qt    <- as.numeric(Fpara) / Qse
  dferes <- cbind(estimate = as.numeric(Fpara), se = Qse, t = Qt)

  # Robust SE
  meat <- matrix(0, dimind, dimind)
  tini <- 1L
  for (i in seq_len(n)) {
    diffT <- if (global_mxlag == 0L) Ti[i] - 1L else Ti[i] - global_mxlag
    tfini <- tini + diffT - 1L
    deptmp <- if (global_mxlag == 0L) QY0[tini:tfini] else QDY[tini:tfini]
    indtmp <- Qreg[tini:tfini, , drop = FALSE]
    Qu     <- deptmp - indtmp %*% Fpara
    meat   <- meat + crossprod(indtmp, Qu) %*% crossprod(Qu, indtmp)
    tini   <- tfini + 1L
  }
  RQvar   <- Qtemp %*% meat %*% Qtemp
  RQse    <- sqrt(diag(RQvar))
  RQt     <- as.numeric(Fpara) / RQse
  rdferes <- cbind(estimate = as.numeric(Fpara), robust_se = RQse, robust_t = RQt)

  # Fixed effects
  Qind   <- numeric(n)
  Qindse <- numeric(n)
  for (i in seq_len(n)) {
    Qind[i]   <- meandep[i] - sum(meanind[i, ] * as.numeric(Fpara))
    Qindse[i] <- sqrt(Qsig2 / Ti[i])
  }
  Qindt  <- Qind / Qindse
  dfeind <- cbind(alpha = Qind, se = Qindse, t = Qindt)

  # Long-run parameters
  ktot <- k1 + k2
  if (pp == 0L) {
    beta1    <- as.numeric(Fpara[1:k1])
    beta2    <- as.numeric(Fpara[(k1 + 1):ktot])
    Qtheta1  <- beta1
    Qtheta2  <- beta2
    Qtheta   <- c(Qtheta1, Qtheta2)
    setheta  <- Qse[1:ktot]
    ttheta   <- Qtheta / setheta
    thres    <- cbind(theta = Qtheta, se = setheta, t = ttheta)

    Rsetheta <- RQse[1:ktot]
    Rttheta  <- Qtheta / Rsetheta
    rthres   <- cbind(theta = Qtheta, robust_se = Rsetheta, robust_t = Rttheta)

  } else {
    phi     <- as.numeric(Fpara[1])
    beta1   <- as.numeric(Fpara[2:(k1 + 1)])
    beta2   <- as.numeric(Fpara[(k1 + 2):(ktot + 1)])
    Qtheta1 <- -beta1 / phi
    Qtheta2 <- -beta2 / phi
    beta    <- c(beta1, beta2)
    Qtheta  <- c(Qtheta1, Qtheta2)

    dtheta1 <- beta / phi^2
    dtheta2 <- diag(-1 / phi, nrow = ktot)
    dtheta  <- cbind(dtheta1, dtheta2)

    vtheta   <- dtheta %*% Qvar[1:(ktot + 1), 1:(ktot + 1)] %*% t(dtheta)
    setheta  <- sqrt(diag(vtheta))
    ttheta   <- Qtheta / setheta
    thres    <- cbind(theta = Qtheta, se = setheta, t = ttheta)

    Rvtheta  <- dtheta %*% RQvar[1:(ktot + 1), 1:(ktot + 1)] %*% t(dtheta)
    Rsetheta <- sqrt(diag(Rvtheta))
    Rttheta  <- Qtheta / Rsetheta
    rthres   <- cbind(theta = Qtheta, robust_se = Rsetheta, robust_t = Rttheta)
  }

  # Diagnostics
  reg_diag   <- if (is.null(QWW)) cbind(QY1, QX10, QX20) else cbind(QY1, QX10, QX20, QWW)
  dep_diag   <- QDY
  Qpara_diag <- if (pp == 0L) c(-1, as.numeric(Fpara)) else as.numeric(Fpara)
  yhat       <- reg_diag %*% Qpara_diag
  dfedres    <- diagnos(dep_diag, reg_diag, as.numeric(dep_diag - yhat), yhat, nreg)

  dferes_full  <- rbind(thres, dferes)
  rdferes_full <- rbind(rthres, rdferes)

  list(dferes  = dferes_full,
       dfeind  = dfeind,
       rdferes = rdferes_full,
       dfedres = dfedres,
       qobs    = Qobs,
       qtheta1 = Qtheta1,
       qtheta2 = Qtheta2)
}
