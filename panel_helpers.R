# ==============================================================================
# Panel Data Helpers
#
# Convenience layer between standard R panel data.frames and the
# stacked-vector API used by the estimators in panel_ecm_fe.R / pmg_mg.R.
#
# The estimators require:
#   Y  - numeric vector, stacked by group
#   X  - numeric matrix, stacked by group (sum(Ti) x k)
#   Z  - numeric matrix, fixed regressors (intercept last column)
#   Ti - integer vector, number of time periods per group
#   n  - number of groups
#
# This file provides:
#   prepare_panel()   - convert a data.frame to the stacked format
#   pmg_panel()       - run PMG directly from a data.frame
#   mg_panel()        - run MG directly from a data.frame
#   dfe_panel()       - run DFE directly from a data.frame
#   sfe_panel()       - run SFE directly from a data.frame
# ==============================================================================

source("panel_ecm_fe.R")
source("pmg_mg.R")


# ==============================================================================
# prepare_panel: Convert a data.frame to stacked estimation format
# ==============================================================================
#'
#' @param data      A data.frame (also accepts data.table or pdata.frame).
#' @param y         Name of the dependent variable (string).
#' @param x         Names of the X regressors (character vector).
#' @param unit      Name of the unit/group column (string).
#'                  If NULL, tries to extract from pdata.frame index.
#' @param time      Name of the time column (string).
#'                  If NULL, tries to extract from pdata.frame index.
#' @param z         Names of additional fixed regressors to include in Z
#'                  (character vector, default NULL = intercept only).
#'                  An intercept is always appended as the last column of Z.
#' @param plag      Lag order for Delta-y. Either a single integer (applied
#'                  to all groups) or a vector of length n. Default 1.
#' @param qlag      Lag orders for Delta-X. Either a single integer (applied
#'                  to all regressors and groups), a vector of length k
#'                  (applied to all groups), or an n x k matrix. Default 1.
#' @param ref_group Name of the reference group for relative variables
#'                  (string, default NULL). When set, the last column of X
#'                  is treated as relative to this group: its qlag is forced
#'                  to 0, and this group is placed last and excluded from
#'                  MG/PMG (NN = n-1). SFE/DFE still use all groups.
#'
#' @return A list suitable for passing to MG(), PMG(), SFE(), DFE1():
#'   Y, X, Z, Ti, n, k, plag, qlag, NN, group_names, x_names
#'
prepare_panel <- function(data, y, x, unit = NULL, time = NULL,
                          z = NULL, plag = 1L, qlag = 1L,
                          ref_group = NULL) {

  # --- Extract panel indices ---
  # Handle pdata.frame (plm package)
  if (is.null(unit) && is.null(time) && !is.null(attr(data, "index"))) {
    idx <- attr(data, "index")
    unit <- names(idx)[1]
    time <- names(idx)[2]
  }
  if (is.null(unit) || is.null(time)) {
    stop("Must specify 'unit' and 'time' column names ",
         "(or pass a pdata.frame with index attributes)")
  }

  # Convert to plain data.frame (handles data.table, tibble, pdata.frame)
  df <- as.data.frame(data)

  # Validate columns exist
  all_cols <- c(y, x, unit, time, z)
  missing  <- setdiff(all_cols, names(df))
  if (length(missing) > 0) {
    stop("Columns not found in data: ", paste(missing, collapse = ", "))
  }

  rel_var <- !is.null(ref_group)

  # --- Sort by unit then time ---
  df <- df[order(df[[unit]], df[[time]]), ]

  # --- Drop rows with NA in any estimation variable ---
  est_cols <- c(y, x, z)
  complete <- complete.cases(df[, est_cols, drop = FALSE])
  if (sum(!complete) > 0) {
    message(sprintf("Dropped %d rows with NA values", sum(!complete)))
    df <- df[complete, ]
  }

  # --- Build group structure ---
  groups      <- unique(df[[unit]])
  n           <- length(groups)
  group_names <- as.character(groups)

  # When a reference group is specified, move it to the last position
  # so that NN = n-1 correctly excludes it from MG/PMG
  if (rel_var) {
    if (!(ref_group %in% group_names)) {
      stop("ref_group '", ref_group, "' not found in unit column. ",
           "Available groups: ", paste(head(group_names, 10), collapse = ", "),
           if (n > 10) ", ..." else "")
    }
    ref_idx <- which(group_names == ref_group)
    groups      <- c(groups[-ref_idx], groups[ref_idx])
    group_names <- as.character(groups)
  }

  Ti <- integer(n)
  Y_list <- vector("list", n)
  X_list <- vector("list", n)
  Z_list <- vector("list", n)

  for (i in seq_len(n)) {
    mask <- df[[unit]] == groups[i]
    gi   <- df[mask, ]
    Ti[i]     <- nrow(gi)
    Y_list[[i]] <- gi[[y]]
    X_list[[i]] <- as.matrix(gi[, x, drop = FALSE])
    if (!is.null(z)) {
      Z_list[[i]] <- cbind(as.matrix(gi[, z, drop = FALSE]),
                           rep(1, nrow(gi)))
    } else {
      Z_list[[i]] <- matrix(1, nrow = nrow(gi), ncol = 1)
    }
  }

  Y_out <- do.call(c, Y_list)
  X_out <- do.call(rbind, X_list)
  Z_out <- do.call(rbind, Z_list)
  k     <- length(x)

  # Remove column names from matrices (estimators expect plain matrices)
  x_names <- x
  colnames(X_out) <- NULL

  # --- Build lag specification ---
  plag_vec <- if (length(plag) == 1L) rep(as.integer(plag), n) else as.integer(plag)
  if (length(plag_vec) != n) {
    stop("plag must be length 1 or length n (", n, ")")
  }

  if (is.matrix(qlag)) {
    qlag_mat <- qlag
  } else if (length(qlag) == 1L) {
    qlag_mat <- matrix(as.integer(qlag), nrow = n, ncol = k)
  } else if (length(qlag) == k) {
    qlag_mat <- matrix(rep(as.integer(qlag), each = n), nrow = n, ncol = k)
  } else {
    stop("qlag must be a scalar, a vector of length k (", k, "), or an n x k matrix")
  }

  # Force relative variable lag to 0
  if (rel_var) {
    qlag_mat[, k] <- 0L
  }

  NN <- if (rel_var) n - 1L else n

  list(
    Y = Y_out, X = X_out, Z = Z_out,
    Ti = Ti, n = n, k = k,
    plag = plag_vec, qlag = qlag_mat,
    NN = NN, rel_var = rel_var,
    group_names = group_names,
    x_names = x_names,
    y_name = y
  )
}


# ==============================================================================
# Label helpers (attach names to raw estimator output)
# ==============================================================================

label_theta <- function(theta, theta_se, theta_t, x_names) {
  out <- data.frame(
    coef = theta, se = theta_se, t = theta_t,
    row.names = x_names
  )
  out
}

label_phi <- function(phi_all, group_names) {
  names(phi_all) <- group_names
  phi_all
}


# ==============================================================================
# Convenience wrappers: data.frame in, labelled results out
# ==============================================================================

#' Run PMG from a data.frame
#'
#' @param data   data.frame with unit, time, and variable columns
#' @param y      Dependent variable name
#' @param x      Regressor names (character vector)
#' @param unit   Unit/group column name
#' @param time   Time column name
#' @param ...    Additional arguments passed to prepare_panel() and PMG_NR()
#'
#' @return PMG result list with added $lr (labelled long-run data.frame),
#'   $group_phi (named vector), $summary (printed table)
#'
pmg_panel <- function(data, y, x, unit, time,
                      z = NULL, plag = 1L, qlag = 1L,
                      ref_group = NULL, ini_theta = NULL,
                      maxiter = 1000L, tol = 1e-4, verbose = FALSE) {

  p   <- prepare_panel(data, y, x, unit, time, z, plag, qlag, ref_group)
  res <- PMG_NR(p$Y, p$X, p$Z, p$n, p$Ti, p$k, p$plag, p$qlag,
                NN = p$NN, ini_theta = ini_theta,
                maxiter = maxiter, tol = tol, verbose = verbose)

  # Attach labels
  res$lr        <- label_theta(res$theta, res$theta_se, res$theta_t, p$x_names)
  res$group_phi <- label_phi(res$phi_all, p$group_names[seq_len(p$NN)])
  res$panel     <- p

  class(res) <- "pmg_result"
  res
}

#' @export
print.pmg_result <- function(x, ...) {
  cat("Pooled Mean Group Estimator (Newton-Raphson)\n")
  cat(sprintf("Converged: %s (%d iterations)\n", x$converged, x$iterations))
  cat(sprintf("Log-likelihood: %.4f\n\n", x$loglik))
  cat("Long-Run Coefficients:\n")
  print(round(x$lr, 4))
  cat(sprintf("\nError Correction (phi): %.4f (se=%.4f)\n", x$phi, x$phi_se))
  cat(sprintf("Groups: %d\n", length(x$group_phi)))
}


#' Run MG from a data.frame
mg_panel <- function(data, y, x, unit, time,
                     z = NULL, plag = 1L, qlag = 1L, ref_group = NULL) {

  p   <- prepare_panel(data, y, x, unit, time, z, plag, qlag, ref_group)
  res <- MG(p$Y, p$X, p$Z, p$n, p$Ti, p$k, p$plag, p$qlag, NN = p$NN)

  res$lr        <- label_theta(res$theta, res$theta_se, res$theta_t, p$x_names)
  res$group_phi <- label_phi(res$phi_all, p$group_names[seq_len(p$NN)])
  res$panel     <- p

  class(res) <- "mg_result"
  res
}

#' @export
print.mg_result <- function(x, ...) {
  cat("Mean Group Estimator\n\n")
  cat("Long-Run Coefficients:\n")
  print(round(x$lr, 4))
  cat(sprintf("\nError Correction (phi): %.4f (se=%.4f)\n", x$phi, x$phi_se))
  cat(sprintf("Groups: %d\n", length(x$group_phi)))
}


#' Run DFE from a data.frame
dfe_panel <- function(data, y, x, unit, time,
                      z = NULL, plag = 1L, qlag = 1L, ref_group = NULL) {

  p <- prepare_panel(data, y, x, unit, time, z, plag, qlag, ref_group)

  # Compute numot from first group
  ytemp   <- p$Y[1:p$Ti[1]]
  xtemp   <- as.matrix(p$X[1:p$Ti[1], , drop = FALSE])
  ztemp   <- as.matrix(p$Z[1:p$Ti[1], , drop = FALSE])
  dgp_tmp <- FDGP1(ytemp, xtemp, ztemp, p$plag, p$qlag, p$k, p$Ti[1], 1L)
  numot   <- if (is.null(dgp_tmp$ww)) 0L else ncol(dgp_tmp$ww)

  res <- DFE1(p$Y, p$X, p$Z, p$n, p$Ti, p$k, numot, p$plag, p$qlag)

  # Label the long-run rows
  lr_df <- data.frame(
    coef = res$dferes[1:p$k, 1],
    se   = res$dferes[1:p$k, 2],
    t    = res$dferes[1:p$k, 3],
    row.names = p$x_names
  )
  res$lr    <- lr_df
  res$panel <- p

  class(res) <- "dfe_result"
  res
}

#' @export
print.dfe_result <- function(x, ...) {
  cat("Dynamic Fixed Effects Estimator\n")
  cat(sprintf("Observations: %d\n\n", x$qobs))
  cat("Long-Run Coefficients:\n")
  print(round(x$lr, 4))
  k <- nrow(x$lr)
  cat(sprintf("\nError Correction (phi): %.4f (se=%.4f)\n",
              x$dferes[k + 1, 1], x$dferes[k + 1, 2]))
}


#' Run SFE from a data.frame
sfe_panel <- function(data, y, x, unit, time,
                      z = NULL, plag = 1L, qlag = 1L, ref_group = NULL) {

  p   <- prepare_panel(data, y, x, unit, time, z, plag, qlag, ref_group)
  res <- SFE(p$Y, p$X, p$Z, p$n, p$Ti, p$k)

  lr_df <- data.frame(
    coef = res$sferes[1:p$k, 1],
    se   = res$sferes[1:p$k, 2],
    t    = res$sferes[1:p$k, 3],
    row.names = p$x_names
  )
  res$lr    <- lr_df
  res$panel <- p

  class(res) <- "sfe_result"
  res
}

#' @export
print.sfe_result <- function(x, ...) {
  cat("Static Fixed Effects Estimator\n")
  cat(sprintf("Observations: %d\n\n", x$qobs))
  cat("Long-Run Coefficients:\n")
  print(round(x$lr, 4))
}
