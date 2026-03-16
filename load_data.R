# ==============================================================================
# Data Loading Functions for Pesaran, Shin & Smith (1999) Replication
#
# Loads the original GAUSS data files for the two empirical examples:
#   Example 1: OECD Consumption Function (24 countries, 1960-1993)
#   Example 2: Asian Energy Demand (10 countries, 1974-1990)
# ==============================================================================


#' Load Example 1: OECD Consumption Function Data
#'
#' Variables:
#'   Y = log real per-capita private consumption (LPC)
#'   X = [log real per-capita national disposable income (LNDI),
#'        inflation rate (DP)]
#'   Z = intercept
#'
#' The last X variable (DP/inflation) is treated as a "relative variable"
#' in the original GAUSS code, meaning its lag order is set to 0 and
#' the last group (U.K.) serves as reference for MG/PMG estimation.
#'
#' @param data_dir Path to directory containing .DAT files
#' @return List with Y, X, Z, Ti, n, country_names, k, plag, qlag, NN
#'
load_example1 <- function(data_dir = "prog") {
  MAXT <- 34L
  MAXN <- 24L
  MISSING <- 8934567

  n_elem <- MAXT * MAXN  # 816

  # Load the three variable files (one value per line, stacked by country)
  lpc  <- scan(file.path(data_dir, "LPC.DAT"),  what = numeric(),
               n = n_elem, quiet = TRUE)
  lndi <- scan(file.path(data_dir, "LNDI.DAT"), what = numeric(),
               n = n_elem, quiet = TRUE)
  dp   <- scan(file.path(data_dir, "DP.DAT"),   what = numeric(),
               n = n_elem, quiet = TRUE)

  country_names <- c(
    "Australi", "Austria",  "Belgium",  "Canada",   "Denmark",
    "Finland",  "France",   "Germany",  "Greece",   "Iceland",
    "Ireland",  "Italy",    "Japan",    "Luxembou", "Netherla",
    "NewZeala", "Norway",   "Portugal", "Spain",    "Sweden",
    "Switzerl", "Turkey",   "U.S.",     "U.K."
  )

  # Delete missing values and count Ti for each group
  Y_list <- vector("list", MAXN)
  X_list <- vector("list", MAXN)
  Ti     <- integer(MAXN)

  for (i in seq_len(MAXN)) {
    idx <- ((i - 1L) * MAXT + 1L):(i * MAXT)
    yi  <- lpc[idx]
    x1i <- lndi[idx]
    x2i <- dp[idx]

    # Drop rows where ANY variable equals the missing value
    valid <- (yi != MISSING) & (x1i != MISSING) & (x2i != MISSING)
    Ti[i]     <- sum(valid)
    Y_list[[i]] <- yi[valid]
    X_list[[i]] <- cbind(x1i[valid], x2i[valid])
  }

  Y <- do.call(c, Y_list)
  X <- do.call(rbind, X_list)
  Z <- matrix(1, nrow = length(Y), ncol = 1)

  # Fixed lag specification matching the OUTPUT.OUT run:
  #   plag = 1 (ARDL lag order for y)
  #   qlag[,1] = 2 (lag order for Δlndi)
  #   qlag[,2] = 0 (lag order for Δdp, the relative variable)
  plag <- rep(1L, MAXN)
  qlag <- cbind(rep(2L, MAXN), rep(0L, MAXN))

  list(
    Y = Y, X = X, Z = Z,
    Ti = Ti, n = MAXN,
    country_names = country_names,
    k = 2L,
    rel_var = TRUE,
    NN = MAXN - 1L,    # MG/PMG use 23 groups (exclude U.K. reference)
    plag = plag,
    qlag = qlag
  )
}


#' Load Example 2: Asian Energy Demand Data
#'
#' Variables:
#'   Y = ln(total energy consumption / population)
#'   X = [ln(real GDP per capita), ln(average energy real price)]
#'   Z = intercept
#'
#' @param data_dir Path to directory containing .DAT files
#' @return List with Y, X, Z, Ti, n, country_names, k, plag, qlag, NN
#'
load_example2 <- function(data_dir = "prog") {
  MAXT <- 18L
  MAXN <- 10L
  MISSING <- 8934567

  n_elem <- MAXT * MAXN  # 180

  # EDA.DAT has 4 columns: AERP, POPUL, TOTC, RGDPL
  eda <- matrix(scan(file.path(data_dir, "EDA.DAT"), what = numeric(),
                     n = n_elem * 4, quiet = TRUE),
                ncol = 4, byrow = TRUE)

  AERP  <- eda[, 1]
  POPUL <- eda[, 2]
  TOTC  <- eda[, 3]
  RGDPL <- eda[, 4]

  # Construct variables as in jasa2.prg
  y_all  <- log(TOTC / POPUL)
  x1_all <- log(RGDPL)
  x2_all <- log(AERP)

  # Try reading country names
  j2name_path <- file.path(data_dir, "J2NAME.DAT")
  if (file.exists(j2name_path)) {
    country_names <- tryCatch(
      trimws(readLines(j2name_path, warn = FALSE)),
      error = function(e) paste0("Country", seq_len(MAXN))
    )
    country_names <- country_names[nchar(country_names) > 0]
    if (length(country_names) < MAXN) {
      country_names <- c(country_names,
                         paste0("Country", seq(length(country_names) + 1, MAXN)))
    }
  } else {
    country_names <- paste0("Country", seq_len(MAXN))
  }

  # Delete missing values
  Y_list <- vector("list", MAXN)
  X_list <- vector("list", MAXN)
  Ti     <- integer(MAXN)

  for (i in seq_len(MAXN)) {
    idx <- ((i - 1L) * MAXT + 1L):(i * MAXT)
    yi  <- y_all[idx]
    x1i <- x1_all[idx]
    x2i <- x2_all[idx]

    valid <- is.finite(yi) & is.finite(x1i) & is.finite(x2i) &
             (TOTC[idx] != MISSING) & (POPUL[idx] != MISSING) &
             (RGDPL[idx] != MISSING) & (AERP[idx] != MISSING)
    Ti[i]     <- sum(valid)
    Y_list[[i]] <- yi[valid]
    X_list[[i]] <- cbind(x1i[valid], x2i[valid])
  }

  Y <- do.call(c, Y_list)
  X <- do.call(rbind, X_list)
  Z <- matrix(1, nrow = length(Y), ncol = 1)

  # Default lag specification: ARDL(1,1,1)
  plag <- rep(1L, MAXN)
  qlag <- cbind(rep(1L, MAXN), rep(1L, MAXN))

  list(
    Y = Y, X = X, Z = Z,
    Ti = Ti, n = MAXN,
    country_names = country_names,
    k = 2L,
    rel_var = FALSE,
    NN = MAXN,
    plag = plag,
    qlag = qlag
  )
}
