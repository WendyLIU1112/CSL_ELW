requested_steps <- "testing"
run_step <- function(step) step %in% requested_steps

if (run_step("testing")) local({
  ## ============================================================
  ## CSL-informed testing example for manuscript Section 3.4
  ## ============================================================
  
  
  get_script_dir <- function() {
    args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed = TRUE)
    file_arg <- grep("^--file=", args, value = TRUE)
    if (length(file_arg) > 0L) {
      script_path <- sub("^--file=", "", file_arg[1])
      script_path <- gsub("~\\+~", " ", script_path)
      return(dirname(normalizePath(script_path)))
    }
    if (!is.null(sys.frames()[[1]]$ofile)) {
      return(dirname(normalizePath(sys.frames()[[1]]$ofile)))
    }
    getwd()
  }
  
  main_dir <- get_script_dir()
  save_dir <- Sys.getenv(
    "CSL_TEST_SAVE_DIR",
    file.path(main_dir, "results", "bootstrap-testing")
  )
  manuscript_dir <- Sys.getenv(
    "CSL_TEST_MANUSCRIPT_DIR",
    file.path(main_dir, "results", "testing-figures")
  )
  dir.create(save_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
  
  setwd(main_dir)
  source(file.path(main_dir, "csl_smoking_all_functions.R"))
  
  input_file <- Sys.getenv(
    "CSL_SIM_DATA",
    file.path(main_dir, "data", "csl_smoking_birthweight_complete_case.csv")
  )
  if (!file.exists(input_file)) {
    stop("Cannot find complete-case CSL file: ", input_file)
  }
  
  set.seed(20260704)
  dat <- read.csv(input_file, check.names = FALSE)
  base <- make_csl_simulation_base(dat)
  
  ## The testing example uses the CSL covariate and outcome distributions but
  ## allows a regularized treatment-assignment design.  A scale below one
  ## preserves the ordering of the CSL logit while moderating limited overlap.
  treatment_prevalence <- as.numeric(
    Sys.getenv("CSL_TEST_TREATMENT_PREVALENCE", "0.15")
  )[1]
  treatment_logit_scale <- as.numeric(
    Sys.getenv("CSL_TEST_TREATMENT_LOGIT_SCALE", "0.75")
  )[1]
  supplied_treatment_logit_offset <- suppressWarnings(as.numeric(
    Sys.getenv("CSL_TEST_TREATMENT_LOGIT_OFFSET", "")
  ))[1]
  if (!is.finite(treatment_prevalence) || treatment_prevalence <= 0 ||
      treatment_prevalence >= 1) {
    stop("CSL_TEST_TREATMENT_PREVALENCE must lie strictly between zero and one.")
  }
  if (!is.finite(treatment_logit_scale) || treatment_logit_scale <= 0) {
    stop("CSL_TEST_TREATMENT_LOGIT_SCALE must be positive.")
  }
  treatment_logit_offset <- if (is.finite(supplied_treatment_logit_offset)) {
    supplied_treatment_logit_offset
  } else {
    stats::uniroot(
      function(offset) {
        mean(stats::plogis(offset + treatment_logit_scale * base$ell)) -
          treatment_prevalence
      },
      interval = c(-20, 20)
    )$root
  }
  realized_treatment_prevalence <- mean(stats::plogis(
    treatment_logit_offset + treatment_logit_scale * base$ell
  ))
  
  mc_reps <- as.integer(Sys.getenv("CSL_TEST_REPS", "1000"))
  mc_reps_equality <- as.integer(Sys.getenv("CSL_TEST_REPS_EQUALITY", "1000"))
  boot_reps <- as.integer(Sys.getenv("CSL_TEST_BOOT", "500"))
  n_values <- as.integer(strsplit(Sys.getenv("CSL_TEST_N", "2000,10000"), ",")[[1]])
  scenario_values <- as.integer(strsplit(
    Sys.getenv("CSL_TEST_SCENARIOS", "2,3"), ","
  )[[1]])
  if (length(scenario_values) == 0L || any(!scenario_values %in% c(2L, 3L))) {
    stop("CSL_TEST_SCENARIOS must contain a subset of 2 and 3.")
  }
  no_cores <- as.integer(Sys.getenv(
    "CSL_TEST_CORES",
    as.character(max(1L, parallel::detectCores(logical = FALSE) - 1L))
  ))
  taus <- c(0.10, 0.25, 0.50, 0.75)
  fsd_alpha <- as.numeric(strsplit(Sys.getenv("CSL_TEST_ALPHA", "0.01,0.05"), ",")[[1]])
  fsd_alpha <- sort(unique(fsd_alpha[is.finite(fsd_alpha) & fsd_alpha > 0 & fsd_alpha < 1]))
  if (length(fsd_alpha) == 0L) fsd_alpha <- c(0.01, 0.05)
  ps_method <- Sys.getenv("CSL_TEST_PS_METHOD", "logistic_linear")
  kernel_bandwidth_constant <- as.numeric(
    Sys.getenv("CSL_TEST_KERNEL_BW_CONSTANT", "1.06")
  )
  if (!is.finite(kernel_bandwidth_constant) || kernel_bandwidth_constant <= 0) {
    kernel_bandwidth_constant <- 1.06
  }
  grid_lower <- as.numeric(Sys.getenv("CSL_TEST_GRID_LOWER", "0"))[1]
  grid_upper <- as.numeric(Sys.getenv("CSL_TEST_GRID_UPPER", "1"))[1]
  if (!is.finite(grid_lower) || grid_lower < 0 || grid_lower >= 1) grid_lower <- 0
  if (!is.finite(grid_upper) || grid_upper <= 0 || grid_upper > 1) grid_upper <- 1
  if (grid_lower >= grid_upper) {
    grid_lower <- 0
    grid_upper <- 1
  }
  qte_alpha <- as.numeric(Sys.getenv("CSL_TEST_QTE_ALPHA", "0.05"))[1]
  negative_shift_intercept <- as.numeric(
    Sys.getenv("CSL_TEST_NEGATIVE_SHIFT_INTERCEPT", "-120")
  )[1]
  negative_shift_slope <- as.numeric(
    Sys.getenv("CSL_TEST_NEGATIVE_SHIFT_SLOPE", "-50")
  )[1]
  if (!is.finite(negative_shift_intercept) || !is.finite(negative_shift_slope) ||
      negative_shift_intercept >= 0 || negative_shift_slope > 0) {
    stop("The negative-shift intercept must be negative and its slope nonpositive.")
  }
  test_mode <- Sys.getenv("CSL_TEST_MODE", "auto")
  cache_file <- Sys.getenv(
    "CSL_TEST_CACHE",
    file.path(save_dir, paste0("csl_testing_bootstrap_cache_", ps_method, ".rds"))
  )
  available_methods <- c("Naive", "SIPW", "ELW")
  methods <- trimws(strsplit(
    Sys.getenv("CSL_TEST_METHODS", "SIPW,ELW"), ",",
    fixed = TRUE
  )[[1]])
  if (length(methods) == 0L || any(!methods %in% available_methods)) {
    stop("CSL_TEST_METHODS must contain a subset of: ",
         paste(available_methods, collapse = ", "))
  }
  
  scenario_effect <- function(r, scenario) {
    switch(
      as.character(scenario),
      "2" = rep(0, length(r)),
      "3" = negative_shift_intercept + negative_shift_slope * r,
      stop("Unknown scenario: ", scenario)
    )
  }
  
  generate_testing_data <- function(base, n, scenario) {
    id <- sample.int(nrow(base$dat), n, replace = TRUE)
    ell <- base$ell[id]
    m0 <- base$m0[id]
    r <- base$risk.percentile[id]
    eps <- sample(base$eps.pool, n, replace = TRUE)
    Y0 <- m0 + base$sigma0 * eps
    delta <- scenario_effect(r, scenario)
    Y1 <- Y0 + delta
    ell <- treatment_logit_offset + treatment_logit_scale * ell
    pi_true <- plogis(ell)
    T <- rbinom(n, 1, pi_true)
    Y <- ifelse(T == 1, Y1, Y0)
    data.frame(
      T = T, Y = Y, Y0 = Y0, Y1 = Y1, pi_true = pi_true, ell = ell,
      risk.percentile = r, delta = delta, check.names = FALSE
    )
  }
  
  weighted_cdf_at_grid <- function(y, weight, grid) {
    ord <- order(y)
    y_ord <- y[ord]
    cw <- cumsum(weight[ord])
    idx <- findInterval(grid, y_ord)
    out <- numeric(length(grid))
    keep <- idx > 0L
    out[keep] <- cw[idx[keep]]
    pmin(pmax(out, 0), 1)
  }
  
  make_distribution_cache <- function(y, grid) {
    ord <- order(y)
    y.ord <- y[ord]
    list(
      ord = ord,
      y.ord = y.ord,
      grid.index = findInterval(grid, y.ord),
      grid.length = length(grid)
    )
  }
  
  distribution_from_weight <- function(weight, cache, probs) {
    weight <- weight[cache$ord]
    weight <- weight / sum(weight)
    cw <- cumsum(weight)
    quantile <- vapply(probs, function(p) {
      cache$y.ord[which(cw >= p)[1L]]
    }, numeric(1))
    cdf <- numeric(cache$grid.length)
    keep <- cache$grid.index > 0L
    cdf[keep] <- cw[cache$grid.index[keep]]
    list(quantile = quantile, cdf = pmin(pmax(cdf, 0), 1))
  }
  
  distribution_matrix_from_weights <- function(weight, cache, probs) {
    weight <- weight[cache$ord, , drop = FALSE]
    weight <- sweep(weight, 2, colSums(weight), "/")
    cw <- apply(weight, 2, cumsum)
    if (is.null(dim(cw))) cw <- matrix(cw, ncol = 1L)
    quantile <- t(vapply(seq_len(ncol(cw)), function(b) {
      vapply(probs, function(p) cache$y.ord[which(cw[, b] >= p)[1L]], numeric(1))
    }, numeric(length(probs))))
    cdf <- matrix(0, nrow = ncol(cw), ncol = cache$grid.length)
    keep <- cache$grid.index > 0L
    cdf[, keep] <- t(cw[cache$grid.index[keep], , drop = FALSE])
    list(quantile = quantile, cdf = pmin(pmax(cdf, 0), 1))
  }
  
  weighted_bin_sum <- function(bin, weight, nb) {
    out <- numeric(nb)
    ok <- is.finite(bin) & !is.na(bin)
    if (any(ok)) {
      z <- rowsum(weight[ok], group = bin[ok], reorder = FALSE)
      out[as.integer(rownames(z))] <- as.numeric(z[, 1])
    }
    out
  }
  
  make_test_grid <- function(y) {
    if (grid_lower <= 0 && grid_upper >= 1) return(sort(unique(y)))
    limits <- stats::quantile(
      y, probs = c(grid_lower, grid_upper), names = FALSE, type = 8
    )
    sort(unique(y[y >= limits[1] & y <= limits[2]]))
  }
  
  logistic_linear_propensity <- function(ell, T, weights = NULL,
                                         max_iter = 30L, tol = 1e-9) {
    mat <- if (is.null(weights)) matrix(1, nrow = length(T), ncol = 1L) else {
      if (is.null(dim(weights))) matrix(as.numeric(weights), ncol = 1L) else weights
    }
    n <- length(T)
    B <- ncol(mat)
    ell <- as.numeric(ell)
    T <- as.numeric(T)
    mat[!is.finite(mat) | mat < 0] <- 0
    empty <- colSums(mat) <= 0
    if (any(empty)) mat[, empty] <- 1
  
    beta0 <- rep(qlogis(weighted.mean(T, rep(1, n))), B)
    beta1 <- rep(0, B)
    ell2 <- ell^2
    for (iter in seq_len(max_iter)) {
      eta <- outer(ell, beta1, "*")
      eta <- sweep(eta, 2, beta0, "+")
      mu <- stats::plogis(eta)
      resid <- mat * (T - mu)
      v <- mat * mu * (1 - mu)
      s0 <- colSums(resid)
      s1 <- colSums(resid * ell)
      h00 <- colSums(v)
      h01 <- colSums(v * ell)
      h11 <- colSums(v * ell2)
      det <- h00 * h11 - h01^2
      ok <- is.finite(det) & det > 1e-12
      step0 <- rep(0, B)
      step1 <- rep(0, B)
      step0[ok] <- (h11[ok] * s0[ok] - h01[ok] * s1[ok]) / det[ok]
      step1[ok] <- (-h01[ok] * s0[ok] + h00[ok] * s1[ok]) / det[ok]
      beta0 <- beta0 + step0
      beta1 <- beta1 + step1
      if (max(abs(step0) + abs(step1), na.rm = TRUE) < tol) break
    }
    eta <- outer(ell, beta1, "*")
    eta <- sweep(eta, 2, beta0, "+")
    out <- stats::plogis(eta)
    if (B == 1L) as.numeric(out[, 1]) else out
  }
  
  
  make_kernel_propensity_cache <- function(ell, bandwidth = NULL) {
    u <- as.numeric(scale(ell))
    n <- length(u)
    if (is.null(bandwidth)) {
      bandwidth <- kernel_bandwidth_constant * stats::sd(u) * n^(-1 / 5)
    }
    if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 0.10
    nbin <- min(500L, max(80L, ceiling(sqrt(n) * 4)))
    breaks <- seq(min(u), max(u), length.out = nbin + 1L)
    if (length(unique(breaks)) < 2L) {
      breaks <- seq(min(u) - 0.5, max(u) + 0.5, length.out = nbin + 1L)
    }
    bin <- findInterval(u, breaks, all.inside = TRUE)
    centers <- (breaks[-1L] + breaks[-length(breaks)]) / 2
    kern <- exp(-0.5 * (outer(u, centers, "-") / bandwidth)^2)
    list(u = u, nbin = nbin, bin = bin, kern = kern)
  }
  
  kernel_propensity <- function(ell, T, weights = NULL, bandwidth = NULL,
                                cache = NULL) {
    n <- length(ell)
    if (is.null(weights)) weights <- rep(1, n)
    weights <- as.numeric(weights)
    weights[!is.finite(weights) | weights < 0] <- 0
    if (!is.finite(sum(weights)) || sum(weights) <= 0) weights <- rep(1, n)
    if (is.null(cache)) cache <- make_kernel_propensity_cache(ell, bandwidth)
    fallback <- stats::weighted.mean(T, weights)
    bin_weight <- weighted_bin_sum(cache$bin, weights, cache$nbin)
    bin_sum <- weighted_bin_sum(cache$bin, weights * T, cache$nbin)
    denom <- as.vector(cache$kern %*% bin_weight)
    numer <- as.vector(cache$kern %*% bin_sum)
    pi_hat <- numer / denom
    pi_hat[!is.finite(pi_hat)] <- fallback
    pi_hat
  }
  
  kernel_propensity_bootstrap <- function(T, multipliers, cache) {
    n <- length(T)
    B <- ncol(multipliers)
    bin_weight <- matrix(0, nrow = cache$nbin, ncol = B)
    bin_sum <- matrix(0, nrow = cache$nbin, ncol = B)
    for (b in seq_len(B)) {
      w <- multipliers[, b]
      bin_weight[, b] <- weighted_bin_sum(cache$bin, w, cache$nbin)
      bin_sum[, b] <- weighted_bin_sum(cache$bin, w * T, cache$nbin)
    }
    denom <- cache$kern %*% bin_weight
    numer <- cache$kern %*% bin_sum
    pi_hat <- numer / denom
    fallback <- colSums(multipliers * T) / colSums(multipliers)
    bad <- !is.finite(pi_hat)
    if (any(bad)) {
      pi_hat[bad] <- fallback[col(pi_hat)[bad]]
    }
    pi_hat
  }
  
  estimate_propensity <- function(dat, weights = NULL, kernel_cache = NULL) {
    if (ps_method == "logistic_linear") {
      return(logistic_linear_propensity(dat$ell, dat$T, weights = weights))
    }
    if (ps_method == "kernel") {
      return(kernel_propensity(dat$ell, dat$T, weights = weights, cache = kernel_cache))
    }
    stop("Unknown CSL_TEST_PS_METHOD: ", ps_method)
  }
  
  estimate_bootstrap_propensity <- function(dat, multipliers, kernel_cache = NULL) {
    if (ps_method == "logistic_linear") {
      return(logistic_linear_propensity(dat$ell, dat$T, weights = multipliers))
    }
    if (ps_method == "kernel") {
      return(kernel_propensity_bootstrap(dat$T, multipliers, kernel_cache))
    }
    stop("Unknown CSL_TEST_PS_METHOD: ", ps_method)
  }
  
  make_g_est_testing_cache <- function(dat) {
    index <- if ("pi_hat" %in% names(dat)) {
      stats::qlogis(dat$pi_hat)
    } else {
      dat$ell
    }
    U <- as.numeric(scale(index))
    U[!is.finite(U)] <- 0
    n <- nrow(dat)
    nbin <- min(240L, max(180L, ceiling(2 * sqrt(n))))
    breaks <- seq(min(U), max(U), length.out = nbin + 1L)
    if (length(unique(breaks)) < 2L) {
      breaks <- seq(min(U) - 0.5, max(U) + 0.5, length.out = nbin + 1L)
    }
    bin <- findInterval(U, breaks, all.inside = TRUE)
    centers <- (breaks[-1L] + breaks[-length(breaks)]) / 2
  
    make_group_cache <- function(group) {
      idx <- dat$T == group
      if (sum(idx) < 5L) stop("Too few observations in treatment group ", group)
      bandwidth <- g_est_bandwidth_constant * stats::sd(U[idx]) *
        sum(idx)^(-1 / 5)
      if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 0.25
      K <- matrix(
        gaussian_kernel(outer(centers, centers, function(a, b) (b - a) / bandwidth)),
        nrow = nbin
      )
      occupied <- tabulate(bin[idx], nbins = nbin) > 0L
      empty.rows <- rowSums(K[, occupied, drop = FALSE]) <= 1e-12
      if (any(empty.rows)) {
        occupied.index <- which(occupied)
        for (j in which(empty.rows)) {
          nearest <- occupied.index[which.min(abs(centers[occupied.index] - centers[j]))]
          K[j, ] <- 0
          K[j, nearest] <- 1
        }
      }
      K <- Matrix::Matrix(K, sparse = TRUE)
      list(idx = idx, K = K, bandwidth = bandwidth)
    }
  
    list(
      nbin = nbin, bin = bin, centers = centers,
      groups = list(`0` = make_group_cache(0L), `1` = make_group_cache(1L))
    )
  }
  
  g_est_group_weights <- function(dat, mult = NULL, group = 1L) {
    if (is.null(mult)) mult <- rep(1, nrow(dat))
    cache <- attr(dat, "g_est_cache", exact = TRUE)
    if (is.null(cache)) cache <- make_g_est_testing_cache(dat)
    gc <- cache$groups[[as.character(group)]]
    obs.mass <- weighted_bin_sum(cache$bin[gc$idx], mult[gc$idx], cache$nbin)
    target.mass <- weighted_bin_sum(cache$bin, mult, cache$nbin)
    denom <- as.vector(gc$K %*% obs.mass)
    coef.target <- (target.mass / sum(mult)) / denom
    coef.target[!is.finite(coef.target)] <- 0
    coef.obs <- as.vector(crossprod(gc$K, coef.target))
    weight <- mult[gc$idx] * coef.obs[cache$bin[gc$idx]]
    weight / sum(weight)
  }
  
  g_est_group_weights_bootstrap <- function(dat, multipliers, group = 1L) {
    cache <- attr(dat, "g_est_cache", exact = TRUE)
    if (is.null(cache)) cache <- make_g_est_testing_cache(dat)
    gc <- cache$groups[[as.character(group)]]
    B <- ncol(multipliers)
    obs.mass <- matrix(0, nrow = cache$nbin, ncol = B)
    target.mass <- matrix(0, nrow = cache$nbin, ncol = B)
    for (b in seq_len(B)) {
      obs.mass[, b] <- weighted_bin_sum(
        cache$bin[gc$idx], multipliers[gc$idx, b], cache$nbin
      )
      target.mass[, b] <- weighted_bin_sum(
        cache$bin, multipliers[, b], cache$nbin
      )
    }
    denom <- as.matrix(gc$K %*% obs.mass)
    coef.target <- sweep(target.mass, 2, colSums(multipliers), "/") / denom
    coef.target[!is.finite(coef.target)] <- 0
    coef.obs <- as.matrix(crossprod(gc$K, coef.target))
    weight <- multipliers[gc$idx, , drop = FALSE] *
      coef.obs[cache$bin[gc$idx], , drop = FALSE]
    sweep(weight, 2, colSums(weight), "/")
  }
  
  elw_group_multiplier <- function(A, pi_hat, mult, group = c("treated", "control"),
                                   tol = 1e-8) {
    group <- match.arg(group)
    mult <- as.numeric(mult)
    if (any(!is.finite(mult)) || any(mult < 0)) {
      stop("Exchangeable multipliers must be finite and nonnegative.")
    }
    n.b <- sum(mult)
    n1.b <- sum(mult * (A == 1))
    r1 <- n1.b / n.b
    if (!is.finite(r1) || r1 <= 0 || r1 >= 1) {
      warning("Degenerate weighted treatment share; using multiplier-normalized weights.")
      idx <- if (group == "treated") A == 1 else A == 0
      w <- mult[idx]
      return(list(alpha = NA_real_, weight = w / sum(w)))
    }
  
    if (group == "treated") {
      idx <- A == 1
      pw <- pi_hat[idx]
      mw <- mult[idx]
      xi <- r1 + (1 - r1) * pw
      low <- min(pw)
      up <- min(xi) - .Machine$double.eps
      root_fun <- function(alpha) {
        denom <- 1 + (1 / r1 - 1) * (pw - alpha) / (1 - alpha)
        sum(mw * (pw - alpha) / denom)
      }
      denom_fun <- function(alpha) {
        1 + (1 / r1 - 1) * (pw - alpha) / (1 - alpha)
      }
    } else {
      idx <- A == 0
      pw <- pi_hat[idx]
      mw <- mult[idx]
      xi <- 1 - r1 * pw
      low <- min(1 - pw)
      up <- min(xi) - .Machine$double.eps
      root_fun <- function(beta) {
        denom <- 1 + (1 - pw - beta) / ((1 / r1 - 1) * (1 - beta))
        sum(mw * (1 - pw - beta) / denom)
      }
      denom_fun <- function(beta) {
        1 + (1 - pw - beta) / ((1 / r1 - 1) * (1 - beta))
      }
    }
  
    f_low <- root_fun(low)
    f_up <- root_fun(up)
    if (!is.finite(f_low) || !is.finite(f_up) || f_low * f_up > 0) {
      warning("Weighted ELW root finding failed for ", group,
              "; using multiplier-normalized weights.")
      w <- mult[idx]
      return(list(alpha = NA_real_, weight = w / sum(w)))
    }
  
    alpha <- stats::uniroot(
      root_fun, c(low, up), f.lower = f_low, f.upper = f_up, tol = tol
    )$root
    weight <- mult[idx] / denom_fun(alpha)
    weight <- weight / sum(weight)
    list(alpha = alpha, weight = weight)
  }
  
  method_weights <- function(T, pi_hat, method, mult = NULL) {
    if (is.null(mult)) mult <- rep(1, length(T))
    idx1 <- T == 1
    idx0 <- T == 0
    if (method == "Naive") {
      w1 <- mult[idx1]
      w0 <- mult[idx0]
      return(list(w1 = w1 / sum(w1), w0 = w0 / sum(w0)))
    }
    if (method == "SIPW") {
      w1 <- mult[idx1] / pi_hat[idx1]
      w0 <- mult[idx0] / (1 - pi_hat[idx0])
      if (any(!is.finite(w1)) || any(!is.finite(w0)) ||
          !is.finite(sum(w1)) || !is.finite(sum(w0)) ||
          sum(w1) <= 0 || sum(w0) <= 0) {
        return(list(w1 = rep(NA_real_, sum(idx1)), w0 = rep(NA_real_, sum(idx0))))
      }
      return(list(w1 = w1 / sum(w1), w0 = w0 / sum(w0)))
    }
    if (method == "ELW") {
      if (all(abs(mult - 1) < 1e-12)) {
        fit1 <- elw_group(T, pi_hat, "treated")
        fit0 <- elw_group(T, pi_hat, "control")
      } else {
        fit1 <- elw_group_multiplier(T, pi_hat, mult, "treated")
        fit0 <- elw_group_multiplier(T, pi_hat, mult, "control")
      }
      return(list(w1 = fit1$weight, w0 = fit0$weight))
    }
    stop("Unknown method: ", method)
  }
  
  estimate_method <- function(dat, method, taus, grid, mult = NULL, pi_hat = NULL) {
    if (method == "G-est") {
      wt <- list(
        w1 = g_est_group_weights(dat, mult = mult, group = 1L),
        w0 = g_est_group_weights(dat, mult = mult, group = 0L)
      )
    } else {
      if (is.null(pi_hat)) pi_hat <- dat$pi_hat
      wt <- method_weights(dat$T, pi_hat, method, mult = mult)
    }
    y1 <- dat$Y[dat$T == 1]
    y0 <- dat$Y[dat$T == 0]
    if (any(!is.finite(wt$w1)) || any(!is.finite(wt$w0))) {
      return(list(
        qte = rep(NA_real_, length(taus)),
        F1 = rep(NA_real_, length(grid)),
        F0 = rep(NA_real_, length(grid))
      ))
    }
    q1 <- weighted_quantile(y1, wt$w1, taus)
    q0 <- weighted_quantile(y0, wt$w0, taus)
    F1 <- weighted_cdf_at_grid(y1, wt$w1, grid)
    F0 <- weighted_cdf_at_grid(y0, wt$w0, grid)
    list(qte = q1 - q0, F1 = F1, F0 = F0)
  }
  
  bootstrap_method <- function(dat, method, taus, grid, multipliers, pi_boot = NULL) {
    B <- ncol(multipliers)
    qte <- matrix(NA_real_, nrow = B, ncol = length(taus))
    cdf_diff <- matrix(NA_real_, nrow = B, ncol = length(grid))
    idx1 <- dat$T == 1
    idx0 <- dat$T == 0
    cache1 <- make_distribution_cache(dat$Y[idx1], grid)
    cache0 <- make_distribution_cache(dat$Y[idx0], grid)
    if (method == "G-est") {
      W1 <- g_est_group_weights_bootstrap(dat, multipliers, group = 1L)
      W0 <- g_est_group_weights_bootstrap(dat, multipliers, group = 0L)
      d1 <- distribution_matrix_from_weights(W1, cache1, taus)
      d0 <- distribution_matrix_from_weights(W0, cache0, taus)
      return(list(qte = d1$quantile - d0$quantile,
                  cdf_diff = d1$cdf - d0$cdf))
    }
    if (method == "Naive") {
      d1 <- distribution_matrix_from_weights(
        multipliers[idx1, , drop = FALSE], cache1, taus
      )
      d0 <- distribution_matrix_from_weights(
        multipliers[idx0, , drop = FALSE], cache0, taus
      )
      return(list(qte = d1$quantile - d0$quantile,
                  cdf_diff = d1$cdf - d0$cdf))
    }
    for (b in seq_len(B)) {
      mult <- multipliers[, b]
      pi_b <- if (is.null(pi_boot)) NULL else pi_boot[, b]
      wt <- tryCatch(
        method_weights(dat$T, pi_b, method, mult = mult),
        error = function(e) {
          stop(
            sprintf(
              "%s bootstrap replicate %d failed while computing weights: %s",
              method, b, conditionMessage(e)
            ),
            call. = FALSE
          )
        }
      )
      if (!all(is.finite(wt$w1)) || !all(is.finite(wt$w0))) {
        stop(
          sprintf("%s bootstrap replicate %d produced nonfinite weights.", method, b),
          call. = FALSE
        )
      }
      d1 <- distribution_from_weight(wt$w1, cache1, taus)
      d0 <- distribution_from_weight(wt$w0, cache0, taus)
      qte[b, ] <- d1$quantile - d0$quantile
      cdf_diff[b, ] <- d1$cdf - d0$cdf
    }
    list(qte = qte, cdf_diff = cdf_diff)
  }
  
  mean_or_na <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    mean(x, na.rm = TRUE)
  }
  
  max_or_na <- function(x) {
    if (all(is.na(x))) return(NA_real_)
    max(x, na.rm = TRUE)
  }
  
  run_testing_replication <- function(n, scenario, iter) {
    ## Keep the Monte Carlo stream attached to the effect configuration after
    ## relabeling Scenarios 1 and 3.
    seed_scenario <- 4L - scenario
    set.seed(20260704 + 100000L * seed_scenario + 1000L * n + iter)
    dat <- generate_testing_data(base, n = n, scenario = scenario)
    kernel_cache <- if (ps_method == "kernel") make_kernel_propensity_cache(dat$ell) else NULL
    dat$pi_hat <- estimate_propensity(dat, kernel_cache = kernel_cache)
    grid <- make_test_grid(dat$Y)
    boot_mult <- matrix(stats::rexp(n * boot_reps, rate = 1), nrow = n)
    needs_ps_boot <- any(methods %in% c("SIPW", "ELW"))
    pi_boot <- NULL
    if (needs_ps_boot) {
      pi_boot <- estimate_bootstrap_propensity(dat, boot_mult, kernel_cache = kernel_cache)
    }
  
    out_fsd <- list()
    out_qte <- list()
    for (method in methods) {
      est <- estimate_method(dat, method, taus, grid)
      boot <- bootstrap_method(
        dat, method, taus, grid, boot_mult,
        pi_boot = if (method %in% c("SIPW", "ELW")) pi_boot else NULL
      )
      if (nrow(boot$qte) != boot_reps || nrow(boot$cdf_diff) != boot_reps ||
          any(!is.finite(boot$qte)) || any(!is.finite(boot$cdf_diff))) {
        stop(
          sprintf(
            "%s produced fewer than %d valid bootstrap replications in scenario %d, n=%d, iteration %d.",
            method, boot_reps, scenario, n, iter
          ),
          call. = FALSE
        )
      }
  
      cdf_diff <- est$F1 - est$F0
      t_fsd <- sqrt(n) * max_or_na(cdf_diff)
      centered_fsd <- sqrt(n) * sweep(boot$cdf_diff, 2, cdf_diff, "-")
      t_fsd_boot <- apply(centered_fsd, 1, max_or_na)
      if (!is.finite(t_fsd) || length(t_fsd_boot) != boot_reps ||
          any(!is.finite(t_fsd_boot))) {
        stop(
          sprintf(
            "%s produced an invalid FSD statistic in scenario %d, n=%d, iteration %d.",
            method, scenario, n, iter
          ),
          call. = FALSE
        )
      }
      p_fsd <- mean(t_fsd_boot >= t_fsd)
      fsd_row <- data.frame(
        scenario = scenario, n = n, iter = iter, method = method,
        t_fsd = t_fsd, p_fsd = p_fsd,
        valid_boot = sum(!is.na(t_fsd_boot)),
        stringsAsFactors = FALSE
      )
      fsd_row$t_fsd_boot <- I(list(t_fsd_boot))
      out_fsd[[method]] <- fsd_row
  
      centered_qte <- sweep(boot$qte, 2, est$qte, "-")
      qte_rows <- do.call(rbind, lapply(seq_along(taus), function(j) {
        qte_boot_centered <- centered_qte[, j]
        qte_boot <- boot$qte[, j]
        if (!is.finite(est$qte[j]) || length(qte_boot) != boot_reps ||
            any(!is.finite(qte_boot))) {
          stop(
            sprintf(
              "%s produced an invalid QTE bootstrap sample at tau=%.2f in scenario %d, n=%d, iteration %d.",
              method, taus[j], scenario, n, iter
            ),
            call. = FALSE
          )
        }
        t_qte <- sqrt(n) * abs(est$qte[j])
        t_qte_boot <- sqrt(n) * abs(qte_boot_centered)
        p_qte <- mean(t_qte_boot >= t_qte)
        qte_row <- data.frame(
          scenario = scenario, n = n, iter = iter, method = method,
          tau = taus[j], qte_hat = est$qte[j], p_qte = p_qte,
          valid_boot = sum(!is.na(qte_boot_centered)),
          stringsAsFactors = FALSE
        )
        qte_row$qte_boot_centered <- I(list(qte_boot_centered))
        qte_row
      }))
      out_qte[[method]] <- qte_rows
    }
    list(
      fsd = do.call(rbind, out_fsd),
      qte = do.call(rbind, out_qte)
    )
  }
  
  run_testing_simulation <- function() {
    tasks <- do.call(rbind, lapply(scenario_values, function(scenario) {
      reps <- if (scenario == 2L) mc_reps_equality else mc_reps
      expand.grid(
        scenario = scenario, n = n_values, iter = seq_len(reps),
        KEEP.OUT.ATTRS = FALSE
      )
    }))
    cat("Running", nrow(tasks), "testing replications on", no_cores, "core(s).\n")
    run_task <- function(k) {
      task <- tasks[k, ]
      run_testing_replication(n = task$n, scenario = task$scenario, iter = task$iter)
    }
    if (.Platform$OS.type == "unix" && no_cores > 1L) {
      results <- parallel::mclapply(
        seq_len(nrow(tasks)), run_task,
        mc.cores = no_cores, mc.preschedule = TRUE
      )
    } else if (no_cores > 1L) {
      if (!requireNamespace("future.apply", quietly = TRUE)) {
        stop("Install package 'future.apply' for parallel execution on Windows.")
      }
      future::plan(future::multisession, workers = no_cores)
      on.exit(future::plan(future::sequential), add = TRUE)
      results <- future.apply::future_lapply(
        seq_len(nrow(tasks)), run_task, future.seed = FALSE
      )
    } else {
      results <- lapply(seq_len(nrow(tasks)), run_task)
    }
  
    fsd_raw <- do.call(rbind, lapply(results, `[[`, "fsd"))
    qte_raw <- do.call(rbind, lapply(results, `[[`, "qte"))
    cache <- list(
      config = list(
        cache_version = 18L,
        outcome_calibration = "nonsmoker_expanded_stepwise_aic",
        one_sided_direction = "harmful_F1_minus_F0",
        fsd_rejection_rule = "bootstrap_p_value",
        qte_rejection_rule = "centered_bootstrap_tail_probability",
        scenario_roles = c(`2` = "equality_boundary",
                           `3` = "fixed_negative_shift"),
        scenarios = scenario_values,
        treatment_prevalence = treatment_prevalence,
        treatment_logit_scale = treatment_logit_scale,
        treatment_logit_offset = treatment_logit_offset,
        realized_treatment_prevalence = realized_treatment_prevalence,
        negative_shift_intercept = negative_shift_intercept,
        negative_shift_slope = negative_shift_slope,
        bootstrap_refit_propensity = TRUE,
        ps_method = ps_method,
        mc_reps = mc_reps, mc_reps_equality = mc_reps_equality,
        boot_reps = boot_reps, n_values = n_values,
        taus = taus, methods = methods,
        kernel_bandwidth_constant = kernel_bandwidth_constant,
        grid_lower = grid_lower, grid_upper = grid_upper
      ),
      fsd_raw = fsd_raw,
      qte_raw = qte_raw
    )
    saveRDS(cache, cache_file)
    cache
  }
  
  validate_cache <- function(cache) {
    if (!identical(cache$config$cache_version, 18L) ||
        !identical(cache$config$outcome_calibration,
                   "nonsmoker_expanded_stepwise_aic") ||
        !identical(cache$config$one_sided_direction, "harmful_F1_minus_F0")) {
      stop(
        "Bootstrap cache uses an earlier testing specification. ",
        "Run with CSL_TEST_MODE=run to regenerate the cache."
      )
    }
    if (!identical(cache$config$fsd_rejection_rule, "bootstrap_p_value") ||
        !identical(cache$config$qte_rejection_rule,
                   "centered_bootstrap_tail_probability")) {
      stop(
        "Bootstrap cache uses different testing decision rules. ",
        "Run with CSL_TEST_MODE=run to regenerate the cache."
      )
    }
    if (!isTRUE(cache$config$bootstrap_refit_propensity)) {
      stop(
        "Bootstrap cache was created without refitting the propensity score. ",
        "Run with CSL_TEST_MODE=run to regenerate the cache."
      )
    }
    if (!identical(cache$config$ps_method, ps_method)) {
      stop(
        "Bootstrap cache was created with CSL_TEST_PS_METHOD=",
        cache$config$ps_method, ", but current CSL_TEST_PS_METHOD=", ps_method,
        ". Run with CSL_TEST_MODE=run or choose the matching cache."
      )
    }
    if (!isTRUE(all.equal(cache$config$treatment_prevalence,
                          treatment_prevalence)) ||
        !isTRUE(all.equal(cache$config$treatment_logit_scale,
                          treatment_logit_scale)) ||
        !isTRUE(all.equal(cache$config$treatment_logit_offset,
                          treatment_logit_offset)) ||
        !isTRUE(all.equal(cache$config$negative_shift_intercept,
                          negative_shift_intercept)) ||
        !isTRUE(all.equal(cache$config$negative_shift_slope,
                          negative_shift_slope)) ||
        !identical(cache$config$scenarios, scenario_values)) {
      stop(
        "Bootstrap cache uses a different treatment-assignment design. ",
        "Run with CSL_TEST_MODE=run to regenerate the cache."
      )
    }
    if (!identical(cache$config$boot_reps, boot_reps)) {
      stop(
        "Bootstrap cache was created with ", cache$config$boot_reps,
        " replications, but the current setting uses ", boot_reps,
        ". Run with CSL_TEST_MODE=run to regenerate the cache."
      )
    }
    if (!identical(cache$config$mc_reps, mc_reps) ||
        !identical(cache$config$mc_reps_equality, mc_reps_equality)) {
      stop(
        "Bootstrap cache uses different Monte Carlo replication counts. ",
        "Run with CSL_TEST_MODE=run or use matching CSL_TEST_REPS settings."
      )
    }
    if (!all(methods %in% cache$config$methods)) {
      stop(
        "Bootstrap cache does not contain all requested methods. ",
        "Run with CSL_TEST_MODE=run or use matching CSL_TEST_METHODS."
      )
    }
    fsd_check <- cache$fsd_raw[cache$fsd_raw$method %in% methods, , drop = FALSE]
    qte_check <- cache$qte_raw[cache$qte_raw$method %in% methods, , drop = FALSE]
    if (any(fsd_check$valid_boot != boot_reps) ||
        any(qte_check$valid_boot != boot_reps)) {
      stop(
        "Bootstrap cache contains incomplete replications. ",
        "Repair or regenerate the cache before summarizing it."
      )
    }
    invisible(TRUE)
  }
  
  summarize_fsd <- function(fsd_raw, alpha = fsd_alpha) {
    rows <- lapply(seq_len(nrow(fsd_raw)), function(i) {
      p_value <- fsd_raw$p_fsd[i]
      reject <- if (!is.finite(p_value)) NA_real_ else as.numeric(p_value < alpha)
      data.frame(
        scenario = fsd_raw$scenario[i],
        n = fsd_raw$n[i],
        iter = fsd_raw$iter[i],
        method = fsd_raw$method[i],
        alpha = alpha,
        reject = reject,
        p_value = p_value,
        stringsAsFactors = FALSE
      )
    })
    fsd_results <- do.call(rbind, rows)
    fsd_summary <- aggregate(
      reject ~ alpha + scenario + n + method,
      data = fsd_results,
      FUN = mean_or_na,
      na.action = NULL
    )
    list(results = fsd_results, summary = fsd_summary)
  }
  
  summarize_qte <- function(qte_raw, alpha = qte_alpha) {
    p_qte <- vapply(seq_len(nrow(qte_raw)), function(i) {
      qte_boot_centered <- qte_raw$qte_boot_centered[[i]]
      qte_hat <- qte_raw$qte_hat[i]
      if (!is.finite(qte_hat) || length(qte_boot_centered) != boot_reps ||
          any(!is.finite(qte_boot_centered))) {
        stop("Cached QTE results contain an incomplete bootstrap sample.")
      }
      mean(abs(qte_boot_centered) >= abs(qte_hat))
    }, numeric(1))
    reject <- ifelse(is.finite(p_qte), as.numeric(p_qte < alpha), NA_real_)
    qte_results <- qte_raw[, c("scenario", "n", "iter", "method", "tau", "p_qte")]
    qte_results$p_qte <- p_qte
    qte_results$alpha <- alpha
    qte_results$reject <- reject
    qte_summary <- aggregate(
      reject ~ scenario + n + tau + method,
      data = subset(qte_results, scenario %in% c(2, 3)),
      FUN = mean_or_na,
      na.action = NULL
    )
    list(results = qte_results, summary = qte_summary)
  }
  
  write_testing_outputs <- function(cache) {
    if (!is.null(cache$config$n_values)) n_values <<- cache$config$n_values
    if (!is.null(cache$config$taus)) taus <<- cache$config$taus
  
    fsd_raw <- cache$fsd_raw[cache$fsd_raw$method %in% methods, , drop = FALSE]
    qte_raw <- cache$qte_raw[cache$qte_raw$method %in% methods, , drop = FALSE]
    fsd <- summarize_fsd(fsd_raw, fsd_alpha)
    qte <- summarize_qte(qte_raw, qte_alpha)
  
    write.csv(fsd$results, file.path(save_dir, "testing_fsd_results.csv"), row.names = FALSE)
    write.csv(qte$results, file.path(save_dir, "testing_qte_results.csv"), row.names = FALSE)
    write.csv(fsd$summary, file.path(save_dir, "testing_fsd_summary.csv"), row.names = FALSE)
    write.csv(qte$summary, file.path(save_dir, "testing_qte_summary.csv"), row.names = FALSE)
    saveRDS(
      list(fsd = fsd$results, qte = qte$results,
           fsd_summary = fsd$summary, qte_summary = qte$summary,
           cache_file = cache_file),
      file.path(save_dir, "csl_testing_example_results.rds")
    )
    list(fsd_summary = fsd$summary, qte_summary = qte$summary)
  }
  
  fmt <- function(x) {
    if (length(x) == 0L || is.na(x)) return("\\cellcolor{gray!20}{NA}")
    sprintf("%.3f", x)
  }
  
  fmt_n <- function(x) formatC(as.integer(x), format = "d", big.mark = ",")
  
  write_fsd_latex <- function(tab, path) {
    lookup <- function(scenario, alpha, n, method) {
      z <- tab[tab$scenario == scenario & tab$n == n &
                 abs(tab$alpha - alpha) < 1e-10 & tab$method == method, "reject"]
      if (length(z) == 0L) return(NA_real_)
      z[1]
    }
    lines <- c(
      "\\begin{tabularx}{\\linewidth}{lccYY}",
      "\\toprule",
      "\\textbf{$\\alpha$} & \\textbf{Scenario} & \\textbf{$n$} & \\textbf{SIPW} & \\textbf{ELW} \\\\",
      "\\midrule"
    )
    retained_scenarios <- c(2L, 3L)
    alpha_span <- length(retained_scenarios) * length(n_values)
    scenario_span <- length(n_values)
    for (alpha_i in seq_along(fsd_alpha)) {
      alpha <- fsd_alpha[alpha_i]
      if (alpha_i > 1L) lines <- c(lines, "\\midrule")
      for (scenario in retained_scenarios) {
        for (n_i in seq_along(n_values)) {
          n <- n_values[n_i]
          alpha_lab <- if (scenario == retained_scenarios[1L] && n_i == 1L && alpha_span > 1L) {
            sprintf("\\multirow{%d}{*}{%.2f}", alpha_span, alpha)
          } else if (scenario == retained_scenarios[1L] && n_i == 1L) {
            sprintf("%.2f", alpha)
          } else {
            ""
          }
          scenario_label <- if (scenario == 2L) "Equality" else "Negative shift"
          scen <- if (n_i == 1L && scenario_span > 1L) {
            sprintf("\\multirow{%d}{*}{%s}", scenario_span, scenario_label)
          } else if (n_i == 1L) {
            scenario_label
          } else {
            ""
          }
          vals <- vapply(c("SIPW", "ELW"), function(m) {
            fmt(lookup(scenario, alpha, n, m))
          }, character(1))
          lines <- c(lines, sprintf("%s & %s & %s & %s & %s \\\\",
                                    alpha_lab, scen, fmt_n(n), vals[1], vals[2]))
        }
      }
    }
    writeLines(c(lines, "\\bottomrule", "\\end{tabularx}"), path)
  }
  
  write_qte_latex <- function(tab, path) {
    lookup <- function(scenario, n, tau, method) {
      z <- tab[tab$scenario == scenario & tab$n == n &
                 abs(tab$tau - tau) < 1e-8 & tab$method == method, "reject"]
      if (length(z) == 0L) return(NA_real_)
      z[1]
    }
    lines <- c(
      "\\begin{tabularx}{\\linewidth}{llcYY}",
      "\\toprule",
      "\\textbf{Scenario} & \\textbf{$n$} & \\textbf{$\\tau$} & \\textbf{SIPW} & \\textbf{ELW} \\\\",
      "\\midrule"
    )
    for (s_i in seq_along(c(2L, 3L))) {
      scenario <- c(2L, 3L)[s_i]
      scenario_label <- if (scenario == 2L) "Equality" else "Negative shift"
      lines <- c(lines, sprintf("\\multirow{%d}{*}{%s}",
                                length(n_values) * length(taus), scenario_label))
      for (n_i in seq_along(n_values)) {
        n <- n_values[n_i]
        if (n_i > 1L) lines <- c(lines, "\\addlinespace")
        for (tau_i in seq_along(taus)) {
          tau <- taus[tau_i]
          nlab <- if (tau_i == 1L && length(taus) > 1L) {
            sprintf("\\multirow{%d}{*}{%s}", length(taus), fmt_n(n))
          } else if (tau_i == 1L) {
            fmt_n(n)
          } else {
            ""
          }
          vals <- vapply(c("SIPW", "ELW"), function(m) {
            fmt(lookup(scenario, n, tau, m))
          }, character(1))
          lines <- c(lines, sprintf("& %s & 0.%02d & %s & %s \\\\",
                                    nlab, as.integer(round(100 * tau)),
                                    vals[1], vals[2]))
        }
      }
      if (s_i == 1L) lines <- c(lines, "\\midrule")
    }
    writeLines(c(lines, "\\bottomrule", "\\end{tabularx}"), path)
  }
  
  if (test_mode == "run" || (test_mode == "auto" && !file.exists(cache_file))) {
    cache <- run_testing_simulation()
  } else if (test_mode %in% c("summarize", "auto")) {
    if (!file.exists(cache_file)) {
      stop("Bootstrap cache not found: ", cache_file,
           ". Run with CSL_TEST_MODE=run to create it.")
    }
    cache <- readRDS(cache_file)
    validate_cache(cache)
    cat("Loaded bootstrap cache from ", cache_file, ".\n", sep = "")
  } else {
    stop("Unknown CSL_TEST_MODE: ", test_mode)
  }
  
  summaries <- write_testing_outputs(cache)
  write_fsd_latex(
    summaries$fsd_summary,
    file.path(manuscript_dir, "csl_simulation_fsd_test_table.tex")
  )
  write_qte_latex(
    summaries$qte_summary,
    file.path(manuscript_dir, "csl_simulation_point_test_table.tex")
  )
  
  cat("\nTesting example complete in ", test_mode, " mode.\n", sep = "")
})
