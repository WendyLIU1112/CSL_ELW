# CSL application workflow in one file. Restricted data are not distributed.
# Usage: Rscript application.R STEP CLEAN_CSV OUTPUT_DIR --run
# Without --run, only the requested plan is printed.
# STEP: propensity, descriptive, distribution, heterogeneity, ps_robustness,
#       weight_stress, msm, lmz, lmz_gamma
# Scientific settings are retained from the separate application scripts.

.application_load_base <- function(envir = parent.frame()) {
evalq({
## ============================================================
## General numeric utilities
## Basic matrix inversion and standardization helpers.
## ============================================================

safe_inverse <- function(M, ridge = 1e-8) {
  M <- as.matrix(M)
  out <- tryCatch(
    solve(M + diag(ridge, ncol(M))),
    error = function(e) NULL
  )
  if (!is.null(out)) return(out)
  MASS::ginv(M + diag(sqrt(ridge), ncol(M)))
}

## ============================================================
## Manuscript table and display helpers
## Format descriptive, quantile, DTE, and sensitivity tables.
## ============================================================

## Grouped LaTeX table writers suppress repeated method or Lambda labels.
write_quantile_latex <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccc}",
    "\\toprule",
    "\\textbf{Method} & \\textbf{Quantile} & \\textbf{$\\hat q_1$} & \\textbf{$\\hat q_0$} & \\textbf{Difference} \\\\",
    "\\midrule"
  )
  groups <- split(seq_len(nrow(tab)), match(tab$method, unique(tab$method)))
  body <- unlist(lapply(seq_along(groups), function(g) {
    idx <- groups[[g]]
    rows <- vapply(seq_along(idx), function(j) {
      row <- tab[idx[j], ]
      method <- if (j == 1L) latex_escape(row$method) else ""
      sprintf(
        "%s & %s & %.0f & %.0f & %.0f \\\\",
        method, row$quantile, row$q1, row$q0, row$difference
      )
    }, character(1))
    if (g < length(groups)) c(rows, "\\addlinespace[0.3em]") else rows
  }), use.names = FALSE)
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}

write_msm_latex <- function(tab, path) {
  tau_label <- function(x) {
    x <- as.character(x)
    x <- sub("th$", "", x)
    sprintf("$0.%02d$", as.integer(x))
  }
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccc}",
    "\\toprule",
    "\\textbf{$\\Lambda$} & \\textbf{$\\tau$} & \\textbf{Lower effect} & \\textbf{Upper effect} & \\textbf{Contains zero} \\\\",
    "\\midrule"
  )
  groups <- split(seq_len(nrow(tab)), match(tab$Lambda, unique(tab$Lambda)))
  body <- unlist(lapply(seq_along(groups), function(g) {
    idx <- groups[[g]]
    rows <- vapply(seq_along(idx), function(j) {
      row <- tab[idx[j], ]
      Lambda <- if (j == 1L) sprintf("%.2f", row$Lambda) else ""
      sprintf(
        "%s & %s & %.0f & %.0f & %s \\\\",
        Lambda, tau_label(row$quantile), row$effect_lower, row$effect_upper,
        ifelse(as.logical(row$contains_zero), "Yes", "No")
      )
    }, character(1))
    if (g < length(groups)) c(rows, "\\addlinespace[0.3em]") else rows
  }), use.names = FALSE)
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}

standardize_matrix <- function(x) {
  means <- colMeans(x)
  sds <- apply(x, 2, stats::sd)
  keep <- is.finite(sds) & sds > 0
  sweep(
    sweep(x[, keep, drop = FALSE], 2, means[keep], "-"),
    2, sds[keep], "/"
  )
}

## ============================================================
## ELW core and weighted empirical summaries
## Compute ELW weights, weighted quantiles, and probabilities.
## ============================================================

elw_group <- function(A, pi_hat, group = c("treated", "control"),
                      tol = 1e-8) {
  group <- match.arg(group)
  if (length(A) != length(pi_hat) || !length(A) || anyNA(A) ||
      any(!A %in% c(0, 1)) || any(!is.finite(pi_hat)) ||
      any(pi_hat < 0 | pi_hat > 1)) stop("Invalid treatment/probability inputs.")
  r1 <- mean(A == 1)
  if (r1 == 0 || r1 == 1) stop("ELW requires two nonempty treatment groups.")
  selected <- A == if (group == "treated") 1 else 0
  supplied <- if (group == "treated") pi_hat[selected] else 1 - pi_hat[selected]
  if (all(supplied == supplied[1])) {
    if (supplied[1] == 1) stop("All-one supplied selection probabilities are incompatible with n_t < n.")
    return(list(alpha = supplied[1], weight = rep(1 / sum(selected), sum(selected))))
  }

  if (group == "treated") {
    idx <- A == 1
    pw <- pi_hat[idx]
    xi <- r1 + (1 - r1) * pw
    low <- min(pw)
    up <- min(xi) - .Machine$double.eps
    root_fun <- function(alpha) sum((pw - alpha) / (xi - alpha))
    denom_fun <- function(alpha) {
      1 + (1 / r1 - 1) * (pw - alpha) / (1 - alpha)
    }
  } else {
    idx <- A == 0
    pw <- pi_hat[idx]
    xi <- 1 - r1 * pw
    low <- min(1 - pw)
    up <- min(xi) - .Machine$double.eps
    root_fun <- function(beta) sum((1 - pw - beta) / (xi - beta))
    denom_fun <- function(beta) {
      1 + (1 - pw - beta) / ((1 / r1 - 1) * (1 - beta))
    }
  }

  f_low <- root_fun(low)
  f_up <- root_fun(up)
  if (!is.finite(f_low) || !is.finite(f_up) || f_low * f_up > 0) {
    warning("ELW root finding failed for ", group, "; using uniform weights.")
    return(list(alpha = NA_real_, weight = rep(1 / sum(idx), sum(idx))))
  }

  alpha <- stats::uniroot(
    root_fun, c(low, up), f.lower = f_low, f.upper = f_up, tol = tol
  )$root
  weight <- 1 / denom_fun(alpha)
  weight <- weight / sum(weight)
  list(alpha = alpha, weight = weight)
}

weighted_quantile <- function(y, weight, probs) {
  ord <- order(y)
  y <- y[ord]
  weight <- weight[ord] / sum(weight)
  cdf <- cumsum(weight)
  vapply(probs, function(p) y[which(cdf >= p)[1]], numeric(1))
}

weighted_probability <- function(y, weight, predicate) {
  sum(weight * predicate(y)) / sum(weight)
}

effective_sample_size <- function(weight) {
  sum(weight)^2 / sum(weight^2)
}

weighted_mean_numeric <- function(x, weight) {
  weight <- weight / sum(weight)
  as.numeric(crossprod(weight, x))
}

weighted_var_numeric <- function(x, weight) {
  weight <- weight / sum(weight)
  mu <- weighted_mean_numeric(x, weight)
  as.numeric(crossprod(weight, (x - mu)^2))
}

simulation_balance_matrix <- function(sim.dat) {
  drop.vars <- c(
    "id", "T", "Y", "Y0", "Y1", "pi.true", "rho1.true", "rho0.true",
    "ell", "r", "sigma", "Z", "lbw.risk", "risk.percentile", "delta"
  )
  xdat <- sim.dat[setdiff(names(sim.dat), drop.vars)]
  xdat <- xdat[vapply(xdat, function(x) is.numeric(x) || is.factor(x), logical(1))]
  mm <- stats::model.matrix(~ . - 1, data = xdat)
  mm[, apply(mm, 2, stats::var, na.rm = TRUE) > 0, drop = FALSE]
}

simulation_smd_summary <- function(sim.dat, weight1 = NULL, weight0 = NULL) {
  X <- simulation_balance_matrix(sim.dat)
  if (ncol(X) == 0L) {
    return(c(max.abs.smd = NA_real_, mean.abs.smd = NA_real_))
  }
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  if (is.null(weight1)) weight1 <- rep(1, sum(idx1))
  if (is.null(weight0)) weight0 <- rep(1, sum(idx0))
  if (any(!is.finite(weight1)) || any(!is.finite(weight0)) ||
      sum(weight1) <= 0 || sum(weight0) <= 0) {
    return(c(max.abs.smd = NA_real_, mean.abs.smd = NA_real_))
  }
  smd <- vapply(seq_len(ncol(X)), function(j) {
    x1 <- X[idx1, j]
    x0 <- X[idx0, j]
    m1 <- weighted_mean_numeric(x1, weight1)
    m0 <- weighted_mean_numeric(x0, weight0)
    ## Keep the SMD scale fixed across methods by using the unweighted
    ## within-group variances; only the group means are reweighted.
    v1 <- stats::var(x1)
    v0 <- stats::var(x0)
    den <- sqrt((v1 + v0) / 2)
    if (!is.finite(den) || den <= 0) return(NA_real_)
    (m1 - m0) / den
  }, numeric(1))
  abs.smd <- abs(smd[is.finite(smd)])
  if (length(abs.smd) == 0L) {
    c(max.abs.smd = NA_real_, mean.abs.smd = NA_real_)
  } else {
    c(max.abs.smd = max(abs.smd), mean.abs.smd = mean(abs.smd))
  }
}

format_median_iqr <- function(x, digits = 1) {
  q <- stats::quantile(x, c(0.25, 0.50, 0.75))
  sprintf(
    paste0("%.", digits, "f (%.", digits, "f, %.", digits, "f)"),
    q[2], q[1], q[3]
  )
}

format_binary <- function(x, positive = 2) {
  sprintf("%d (%.1f\\%%)", sum(x == positive), 100 * mean(x == positive))
}

format_p <- function(p) {
  ifelse(p < 0.001, "$<0.001$", sprintf("%.3f", p))
}

latex_escape <- function(x) {
  x <- gsub("&", "\\\\&", x, fixed = TRUE)
  x <- gsub("%", "\\\\%", x, fixed = TRUE)
  x <- gsub("_", "\\\\_", x, fixed = TRUE)
  x
}

write_descriptive_latex <- function(tab, path, n_all, n0, n1) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccc}",
    "\\toprule",
    "\\textbf{Characteristic} & \\textbf{Overall} & \\textbf{No smoking} & \\textbf{Smoking} & \\textbf{$p$-value} \\\\",
    sprintf(
      "& \\textbf{N = %s} & \\textbf{N = %s} & \\textbf{N = %s} & \\\\",
      format(n_all, big.mark = ","), format(n0, big.mark = ","),
      format(n1, big.mark = ",")
    ),
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    paste(
      latex_escape(row[["characteristic"]]), row[["overall"]],
      row[["T0"]], row[["T1"]], row[["p_value"]],
      sep = " & "
    )
  })
  lines <- c(lines, paste0(body, " \\\\"), "\\bottomrule", "\\end{tabular*}")
  writeLines(lines, path)
}

write_quantile_latex_basic <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccc}",
    "\\toprule",
    "\\textbf{Method} & \\textbf{Quantile} & \\textbf{$\\hat q_1$} & \\textbf{$\\hat q_0$} & \\textbf{Difference} \\\\",
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    sprintf(
      "%s & %s & %.0f & %.0f & %.0f \\\\",
      latex_escape(row[["method"]]), row[["quantile"]],
      as.numeric(row[["q1"]]), as.numeric(row[["q0"]]),
      as.numeric(row[["difference"]])
    )
  })
  lines <- c(lines, body, "\\bottomrule", "\\end{tabular*}")
  writeLines(lines, path)
}

write_msm_latex_basic <- function(tab, path) {
  tau_label <- function(x) {
    x <- as.character(x)
    x <- sub("th$", "", x)
    sprintf("$0.%02d$", as.integer(x))
  }
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccc}",
    "\\toprule",
    "\\textbf{$\\Lambda$} & \\textbf{$\\tau$} & \\textbf{Lower effect} & \\textbf{Upper effect} & \\textbf{Contains zero} \\\\",
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    sprintf(
      "%.2f & %s & %.0f & %.0f & %s \\\\",
      as.numeric(row[["Lambda"]]), tau_label(row[["quantile"]]),
      as.numeric(row[["effect_lower"]]),
      as.numeric(row[["effect_upper"]]),
      ifelse(as.logical(row[["contains_zero"]]), "Yes", "No")
    )
  })
  lines <- c(lines, body, "\\bottomrule", "\\end{tabular*}")
  writeLines(lines, path)
}

write_dte_latex <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccc}",
    "\\toprule",
    "\\textbf{Method} & \\textbf{$\\hat F_1(2499)$} & \\textbf{$\\hat F_0(2499)$} & \\textbf{Risk difference} \\\\",
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    sprintf(
      "%s & %.3f & %.3f & %.3f \\\\",
      latex_escape(row[["method"]]),
      as.numeric(row[["risk_smoking"]]),
      as.numeric(row[["risk_no_smoking"]]),
      as.numeric(row[["risk_difference"]])
    )
  })
  lines <- c(lines, body, "\\bottomrule", "\\end{tabular*}")
  writeLines(lines, path)
}

## ============================================================
## Marginal sensitivity model utilities
## Build ELW sensitivity intervals under bounded odds-ratio departures.
## ============================================================

msm_probability_bounds <- function(probability, Lambda) {
  list(
    lower = probability /
      (probability + Lambda * (1 - probability)),
    upper = Lambda * probability /
      (1 - probability + Lambda * probability)
  )
}

elw_from_observation_probability <- function(A, rho, treatment_level) {
  if (treatment_level == 1L) {
    elw_group(A, rho, "treated")
  } else {
    elw_group(A, 1 - rho, "control")
  }
}

elw_msm_side_interval <- function(A, Y, propensity, treatment_level,
                                  Lambda, probs, threshold = 2500) {
  observed <- A == treatment_level
  y <- Y[observed]
  e_t_full <- if (treatment_level == 1L) propensity else 1 - propensity
  e_t <- e_t_full[observed]

  ## Lambda = 1 is the MAR model and must reproduce the fitted ELW analysis.
  ## Bypass the endpoint search so that this identity holds exactly.
  if (abs(Lambda - 1) <= sqrt(.Machine$double.eps)) {
    fit <- suppressWarnings(
      elw_from_observation_probability(A, e_t_full, treatment_level)
    )
    q <- weighted_quantile(y, fit$weight, probs)
    risk <- weighted_probability(y, fit$weight, function(x) x < threshold)
    return(list(
      q_lower = q, q_upper = q,
      risk_lower = risk, risk_upper = risk
    ))
  }

  bounds <- msm_probability_bounds(e_t, Lambda)
  support <- sort(unique(y))

  fit_for_partition <- function(low_outcome, low_upweighted) {
    rho <- if (low_upweighted) {
      ifelse(low_outcome, bounds$lower, bounds$upper)
    } else {
      ifelse(low_outcome, bounds$upper, bounds$lower)
    }
    rho_full <- e_t_full
    rho_full[observed] <- rho
    suppressWarnings(
      elw_from_observation_probability(A, rho_full, treatment_level)
    )
  }

  ## For a fixed y, the two endpoint configurations maximize and minimize
  ## the fitted CDF at y. Inverting these CDF bounds by binary search gives
  ## the same quantile extrema as enumerating every outcome switch point.
  cdf_cache <- new.env(parent = emptyenv())
  cdf_at_index <- function(index, low_upweighted) {
    key <- paste0(if (low_upweighted) "L" else "H", index)
    if (exists(key, envir = cdf_cache, inherits = FALSE)) {
      return(get(key, envir = cdf_cache, inherits = FALSE))
    }
    cut <- support[index]
    fit <- fit_for_partition(y <= cut, low_upweighted)
    value <- weighted_probability(y, fit$weight, function(x) x <= cut)
    assign(key, value, envir = cdf_cache)
    value
  }
  envelope_quantile <- function(prob, low_upweighted) {
    left <- 1L
    right <- length(support)
    while (left < right) {
      middle <- floor((left + right) / 2)
      if (cdf_at_index(middle, low_upweighted) >= prob) {
        right <- middle
      } else {
        left <- middle + 1L
      }
    }
    support[left]
  }

  q_lower <- vapply(
    probs, envelope_quantile, numeric(1), low_upweighted = TRUE
  )
  q_upper <- vapply(
    probs, envelope_quantile, numeric(1), low_upweighted = FALSE
  )
  risk_at_direction <- function(low_upweighted) {
    fit <- fit_for_partition(y < threshold, low_upweighted)
    weighted_probability(y, fit$weight, function(x) x < threshold)
  }
  risks <- c(risk_at_direction(TRUE), risk_at_direction(FALSE))

  list(
    q_lower = q_lower, q_upper = q_upper,
    risk_lower = min(risks), risk_upper = max(risks)
  )
}

elw_msm_sensitivity <- function(A, Y, propensity, Lambda_grid, probs,
                                threshold = 2500) {
  rows <- lapply(Lambda_grid, function(Lambda) {
    side1 <- elw_msm_side_interval(
      A, Y, propensity, 1L, Lambda, probs, threshold
    )
    side0 <- elw_msm_side_interval(
      A, Y, propensity, 0L, Lambda, probs, threshold
    )
    data.frame(
      Lambda = Lambda,
      tau = probs,
      q1_lower = side1$q_lower,
      q1_upper = side1$q_upper,
      q0_lower = side0$q_lower,
      q0_upper = side0$q_upper,
      effect_lower = side1$q_lower - side0$q_upper,
      effect_upper = side1$q_upper - side0$q_lower,
      lbw_rd_lower = side1$risk_lower - side0$risk_upper,
      lbw_rd_upper = side1$risk_upper - side0$risk_lower
    )
  })
  do.call(rbind, rows)
}

## ============================================================
## Shao--Wang single-index MNAR utilities
## Estimate observation probabilities when treatment depends on latent outcomes.
## ============================================================

make_kernel_cache <- function(U, delta, Y, max_bins = 2000L) {
  ord <- order(U)
  out <- list(
    U = U, U_sorted = U[ord], order = ord,
    delta_sorted = delta[ord], Y_sorted = Y[ord]
  )
  if (length(U) > 10000L) {
    bin_sorted <- pmin(
      max_bins,
      ceiling(seq_along(ord) * max_bins / length(ord))
    )
    out$bin <- integer(length(U))
    out$bin[ord] <- bin_sorted
    out$bin_center <- as.numeric(rowsum(U[ord], bin_sorted)) /
      as.numeric(tabulate(bin_sorted, nbins = max_bins))
    out$numerator_bin <- as.numeric(rowsum(1 - delta[ord], bin_sorted))
    out$binned <- TRUE
  } else {
    out$binned <- FALSE
  }
  out
}

kernel_ghat <- function(gamma, cache, bandwidth) {
  if (isTRUE(cache$binned)) {
    denominator_bin <- as.numeric(rowsum(
      cache$delta_sorted * exp(gamma * cache$Y_sorted),
      cache$bin[cache$order]
    ))
    ## Here bandwidth is the standard Gaussian scale in K((u-v)/bandwidth).
    z <- outer(cache$bin_center, cache$bin_center, "-") / bandwidth
    kernel_weight <- exp(-0.5 * z^2)
    numerator_fit <- as.vector(kernel_weight %*% cache$numerator_bin)
    denominator_fit <- as.vector(kernel_weight %*% denominator_bin)
    g_bin <- log(numerator_fit) - log(denominator_fit)
    return(stats::approx(
      cache$bin_center, g_bin, xout = cache$U,
      method = "linear", rule = 2, ties = "ordered"
    )$y)
  }
  numerator <- 1 - cache$delta_sorted
  denominator <- cache$delta_sorted * exp(gamma * cache$Y_sorted)
  numerator_fit_sorted <- stats::ksmooth(
    cache$U_sorted, numerator, kernel = "normal",
    bandwidth = bandwidth / 0.3706506, x.points = cache$U_sorted
  )$y
  denominator_fit_sorted <- stats::ksmooth(
    cache$U_sorted, denominator, kernel = "normal",
    bandwidth = bandwidth / 0.3706506, x.points = cache$U_sorted
  )$y
  numerator_fit <- denominator_fit <- numeric(length(cache$U))
  numerator_fit[cache$order] <- numerator_fit_sorted
  denominator_fit[cache$order] <- denominator_fit_sorted
  log(numerator_fit) - log(denominator_fit)
}

moment_matrix_shao <- function(gamma, cache, delta, Y, Z, bandwidth) {
  g_hat <- kernel_ghat(gamma, cache, bandwidth)
  out <- sapply(levels(Z), function(level) {
    as.numeric(Z == level) *
      (delta * (1 + exp(g_hat + gamma * Y)) - 1)
  })
  colnames(out) <- levels(Z)
  out
}

estimate_shao_side <- function(U, delta, Y, Z,
                               gamma_range = c(-0.01, 0.01),
                               bandwidth = NULL) {
  n <- length(Y)
  if (is.null(bandwidth)) {
    bandwidth <- 1.06 * stats::sd(U) * n^(-1 / 5)
  }
  cache <- make_kernel_cache(U, delta, Y)
  moment_mean <- function(gamma) {
    colMeans(moment_matrix_shao(
      gamma, cache, delta, Y, Z, bandwidth
    ))
  }
  objective <- function(gamma, W = NULL) {
    m <- moment_mean(gamma)
    if (is.null(W)) W <- diag(length(m))
    as.numeric(t(m) %*% W %*% m)
  }

  first <- stats::optimize(objective, gamma_range)
  first_moments <- moment_matrix_shao(
    first$minimum, cache, delta, Y, Z, bandwidth
  )
  W <- safe_inverse(stats::cov(first_moments))
  second <- stats::optimize(function(g) objective(g, W), gamma_range)
  gamma_hat <- second$minimum
  g_hat <- kernel_ghat(gamma_hat, cache, bandwidth)
  rho <- stats::plogis(-(g_hat + gamma_hat * Y))

  list(gamma = gamma_hat, rho = rho, bandwidth = bandwidth)
}

bootstrap_shao_gamma <- function(U_matrix, treatment, Y, Z,
                                 gamma_hat, B = 500,
                                 gamma_range = c(-0.01, 0.01),
                                 seed = 20260615) {
  n <- length(Y)
  Z <- factor(Z)
  z_levels <- levels(Z)
  estimates <- matrix(
    NA_real_, nrow = B, ncol = 2,
    dimnames = list(NULL, c("Y(1)", "Y(0)"))
  )

  set.seed(seed)
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    Ab <- treatment[idx]
    if (length(unique(Ab)) < 2L) next

    Xb <- cbind("(Intercept)" = 1, U_matrix[idx, , drop = FALSE])
    index_fit <- tryCatch(
      suppressWarnings(stats::glm.fit(Xb, Ab, family = stats::binomial())),
      error = function(e) NULL
    )
    if (is.null(index_fit)) next
    coefficients <- index_fit$coefficients
    coefficients[!is.finite(coefficients)] <- 0
    eta <- as.vector(Xb %*% coefficients)
    eta_sd <- stats::sd(eta)
    if (!is.finite(eta_sd) || eta_sd <= 0) next

    score <- (eta - mean(eta)) / eta_sd
    Zb <- factor(Z[idx], levels = z_levels)
    Yb <- Y[idx]
    fits <- tryCatch(
      list(
        estimate_shao_side(
          score, Ab, Yb, Zb, gamma_range = gamma_range
        ),
        estimate_shao_side(
          score, 1 - Ab, Yb, Zb, gamma_range = gamma_range
        )
      ),
      error = function(e) NULL
    )
    if (!is.null(fits)) {
      estimates[b, ] <- c(fits[[1]]$gamma, fits[[2]]$gamma)
    }
  }

  successful <- colSums(is.finite(estimates))
  if (any(successful < max(30, ceiling(0.8 * B)))) {
    warning("Too few successful Shao--Wang bootstrap replicates.")
  }
  se <- apply(estimates, 2, stats::sd, na.rm = TRUE)
  statistic <- gamma_hat / se
  p_value <- 2 * stats::pnorm(-abs(statistic))
  boundary_fraction <- colMeans(
    abs(estimates) >= 0.99 * max(abs(gamma_range)),
    na.rm = TRUE
  )
  ci <- apply(
    estimates, 2, stats::quantile, probs = c(0.025, 0.975),
    na.rm = TRUE, names = FALSE
  )

  list(
    estimates = estimates,
    summary = data.frame(
      potential_outcome = c("Y(1)", "Y(0)"),
      gamma = gamma_hat,
      bootstrap_se = se,
      statistic = statistic,
      p_value = p_value,
      ci_lower = ci[1, ],
      ci_upper = ci[2, ],
      boundary_fraction = boundary_fraction,
      successful_replicates = successful,
      requested_replicates = B
    )
  )
}


## ============================================================
## Li--Qin--Liu conditional-likelihood MNAR utilities
## Estimate outcome-dependent selection probabilities for ELW sensitivity.
## ============================================================

stable_log_mean_exp <- function(x) {
  x <- as.numeric(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_real_)
  m <- max(x)
  m + log(mean(exp(x - m)))
}

lql_mnar_fit_side <- function(Y, X1, X2 = NULL, R) {
  Y <- as.numeric(Y)
  R <- as.integer(R)
  X1 <- as.matrix(X1)
  if (is.null(X2)) {
    X2 <- matrix(numeric(0), nrow = length(Y), ncol = 0)
  } else {
    X2 <- as.matrix(X2)
  }
  if (nrow(X1) != length(Y) || nrow(X2) != length(Y) || length(R) != length(Y)) {
    stop("Y, X1, X2, and R must have compatible dimensions.")
  }

  observed <- R == 1L
  if (sum(observed) <= ncol(X1) + ncol(X2) + 2L) {
    stop("Too few observed outcomes for the LQL observed-outcome model.")
  }

  X.outcome <- cbind("(Intercept)" = 1, X1, X2)
  outcome.fit <- stats::lm.fit(
    x = X.outcome[observed, , drop = FALSE],
    y = Y[observed]
  )
  xi <- outcome.fit$coefficients
  xi[!is.finite(xi)] <- 0
  mu <- as.numeric(X.outcome %*% xi)
  eps <- Y - mu

  X.selection <- cbind("(Intercept)" = 1, X1, mu = mu)
  missing.fit <- tryCatch(
    suppressWarnings(stats::glm.fit(
      x = X.selection,
      y = 1L - R,
      family = stats::binomial()
    )),
    error = function(e) NULL
  )
  if (is.null(missing.fit)) {
    stop("LQL missingness logistic regression failed.")
  }
  eta <- missing.fit$coefficients
  eta[!is.finite(eta)] <- 0
  p.missing <- stats::plogis(as.numeric(X.selection %*% eta))
  info <- crossprod(X.selection, X.selection * as.numeric(p.missing * (1 - p.missing)))
  vcov.eta <- safe_inverse(info)
  eta.se <- sqrt(pmax(diag(vcov.eta), 0))

  a <- unname(eta[1])
  gamma <- unname(eta[length(eta)])
  gamma.se <- unname(eta.se[length(eta.se)])
  gamma.z <- gamma / gamma.se
  gamma.p <- 2 * stats::pnorm(-abs(gamma.z))
  beta <- eta[-c(1, length(eta))]
  log.M <- stable_log_mean_exp(gamma * eps[observed])
  a0 <- a - log.M

  linear.observed <- as.numeric(
    a0 + X1 %*% beta + gamma * Y
  )
  rho <- stats::plogis(-linear.observed)

  list(
    alpha0 = a0,
    alpha = a,
    beta = beta,
    gamma = gamma,
    gamma_se = gamma.se,
    gamma_z = gamma.z,
    gamma_p = gamma.p,
    xi = xi,
    mu = mu,
    residual = eps,
    rho = as.numeric(rho),
    observed = observed,
    log_mgf = log.M,
    outcome_rank = outcome.fit$rank,
    selection_converged = isTRUE(missing.fit$converged)
  )
}

lql_mnar_elw_analysis <- function(A, Y, X1, X2 = NULL,
                                  probs = c(0.10, 0.25, 0.50, 0.75),
                                  threshold = 2500) {
  A <- as.integer(A)
  Y <- as.numeric(Y)
  fit1 <- lql_mnar_fit_side(Y, X1, X2, R = as.integer(A == 1L))
  fit0 <- lql_mnar_fit_side(Y, X1, X2, R = as.integer(A == 0L))

  elw1 <- elw_from_observation_probability(A, fit1$rho, 1L)
  elw0 <- elw_from_observation_probability(A, fit0$rho, 0L)
  q1 <- weighted_quantile(Y[A == 1L], elw1$weight, probs)
  q0 <- weighted_quantile(Y[A == 0L], elw0$weight, probs)
  lbw1 <- weighted_probability(Y[A == 1L], elw1$weight, function(x) x < threshold)
  lbw0 <- weighted_probability(Y[A == 0L], elw0$weight, function(x) x < threshold)

  diagnostics <- data.frame(
    side = c("Y(1)", "Y(0)"),
    gamma = c(fit1$gamma, fit0$gamma),
    gamma_se = c(fit1$gamma_se, fit0$gamma_se),
    gamma_z = c(fit1$gamma_z, fit0$gamma_z),
    gamma_p = c(fit1$gamma_p, fit0$gamma_p),
    alpha0 = c(fit1$alpha0, fit0$alpha0),
    log_mgf = c(fit1$log_mgf, fit0$log_mgf),
    rho_min = c(min(fit1$rho[A == 1L]), min(fit0$rho[A == 0L])),
    rho_q25 = c(
      unname(stats::quantile(fit1$rho[A == 1L], 0.25)),
      unname(stats::quantile(fit0$rho[A == 0L], 0.25))
    ),
    rho_median = c(
      stats::median(fit1$rho[A == 1L]),
      stats::median(fit0$rho[A == 0L])
    ),
    rho_q75 = c(
      unname(stats::quantile(fit1$rho[A == 1L], 0.75)),
      unname(stats::quantile(fit0$rho[A == 0L], 0.75))
    ),
    rho_max = c(max(fit1$rho[A == 1L]), max(fit0$rho[A == 0L])),
    max_weight = c(max(elw1$weight), max(elw0$weight)),
    outcome_rank = c(fit1$outcome_rank, fit0$outcome_rank),
    selection_converged = c(fit1$selection_converged, fit0$selection_converged)
  )

  list(
    fit1 = fit1,
    fit0 = fit0,
    elw1 = elw1,
    elw0 = elw0,
    q1 = q1,
    q0 = q0,
    qte = q1 - q0,
    lbw1 = lbw1,
    lbw0 = lbw0,
    lbw_rd = lbw1 - lbw0,
    diagnostics = diagnostics
  )
}

weighted_interval_probability <- function(y, weight, lower, upper,
                                          lower_open = TRUE,
                                          upper_closed = TRUE) {
  left <- if (lower_open) y > lower else y >= lower
  right <- if (upper_closed) y <= upper else y < upper
  weighted_probability(y, weight, function(x) {
    left <- if (lower_open) x > lower else x >= lower
    right <- if (upper_closed) x <= upper else x < upper
    left & right
  })
}

elw_pte_intervals <- function(
    Y, A, weight1, weight0,
    lower = c(-Inf, 1500, 2500, 4000),
    upper = c(1500, 2500, 4000, Inf),
    labels = c("<1500", "1500--2499", "2500--3999", ">=4000")) {
  rows <- lapply(seq_along(lower), function(j) {
    p1 <- weighted_interval_probability(
      Y[A == 1L], weight1, lower[j], upper[j],
      lower_open = FALSE, upper_closed = FALSE
    )
    p0 <- weighted_interval_probability(
      Y[A == 0L], weight0, lower[j], upper[j],
      lower_open = FALSE, upper_closed = FALSE
    )
    data.frame(
      interval = labels[j],
      P1 = p1,
      P0 = p0,
      PTE = p1 - p0
    )
  })
  do.call(rbind, rows)
}

write_lql_mnar_estimate_latex <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrrrrr}",
    "\\toprule",
    "\\textbf{Method} & \\multicolumn{3}{c}{\\textbf{DTE}} & \\multicolumn{4}{c}{\\textbf{QTE}} \\\\",
    "\\cmidrule(lr){2-4}\\cmidrule(lr){5-8}",
    paste(
      " & \\textbf{$\\hat F_1(2499)$} & \\textbf{$\\hat F_0(2499)$}",
      "& \\textbf{Risk difference} & \\textbf{$\\tau=0.10$} & \\textbf{$\\tau=0.25$}",
      "& \\textbf{$\\tau=0.50$} & \\textbf{$\\tau=0.75$} \\\\"
    ),
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    sprintf(
      "%s & %.3f & %.3f & %.3f & %.0f & %.0f & %.0f & %.0f \\\\",
      latex_escape(row[["method"]]),
      as.numeric(row[["F1_2500"]]),
      as.numeric(row[["F0_2500"]]),
      as.numeric(row[["RD_2500"]]),
      as.numeric(row[["qte_10"]]),
      as.numeric(row[["qte_25"]]),
      as.numeric(row[["qte_50"]]),
      as.numeric(row[["qte_75"]])
    )
  })
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}

write_lql_mnar_dte_latex <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrr}",
    "\\toprule",
    "\\textbf{Method} & \\textbf{$\\hat F_1(2499)$} & \\textbf{$\\hat F_0(2499)$} & \\textbf{Risk difference} \\\\",
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    sprintf(
      "%s & %.3f & %.3f & %.3f \\",
      latex_escape(row[["method"]]),
      as.numeric(row[["F1_2500"]]),
      as.numeric(row[["F0_2500"]]),
      as.numeric(row[["RD_2500"]])
    )
  })
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}

write_lql_mnar_qte_latex <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lrrrr}",
    "\\toprule",
    "\\textbf{Method} & \\textbf{$\\tau=0.10$} & \\textbf{$\\tau=0.25$} & \\textbf{$\\tau=0.50$} & \\textbf{$\\tau=0.75$} \\\\",
    "\\midrule"
  )
  body <- apply(tab, 1, function(row) {
    sprintf(
      "%s & %.0f & %.0f & %.0f & %.0f \\\\",
      latex_escape(row[["method"]]),
      as.numeric(row[["qte_10"]]),
      as.numeric(row[["qte_25"]]),
      as.numeric(row[["qte_50"]]),
      as.numeric(row[["qte_75"]])
    )
  })
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}

write_lql_mnar_pte_latex <- function(tab, path) {
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}llrrr}",
    "\\toprule",
    "\\textbf{Method} & \\textbf{Interval} & \\textbf{$\\hat P_1$} & \\textbf{$\\hat P_0$} & \\textbf{PTE} \\\\",
    "\\midrule"
  )
  groups <- split(seq_len(nrow(tab)), match(tab$method, unique(tab$method)))
  body <- unlist(lapply(seq_along(groups), function(g) {
    idx <- groups[[g]]
    rows <- vapply(seq_along(idx), function(j) {
      row <- tab[idx[j], ]
      method <- if (j == 1L) latex_escape(row$method) else ""
      sprintf(
        "%s & %s & %.3f & %.3f & %.3f \\\\",
        method, row$interval,
        as.numeric(row$P1), as.numeric(row$P0), as.numeric(row$PTE)
      )
    }, character(1))
    if (g < length(groups)) c(rows, "\\addlinespace[0.25em]") else rows
  }), use.names = FALSE)
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}

write_lql_mnar_diagnostic_latex <- function(tab, path, n_by_side = NULL) {
  if (is.null(n_by_side)) {
    tab$max_relative_weight <- tab$max_weight
    max_weight_label <- "\\textbf{Max. weight}"
    max_weight_digits <- "%.5f"
  } else {
    if (!all(tab$side %in% names(n_by_side))) {
      stop("n_by_side must be named for every side in the diagnostic table.")
    }
    tab$max_relative_weight <- as.numeric(tab$max_weight) * n_by_side[tab$side]
    max_weight_label <- "\\textbf{Max. $n_t w$}"
    max_weight_digits <- "%.2f"
  }
  lines <- c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}llrrrrrrr}",
    "\\toprule",
    paste(
      "\\textbf{Model} & \\textbf{Side} & \\textbf{$\\hat\\gamma$} & \\textbf{$p_\\gamma$}",
      "& \\textbf{Min. $\\hat\\rho$}",
      "& \\textbf{Median $\\hat\\rho$} & \\textbf{Max. $\\hat\\rho$}",
      "&", max_weight_label, "& \\textbf{Converged} \\\\"
    ),
    "\\midrule"
  )
  models <- unique(tab$model)
  body <- unlist(lapply(seq_along(models), function(g) {
    rows <- tab[tab$model == models[g], , drop = FALSE]
    out <- vapply(seq_len(nrow(rows)), function(j) {
      row <- rows[j, , drop = FALSE]
      p_gamma <- as.numeric(row$gamma_p)
      p_gamma_cell <- if (is.finite(p_gamma) && p_gamma < 0.001) {
        "$<0.001$"
      } else {
        sprintf("%.3f", p_gamma)
      }
      model_cell <- if (j == 1L) {
        sprintf("\\multirow{%d}{*}{%s}", nrow(rows), latex_escape(row$model))
      } else {
        ""
      }
      sprintf(
        paste0("%s & %s & %.5f & %s & %.3f & %.3f & %.3f & ", max_weight_digits, " & %s \\\\"),
        model_cell,
        row$side,
        as.numeric(row$gamma),
        p_gamma_cell,
        as.numeric(row$rho_min),
        as.numeric(row$rho_median),
        as.numeric(row$rho_max),
        as.numeric(row$max_relative_weight),
        ifelse(tolower(row$selection_converged) == "true", "Yes", "No")
      )
    }, character(1))
    if (g < length(models)) c(out, "\\addlinespace[0.25em]") else out
  }), use.names = FALSE)
  writeLines(c(lines, body, "\\bottomrule", "\\end{tabular*}"), path)
}


## ============================================================
## Simulation LaTeX table helpers
## Write compact manuscript tables from Monte Carlo summaries.
## ============================================================

format_sim_number <- function(x, digits = 2) {
  out <- rep("--", length(x))
  keep <- is.finite(x)
  out[keep] <- formatC(
    x[keep], format = "f", digits = digits, big.mark = ",", preserve.width = "none"
  )
  out
}

format_sim_percent <- function(x, digits = 1) {
  ifelse(is.finite(x), sprintf(paste0("%.", digits, "f"), 100 * x), "--")
}

format_sim_best <- function(x, best, digits = 2, percent = FALSE) {
  out <- if (percent) format_sim_percent(x, digits) else format_sim_number(x, digits)
  best.out <- if (percent) {
    format_sim_percent(best, digits)
  } else {
    format_sim_number(best, digits)
  }
  ## Treat values that are tied at the reported precision as tied in the table.
  if (is.finite(x) && is.finite(best) && identical(out, best.out)) {
    paste0("\\textbf{", out, "}")
  } else {
    out
  }
}

simulation_best <- function(x, direction = c("min", "max")) {
  direction <- match.arg(direction)
  x <- x[is.finite(x)]
  if (length(x) == 0L) return(NA_real_)
  if (direction == "min") min(x) else max(x)
}

simulation_complete_rows <- function(tab) {
  if (!("failure.rate" %in% names(tab))) return(tab)
  tab[!is.finite(tab$failure.rate) | tab$failure.rate <= 0, , drop = FALSE]
}

simulation_benchmark_rows <- function(tab) {
  simulation_complete_rows(tab)
}

simulation_row_is_gray <- function(rows) {
  "failure.rate" %in% names(rows) && any(is.finite(rows$failure.rate) & rows$failure.rate > 0)
}

simulation_row_prefix <- function(rows) {
  if (simulation_row_is_gray(rows)) "\\rowcolor{gray!12} " else ""
}

add_simulation_failure_info <- function(tab, diagnostics = NULL) {
  if (is.null(diagnostics) || !("failure.rate" %in% names(diagnostics))) {
    tab$failure.rate <- 0
    return(tab)
  }
  diag <- prepare_simulation_display_table(diagnostics)
  diag <- diag[, c("scenario_label", "method_label", "failure.rate"), drop = FALSE]
  diag <- diag[!duplicated(diag[c("scenario_label", "method_label")]), , drop = FALSE]
  out <- merge(tab, diag, by = c("scenario_label", "method_label"), all.x = TRUE, sort = FALSE)
  out$failure.rate[is.na(out$failure.rate)] <- 0
  out
}

simulation_ps_label <- function(ps.method) {
  ifelse(ps.method == "PLE", "PLE", ifelse(ps.method == "SLE", "SLE", "Kernel"))
}

simulation_design_label <- function(gamma0, gamma1, b = 1) {
  ifelse(
    gamma0 == 0 & gamma1 == 0,
    ifelse(is.finite(b) & b > 1.0001, "MAR-EWO", "MAR"),
    "MNAR"
  )
}

simulation_display_method <- function(method, ps.method = NULL) {
  if (method %in% c("G-est", "Initial")) return("CQR--G")
  if (method == "SIPW") return(paste0("SIPW-", simulation_ps_label(ps.method)))
  if (method == "SIPW-MNAR") return("SIPW-IV")
  if (method == "ELW-MAR") return(paste0("ELW-", simulation_ps_label(ps.method)))
  if (method == "ELW-MNAR") return("ELW-IV")
  method
}

scenario_label <- function(gamma0, gamma1, ps.method = NULL, b = 1) {
  simulation_design_label(gamma0, gamma1, b)
}

simulation_scenario_levels <- function(labels) {
  lev <- c("MAR", "MAR-EWO", "MNAR")
  lev[lev %in% labels]
}

simulation_method_levels <- function(labels) {
  lev <- c(
    "ELW-PLE", "ELW-SLE", "ELW-Kernel", "ELW-IV",
    "SIPW-PLE", "SIPW-SLE", "SIPW-Kernel", "SIPW-IV",
    "CQR--G", "Naive"
  )
  lev[lev %in% labels]
}

prepare_simulation_display_table <- function(tab) {
  tab$scenario_label <- simulation_design_label(tab$gamma0, tab$gamma1, tab$b)
  tab$method_label <- mapply(
    simulation_display_method, tab$method, tab$ps.method,
    USE.NAMES = FALSE
  )
  key <- intersect(
    c("scenario_label", "method_label", "estimand", "tau_label"),
    names(tab)
  )
  tab[!duplicated(tab[key]), , drop = FALSE]
}

write_simulation_qte_latex <- function(summary_estimates, path,
                                       summary_diagnostics = NULL,
                                       taus = c(0.10, 0.25, 0.50, 0.75)) {
  tab <- summary_estimates[
    summary_estimates$estimand == "qte" & summary_estimates$tau %in% taus,
    , drop = FALSE
  ]
  tab$tau_label <- sprintf("%.2f", as.numeric(tab$tau))
  tab <- prepare_simulation_display_table(tab)
  tab <- add_simulation_failure_info(tab, summary_diagnostics)
  methods <- simulation_method_levels(unique(tab$method_label))
  scenarios <- simulation_scenario_levels(unique(tab$scenario_label))
  lines <- c(
    "\\resizebox{\\linewidth}{!}{%",
    "\\begin{tabular}{ll*{12}{r}}",
    "\\toprule",
    paste(
      "\\multirow{2}{*}{\\textbf{Scenario}} & \\multirow{2}{*}{\\textbf{Method}}",
      "& \\multicolumn{3}{c}{\\textbf{$\\tau=0.10$}}",
      "& \\multicolumn{3}{c}{\\textbf{$\\tau=0.25$}}",
      "& \\multicolumn{3}{c}{\\textbf{$\\tau=0.50$}}",
      "& \\multicolumn{3}{c}{\\textbf{$\\tau=0.75$}}",
      "\\\\"
    ),
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}\\cmidrule(lr){12-14}",
    paste(
      " &",
      "& \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE}",
      "& \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE}",
      "& \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE}",
      "& \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE}",
      "\\\\"
    ),
    "\\midrule"
  )
  for (sc in scenarios) {
    first <- TRUE
    for (m in methods) {
      rows <- tab[tab$scenario_label == sc & tab$method_label == m, , drop = FALSE]
      if (nrow(rows) == 0L) next
      eligible <- !simulation_row_is_gray(rows)
      vals <- unlist(lapply(sprintf("%.2f", taus), function(tt) {
        rr <- rows[rows$tau_label == tt, , drop = FALSE]
        sc.tt <- simulation_benchmark_rows(
          tab[tab$scenario_label == sc & tab$tau_label == tt, , drop = FALSE]
        )
        best.bias <- simulation_best(abs(sc.tt$bias.mean), "min")
        best.sd <- simulation_best(sc.tt$bias.sd, "min")
        best.rmse <- simulation_best(sc.tt$rmse, "min")
        if (nrow(rr) == 0L) {
          c("--", "--", "--")
        } else {
          c(
            if (eligible) format_sim_best(abs(rr$bias.mean[1]), best.bias, 1) else
              format_sim_number(abs(rr$bias.mean[1]), 1),
            if (eligible) format_sim_best(rr$bias.sd[1], best.sd, 1) else
              format_sim_number(rr$bias.sd[1], 1),
            if (eligible) format_sim_best(rr$rmse[1], best.rmse, 1) else
              format_sim_number(rr$rmse[1], 1)
          )
        }
      }), use.names = FALSE)
      lines <- c(lines, sprintf(
        paste(
          "%s%s & %s",
          "& %s & %s & %s",
          "& %s & %s & %s",
          "& %s & %s & %s",
          "& %s & %s & %s",
          "\\\\"
        ),
        simulation_row_prefix(rows),
        if (first) sc else "", latex_escape(m),
        vals[1], vals[2], vals[3],
        vals[4], vals[5], vals[6],
        vals[7], vals[8], vals[9],
        vals[10], vals[11], vals[12]
      ))
      first <- FALSE
    }
    lines <- c(lines, "\\addlinespace[0.25em]")
  }
  lines <- c(
    lines[-length(lines)],
    "\\bottomrule",
    "\\end{tabular}%",
    "}"
  )
  writeLines(lines, path)
}

write_simulation_lbw_latex <- function(summary_estimates, path,
                                       summary_diagnostics = NULL) {
  tab <- summary_estimates[
    summary_estimates$estimand %in% c("lbw1", "lbw0", "lbw_rd"),
    , drop = FALSE
  ]
  tab <- prepare_simulation_display_table(tab)
  tab <- add_simulation_failure_info(tab, summary_diagnostics)
  methods <- simulation_method_levels(unique(tab$method_label))
  scenarios <- simulation_scenario_levels(unique(tab$scenario_label))
  get_row <- function(rows, estimand) {
    rows[rows$estimand == estimand, , drop = FALSE]
  }
  get_sd <- function(row) {
    if (nrow(row) == 0L || !is.finite(row$bias.sd[1])) NA_real_ else row$bias.sd[1]
  }
  lines <- c(
    "\\begin{tabularx}{\\linewidth}{ll*{9}{Y}}",
    "\\toprule",
    paste(
      "\\multirow{2}{*}{\\textbf{Scenario}} & \\multirow{2}{*}{\\textbf{Method}} & \\multicolumn{3}{c}{\\textbf{$\\Pr\\{Y(1)<2{,}500\\}$}} & \\multicolumn{3}{c}{\\textbf{$\\Pr\\{Y(0)<2{,}500\\}$}} & \\multicolumn{3}{c}{\\textbf{Risk difference}}",
      "\\\\"
    ),
    "\\cmidrule(lr){3-5}\\cmidrule(lr){6-8}\\cmidrule(lr){9-11}",
    paste(
      " & & \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE} & \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE} & \\textbf{AB} & \\textbf{SD} & \\textbf{RMSE}",
      "\\\\"
    ),
    "\\midrule"
  )
  for (sc in scenarios) {
    first <- TRUE
    sc.rows <- tab[tab$scenario_label == sc, , drop = FALSE]
    sc.adj <- simulation_benchmark_rows(sc.rows)
    best.bias1 <- simulation_best(abs(sc.adj$bias.mean[sc.adj$estimand == "lbw1"]), "min")
    best.bias0 <- simulation_best(abs(sc.adj$bias.mean[sc.adj$estimand == "lbw0"]), "min")
    best.bias <- simulation_best(abs(sc.adj$bias.mean[sc.adj$estimand == "lbw_rd"]), "min")
    best.sd1 <- simulation_best(sc.adj$bias.sd[sc.adj$estimand == "lbw1"], "min")
    best.sd0 <- simulation_best(sc.adj$bias.sd[sc.adj$estimand == "lbw0"], "min")
    best.sd <- simulation_best(sc.adj$bias.sd[sc.adj$estimand == "lbw_rd"], "min")
    best.rmse1 <- simulation_best(sc.adj$rmse[sc.adj$estimand == "lbw1"], "min")
    best.rmse0 <- simulation_best(sc.adj$rmse[sc.adj$estimand == "lbw0"], "min")
    best.rmse <- simulation_best(sc.adj$rmse[sc.adj$estimand == "lbw_rd"], "min")
    for (m in methods) {
      rows <- tab[tab$scenario_label == sc & tab$method_label == m, , drop = FALSE]
      if (nrow(rows) == 0L) next
      eligible <- !simulation_row_is_gray(rows)
      row1 <- get_row(rows, "lbw1")
      row0 <- get_row(rows, "lbw0")
      rowd <- get_row(rows, "lbw_rd")
      lines <- c(lines, sprintf(
        paste(
          "%s%s & %s",
          "& %s & %s & %s",
          "& %s & %s & %s",
          "& %s & %s & %s",
          "\\\\"
        ),
        simulation_row_prefix(rows),
        if (first) sc else "", latex_escape(m),
        if (nrow(row1) == 0L) "--" else if (eligible) {
          format_sim_best(abs(row1$bias.mean[1]), best.bias1, 3)
        } else format_sim_number(abs(row1$bias.mean[1]), 3),
        if (eligible) format_sim_best(get_sd(row1), best.sd1, 3) else
          format_sim_number(get_sd(row1), 3),
        if (nrow(row1) == 0L) "--" else if (eligible) {
          format_sim_best(row1$rmse[1], best.rmse1, 3)
        } else format_sim_number(row1$rmse[1], 3),
        if (nrow(row0) == 0L) "--" else if (eligible) {
          format_sim_best(abs(row0$bias.mean[1]), best.bias0, 3)
        } else format_sim_number(abs(row0$bias.mean[1]), 3),
        if (eligible) format_sim_best(get_sd(row0), best.sd0, 3) else
          format_sim_number(get_sd(row0), 3),
        if (nrow(row0) == 0L) "--" else if (eligible) {
          format_sim_best(row0$rmse[1], best.rmse0, 3)
        } else format_sim_number(row0$rmse[1], 3),
        if (nrow(rowd) == 0L) "--" else if (eligible) {
          format_sim_best(abs(rowd$bias.mean[1]), best.bias, 3)
        } else format_sim_number(abs(rowd$bias.mean[1]), 3),
        if (eligible) format_sim_best(get_sd(rowd), best.sd, 3) else
          format_sim_number(get_sd(rowd), 3),
        if (nrow(rowd) == 0L) "--" else if (eligible) {
          format_sim_best(rowd$rmse[1], best.rmse, 3)
        } else format_sim_number(rowd$rmse[1], 3)
      ))
      first <- FALSE
    }
    lines <- c(lines, "\\addlinespace[0.25em]")
  }
  lines <- c(
    lines[-length(lines)], "\\bottomrule", "\\end{tabularx}"
  )
  writeLines(lines, path)
}

write_simulation_diagnostics_latex <- function(summary_diagnostics, path) {
  tab <- summary_diagnostics
  tab <- prepare_simulation_display_table(tab)
  methods <- simulation_method_levels(unique(tab$method_label))
  scenarios <- simulation_scenario_levels(unique(tab$scenario_label))
  lines <- c(
    "\\begin{tabularx}{\\linewidth}{ll*{5}{Y}}",
    "\\toprule",
    paste("\\textbf{Scenario} & \\textbf{Method} & \\textbf{max $|$SMD$|$} & \\textbf{mean $|$SMD$|$} & \\textbf{min $\\hat\\rho$} & \\textbf{ISE} & \\textbf{NA\\%}", "\\\\"),
    "\\midrule"
  )
  for (sc in scenarios) {
    first <- TRUE
    sc.rows <- tab[tab$scenario_label == sc, , drop = FALSE]
    sc.adj <- simulation_benchmark_rows(sc.rows)
    best.max.smd <- simulation_best(sc.adj$max.abs.smd, "min")
    best.mean.smd <- simulation_best(sc.adj$mean.abs.smd, "min")
    best.ise <- simulation_best(sc.adj$ise, "min")
    for (m in methods) {
      row <- tab[tab$scenario_label == sc & tab$method_label == m, , drop = FALSE]
      if (nrow(row) == 0L) next
      eligible <- !simulation_row_is_gray(row)
      ps.value <- if (row$method[1] %in% c("SIPW", "SIPW-MNAR", "ELW-MAR", "ELW-MNAR")) {
        format_sim_number(row$ps.min[1], 4)
      } else {
        "--"
      }
      lines <- c(lines, sprintf(
        paste("%s%s & %s & %s & %s & %s & %s & %s", "\\\\"),
        simulation_row_prefix(row),
        if (first) sc else "", latex_escape(m),
        if (eligible) format_sim_best(row$max.abs.smd[1], best.max.smd, 3) else
          format_sim_number(row$max.abs.smd[1], 3),
        if (eligible) format_sim_best(row$mean.abs.smd[1], best.mean.smd, 3) else
          format_sim_number(row$mean.abs.smd[1], 3),
        ps.value,
        if (eligible) format_sim_best(row$ise[1], best.ise, 4) else
          format_sim_number(row$ise[1], 4),
        format_sim_percent(row$failure.rate[1], 1)
      ))
      first <- FALSE
    }
    lines <- c(lines, "\\addlinespace[0.25em]")
  }
  lines <- c(
    lines[-length(lines)],
    "\\bottomrule",
    "\\end{tabularx}"
  )
  writeLines(lines, path)
}

write_simulation_latex_outputs <- function(summary, save_dir, manuscript_dir = NULL) {
  out_dirs <- save_dir
  if (!is.null(manuscript_dir)) {
    manuscript_dir <- if (basename(manuscript_dir) == "FIGURES") {
      manuscript_dir
    } else {
      file.path(manuscript_dir, "FIGURES")
    }
    dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
    out_dirs <- unique(c(out_dirs, manuscript_dir))
  }
  for (d in out_dirs) {
    write_simulation_qte_latex(
      summary$estimates, file.path(d, "csl_simulation_qte_table.tex"),
      summary$diagnostics
    )
    write_simulation_lbw_latex(
      summary$estimates, file.path(d, "csl_simulation_lbw_table.tex"),
      summary$diagnostics
    )
    write_simulation_diagnostics_latex(
      summary$diagnostics, file.path(d, "csl_simulation_diagnostics_table.tex")
    )
  }
  invisible(out_dirs)
}

## ============================================================
## CSL-calibrated simulation functions
## Generate DGP data, run competing estimators, and summarize Monte Carlo output.
## ============================================================

simulation_validate_probability <- function(p, name = "probability",
                                            allow.boundary = FALSE) {
  p <- as.numeric(p)
  if (any(!is.finite(p))) {
    stop(name, " contains nonfinite values.")
  }
  invalid <- if (allow.boundary) p < 0 | p > 1 else p <= 0 | p >= 1
  if (any(invalid)) {
    qualifier <- if (allow.boundary) "outside [0, 1]" else "at or outside (0, 1)"
    stop(name, " contains values ", qualifier, ".")
  }
  p
}

simulation_bind_rows <- function(x) {
  cols <- unique(unlist(lapply(x, names), use.names = FALSE))
  x <- lapply(x, function(d) {
    missing <- setdiff(cols, names(d))
    for (m in missing) d[[m]] <- NA
    d[cols]
  })
  do.call(rbind, x)
}

simulation_auc_rank <- function(y, p) {
  n1 <- sum(y == 1)
  n0 <- sum(y == 0)
  if (n1 == 0 || n0 == 0) return(NA_real_)
  (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) /
    (n1 * n0)
}

simulation_standardize <- function(x) {
  x <- as.matrix(x)
  means <- colMeans(x)
  sds <- apply(x, 2, stats::sd)
  keep <- is.finite(sds) & sds > 0
  x <- x[, keep, drop = FALSE]
  sweep(sweep(x, 2, means[keep], "-"), 2, sds[keep], "/")
}

simulation_reference_last_matrix <- function(xdat) {
  positive_level <- c(
    Parity = "1", Hxpreterm = "2", Diabetes = "2", Hxasthma = "2",
    HIV = "2", Alcohol = "2", Sex = "2"
  )
  blocks <- lapply(names(xdat), function(v) {
    x <- xdat[[v]]
    if (is.factor(x) || is.character(x)) {
      x <- factor(x)
      lev <- levels(x)
      if (length(lev) <= 1L) return(NULL)
      if (length(lev) == 2L && v %in% names(positive_level)) {
        out <- matrix(as.numeric(as.character(x) == positive_level[[v]]), ncol = 1)
        colnames(out) <- make.names(v)
        return(out)
      }
      out <- sapply(lev[-length(lev)], function(a) as.numeric(x == a))
      out <- as.matrix(out)
      colnames(out) <- make.names(paste0(v, lev[-length(lev)]))
      out
    } else {
      out <- matrix(as.numeric(x), ncol = 1)
      colnames(out) <- make.names(v)
      out
    }
  })
  blocks <- blocks[!vapply(blocks, is.null, logical(1))]
  out <- do.call(cbind, blocks)
  storage.mode(out) <- "double"
  out
}

make_csl_simulation_base <- function(dat,
                                     treatment.var = "treatment_smoking",
                                     outcome.var = "outcome_birthweight",
                                     exclude.vars = c(
                                       "row_id", "Smoke", "Birthweight",
                                       "treatment_smoking", "outcome_birthweight"
                                     ),
                                     factor.vars = c(
                                       "Momrace", "Marital", "Education", "Parity",
                                       "Hxpreterm", "Hxstillbirth", "Diabetes",
                                       "Hxasthma", "HIV", "Alcohol", "Preeclampsia",
                                       "AnteGDM", "AnteSGA", "AnteLGA", "Sex"
                                     )) {
  A <- as.integer(dat[[treatment.var]])
  Y <- as.numeric(dat[[outcome.var]])
  ## Prespecified pretreatment adjustment set used in the CSL application.
  covariates <- intersect(c(
    "Momage", "Height", "PreBMI", "Momrace", "Marital", "Parity",
    "Hxpreterm", "Diabetes", "Hxasthma", "HIV", "Alcohol", "Sex"
  ), names(dat))
  xdat <- dat[covariates]
  factor.vars <- intersect(factor.vars, names(xdat))
  xdat[factor.vars] <- lapply(xdat[factor.vars], factor)
  varying <- vapply(xdat, function(x) length(unique(x)) > 1L, logical(1))
  xdat <- xdat[varying]

  X <- simulation_reference_last_matrix(xdat)
  X.std <- simulation_standardize(X)

  continuous <- intersect(c("Momage", "Height", "PreBMI"), colnames(X.std))
  X.expanded <- X.std
  for (v in continuous) {
    X.expanded <- cbind(X.expanded, X.std[, v]^2)
    colnames(X.expanded)[ncol(X.expanded)] <- paste0(v, "_sq")
  }
  if (all(c("Momage", "PreBMI") %in% continuous)) {
    X.expanded <- cbind(X.expanded, X.std[, "Momage"] * X.std[, "PreBMI"])
    colnames(X.expanded)[ncol(X.expanded)] <- "Momage_x_PreBMI"
  }
  X.expanded <- X.expanded[, !duplicated(colnames(X.expanded)), drop = FALSE]
  X.ps <- as.data.frame(X.expanded, check.names = TRUE)
  ## The DGP propensity calibration is a fixed full-main-effect PLE. No
  ## data-driven selection is performed either here or in a replication.
  X.main <- as.data.frame(X.std, check.names = TRUE)
  ps.full.formula <- stats::as.formula(paste("A ~", paste(names(X.main), collapse = " + ")))
  ps.full.fit <- suppressWarnings(stats::glm(
    ps.full.formula, data = data.frame(A = A, X.main), family = stats::binomial()
  ))
  ps.selected.terms <- attr(stats::terms(ps.full.fit), "term.labels")
  ell <- as.vector(stats::predict(ps.full.fit, type = "link"))

  outcome.data <- data.frame(Y = Y[A == 0], X.ps[A == 0, , drop = FALSE])
  outcome.full.formula <- stats::as.formula(
    paste("Y ~", paste(names(X.ps), collapse = " + "))
  )
  outcome.full.fit <- suppressWarnings(stats::lm(
    outcome.full.formula, data = outcome.data
  ))
  outcome.selected.fit <- tryCatch(
    suppressWarnings(MASS::stepAIC(
      outcome.full.fit,
      scope = list(
        lower = stats::as.formula("Y ~ 1"),
        upper = outcome.full.formula
      ),
      direction = "both", trace = FALSE
    )),
    error = function(e) outcome.full.fit
  )
  outcome.selected.terms <- attr(stats::terms(outcome.selected.fit), "term.labels")
  m0 <- as.numeric(stats::predict(
    outcome.selected.fit,
    newdata = data.frame(X.ps, check.names = TRUE)
  ))
  resid0 <- as.numeric(Y[A == 0] - m0[A == 0])
  sigma0 <- stats::sd(resid0, na.rm = TRUE)
  eps.pool <- (resid0 - mean(resid0, na.rm = TRUE)) / sigma0
  eps.pool <- eps.pool[is.finite(eps.pool)]

  lbw.fit <- suppressWarnings(stats::glm.fit(
    cbind("(Intercept)" = 1, X.std), as.integer(Y < 2500), family = stats::binomial()
  ))
  lbw.coef <- lbw.fit$coefficients
  lbw.coef[!is.finite(lbw.coef)] <- 0
  lbw.risk <- stats::plogis(as.vector(cbind("(Intercept)" = 1, X.std) %*% lbw.coef))
  risk.percentile <- rank(lbw.risk, ties.method = "average") / length(lbw.risk)

  list(
    dat = dat,
    A = A,
    Y = Y,
    X = X.std,
    X.expanded = X.expanded,
    ell = ell,
    m0 = m0,
    sigma0 = sigma0,
    eps.pool = eps.pool,
    lbw.risk = lbw.risk,
    risk.percentile = risk.percentile,
    smoking.prevalence = mean(A == 1),
    ps.selected.terms = ps.selected.terms,
    ps.selected.formula = stats::formula(ps.full.fit),
    outcome.selected.terms = outcome.selected.terms,
    outcome.selected.formula = stats::formula(outcome.selected.fit),
    outcome.selected.coefficients = stats::coef(outcome.selected.fit),
    covariates = colnames(X.std),
    Z = if ("Sex" %in% names(xdat)) factor(xdat$Sex) else factor(rep(1, length(A)))
  )
}

calibrate_logit_intercept <- function(lp, target) {
  target <- min(max(target, 1e-5), 1 - 1e-5)
  f <- function(a) mean(stats::plogis(a + lp)) - target
  stats::uniroot(f, interval = c(-40, 40), tol = 1e-10)$root
}

generate_csl_simulation_data <- function(base, n = nrow(base$dat),
                                         b = 1,
                                         gamma0 = 0,
                                         gamma1 = 0,
                                         effect.low = -150,
                                         effect.high = -300,
                                         eta.sd = 50,
                                         instrument.effect = 300,
                                         heteroskedastic = FALSE,
                                         outcome.nonlinearity = 0,
                                         effect.nonlinearity = 0,
                                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  is.mnar <- gamma0 != 0 || gamma1 != 0
  id <- sample.int(nrow(base$dat), n, replace = TRUE)
  ell <- base$ell[id]
  X <- base$X[id, , drop = FALSE]
  X.expanded <- base$X.expanded[id, , drop = FALSE]
  m0 <- base$m0[id]
  r <- base$risk.percentile[id]
  Z <- if (is.mnar) {
    factor(stats::rbinom(n, 1, 0.5), levels = 0:1)
  } else {
    factor(base$Z[id], levels = levels(base$Z))
  }

  eps <- sample(base$eps.pool, n, replace = TRUE)
  if (heteroskedastic) {
    sigma <- base$sigma0 * (0.75 + 0.50 * r)
  } else {
    sigma <- rep(base$sigma0, n)
  }
  Y0 <- m0 + sigma * eps
  ## Outcome-model stress component deliberately omitted from the expanded
  ## linear CQR--G basis. Centering preserves the approximate marginal scale.
  nonlinear.shape <- (r - 0.5)^2 + 0.75 * as.numeric(r > 0.80)
  nonlinear.shape <- nonlinear.shape - mean(nonlinear.shape)
  Y0 <- Y0 + outcome.nonlinearity * nonlinear.shape
  if (is.mnar) {
    Y0 <- Y0 + instrument.effect * (as.numeric(Z) - 1.5)
  }
  delta <- effect.low + (effect.high - effect.low) * r
  ## A negative value adds a sharp adverse effect in the highest baseline-risk
  ## quintile. This threshold term is omitted from the linear CQR--G basis.
  effect.stress <- as.numeric(r > 0.80)
  Y1 <- Y0 + delta + effect.nonlinearity * effect.stress +
    stats::rnorm(n, mean = 0, sd = eta.sd)

  z0 <- as.numeric(scale(Y0))
  z1 <- as.numeric(scale(Y1))
  z0[!is.finite(z0)] <- 0
  z1[!is.finite(z1)] <- 0
  treatment.lp <- b * ell
  if (isTRUE(all.equal(b, 1)) && gamma0 == 0 && gamma1 == 0) {
    a <- 0
  } else {
    a <- calibrate_logit_intercept(
      treatment.lp + gamma1 * z1 + gamma0 * z0,
      base$smoking.prevalence
    )
  }
  lp <- a + treatment.lp + gamma1 * z1 + gamma0 * z0
  pi.true <- stats::plogis(lp)
  T <- stats::rbinom(n, 1, pi.true)
  Y <- ifelse(T == 1, Y1, Y0)

  data.frame(
    id = id,
    Z = Z,
    T = T,
    Y = Y,
    Y0 = Y0,
    Y1 = Y1,
    pi.true = pi.true,
    rho1.true = pi.true,
    rho0.true = 1 - pi.true,
    ell = ell,
    risk.percentile = r,
    delta = delta,
    check.names = FALSE
  ) |>
    cbind(as.data.frame(X, check.names = FALSE)) |>
    structure(X.expanded = X.expanded, X.main = X,
              ps.selected.terms = base$ps.selected.terms,
              intercept.a = a, b = b, gamma0 = gamma0, gamma1 = gamma1,
              outcome.nonlinearity = outcome.nonlinearity,
              effect.nonlinearity = effect.nonlinearity)
}

fit_simulation_propensity <- function(sim.dat, use.expanded = TRUE) {
  X <- if (!is.null(attr(sim.dat, "X.main"))) {
    attr(sim.dat, "X.main")
  } else {
    as.matrix(sim.dat[, setdiff(names(sim.dat), c(
      "id", "Z", "T", "Y", "Y0", "Y1", "pi.true", "rho1.true", "rho0.true",
      "ell", "risk.percentile", "delta"
    )), drop = FALSE])
  }
  X.ps <- as.data.frame(X, check.names = TRUE)
  form <- stats::as.formula(paste("T ~", paste(names(X.ps), collapse = " + ")))
  fit <- suppressWarnings(stats::glm(form, data = data.frame(T = sim.dat$T, X.ps),
                                     family = stats::binomial()))
  p <- as.numeric(stats::predict(fit, type = "response"))
  eta <- as.numeric(stats::predict(fit, type = "link"))
  list(
    probability = p,
    linear_predictor = eta,
    coefficient = stats::coef(fit),
    auc = simulation_auc_rank(sim.dat$T, p)
  )
}


fit_simulation_propensity_sle <- function(sim.dat) {
  if (is.null(attr(sim.dat, "X.expanded"))) {
    stop("SLE propensity estimation requires expanded covariates.")
  }
  X <- as.data.frame(attr(sim.dat, "X.expanded"), check.names = TRUE)
  form <- stats::as.formula(paste("T ~", paste(names(X), collapse = " + ")))
  fit <- tryCatch(
    suppressWarnings(stats::glm(form, data = data.frame(T = sim.dat$T, X),
                                family = stats::binomial(),
                                control = stats::glm.control(maxit = 50))),
    error = function(e) NULL
  )
  if (is.null(fit)) {
    return(fit_simulation_propensity(sim.dat))
  }
  p <- as.numeric(stats::predict(fit, type = "response"))
  eta <- as.numeric(stats::predict(fit, type = "link"))
  list(
    probability = p,
    linear_predictor = eta,
    coefficient = stats::coef(fit),
    auc = simulation_auc_rank(sim.dat$T, p),
    min.prob = min(p, na.rm = TRUE),
    p01.prob = as.numeric(stats::quantile(p, 0.01, names = FALSE, na.rm = TRUE)),
    p99.prob = as.numeric(stats::quantile(p, 0.99, names = FALSE, na.rm = TRUE)),
    max.prob = max(p, na.rm = TRUE)
  )
}



fit_simulation_propensity_kernel <- function(sim.dat, pilot.fit = NULL,
                                             bandwidth = NULL) {
  X <- attr(sim.dat, "X.main")
  if (is.null(X)) stop("PCA--Kernel requires the fixed main-effect design matrix.")
  fit <- pca_kernel_predict(
    X_train = X, treatment_train = sim.dat$T, X_test = X,
    n_components = 3L, grid_bins = 35L
  )
  p <- fit$probability
  list(
    probability = p,
    linear_predictor = rep(NA_real_, length(p)),
    coefficient = NULL,
    bandwidth = fit$bandwidth,
    explained = fit$explained,
    auc = simulation_auc_rank(sim.dat$T, p),
    min.prob = min(p, na.rm = TRUE),
    p01.prob = as.numeric(stats::quantile(p, 0.01, names = FALSE, na.rm = TRUE)),
    p99.prob = as.numeric(stats::quantile(p, 0.99, names = FALSE, na.rm = TRUE)),
    max.prob = max(p, na.rm = TRUE)
  )
}

simulation_true_estimands <- function(sim.dat, probs = c(0.10, 0.25, 0.50, 0.75),
                                      threshold = 2500,
                                      grid = NULL) {
  if (is.null(grid)) {
    grid <- seq(
      stats::quantile(c(sim.dat$Y0, sim.dat$Y1), 0.01, names = FALSE),
      stats::quantile(c(sim.dat$Y0, sim.dat$Y1), 0.99, names = FALSE),
      length.out = 101
    )
  }
  q0 <- as.numeric(stats::quantile(sim.dat$Y0, probs, names = FALSE, type = 1))
  q1 <- as.numeric(stats::quantile(sim.dat$Y1, probs, names = FALSE, type = 1))
  list(
    probs = probs,
    grid = grid,
    q0 = q0,
    q1 = q1,
    qte = q1 - q0,
    lbw0 = mean(sim.dat$Y0 < threshold),
    lbw1 = mean(sim.dat$Y1 < threshold),
    lbw.rd = mean(sim.dat$Y1 < threshold) - mean(sim.dat$Y0 < threshold),
    F0 = vapply(grid, function(y) mean(sim.dat$Y0 <= y), numeric(1)),
    F1 = vapply(grid, function(y) mean(sim.dat$Y1 <= y), numeric(1))
  )
}

estimate_weighted_distribution <- function(y, weight, probs, threshold, grid) {
  weight <- as.numeric(weight)
  weight <- weight / sum(weight)
  list(
    quantile = weighted_quantile(y, weight, probs),
    lbw = weighted_probability(y, weight, function(x) x < threshold),
    F = vapply(grid, function(g) sum(weight * (y <= g)), numeric(1)),
    ess = effective_sample_size(weight),
    max.weight = max(weight),
    p99.weight = as.numeric(stats::quantile(weight, 0.99, names = FALSE))
  )
}


gaussian_kernel <- function(u) {
  out <- numeric(length(u))
  keep <- is.finite(u)
  out[keep] <- exp(-0.5 * u[keep]^2)
  out
}

g_est_average_nw_weights <- function(u.obs, u.target, bandwidth = NULL,
                                     obs.mult = NULL, target.mult = NULL,
                                     bandwidth.constant = 1.06,
                                     chunk.size = 500L) {
  u.obs <- as.numeric(u.obs)
  u.target <- as.numeric(u.target)
  keep <- is.finite(u.obs)
  u.obs <- u.obs[keep]
  if (length(u.obs) < 5L) {
    stop("Single-index g-estimation requires at least five observations per group.")
  }
  if (is.null(obs.mult)) obs.mult <- rep(1, length(keep))
  obs.mult <- as.numeric(obs.mult)[keep]
  if (is.null(target.mult)) target.mult <- rep(1, length(u.target))
  target.mult <- as.numeric(target.mult)
  obs.mult[!is.finite(obs.mult) | obs.mult < 0] <- 0
  target.mult[!is.finite(target.mult) | target.mult < 0] <- 0
  if (sum(obs.mult) <= 0 || sum(target.mult) <= 0) {
    stop("Degenerate multiplier weights in single-index g-estimation.")
  }
  if (is.null(bandwidth)) {
    bandwidth <- bandwidth.constant * stats::sd(u.obs) *
      length(u.obs)^(-1 / 5)
  }
  if (!is.finite(bandwidth) || bandwidth <= 0) bandwidth <- 0.25

  acc <- numeric(length(u.obs))
  target.mass <- target.mult / sum(target.mult)
  starts <- seq.int(1L, length(u.target), by = chunk.size)
  for (st in starts) {
    en <- min(length(u.target), st + chunk.size - 1L)
    ut <- u.target[st:en]
    D <- outer(ut, u.obs, function(a, b) (b - a) / bandwidth)
    K <- matrix(gaussian_kernel(D), nrow = length(ut))
    weighted.K <- sweep(K, 2, obs.mult, "*")
    denom <- rowSums(weighted.K)
    bad <- !is.finite(denom) | denom <= 1e-12
    if (any(bad)) {
      distance <- abs(D[bad, , drop = FALSE])
      distance[, obs.mult <= 0] <- Inf
      nearest <- max.col(-distance, ties.method = "first")
      weighted.K[bad, ] <- 0
      weighted.K[cbind(which(bad), nearest)] <- obs.mult[nearest]
      denom[bad] <- obs.mult[nearest]
    }
    W <- weighted.K / denom
    acc <- acc + colSums(sweep(W, 1, target.mass[st:en], "*"))
  }
  out <- numeric(length(keep))
  out[keep] <- acc
  out <- out / sum(out)
  out
}

g_est_side_distribution <- function(y.obs, u.obs, u.target, probs, threshold, grid,
                                    bandwidth = NULL, bandwidth.constant = 1.06) {
  y.obs <- as.numeric(y.obs)
  w <- g_est_average_nw_weights(
    u.obs, u.target, bandwidth = bandwidth,
    bandwidth.constant = bandwidth.constant
  )
  dist <- estimate_weighted_distribution(y.obs, w, probs, threshold, grid)
  list(
    quantile = dist$quantile,
    lbw = dist$lbw,
    F = dist$F,
    ess = dist$ess,
    max.weight = dist$max.weight,
    p99.weight = dist$p99.weight,
    g.est.violations = 0,
    g.est.range.error = 0,
    bandwidth = ifelse(is.null(bandwidth), NA_real_, bandwidth)
  )
}

estimate_g_est_simulation <- function(sim.dat, index, truth, threshold = 2500,
                                      method = "G-est") {
  if (length(index) != nrow(sim.dat)) {
    stop("Single-index g-estimation requires the treatment-assignment index.")
  }
  U <- as.numeric(scale(index))
  U[!is.finite(U)] <- 0
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  side1 <- g_est_side_distribution(
    sim.dat$Y[idx1], U[idx1], U, truth$probs, threshold, truth$grid
  )
  side0 <- g_est_side_distribution(
    sim.dat$Y[idx0], U[idx0], U, truth$probs, threshold, truth$grid
  )
  out <- make_simulation_result_rows(method, side1, side0, truth)
  out$diagnostics$g.est.violations1 <- side1$g.est.violations
  out$diagnostics$g.est.violations0 <- side0$g.est.violations
  out$diagnostics$g.est.range.error1 <- side1$g.est.range.error
  out$diagnostics$g.est.range.error0 <- side0$g.est.range.error
  out
}


initial_qr_side_distribution <- function(y.obs, x.obs, x.target, probs,
                                         threshold, grid,
                                         tau.grid = seq(0.01, 0.99, by = 0.01)) {
  if (!requireNamespace("quantreg", quietly = TRUE)) {
    stop("The G-est estimator requires the quantreg package.")
  }
  y.obs <- as.numeric(y.obs)
  x.obs <- as.matrix(x.obs)
  x.target <- as.matrix(x.target)
  if (nrow(x.obs) != length(y.obs) || ncol(x.obs) != ncol(x.target)) {
    stop("Incompatible outcome and covariate dimensions for G-est.")
  }
  tau.grid <- sort(unique(as.numeric(tau.grid)))
  if (length(tau.grid) < 3L || tau.grid[1] <= 0 ||
      tau.grid[length(tau.grid)] >= 1) {
    stop("G-est requires an interior grid of at least three quantile levels.")
  }
  step <- diff(tau.grid)
  if (max(abs(step - step[1])) > 1e-10) {
    stop("G-est currently requires an equally spaced quantile grid.")
  }
  epsilon <- tau.grid[1]
  if (abs(tau.grid[length(tau.grid)] - (1 - epsilon)) > 1e-10) {
    stop("The G-est quantile grid must run from epsilon to 1 - epsilon.")
  }

  x.obs <- cbind("(Intercept)" = 1, x.obs)
  x.target <- cbind("(Intercept)" = 1, x.target)
  svd.x <- svd(x.obs, nu = 0)
  keep <- svd.x$d > max(svd.x$d) * sqrt(.Machine$double.eps)
  transform <- sweep(
    svd.x$v[, keep, drop = FALSE], 2,
    svd.x$d[keep] / sqrt(nrow(x.obs)), "/"
  )
  x.obs <- x.obs %*% transform
  x.target <- x.target %*% transform

  coefficients <- suppressWarnings(vapply(
    tau.grid,
    function(tau) quantreg::rq.fit.fnb(
      x.obs, y.obs, tau = tau
    )$coefficients,
    numeric(ncol(x.obs))
  ))
  if (any(!is.finite(coefficients))) {
    stop("G-est conditional quantile regression did not converge.")
  }
  fitted.quantiles <- x.target %*% coefficients

  ## Cai et al. (2024), equation (14): use right endpoints of the
  ## equally spaced interior quantile mesh and add the trimmed lower-tail mass.
  support <- sort(as.numeric(fitted.quantiles[, -1L, drop = FALSE]))
  interior.mass <- 1 - 2 * epsilon
  cdf_at <- function(y, strict = FALSE) {
    if (strict) {
      epsilon + interior.mass * (findInterval(y, support, left.open = TRUE) /
                                   length(support))
    } else {
      epsilon + interior.mass * (findInterval(y, support) / length(support))
    }
  }
  adjusted.probs <- (probs - epsilon) / interior.mass
  if (any(adjusted.probs <= 0 | adjusted.probs >= 1)) {
    stop("Requested quantiles fall outside the trimmed G-est distribution.")
  }
  list(
    quantile = as.numeric(stats::quantile(
      support, adjusted.probs, names = FALSE, type = 1
    )),
    lbw = cdf_at(threshold, strict = TRUE),
    F = vapply(grid, cdf_at, numeric(1)),
    ess = NA_real_,
    max.weight = NA_real_,
    p99.weight = NA_real_,
    epsilon = epsilon,
    tau.count = length(tau.grid),
    rank = ncol(x.obs)
  )
}

estimate_initial_qr_simulation <- function(sim.dat, truth, threshold = 2500,
                                           method = "G-est",
                                           tau.grid = seq(0.01, 0.99, by = 0.01)) {
  X <- attr(sim.dat, "X.expanded")
  if (is.null(X)) {
    stop("G-est requires the expanded covariate design matrix.")
  }
  X <- as.matrix(X)
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  side1 <- initial_qr_side_distribution(
    sim.dat$Y[idx1], X[idx1, , drop = FALSE], X,
    truth$probs, threshold, truth$grid, tau.grid = tau.grid
  )
  side0 <- initial_qr_side_distribution(
    sim.dat$Y[idx0], X[idx0, , drop = FALSE], X,
    truth$probs, threshold, truth$grid, tau.grid = tau.grid
  )
  out <- make_simulation_result_rows(method, side1, side0, truth)
  out$diagnostics$raw.endpoint1 <- 1 - side1$epsilon
  out$diagnostics$raw.endpoint0 <- 1 - side0$epsilon
  out$diagnostics$raw.endpoint.error <- max(side1$epsilon, side0$epsilon)
  out
}


make_simulation_result_rows <- function(method, side1, side0, truth) {
  q.rows <- data.frame(
    method = method,
    estimand = "qte",
    tau = truth$probs,
    estimate = as.numeric(side1$quantile - side0$quantile),
    truth = as.numeric(truth$qte),
    q1 = as.numeric(side1$quantile),
    q0 = as.numeric(side0$quantile),
    stringsAsFactors = FALSE
  )
  lbw.rows <- rbind(
    data.frame(
      method = method,
      estimand = "lbw1",
      tau = NA_real_,
      estimate = side1$lbw,
      truth = truth$lbw1,
      q1 = side1$lbw,
      q0 = side0$lbw,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = method,
      estimand = "lbw0",
      tau = NA_real_,
      estimate = side0$lbw,
      truth = truth$lbw0,
      q1 = side1$lbw,
      q0 = side0$lbw,
      stringsAsFactors = FALSE
    ),
    data.frame(
      method = method,
      estimand = "lbw_rd",
      tau = NA_real_,
      estimate = side1$lbw - side0$lbw,
      truth = truth$lbw.rd,
      q1 = side1$lbw,
      q0 = side0$lbw,
      stringsAsFactors = FALSE
    )
  )
  ise <- mean((side1$F - truth$F1)^2) + mean((side0$F - truth$F0)^2)
  diag.row <- data.frame(
    method = method,
    ess1 = side1$ess,
    ess0 = side0$ess,
    max.weight1 = side1$max.weight,
    max.weight0 = side0$max.weight,
    p99.weight1 = side1$p99.weight,
    p99.weight0 = side0$p99.weight,
    ise = ise,
    raw.endpoint1 = 1,
    raw.endpoint0 = 1,
    raw.endpoint.error = 0,
    raw.ise = ise,
    g.est.violations1 = NA_real_,
    g.est.violations0 = NA_real_,
    g.est.range.error1 = NA_real_,
    g.est.range.error0 = NA_real_,
    max.abs.smd = NA_real_,
    mean.abs.smd = NA_real_,
    stringsAsFactors = FALSE
  )
  list(estimates = rbind(q.rows, lbw.rows), diagnostics = diag.row)
}

estimate_naive_simulation <- function(sim.dat, truth, threshold = 2500) {
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  side1 <- estimate_weighted_distribution(
    sim.dat$Y[idx1], rep(1, sum(idx1)), truth$probs, threshold, truth$grid
  )
  side0 <- estimate_weighted_distribution(
    sim.dat$Y[idx0], rep(1, sum(idx0)), truth$probs, threshold, truth$grid
  )
  out <- make_simulation_result_rows("Naive", side1, side0, truth)
  smd <- simulation_smd_summary(sim.dat)
  out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
  out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
  out
}

estimate_sipw_simulation <- function(sim.dat, ps, truth, threshold = 2500,
                                     method = "SIPW") {
  ps <- simulation_validate_probability(ps, "propensity score")
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  w1 <- mean(sim.dat$T == 1) / ps[idx1]
  w0 <- mean(sim.dat$T == 0) / (1 - ps[idx0])
  side1 <- estimate_weighted_distribution(
    sim.dat$Y[idx1], w1, truth$probs, threshold, truth$grid
  )
  side0 <- estimate_weighted_distribution(
    sim.dat$Y[idx0], w0, truth$probs, threshold, truth$grid
  )
  out <- make_simulation_result_rows(method, side1, side0, truth)
  smd <- simulation_smd_summary(sim.dat, w1, w0)
  out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
  out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
  out
}

estimate_elw_simulation <- function(sim.dat, rho, truth, threshold = 2500,
                                    method = "ELW") {
  rho <- simulation_validate_probability(
    rho, "observation probability", allow.boundary = TRUE
  )
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  fit1 <- elw_group(sim.dat$T, rho, "treated")
  fit0 <- elw_group(sim.dat$T, rho, "control")
  side1 <- estimate_weighted_distribution(
    sim.dat$Y[idx1], fit1$weight, truth$probs, threshold, truth$grid
  )
  side0 <- estimate_weighted_distribution(
    sim.dat$Y[idx0], fit0$weight, truth$probs, threshold, truth$grid
  )
  out <- make_simulation_result_rows(method, side1, side0, truth)
  smd <- simulation_smd_summary(sim.dat, fit1$weight, fit0$weight)
  out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
  out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
  out$diagnostics$root1 <- fit1$alpha
  out$diagnostics$root0 <- fit0$alpha
  out
}

estimate_shao_mnar_probabilities <- function(sim.dat, ps.fit,
                                             gamma_range1 = c(0, 0.01),
                                             gamma_range0 = c(-0.01, 0)) {
  if (is.null(ps.fit$linear_predictor)) {
    stop("Shao--Wang estimation requires a fitted single-index score.")
  }
  if (!"Z" %in% names(sim.dat)) {
    stop("Shao--Wang estimation requires an IV or grouping variable Z.")
  }
  Z <- factor(sim.dat$Z)
  if (length(levels(Z)) < 2L) {
    stop("Shao--Wang estimation requires at least two levels of Z.")
  }
  U <- as.numeric(scale(ps.fit$linear_predictor))
  U[!is.finite(U)] <- 0

  shao1 <- estimate_shao_side(
    U = U, delta = sim.dat$T, Y = sim.dat$Y, Z = Z,
    gamma_range = gamma_range1
  )
  shao0 <- estimate_shao_side(
    U = U, delta = 1 - sim.dat$T, Y = sim.dat$Y, Z = Z,
    gamma_range = gamma_range0
  )

  rho1 <- simulation_validate_probability(
    shao1$rho, "Shao--Wang rho1", allow.boundary = TRUE
  )
  rho0 <- simulation_validate_probability(
    shao0$rho, "Shao--Wang rho0", allow.boundary = TRUE
  )
  list(rho1 = rho1, rho0 = rho0, shao1 = shao1, shao0 = shao0)
}

estimate_sipw_mnar_simulation <- function(sim.dat, ps.fit, truth, threshold = 2500,
                                          method = "SIPW-MNAR",
                                          gamma_range1 = c(0, 0.01),
                                          gamma_range0 = c(-0.01, 0)) {
  shao <- estimate_shao_mnar_probabilities(
    sim.dat, ps.fit, gamma_range1 = gamma_range1, gamma_range0 = gamma_range0
  )
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  observed.rho <- c(shao$rho1[idx1], shao$rho0[idx0])
  simulation_validate_probability(
    observed.rho, "observed SIPW selection probability"
  )
  w1 <- mean(idx1) / shao$rho1[idx1]
  w0 <- mean(idx0) / shao$rho0[idx0]
  side1 <- estimate_weighted_distribution(
    sim.dat$Y[idx1], w1, truth$probs, threshold, truth$grid
  )
  side0 <- estimate_weighted_distribution(
    sim.dat$Y[idx0], w0, truth$probs, threshold, truth$grid
  )
  out <- make_simulation_result_rows(method, side1, side0, truth)
  smd <- simulation_smd_summary(sim.dat, w1, w0)
  out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
  out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
  out$diagnostics$ps.min <- min(observed.rho)
  out$diagnostics$ps.p01 <- as.numeric(stats::quantile(observed.rho, 0.01, names = FALSE))
  out$diagnostics$ps.p99 <- as.numeric(stats::quantile(observed.rho, 0.99, names = FALSE))
  out$diagnostics$ps.max <- max(observed.rho)
  out$diagnostics$shao_gamma1 <- shao$shao1$gamma
  out$diagnostics$shao_gamma0 <- shao$shao0$gamma
  out$diagnostics$shao_bandwidth1 <- shao$shao1$bandwidth
  out$diagnostics$shao_bandwidth0 <- shao$shao0$bandwidth
  out
}

estimate_elw_mnar_simulation <- function(sim.dat, ps.fit, truth, threshold = 2500,
                                         method = "ELW-MNAR",
                                         gamma_range1 = c(0, 0.01),
                                         gamma_range0 = c(-0.01, 0)) {
  shao <- estimate_shao_mnar_probabilities(
    sim.dat, ps.fit, gamma_range1 = gamma_range1, gamma_range0 = gamma_range0
  )
  idx1 <- sim.dat$T == 1
  idx0 <- sim.dat$T == 0
  fit1 <- elw_group(sim.dat$T, shao$rho1, "treated")
  fit0 <- elw_group(sim.dat$T, 1 - shao$rho0, "control")
  side1 <- estimate_weighted_distribution(
    sim.dat$Y[idx1], fit1$weight, truth$probs, threshold, truth$grid
  )
  side0 <- estimate_weighted_distribution(
    sim.dat$Y[idx0], fit0$weight, truth$probs, threshold, truth$grid
  )
  out <- make_simulation_result_rows(method, side1, side0, truth)
  smd <- simulation_smd_summary(sim.dat, fit1$weight, fit0$weight)
  out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
  out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
  out$diagnostics$root1 <- fit1$alpha
  out$diagnostics$root0 <- fit0$alpha
  observed.rho <- c(shao$rho1[idx1], shao$rho0[idx0])
  out$diagnostics$ps.min <- min(observed.rho)
  out$diagnostics$ps.p01 <- as.numeric(stats::quantile(observed.rho, 0.01, names = FALSE))
  out$diagnostics$ps.p99 <- as.numeric(stats::quantile(observed.rho, 0.99, names = FALSE))
  out$diagnostics$ps.max <- max(observed.rho)
  out$diagnostics$shao_gamma1 <- shao$shao1$gamma
  out$diagnostics$shao_gamma0 <- shao$shao0$gamma
  out$diagnostics$shao_bandwidth1 <- shao$shao1$bandwidth
  out$diagnostics$shao_bandwidth0 <- shao$shao0$bandwidth
  out
}

safe_simulation_method <- function(expr, method, truth,
                                   estimands = c("qte", "lbw1", "lbw0", "lbw_rd")) {
  tryCatch(
    expr,
    error = function(e) {
      include.qte <- "qte" %in% estimands
      include.lbw <- any(c("lbw1", "lbw0", "lbw_rd") %in% estimands)
      est <- data.frame(
        method = character(0), estimand = character(0), tau = numeric(0),
        estimate = numeric(0), truth = numeric(0), q1 = numeric(0),
        q0 = numeric(0), error = character(0), stringsAsFactors = FALSE
      )
      if (include.qte) {
        est <- rbind(est, data.frame(
          method = method,
          estimand = "qte",
          tau = truth$probs,
          estimate = NA_real_,
          truth = truth$qte,
          q1 = NA_real_,
          q0 = NA_real_,
          error = conditionMessage(e),
          stringsAsFactors = FALSE
        ))
      }
      if (include.lbw) {
        est <- rbind(
          est,
          data.frame(
            method = method,
            estimand = "lbw1",
            tau = NA_real_,
            estimate = NA_real_,
            truth = truth$lbw1,
            q1 = NA_real_,
            q0 = NA_real_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          ),
          data.frame(
            method = method,
            estimand = "lbw0",
            tau = NA_real_,
            estimate = NA_real_,
            truth = truth$lbw0,
            q1 = NA_real_,
            q0 = NA_real_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          ),
          data.frame(
            method = method,
            estimand = "lbw_rd",
            tau = NA_real_,
            estimate = NA_real_,
            truth = truth$lbw.rd,
            q1 = NA_real_,
            q0 = NA_real_,
            error = conditionMessage(e),
            stringsAsFactors = FALSE
          )
        )
      }
      list(
        estimates = est,
        diagnostics = data.frame(
          method = method,
          ess1 = NA_real_, ess0 = NA_real_,
          max.weight1 = NA_real_, max.weight0 = NA_real_,
          p99.weight1 = NA_real_, p99.weight0 = NA_real_,
          ise = NA_real_, raw.endpoint1 = NA_real_, raw.endpoint0 = NA_real_,
          raw.endpoint.error = NA_real_, raw.ise = NA_real_,
          g.est.violations1 = NA_real_, g.est.violations0 = NA_real_,
          g.est.range.error1 = NA_real_, g.est.range.error0 = NA_real_,
          max.abs.smd = NA_real_, mean.abs.smd = NA_real_,
          error = conditionMessage(e), stringsAsFactors = FALSE
        )
      )
    }
  )
}

run_one_csl_simulation <- function(base, iter, n = 5000, b = 1,
                                   gamma0 = 0, gamma1 = 0,
                                   probs = c(0.10, 0.25, 0.50, 0.75),
                                   threshold = 2500,
                                   seed.offset = 20260623,
                                   heteroskedastic = TRUE,
                                   outcome.nonlinearity = 0,
                                   effect.nonlinearity = 0,
                                   include_mnar_sipw = TRUE,
                                   include_mnar_elw = TRUE,
                                   ps.method = c("PLE", "SLE", "kernel")) {
  ps.method <- match.arg(ps.method)
  sim.dat <- generate_csl_simulation_data(
    base = base, n = n, b = b, gamma0 = gamma0, gamma1 = gamma1,
    eta.sd = if (gamma0 != 0 || gamma1 != 0) 0 else 50,
    heteroskedastic = heteroskedastic,
    outcome.nonlinearity = outcome.nonlinearity,
    effect.nonlinearity = effect.nonlinearity,
    seed = seed.offset + iter
  )
  truth <- simulation_true_estimands(sim.dat, probs = probs, threshold = threshold)
  ps.pilot <- fit_simulation_propensity(sim.dat)
  ps.fit <- switch(
    ps.method,
    PLE = ps.pilot,
    SLE = fit_simulation_propensity_sle(sim.dat),
    kernel = fit_simulation_propensity_kernel(sim.dat, pilot.fit = ps.pilot)
  )
  ps <- ps.fit$probability

  results <- list(
    safe_simulation_method(estimate_naive_simulation(sim.dat, truth, threshold),
                           "Naive", truth),
    safe_simulation_method(estimate_sipw_simulation(sim.dat, ps, truth, threshold, "SIPW"),
                           "SIPW", truth),
    safe_simulation_method(estimate_elw_simulation(sim.dat, ps, truth, threshold, "ELW-MAR"),
                           "ELW-MAR", truth)
  )
  if (ps.method == "PLE") {
    results[[length(results) + 1L]] <- safe_simulation_method(
      estimate_initial_qr_simulation(sim.dat, truth, threshold),
      "G-est", truth
    )
  }
  if (include_mnar_sipw && (gamma0 != 0 || gamma1 != 0)) {
    results[[length(results) + 1L]] <- safe_simulation_method(
      estimate_sipw_mnar_simulation(sim.dat, ps.fit, truth, threshold),
      "SIPW-MNAR", truth
    )
  }
  if (include_mnar_elw && (gamma0 != 0 || gamma1 != 0)) {
    results[[length(results) + 1L]] <- safe_simulation_method(
      estimate_elw_mnar_simulation(sim.dat, ps.fit, truth, threshold),
      "ELW-MNAR", truth
    )
  }

  estimates <- simulation_bind_rows(lapply(results, `[[`, "estimates"))
  diagnostics <- simulation_bind_rows(lapply(results, `[[`, "diagnostics"))
  estimates$iter <- iter
  diagnostics$iter <- iter
  estimates$b <- b
  estimates$gamma0 <- gamma0
  estimates$gamma1 <- gamma1
  estimates$ps.method <- ps.method
  diagnostics$b <- b
  diagnostics$gamma0 <- gamma0
  diagnostics$gamma1 <- gamma1
  diagnostics$ps.method <- ps.method
  diagnostics$ps.auc <- ps.fit$auc
  observed.ps <- c(ps[sim.dat$T == 1], 1 - ps[sim.dat$T == 0])
  default.ps.min <- min(observed.ps)
  default.ps.p01 <- as.numeric(stats::quantile(observed.ps, 0.01, names = FALSE))
  default.ps.p99 <- as.numeric(stats::quantile(observed.ps, 0.99, names = FALSE))
  default.ps.max <- max(observed.ps)
  for (v in c("ps.min", "ps.p01", "ps.p99", "ps.max")) {
    if (!v %in% names(diagnostics)) diagnostics[[v]] <- NA_real_
  }
  mar.weighting <- diagnostics$method %in% c("SIPW", "ELW-MAR")
  diagnostics$ps.min[mar.weighting & !is.finite(diagnostics$ps.min)] <- default.ps.min
  diagnostics$ps.p01[mar.weighting & !is.finite(diagnostics$ps.p01)] <- default.ps.p01
  diagnostics$ps.p99[mar.weighting & !is.finite(diagnostics$ps.p99)] <- default.ps.p99
  diagnostics$ps.max[mar.weighting & !is.finite(diagnostics$ps.max)] <- default.ps.max
  diagnostics$smoking.prevalence <- mean(sim.dat$T == 1)
  list(estimates = estimates, diagnostics = diagnostics)
}

summarize_csl_simulation <- function(estimates, diagnostics) {
  estimates$bias <- estimates$estimate - estimates$truth
  estimates$tau_label <- ifelse(
    is.na(estimates$tau), estimates$estimand, sprintf("%.2f", estimates$tau)
  )
  key <- c("b", "gamma0", "gamma1", "ps.method", "method", "estimand", "tau_label")
  est.summary <- stats::aggregate(
    cbind(estimate, truth, bias) ~ b + gamma0 + gamma1 + ps.method + method + estimand + tau_label,
    data = estimates,
    FUN = function(x) c(mean = mean(x, na.rm = TRUE), sd = stats::sd(x, na.rm = TRUE)),
    na.action = stats::na.pass
  )
  unpack <- function(x, name) {
    m <- if (is.matrix(x)) x else do.call(rbind, x)
    colnames(m) <- paste(name, colnames(m), sep = ".")
    as.data.frame(m)
  }
  est.out <- cbind(
    est.summary[key],
    unpack(est.summary$estimate, "estimate"),
    unpack(est.summary$truth, "truth"),
    unpack(est.summary$bias, "bias")
  )
  est.out$tau <- suppressWarnings(as.numeric(est.out$tau_label))
  rmse <- stats::aggregate(
    bias ~ b + gamma0 + gamma1 + ps.method + method + estimand + tau_label,
    data = estimates,
    FUN = function(x) sqrt(mean(x^2, na.rm = TRUE)),
    na.action = stats::na.pass
  )
  names(rmse)[names(rmse) == "bias"] <- "rmse"
  est.out <- merge(est.out, rmse, by = key, all.x = TRUE)

  diag.summary <- stats::aggregate(
    cbind(ess1, ess0, max.weight1, max.weight0, p99.weight1, p99.weight0, ise,
          raw.endpoint1, raw.endpoint0, raw.endpoint.error, raw.ise,
          g.est.violations1, g.est.violations0, g.est.range.error1, g.est.range.error0,
          max.abs.smd, mean.abs.smd,
          ps.auc, ps.min, ps.p01, ps.p99, ps.max, smoking.prevalence) ~
      b + gamma0 + gamma1 + ps.method + method,
    data = diagnostics,
    FUN = function(x) mean(x, na.rm = TRUE),
    na.action = stats::na.pass
  )
  failure.summary <- stats::aggregate(
    is.na(estimate) ~ b + gamma0 + gamma1 + ps.method + method,
    data = estimates,
    FUN = mean,
    na.action = stats::na.pass
  )
  names(failure.summary)[names(failure.summary) == "is.na(estimate)"] <- "failure.rate"
  diag.summary <- merge(
    diag.summary, failure.summary,
    by = c("b", "gamma0", "gamma1", "ps.method", "method"), all.x = TRUE
  )
  list(estimates = est.out, diagnostics = diag.summary)
}
## Fixed three-component PCA Gaussian-kernel propensity-score utilities.

weighted_pca_transform <- function(X_train, X_test = X_train, weights = NULL,
                                   n_components = 3L) {
  X_train <- as.matrix(X_train)
  X_test <- as.matrix(X_test)
  if (is.null(weights)) weights <- rep(1, nrow(X_train))
  weights <- as.numeric(weights)
  weights <- weights / sum(weights)
  center <- colSums(X_train * weights)
  train_centered <- sweep(X_train, 2, center, "-")
  scale <- sqrt(colSums(train_centered^2 * weights))
  keep <- is.finite(scale) & scale > 1e-10
  train_z <- sweep(train_centered[, keep, drop = FALSE], 2, scale[keep], "/")
  test_z <- sweep(
    sweep(X_test[, keep, drop = FALSE], 2, center[keep], "-"),
    2, scale[keep], "/"
  )
  covariance <- crossprod(train_z, train_z * weights)
  eig <- eigen(covariance, symmetric = TRUE)
  rotation <- eig$vectors[, seq_len(n_components), drop = FALSE]
  train_scores <- train_z %*% rotation
  test_scores <- test_z %*% rotation
  score_scale <- sqrt(colSums(train_scores^2 * weights))
  train_scores <- sweep(train_scores, 2, score_scale, "/")
  test_scores <- sweep(test_scores, 2, score_scale, "/")
  list(
    train = train_scores, test = test_scores,
    explained = eig$values / sum(eig$values),
    center = center, scale = scale, rotation = rotation,
    score_scale = score_scale
  )
}

normal_reference_bandwidth <- function(n, dimension = 3L) {
  (4 / (dimension + 2))^(1 / (dimension + 4)) * n^(-1 / (dimension + 4))
}

smooth_array_dimension <- function(x, kernel_matrix, dimension) {
  permutation <- c(dimension, setdiff(seq_along(dim(x)), dimension))
  xp <- aperm(x, permutation)
  smoothed <- kernel_matrix %*% matrix(xp, nrow = dim(xp)[1])
  yp <- array(smoothed, dim = dim(xp))
  aperm(yp, order(permutation))
}

pca_kernel_predict <- function(X_train, treatment_train, X_test = X_train,
                               weights = NULL, n_components = 3L,
                               grid_bins = 45L) {
  if (is.null(weights)) weights <- rep(1, nrow(X_train))
  pca <- weighted_pca_transform(
    X_train, X_test, weights = weights, n_components = n_components
  )
  z_train <- pca$train
  z_test <- pca$test
  dimension <- ncol(z_train)
  ## The bandwidth rule is fixed by the nominal sample size; multiplier sums
  ## fluctuate across draws but do not redefine the bootstrap sample size.
  bandwidth <- normal_reference_bandwidth(nrow(X_train), dimension)

  ## A fixed standardized-PC grid makes the three-dimensional product kernel
  ## computationally feasible for the full CSL sample and multiplier bootstrap.
  lower <- -4
  upper <- 4
  centers <- seq(lower, upper, length.out = grid_bins)
  step <- centers[2] - centers[1]
  to_bin <- function(z) {
    out <- round((pmin(upper, pmax(lower, z)) - lower) / step) + 1L
    pmin(grid_bins, pmax(1L, out))
  }
  train_bin <- apply(z_train, 2, to_bin)
  linear_bin <- train_bin[, 1] +
    grid_bins * (train_bin[, 2] - 1L) +
    grid_bins^2 * (train_bin[, 3] - 1L)
  weighted_tabulate <- function(values) {
    total <- numeric(grid_bins^3)
    aggregated <- rowsum(values, linear_bin, reorder = FALSE)
    total[as.integer(rownames(aggregated))] <- aggregated[, 1]
    total
  }
  denominator <- array(weighted_tabulate(weights), dim = rep(grid_bins, 3))
  numerator <- array(
    weighted_tabulate(weights * as.numeric(treatment_train)),
    dim = rep(grid_bins, 3)
  )
  kernel_matrix <- exp(-0.5 * (outer(centers, centers, "-") / bandwidth)^2)
  for (j in seq_len(3L)) {
    denominator <- smooth_array_dimension(denominator, kernel_matrix, j)
    numerator <- smooth_array_dimension(numerator, kernel_matrix, j)
  }
  fitted_grid <- numerator / denominator
  fallback <- stats::weighted.mean(treatment_train, weights)
  fitted_grid[!is.finite(fitted_grid)] <- fallback

  test_bin <- apply(z_test, 2, to_bin)
  test_linear <- test_bin[, 1] +
    grid_bins * (test_bin[, 2] - 1L) +
    grid_bins^2 * (test_bin[, 3] - 1L)
  probability <- as.numeric(fitted_grid[test_linear])
  list(
    probability = probability,
    bandwidth = bandwidth,
    explained = pca$explained,
    grid_bins = grid_bins
  )
}
}, envir = envir)
}

.application_load_lmz <- function(envir = parent.frame()) {
.application_load_base(envir)
evalq({
elw_group <- function(A, pi_hat, group = c("treated", "control"), tol = 1e-12, check_tol = 1e-07) {
    group <- match.arg(group)
    r1 <- mean(A == 1)
    if (group == "treated") {
        idx <- A == 1
        pw <- pi_hat[idx]
        xi <- r1 + (1 - r1) * pw
        low <- min(pw)
        up <- min(xi) - .Machine$double.eps
        root_fun <- function(alpha) sum((pw - alpha)/(xi - alpha))
        denom_fun <- function(alpha) {
            1 + (1/r1 - 1) * (pw - alpha)/(1 - alpha)
        }
    }
    else {
        idx <- A == 0
        pw <- pi_hat[idx]
        xi <- 1 - r1 * pw
        low <- min(1 - pw)
        up <- min(xi) - .Machine$double.eps
        root_fun <- function(beta) sum((1 - pw - beta)/(xi - beta))
        denom_fun <- function(beta) {
            1 + (1 - pw - beta)/((1/r1 - 1) * (1 - beta))
        }
    }
    f_low <- root_fun(low)
    f_up <- root_fun(up)
    if (!is.finite(f_low) || !is.finite(f_up) || f_low * f_up > 0) {
        warning("ELW root finding failed for ", group, ".")
        return(list(alpha = NA_real_, weight = rep(NA_real_, sum(idx)), converged = FALSE, failure_reason = "invalid root bracket", 
            normalization_error = NA_real_, calibration_error = NA_real_, boundary_gap = NA_real_, lower_bound = low, 
            upper_bound = up, f_lower = f_low, f_upper = f_up))
    }
    alpha <- stats::uniroot(root_fun, c(low, up), f.lower = f_low, f.upper = f_up, tol = tol)$root
    weight <- 1/denom_fun(alpha)
    weight <- weight/sum(weight)
    rho_observed <- if (group == "treated") 
        pw
    else 1 - pw
    normalization_error <- abs(sum(weight) - 1)
    calibration_error <- abs(sum(weight * (rho_observed - alpha)))
    boundary_gap <- up - alpha
    valid <- all(is.finite(weight)) && all(weight > 0) && is.finite(normalization_error) && normalization_error <= 
        check_tol && is.finite(calibration_error) && calibration_error <= check_tol && is.finite(boundary_gap) && 
        boundary_gap > 0
    if (!valid) {
        warning("ELW post-solution checks failed for ", group, ".")
        weight[] <- NA_real_
    }
    list(alpha = alpha, weight = weight, converged = valid, failure_reason = if (valid) NA_character_ else "post-solution check failed", 
        normalization_error = normalization_error, calibration_error = calibration_error, boundary_gap = boundary_gap, 
        lower_bound = low, upper_bound = up, f_lower = f_low, f_upper = f_up)
}
msm_probability_bounds <- function(probability, Lambda) {
    probability <- pmin(pmax(probability, 1e-06), 1 - 1e-06)
    list(lower = probability/(probability + Lambda * (1 - probability)), upper = Lambda * probability/(1 - 
        probability + Lambda * probability))
}
kernel_ghat <- function(gamma, cache, bandwidth) {
    if (isTRUE(cache$binned)) {
        denominator_bin <- as.numeric(rowsum(cache$delta_sorted * exp(gamma * cache$Y_sorted), cache$bin[cache$order]))
        z <- outer(cache$bin_center, cache$bin_center, "-")/bandwidth
        kernel_weight <- exp(-0.5 * z^2)
        numerator_fit <- as.vector(kernel_weight %*% cache$numerator_bin)
        denominator_fit <- as.vector(kernel_weight %*% denominator_bin)
        g_bin <- log(pmax(numerator_fit, 1e-10)) - log(pmax(denominator_fit, 1e-10))
        return(stats::approx(cache$bin_center, g_bin, xout = cache$U, method = "linear", rule = 2, ties = "ordered")$y)
    }
    numerator <- 1 - cache$delta_sorted
    denominator <- cache$delta_sorted * exp(gamma * cache$Y_sorted)
    numerator_fit_sorted <- stats::ksmooth(cache$U_sorted, numerator, kernel = "normal", bandwidth = bandwidth/0.3706506, 
        x.points = cache$U_sorted)$y
    denominator_fit_sorted <- stats::ksmooth(cache$U_sorted, denominator, kernel = "normal", bandwidth = bandwidth/0.3706506, 
        x.points = cache$U_sorted)$y
    numerator_fit <- denominator_fit <- numeric(length(cache$U))
    numerator_fit[cache$order] <- numerator_fit_sorted
    denominator_fit[cache$order] <- denominator_fit_sorted
    log(pmax(numerator_fit, 1e-10)) - log(pmax(denominator_fit, 1e-10))
}
simulation_display_method <- function(method, ps.method = NULL) {
    if (method %in% c("G-est", "Initial")) 
        return("CQR--G")
    if (method == "SIPW") 
        return(paste0("SIPW-", simulation_ps_label(ps.method)))
    if (method == "SIPW-MNAR") 
        return("SIPW-LMZ")
    if (method == "ELW-MAR") 
        return(paste0("ELW-", simulation_ps_label(ps.method)))
    if (method == "ELW-MNAR") 
        return("ELW-LMZ")
    method
}
make_csl_simulation_base <- function(dat, treatment.var = "treatment_smoking", outcome.var = "outcome_birthweight", 
    exclude.vars = c("row_id", "Smoke", "Birthweight", "treatment_smoking", "outcome_birthweight"), factor.vars = c("Momrace", 
        "Marital", "Education", "Parity", "Hxpreterm", "Hxstillbirth", "Diabetes", "Hxasthma", "HIV", 
        "Alcohol", "Preeclampsia", "AnteGDM", "AnteSGA", "AnteLGA", "Sex")) {
    A <- as.integer(dat[[treatment.var]])
    Y <- as.numeric(dat[[outcome.var]])
    covariates <- intersect(c("Momage", "Height", "PreBMI", "Momrace", "Marital", "Parity", "Hxpreterm", 
        "Diabetes", "Hxasthma", "HIV", "Alcohol", "Sex"), names(dat))
    xdat <- dat[covariates]
    factor.vars <- intersect(factor.vars, names(xdat))
    xdat[factor.vars] <- lapply(xdat[factor.vars], factor)
    varying <- vapply(xdat, function(x) length(unique(x)) > 1L, logical(1))
    xdat <- xdat[varying]
    X <- simulation_reference_last_matrix(xdat)
    X.std <- simulation_standardize(X)
    continuous <- intersect(c("Momage", "Height", "PreBMI"), colnames(X.std))
    X.expanded <- X.std
    for (v in continuous) {
        X.expanded <- cbind(X.expanded, X.std[, v]^2)
        colnames(X.expanded)[ncol(X.expanded)] <- paste0(v, "_sq")
    }
    if (all(c("Momage", "PreBMI") %in% continuous)) {
        X.expanded <- cbind(X.expanded, X.std[, "Momage"] * X.std[, "PreBMI"])
        colnames(X.expanded)[ncol(X.expanded)] <- "Momage_x_PreBMI"
    }
    X.expanded <- X.expanded[, !duplicated(colnames(X.expanded)), drop = FALSE]
    X.ps <- as.data.frame(X.expanded, check.names = TRUE)
    X.main <- as.data.frame(X.std, check.names = TRUE)
    ps.full.formula <- stats::as.formula(paste("A ~", paste(names(X.main), collapse = " + ")))
    ps.full.fit <- suppressWarnings(stats::glm(ps.full.formula, data = data.frame(A = A, X.main), family = stats::binomial()))
    ps.selected.fit <- ps.full.fit
    ps.selected.terms <- attr(stats::terms(ps.selected.fit), "term.labels")
    ell <- as.vector(stats::predict(ps.selected.fit, type = "link"))
    shadow.cols <- grepl("^Sex", colnames(X.std))
    X.no.shadow <- X.std[, !shadow.cols, drop = FALSE]
    ps.no.shadow.fit <- suppressWarnings(stats::glm.fit(cbind(`(Intercept)` = 1, X.no.shadow), A, family = stats::binomial()))
    ell.no.shadow <- as.vector(ps.no.shadow.fit$linear.predictors)
    outcome.data <- data.frame(Y = Y[A == 0], X.main[A == 0, , drop = FALSE])
    outcome.full.formula <- stats::as.formula(paste("Y ~", paste(names(X.main), collapse = " + ")))
    outcome.full.fit <- suppressWarnings(stats::lm(outcome.full.formula, data = outcome.data))
    outcome.selected.fit <- outcome.full.fit
    outcome.selected.terms <- attr(stats::terms(outcome.selected.fit), "term.labels")
    m0 <- as.numeric(stats::predict(outcome.selected.fit, newdata = data.frame(X.main, check.names = TRUE)))
    resid0 <- as.numeric(Y[A == 0] - m0[A == 0])
    sigma0 <- stats::sd(resid0, na.rm = TRUE)
    eps.pool <- (resid0 - mean(resid0, na.rm = TRUE))/sigma0
    eps.pool <- eps.pool[is.finite(eps.pool)]
    lbw.fit <- suppressWarnings(stats::glm.fit(cbind(`(Intercept)` = 1, X.std), as.integer(Y < 2500), 
        family = stats::binomial()))
    lbw.coef <- lbw.fit$coefficients
    lbw.coef[!is.finite(lbw.coef)] <- 0
    lbw.risk <- stats::plogis(as.vector(cbind(`(Intercept)` = 1, X.std) %*% lbw.coef))
    risk.percentile <- rank(lbw.risk, ties.method = "average")/length(lbw.risk)
    list(dat = dat, A = A, Y = Y, X = X.std, X.no.shadow = X.no.shadow, X.expanded = X.expanded, ell = ell, 
        ell.no.shadow = ell.no.shadow, m0 = m0, sigma0 = sigma0, eps.pool = eps.pool, lbw.risk = lbw.risk, 
        risk.percentile = risk.percentile, smoking.prevalence = mean(A == 1), ps.selected.terms = ps.selected.terms, 
        ps.selected.formula = stats::formula(ps.selected.fit), outcome.selected.terms = outcome.selected.terms, 
        outcome.selected.formula = stats::formula(outcome.selected.fit), outcome.selected.coefficients = stats::coef(outcome.selected.fit), 
        covariates = colnames(X.std), Z = if ("Sex" %in% names(xdat)) factor(xdat$Sex) else factor(rep(1, 
            length(A))))
}
generate_csl_simulation_data <- function(base, n = nrow(base$dat), b = 1, gamma0 = 0, gamma1 = 0, effect.low = -150, 
    effect.high = -300, eta.sd = 50, instrument.effect = 300, heteroskedastic = FALSE, outcome.nonlinearity = 0, 
    effect.nonlinearity = 0, treatment.nonlinear = FALSE, overlap.step = FALSE, overlap.low.prob = 1e-04, 
    overlap.high.quantile = 0.9, seed = NULL) {
    if (!is.null(seed)) 
        set.seed(seed)
    is.mnar <- gamma0 != 0 || gamma1 != 0
    id <- sample.int(nrow(base$dat), n, replace = TRUE)
    if (is.mnar) {
        X <- base$X[id, , drop = FALSE]
        X.no.shadow <- base$X.no.shadow[id, , drop = FALSE]
        X.expanded <- base$X.expanded[id, , drop = FALSE]
        Z <- factor(base$Z[id], levels = levels(base$Z))
        ell <- base$ell.no.shadow[id]
        u <- as.numeric(scale(ell))
        u[!is.finite(u)] <- 0
        control <- base$A == 0L
        z.mean <- tapply(base$Y[control], base$Z[control], mean)
        centered <- base$Y[control] - as.numeric(z.mean[as.character(base$Z[control])])
        sigma0 <- stats::sd(centered, na.rm = TRUE)
        z.high <- as.numeric(Z == levels(Z)[length(levels(Z))])
        z.share <- mean(base$Z == levels(base$Z)[length(levels(base$Z))])
        mu0 <- mean(base$Y[control]) + instrument.effect * (z.high - z.share)
        tau <- mean(c(effect.low, effect.high))
        y.center <- mean(base$Y, na.rm = TRUE)
        gamma <- gamma1
        log.m <- gamma * (mu0 + tau - y.center)/sigma0 + 0.5 * gamma^2
        h <- b * u
        a <- calibrate_logit_intercept(h + log.m, base$smoking.prevalence)
        p.treated <- stats::plogis(a + h + log.m)
        T <- stats::rbinom(n, 1L, p.treated)
        Y0 <- stats::rnorm(n, mu0 + T * gamma * sigma0, sigma0)
        Y1 <- Y0 + tau
        Y <- ifelse(T == 1L, Y1, Y0)
        pi.true <- stats::plogis(a + h + gamma * (Y1 - y.center)/sigma0)
        r <- rank(mu0, ties.method = "average")/n
        return(structure(cbind(data.frame(id = id, Z = Z, T = T, Y = Y, Y0 = Y0, Y1 = Y1, pi.true = pi.true, 
            rho1.true = pi.true, rho0.true = 1 - pi.true, ell = ell, overlap_stratum = 0, risk.percentile = r, 
            delta = tau, check.names = FALSE), as.data.frame(X, check.names = FALSE)), X.expanded = X.expanded, 
            X.main = X.no.shadow, ps.selected.terms = colnames(X.no.shadow), intercept.a = a, b = b, 
            gamma0 = gamma0, gamma1 = gamma1, outcome.nonlinearity = 0, effect.nonlinearity = 0, lmz.compatible = TRUE, 
            shadow.variable = "Sex"))
    }
    ell <- base$ell[id]
    X <- base$X[id, , drop = FALSE]
    X.expanded <- base$X.expanded[id, , drop = FALSE]
    m0 <- base$m0[id]
    sampled.lbw.risk <- base$lbw.risk[id]
    r <- rank(sampled.lbw.risk, ties.method = "average")/n
    Z <- if (is.mnar) {
        factor(stats::rbinom(n, 1, 0.5), levels = 0:1)
    }
    else {
        factor(base$Z[id], levels = levels(base$Z))
    }
    eps <- sample(base$eps.pool, n, replace = TRUE)
    if (heteroskedastic) {
        sigma <- base$sigma0 * (0.75 + 0.5 * r)
    }
    else {
        sigma <- rep(base$sigma0, n)
    }
    Y0 <- m0 + sigma * eps
    nonlinear.shape <- (r - 0.5)^2 + 0.75 * as.numeric(r > 0.8)
    nonlinear.shape <- nonlinear.shape - mean(nonlinear.shape)
    Y0 <- Y0 + outcome.nonlinearity * nonlinear.shape
    if (is.mnar) {
        Y0 <- Y0 + instrument.effect * (as.numeric(Z) - 1.5)
    }
    delta <- effect.low + (effect.high - effect.low) * r
    effect.stress <- as.numeric(r > 0.8)
    Y1 <- Y0 + delta + effect.nonlinearity * effect.stress + stats::rnorm(n, mean = 0, sd = eta.sd)
    z0 <- as.numeric(scale(Y0))
    z1 <- as.numeric(scale(Y1))
    z0[!is.finite(z0)] <- 0
    z1[!is.finite(z1)] <- 0
    overlap.stratum <- as.numeric(ell > stats::quantile(base$ell, overlap.high.quantile, names = FALSE))
    if (overlap.step && !is.mnar) {
        share.high <- mean(overlap.stratum)
        high.prob <- (base$smoking.prevalence - (1 - share.high) * overlap.low.prob)/share.high
        high.prob <- min(max(high.prob, overlap.low.prob), 1 - 1e-08)
        a <- stats::qlogis(overlap.low.prob)
        treatment.lp <- (stats::qlogis(high.prob) - a) * overlap.stratum
    }
    else {
        treatment.lp <- b * ell
        if (treatment.nonlinear) {
            age <- if ("Momage" %in% colnames(X)) 
                X[, "Momage"]
            else 0
            bmi <- if ("PreBMI" %in% colnames(X)) 
                X[, "PreBMI"]
            else 0
            height <- if ("Height" %in% colnames(X)) 
                X[, "Height"]
            else 0
            treatment.lp <- treatment.lp + 0.6 * (age^2 - 1) + 0.6 * sin(bmi) + 0.5 * age * bmi - 0.4 * 
                as.numeric(height > 0.75)
        }
    }
    if (!overlap.step && isTRUE(all.equal(b, 1)) && gamma0 == 0 && gamma1 == 0) {
        a <- 0
    }
    else if (!overlap.step || is.mnar) {
        a <- calibrate_logit_intercept(treatment.lp + gamma1 * z1 + gamma0 * z0, base$smoking.prevalence)
    }
    if (overlap.step && !is.mnar) {
        a <- a
    }
    else {
        a <- a
    }
    lp <- a + treatment.lp + gamma1 * z1 + gamma0 * z0
    pi.true <- stats::plogis(lp)
    T <- stats::rbinom(n, 1, pi.true)
    Y <- ifelse(T == 1, Y1, Y0)
    structure(cbind(data.frame(id = id, Z = Z, T = T, Y = Y, Y0 = Y0, Y1 = Y1, pi.true = pi.true, rho1.true = pi.true, 
        rho0.true = 1 - pi.true, ell = ell, overlap_stratum = overlap.stratum, risk.percentile = r, delta = delta, 
        check.names = FALSE), as.data.frame(X, check.names = FALSE)), X.expanded = X.expanded, X.main = X, 
        ps.selected.terms = base$ps.selected.terms, intercept.a = a, b = b, gamma0 = gamma0, gamma1 = gamma1, 
        outcome.nonlinearity = outcome.nonlinearity, effect.nonlinearity = effect.nonlinearity)
}
fit_simulation_propensity <- function(sim.dat, use.expanded = FALSE) {
    X <- if (!use.expanded && !is.null(attr(sim.dat, "X.main"))) {
        attr(sim.dat, "X.main")
    }
    else if (!is.null(attr(sim.dat, "X.expanded"))) {
        attr(sim.dat, "X.expanded")
    }
    else {
        as.matrix(sim.dat[, setdiff(names(sim.dat), c("id", "Z", "T", "Y", "Y0", "Y1", "pi.true", "rho1.true", 
            "rho0.true", "ell", "risk.percentile", "delta")), drop = FALSE])
    }
    X.ps <- as.data.frame(X, check.names = TRUE)
    form <- stats::as.formula(paste("T ~", paste(names(X.ps), collapse = " + ")))
    fit <- suppressWarnings(stats::glm(form, data = data.frame(T = sim.dat$T, X.ps), family = stats::binomial()))
    p <- as.numeric(stats::predict(fit, type = "response"))
    eta <- as.numeric(stats::predict(fit, type = "link"))
    list(probability = p, linear_predictor = eta, coefficient = stats::coef(fit), auc = simulation_auc_rank(sim.dat$T, 
        p))
}
fit_simulation_propensity_gbm <- function(sim.dat, seed = 20260811L, n.trees = 500L, interaction.depth = 2L, 
    min.node.size = 100L, shrinkage = 0.02, bag.fraction = 0.5) {
    if (!requireNamespace("gbm", quietly = TRUE)) {
        stop("GBM propensity estimation requires the gbm package.")
    }
    X <- if (!is.null(attr(sim.dat, "X.expanded"))) {
        as.matrix(attr(sim.dat, "X.expanded"))
    }
    else {
        as.matrix(sim.dat[, setdiff(names(sim.dat), c("id", "Z", "T", "Y", "Y0", "Y1", "pi.true", "rho1.true", 
            "rho0.true", "ell", "risk.percentile", "delta")), drop = FALSE])
    }
    set.seed(seed)
    fit <- gbm::gbm.fit(x = X, y = sim.dat$T, distribution = "bernoulli", n.trees = n.trees, interaction.depth = interaction.depth, 
        n.minobsinnode = min.node.size, shrinkage = shrinkage, bag.fraction = bag.fraction, nTrain = nrow(X), 
        keep.data = FALSE, verbose = FALSE)
    p <- as.numeric(stats::predict(fit, newdata = X, n.trees = n.trees, type = "response"))
    list(probability = p, linear_predictor = stats::qlogis(p), coefficient = NULL, auc = simulation_auc_rank(sim.dat$T, 
        p), min.prob = min(p), p01.prob = as.numeric(stats::quantile(p, 0.01, names = FALSE)), p99.prob = as.numeric(stats::quantile(p, 
        0.99, names = FALSE)), max.prob = max(p))
}
fit_simulation_propensity_honest_tree <- function(sim.dat, seed = 20260811L, minbucket = NULL, maxdepth = 30L) {
    if (!requireNamespace("rpart", quietly = TRUE)) {
        stop("Honest-tree propensity estimation requires the rpart package.")
    }
    X <- attr(sim.dat, "X.main")
    if (is.null(X)) 
        stop("Honest tree requires the fixed main-effect design matrix.")
    X <- as.data.frame(X, check.names = TRUE)
    n <- nrow(X)
    set.seed(seed)
    fold <- integer(n)
    for (tt in 0:1) {
        ii <- which(sim.dat$T == tt)
        fold[ii] <- sample(rep(1:2, length.out = length(ii)))
    }
    p <- rep(NA_real_, n)
    for (k in 1:2) {
        train <- fold != k
        test <- fold == k
        leaf_min <- if (is.null(minbucket)) 
            ceiling(sum(train)^(1/3))
        else minbucket
        fit <- rpart::rpart(factor(T, levels = 0:1) ~ ., data = data.frame(T = sim.dat$T[train], X[train, 
            ]), method = "class", control = rpart::rpart.control(cp = 0, minsplit = 2L * leaf_min, minbucket = leaf_min, 
            maxdepth = maxdepth, xval = 0))
        pred <- stats::predict(fit, newdata = X[test, ], type = "prob")
        p[test] <- if ("1" %in% colnames(pred)) 
            pred[, "1"]
        else 0
    }
    list(probability = p, linear_predictor = rep(NA_real_, n), coefficient = NULL, auc = simulation_auc_rank(sim.dat$T, 
        p), min.prob = min(p), p01.prob = as.numeric(stats::quantile(p, 0.01)), p99.prob = as.numeric(stats::quantile(p, 
        0.99)), max.prob = max(p), boundary.zero = sum(p == 0), boundary.one = sum(p == 1), observed.boundary = sum((sim.dat$T == 
        1 & p == 0) | (sim.dat$T == 0 & p == 1)))
}
fit_simulation_propensity_kernel <- function(sim.dat, pilot.fit = NULL, bandwidth = NULL) {
    if (!exists("pca_kernel_predict", mode = "function")) {
        stop("PCA--Kernel requires pca_kernel_utils.R to be sourced.")
    }
    X <- attr(sim.dat, "X.main")
    if (is.null(X)) 
        stop("PCA--Kernel requires the fixed main-effect design matrix.")
    fit <- pca_kernel_predict(X_train = X, treatment_train = sim.dat$T, X_test = X, n_components = 3L, 
        grid_bins = 45L, kernel = "gaussian", bandwidth_multiplier = 1)
    p <- fit$probability
    list(probability = p, linear_predictor = rep(NA_real_, length(p)), coefficient = NULL, bandwidth = fit$bandwidth, 
        explained = fit$explained, auc = simulation_auc_rank(sim.dat$T, p), min.prob = min(p, na.rm = TRUE), 
        p01.prob = as.numeric(stats::quantile(p, 0.01, names = FALSE, na.rm = TRUE)), p99.prob = as.numeric(stats::quantile(p, 
            0.99, names = FALSE, na.rm = TRUE)), max.prob = max(p, na.rm = TRUE))
}
initial_qr_side_distribution <- function(y.obs, x.obs, x.target, probs, threshold, grid, tau.grid = seq(0.05, 
    0.95, by = 0.05)) {
    if (!requireNamespace("quantreg", quietly = TRUE)) {
        stop("The G-est estimator requires the quantreg package.")
    }
    y.obs <- as.numeric(y.obs)
    x.obs <- as.matrix(x.obs)
    x.target <- as.matrix(x.target)
    if (nrow(x.obs) != length(y.obs) || ncol(x.obs) != ncol(x.target)) {
        stop("Incompatible outcome and covariate dimensions for G-est.")
    }
    tau.grid <- sort(unique(as.numeric(tau.grid)))
    if (length(tau.grid) < 3L || tau.grid[1] <= 0 || tau.grid[length(tau.grid)] >= 1) {
        stop("G-est requires an interior grid of at least three quantile levels.")
    }
    step <- diff(tau.grid)
    if (max(abs(step - step[1])) > 1e-10) {
        stop("G-est currently requires an equally spaced quantile grid.")
    }
    epsilon <- tau.grid[1]
    if (abs(tau.grid[length(tau.grid)] - (1 - epsilon)) > 1e-10) {
        stop("The G-est quantile grid must run from epsilon to 1 - epsilon.")
    }
    x.obs <- cbind(`(Intercept)` = 1, x.obs)
    x.target <- cbind(`(Intercept)` = 1, x.target)
    svd.x <- svd(x.obs, nu = 0)
    keep <- svd.x$d > max(svd.x$d) * sqrt(.Machine$double.eps)
    transform <- sweep(svd.x$v[, keep, drop = FALSE], 2, svd.x$d[keep]/sqrt(nrow(x.obs)), "/")
    x.obs <- x.obs %*% transform
    x.target <- x.target %*% transform
    coefficients <- suppressWarnings(vapply(tau.grid, function(tau) quantreg::rq.fit.fnb(x.obs, y.obs, 
        tau = tau)$coefficients, numeric(ncol(x.obs))))
    if (any(!is.finite(coefficients))) {
        stop("G-est conditional quantile regression did not converge.")
    }
    fitted.quantiles <- x.target %*% coefficients
    support <- sort(as.numeric(fitted.quantiles[, -1L, drop = FALSE]))
    interior.mass <- 1 - 2 * epsilon
    cdf_at <- function(y, strict = FALSE) {
        if (strict) {
            epsilon + interior.mass * (findInterval(y, support, left.open = TRUE)/length(support))
        }
        else {
            epsilon + interior.mass * (findInterval(y, support)/length(support))
        }
    }
    adjusted.probs <- (probs - epsilon)/interior.mass
    if (any(adjusted.probs <= 0 | adjusted.probs >= 1)) {
        stop("Requested quantiles fall outside the trimmed G-est distribution.")
    }
    list(quantile = as.numeric(stats::quantile(support, adjusted.probs, names = FALSE, type = 1)), lbw = cdf_at(threshold, 
        strict = TRUE), F = vapply(grid, cdf_at, numeric(1)), ess = NA_real_, max.weight = NA_real_, 
        p99.weight = NA_real_, epsilon = epsilon, tau.count = length(tau.grid), rank = ncol(x.obs))
}
estimate_initial_qr_simulation <- function(sim.dat, truth, threshold = 2500, method = "G-est", tau.grid = seq(0.05, 
    0.95, by = 0.05)) {
    X <- attr(sim.dat, "X.expanded")
    if (is.null(X)) {
        stop("G-est requires the expanded covariate design matrix.")
    }
    X <- as.matrix(X)
    idx1 <- sim.dat$T == 1
    idx0 <- sim.dat$T == 0
    side1 <- initial_qr_side_distribution(sim.dat$Y[idx1], X[idx1, , drop = FALSE], X, truth$probs, threshold, 
        truth$grid, tau.grid = tau.grid)
    side0 <- initial_qr_side_distribution(sim.dat$Y[idx0], X[idx0, , drop = FALSE], X, truth$probs, threshold, 
        truth$grid, tau.grid = tau.grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    out$diagnostics$raw.endpoint1 <- 1 - side1$epsilon
    out$diagnostics$raw.endpoint0 <- 1 - side0$epsilon
    out$diagnostics$raw.endpoint.error <- max(side1$epsilon, side0$epsilon)
    out
}
estimate_sipw_simulation <- function(sim.dat, ps, truth, threshold = 2500, method = "SIPW") {
    ps <- simulation_validate_probability(ps, "propensity score", allow.boundary = TRUE)
    idx1 <- sim.dat$T == 1
    idx0 <- sim.dat$T == 0
    observed.ps <- c(ps[idx1], 1 - ps[idx0])
    if (any(!is.finite(observed.ps)) || any(observed.ps <= 0 | observed.ps > 1)) {
        stop("observed SIPW selection probability contains a zero or invalid value.")
    }
    w1 <- mean(sim.dat$T == 1)/ps[idx1]
    w0 <- mean(sim.dat$T == 0)/(1 - ps[idx0])
    side1 <- estimate_weighted_distribution(sim.dat$Y[idx1], w1, truth$probs, threshold, truth$grid)
    side0 <- estimate_weighted_distribution(sim.dat$Y[idx0], w0, truth$probs, threshold, truth$grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    smd <- simulation_smd_summary(sim.dat, w1, w0)
    out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
    out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
    out
}
estimate_elw_simulation <- function(sim.dat, rho, truth, threshold = 2500, method = "ELW") {
    rho <- simulation_validate_probability(rho, "observation probability", allow.boundary = TRUE)
    idx1 <- sim.dat$T == 1
    idx0 <- sim.dat$T == 0
    fit1 <- elw_group(sim.dat$T, rho, "treated")
    fit0 <- elw_group(sim.dat$T, rho, "control")
    side1 <- estimate_weighted_distribution(sim.dat$Y[idx1], fit1$weight, truth$probs, threshold, truth$grid)
    side0 <- estimate_weighted_distribution(sim.dat$Y[idx0], fit0$weight, truth$probs, threshold, truth$grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    smd <- simulation_smd_summary(sim.dat, fit1$weight, fit0$weight)
    out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
    out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
    out$diagnostics$root1 <- fit1$alpha
    out$diagnostics$root0 <- fit0$alpha
    out$diagnostics$elw.converged1 <- fit1$converged
    out$diagnostics$elw.converged0 <- fit0$converged
    out$diagnostics$elw.norm.error1 <- fit1$normalization_error
    out$diagnostics$elw.norm.error0 <- fit0$normalization_error
    out$diagnostics$elw.cal.error1 <- fit1$calibration_error
    out$diagnostics$elw.cal.error0 <- fit0$calibration_error
    out$diagnostics$elw.boundary.gap1 <- fit1$boundary_gap
    out$diagnostics$elw.boundary.gap0 <- fit0$boundary_gap
    out$diagnostics$elw.failure1 <- fit1$failure_reason
    out$diagnostics$elw.failure0 <- fit0$failure_reason
    out
}
lmz_gauss_hermite <- function(n = 15L) {
    idx <- seq_len(n - 1L)
    J <- matrix(0, n, n)
    J[cbind(idx, idx + 1L)] <- sqrt(idx/2)
    J[cbind(idx + 1L, idx)] <- sqrt(idx/2)
    ee <- eigen(J, symmetric = TRUE)
    ord <- order(ee$values)
    list(nodes = ee$values[ord], weights = sqrt(pi) * ee$vectors[1L, ord]^2)
}
estimate_lmz_side <- function(sim.dat, R, ps.fit, spline_df = 1L, quadrature_n = 15L) {
    R <- as.integer(R)
    observed <- R == 1L
    if (sum(observed) < 50L || sum(!observed) < 20L) {
        stop("LMZ requires both observed and missing outcomes on each side.")
    }
    if (is.null(ps.fit$linear_predictor) && is.null(sim.dat$ell)) {
        stop("LMZ requires a fitted covariate index for h(U; beta).")
    }
    u.source <- if (!is.null(sim.dat$ell)) 
        sim.dat$ell
    else ps.fit$linear_predictor
    u <- as.numeric(scale(u.source))
    u[!is.finite(u)] <- 0
    z <- stats::model.matrix(~Z, data = sim.dat)
    obs_fit <- stats::lm.fit(z[observed, , drop = FALSE], sim.dat$Y[observed])
    alpha <- obs_fit$coefficients
    alpha[!is.finite(alpha)] <- 0
    mu <- as.numeric(z %*% alpha)
    sigma <- sqrt(mean(obs_fit$residuals^2))
    if (!is.finite(sigma) || sigma <= 0) 
        stop("Invalid LMZ complete-case outcome scale.")
    y_center <- mean(sim.dat$Y[observed])
    y_scale <- stats::sd(sim.dat$Y[observed])
    y_std <- (sim.dat$Y[observed] - y_center)/y_scale
    if (spline_df == 1L) {
        basis_at <- function(y) matrix((y - y_center)/y_scale, ncol = 1L)
    }
    else {
        basis_fit <- splines::ns(y_std, df = spline_df, intercept = FALSE)
        basis_at <- function(y) {
            b <- predict(basis_fit, newx = (y - y_center)/y_scale)
            b0 <- predict(basis_fit, newx = 0)
            sweep(b, 2L, as.numeric(b0), "-")
        }
    }
    gh <- lmz_gauss_hermite(quadrature_n)
    y_quad <- as.vector(t(outer(mu, sqrt(2) * sigma * gh$nodes, "+")))
    B_quad <- basis_at(y_quad)
    q_index <- rep(seq_len(nrow(sim.dat)), each = quadrature_n)
    q_weight <- rep(gh$weights/sqrt(pi), times = nrow(sim.dat))
    Xh <- cbind(1, u)
    init_glm <- suppressWarnings(stats::glm.fit(Xh, R, family = stats::binomial()))
    init_beta <- init_glm$coefficients
    init_beta[!is.finite(init_beta)] <- 0
    init <- c(init_beta, rep(0, ncol(B_quad)))
    objective <- function(par) {
        beta <- par[seq_len(ncol(Xh))]
        gamma <- par[-seq_len(ncol(Xh))]
        log_terms <- -as.numeric(B_quad %*% gamma) + log(q_weight)
        max_by_i <- as.numeric(tapply(log_terms, q_index, max))
        exp_sum <- rowsum(exp(log_terms - max_by_i[q_index]), q_index, reorder = FALSE)
        log_m <- max_by_i + log(as.numeric(exp_sum))
        w <- stats::plogis(as.numeric(Xh %*% beta) - log_m)
        nll <- -sum(R * log(pmax(w, 1e-12)) + (1 - R) * log(pmax(1 - w, 1e-12)))
        penalty <- if (length(gamma) >= 3L) 
            sum(diff(gamma, differences = 2L)^2)
        else 0
        nll + 0.001 * penalty
    }
    opt <- stats::optim(init, objective, method = "BFGS", control = list(maxit = 300L, reltol = 1e-08))
    beta <- opt$par[seq_len(ncol(Xh))]
    gamma <- opt$par[-seq_len(ncol(Xh))]
    g_observed <- as.numeric(basis_at(sim.dat$Y[observed]) %*% gamma)
    rho_observed <- stats::plogis(as.numeric(Xh[observed, , drop = FALSE] %*% beta) + g_observed)
    rho <- rep(NA_real_, nrow(sim.dat))
    rho[observed] <- rho_observed
    list(rho = rho, beta = beta, gamma = gamma, alpha = alpha, sigma = sigma, convergence = opt$convergence, 
        objective = opt$value)
}
estimate_lmz_probabilities <- function(sim.dat, ps.fit) {
    fit1 <- estimate_lmz_side(sim.dat, sim.dat$T == 1L, ps.fit)
    fit0 <- estimate_lmz_side(sim.dat, sim.dat$T == 0L, ps.fit)
    list(fit1 = fit1, fit0 = fit0)
}
estimate_lmz_mnar_simulation <- function(sim.dat, ps.fit, truth, threshold = 2500, method = c("ELW-LMZ", 
    "SIPW-LMZ")) {
    method <- match.arg(method)
    lmz <- estimate_lmz_probabilities(sim.dat, ps.fit)
    idx1 <- sim.dat$T == 1L
    idx0 <- sim.dat$T == 0L
    rho1 <- lmz$fit1$rho[idx1]
    rho0 <- lmz$fit0$rho[idx0]
    simulation_validate_probability(c(rho1, rho0), "observed LMZ probability")
    if (method == "ELW-LMZ") {
        rho_for_1 <- rep(0.5, nrow(sim.dat))
        rho_for_1[idx1] <- rho1
        rho_for_0 <- rep(0.5, nrow(sim.dat))
        rho_for_0[idx0] <- 1 - rho0
        fit1 <- elw_group(sim.dat$T, rho_for_1, "treated")
        fit0 <- elw_group(sim.dat$T, rho_for_0, "control")
        w1 <- fit1$weight
        w0 <- fit0$weight
    }
    else {
        w1 <- mean(idx1)/rho1
        w0 <- mean(idx0)/rho0
    }
    side1 <- estimate_weighted_distribution(sim.dat$Y[idx1], w1, truth$probs, threshold, truth$grid)
    side0 <- estimate_weighted_distribution(sim.dat$Y[idx0], w0, truth$probs, threshold, truth$grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    out$diagnostics$lmz.converged1 <- lmz$fit1$convergence == 0L
    out$diagnostics$lmz.converged0 <- lmz$fit0$convergence == 0L
    out$diagnostics$ps.min <- min(c(rho1, rho0))
    out$diagnostics$ps.max <- max(c(rho1, rho0))
    if (method == "ELW-LMZ") {
        out$diagnostics$elw.norm.error1 <- fit1$normalization_error
        out$diagnostics$elw.norm.error0 <- fit0$normalization_error
        out$diagnostics$elw.cal.error1 <- fit1$calibration_error
        out$diagnostics$elw.cal.error0 <- fit0$calibration_error
    }
    out
}
estimate_elw_mnar_simulation <- function(sim.dat, ps.fit, truth, threshold = 2500, method = "ELW-MNAR", 
    gamma_range1 = c(0, 0.01), gamma_range0 = c(-0.01, 0)) {
    shao <- estimate_shao_mnar_probabilities(sim.dat, ps.fit, gamma_range1 = gamma_range1, gamma_range0 = gamma_range0)
    idx1 <- sim.dat$T == 1
    idx0 <- sim.dat$T == 0
    fit1 <- elw_group(sim.dat$T, shao$rho1, "treated")
    fit0 <- elw_group(sim.dat$T, 1 - shao$rho0, "control")
    side1 <- estimate_weighted_distribution(sim.dat$Y[idx1], fit1$weight, truth$probs, threshold, truth$grid)
    side0 <- estimate_weighted_distribution(sim.dat$Y[idx0], fit0$weight, truth$probs, threshold, truth$grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    smd <- simulation_smd_summary(sim.dat, fit1$weight, fit0$weight)
    out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
    out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
    out$diagnostics$root1 <- fit1$alpha
    out$diagnostics$root0 <- fit0$alpha
    out$diagnostics$elw.converged1 <- fit1$converged
    out$diagnostics$elw.converged0 <- fit0$converged
    out$diagnostics$elw.norm.error1 <- fit1$normalization_error
    out$diagnostics$elw.norm.error0 <- fit0$normalization_error
    out$diagnostics$elw.cal.error1 <- fit1$calibration_error
    out$diagnostics$elw.cal.error0 <- fit0$calibration_error
    out$diagnostics$elw.boundary.gap1 <- fit1$boundary_gap
    out$diagnostics$elw.boundary.gap0 <- fit0$boundary_gap
    out$diagnostics$elw.failure1 <- fit1$failure_reason
    out$diagnostics$elw.failure0 <- fit0$failure_reason
    observed.rho <- c(shao$rho1[idx1], shao$rho0[idx0])
    out$diagnostics$ps.min <- min(observed.rho)
    out$diagnostics$ps.p01 <- as.numeric(stats::quantile(observed.rho, 0.01, names = FALSE))
    out$diagnostics$ps.p99 <- as.numeric(stats::quantile(observed.rho, 0.99, names = FALSE))
    out$diagnostics$ps.max <- max(observed.rho)
    out$diagnostics$shao_gamma1 <- shao$shao1$gamma
    out$diagnostics$shao_gamma0 <- shao$shao0$gamma
    out$diagnostics$shao_bandwidth1 <- shao$shao1$bandwidth
    out$diagnostics$shao_bandwidth0 <- shao$shao0$bandwidth
    out
}
estimate_lql_mnar_simulation <- function(sim.dat, truth, threshold = 2500, method = "ELW-LQL") {
    X1 <- attr(sim.dat, "X.main")
    if (is.null(X1)) 
        stop("LQL simulation requires the main-effect covariate matrix.")
    Z <- stats::model.matrix(~Z - 1, data = sim.dat)
    if (ncol(Z) > 1L) 
        Z <- Z[, -1L, drop = FALSE]
    fit1 <- lql_mnar_fit_side(sim.dat$Y, X1, Z, R = as.integer(sim.dat$T == 1L))
    fit0 <- lql_mnar_fit_side(sim.dat$Y, X1, Z, R = as.integer(sim.dat$T == 0L))
    idx1 <- sim.dat$T == 1L
    idx0 <- sim.dat$T == 0L
    elw1 <- elw_group(sim.dat$T, fit1$rho, "treated")
    elw0 <- elw_group(sim.dat$T, 1 - fit0$rho, "control")
    side1 <- estimate_weighted_distribution(sim.dat$Y[idx1], elw1$weight, truth$probs, threshold, truth$grid)
    side0 <- estimate_weighted_distribution(sim.dat$Y[idx0], elw0$weight, truth$probs, threshold, truth$grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    smd <- simulation_smd_summary(sim.dat, elw1$weight, elw0$weight)
    out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
    out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
    out$diagnostics$lql.gamma1 <- fit1$gamma
    out$diagnostics$lql.gamma0 <- fit0$gamma
    out$diagnostics$lql.converged1 <- fit1$selection_converged
    out$diagnostics$lql.converged0 <- fit0$selection_converged
    observed.prob <- c(fit1$rho[idx1], 1 - fit0$rho[idx0])
    out$diagnostics$ps.min <- min(observed.prob)
    out$diagnostics$ps.p01 <- unname(stats::quantile(observed.prob, 0.01))
    out$diagnostics$ps.p99 <- unname(stats::quantile(observed.prob, 0.99))
    out$diagnostics$ps.max <- max(observed.prob)
    out$diagnostics$elw.norm.error1 <- elw1$normalization_error
    out$diagnostics$elw.norm.error0 <- elw0$normalization_error
    out$diagnostics$elw.cal.error1 <- elw1$calibration_error
    out$diagnostics$elw.cal.error0 <- elw0$calibration_error
    out
}
estimate_sipw_lql_mnar_simulation <- function(sim.dat, truth, threshold = 2500, method = "SIPW-LQL") {
    X1 <- attr(sim.dat, "X.main")
    if (is.null(X1)) 
        stop("LQL simulation requires the main-effect covariate matrix.")
    Z <- stats::model.matrix(~Z - 1, data = sim.dat)
    if (ncol(Z) > 1L) 
        Z <- Z[, -1L, drop = FALSE]
    fit1 <- lql_mnar_fit_side(sim.dat$Y, X1, Z, R = as.integer(sim.dat$T == 1L))
    fit0 <- lql_mnar_fit_side(sim.dat$Y, X1, Z, R = as.integer(sim.dat$T == 0L))
    idx1 <- sim.dat$T == 1L
    idx0 <- sim.dat$T == 0L
    observed.prob <- c(fit1$rho[idx1], 1 - fit0$rho[idx0])
    if (any(!is.finite(observed.prob)) || any(observed.prob <= 0 | observed.prob > 1)) {
        stop("observed SIPW-LQL selection probability contains a zero or invalid value.")
    }
    w1 <- mean(idx1)/fit1$rho[idx1]
    w0 <- mean(idx0)/(1 - fit0$rho[idx0])
    side1 <- estimate_weighted_distribution(sim.dat$Y[idx1], w1, truth$probs, threshold, truth$grid)
    side0 <- estimate_weighted_distribution(sim.dat$Y[idx0], w0, truth$probs, threshold, truth$grid)
    out <- make_simulation_result_rows(method, side1, side0, truth)
    smd <- simulation_smd_summary(sim.dat, w1, w0)
    out$diagnostics$max.abs.smd <- smd[["max.abs.smd"]]
    out$diagnostics$mean.abs.smd <- smd[["mean.abs.smd"]]
    out$diagnostics$lql.gamma1 <- fit1$gamma
    out$diagnostics$lql.gamma0 <- fit0$gamma
    out$diagnostics$lql.converged1 <- fit1$selection_converged
    out$diagnostics$lql.converged0 <- fit0$selection_converged
    out$diagnostics$ps.min <- min(observed.prob)
    out$diagnostics$ps.p01 <- unname(stats::quantile(observed.prob, 0.01))
    out$diagnostics$ps.p99 <- unname(stats::quantile(observed.prob, 0.99))
    out$diagnostics$ps.max <- max(observed.prob)
    out
}
run_one_csl_simulation <- function(base, iter, n = 5000, b = 1, gamma0 = 0, gamma1 = 0, probs = c(0.1, 
    0.25, 0.5, 0.75), threshold = 2500, seed.offset = 20260623, heteroskedastic = TRUE, outcome.nonlinearity = 0, 
    effect.nonlinearity = 0, treatment.nonlinear = FALSE, instrument.effect = 300, overlap.step = FALSE, 
    overlap.low.prob = 1e-04, overlap.high.quantile = 0.9, include_mnar_sipw = TRUE, include_mnar_elw = TRUE, 
    include_mnar_lql = TRUE, gbm.n.trees = 500L, gbm.interaction.depth = 2L, gbm.min.node.size = 100L, 
    gbm.shrinkage = 0.02, gbm.bag.fraction = 0.5, ps.method = c("PLE", "GBM", "SLE", "kernel", "honest_tree")) {
    ps.method <- match.arg(ps.method)
    sim.dat <- generate_csl_simulation_data(base = base, n = n, b = b, gamma0 = gamma0, gamma1 = gamma1, 
        eta.sd = if (gamma0 != 0 || gamma1 != 0) 
            0
        else 50, instrument.effect = instrument.effect, heteroskedastic = heteroskedastic, outcome.nonlinearity = outcome.nonlinearity, 
        effect.nonlinearity = effect.nonlinearity, treatment.nonlinear = treatment.nonlinear, overlap.step = overlap.step, 
        overlap.low.prob = overlap.low.prob, overlap.high.quantile = overlap.high.quantile, seed = seed.offset + 
            iter)
    truth <- simulation_true_estimands(sim.dat, probs = probs, threshold = threshold)
    ps.pilot <- fit_simulation_propensity(sim.dat)
    ps.fit <- switch(ps.method, PLE = ps.pilot, GBM = fit_simulation_propensity_gbm(sim.dat, seed = seed.offset + 
        100000L + iter, n.trees = gbm.n.trees, interaction.depth = gbm.interaction.depth, min.node.size = gbm.min.node.size, 
        shrinkage = gbm.shrinkage, bag.fraction = gbm.bag.fraction), honest_tree = fit_simulation_propensity_honest_tree(sim.dat, 
        seed = seed.offset + 100000L + iter), SLE = fit_simulation_propensity_sle(sim.dat), kernel = fit_simulation_propensity_kernel(sim.dat, 
        pilot.fit = ps.pilot))
    ps <- ps.fit$probability
    results <- list(safe_simulation_method(estimate_naive_simulation(sim.dat, truth, threshold), "Naive", 
        truth), safe_simulation_method(estimate_sipw_simulation(sim.dat, ps, truth, threshold, "SIPW"), 
        "SIPW", truth), safe_simulation_method(estimate_elw_simulation(sim.dat, ps, truth, threshold, 
        "ELW-MAR"), "ELW-MAR", truth))
    if (ps.method == "PLE") {
        results[[length(results) + 1L]] <- safe_simulation_method(estimate_initial_qr_simulation(sim.dat, 
            truth, threshold), "G-est", truth)
    }
    if (include_mnar_sipw && (gamma0 != 0 || gamma1 != 0)) {
        results[[length(results) + 1L]] <- safe_simulation_method(estimate_lmz_mnar_simulation(sim.dat, 
            ps.fit, truth, threshold, "SIPW-LMZ"), "SIPW-MNAR", truth)
    }
    if (include_mnar_elw && (gamma0 != 0 || gamma1 != 0)) {
        results[[length(results) + 1L]] <- safe_simulation_method(estimate_lmz_mnar_simulation(sim.dat, 
            ps.fit, truth, threshold, "ELW-LMZ"), "ELW-MNAR", truth)
        if (include_mnar_lql) {
            results[[length(results) + 1L]] <- safe_simulation_method(estimate_lql_mnar_simulation(sim.dat, 
                truth, threshold), "ELW-LQL", truth)
            results[[length(results) + 1L]] <- safe_simulation_method(estimate_sipw_lql_mnar_simulation(sim.dat, 
                truth, threshold), "SIPW-LQL", truth)
        }
    }
    estimates <- simulation_bind_rows(lapply(results, `[[`, "estimates"))
    diagnostics <- simulation_bind_rows(lapply(results, `[[`, "diagnostics"))
    estimates$iter <- iter
    diagnostics$iter <- iter
    estimates$b <- b
    estimates$gamma0 <- gamma0
    estimates$gamma1 <- gamma1
    estimates$ps.method <- ps.method
    diagnostics$b <- b
    diagnostics$gamma0 <- gamma0
    diagnostics$gamma1 <- gamma1
    diagnostics$ps.method <- ps.method
    diagnostics$ps.auc <- ps.fit$auc
    diagnostics$ps.boundary.count <- sum(ps == 0 | ps == 1)
    diagnostics$ps.boundary.any <- as.numeric(any(ps == 0 | ps == 1))
    diagnostics$ps.fatal.boundary.count <- sum((sim.dat$T == 1 & ps == 0) | (sim.dat$T == 0 & ps == 1))
    diagnostics$ps.fatal.boundary.any <- as.numeric(any((sim.dat$T == 1 & ps == 0) | (sim.dat$T == 0 & 
        ps == 1)))
    observed.ps <- c(ps[sim.dat$T == 1], 1 - ps[sim.dat$T == 0])
    default.ps.min <- min(observed.ps)
    default.ps.p01 <- as.numeric(stats::quantile(observed.ps, 0.01, names = FALSE))
    default.ps.p99 <- as.numeric(stats::quantile(observed.ps, 0.99, names = FALSE))
    default.ps.max <- max(observed.ps)
    for (v in c("ps.min", "ps.p01", "ps.p99", "ps.max")) {
        if (!v %in% names(diagnostics)) 
            diagnostics[[v]] <- NA_real_
    }
    mar.weighting <- diagnostics$method %in% c("SIPW", "ELW-MAR")
    diagnostics$ps.min[mar.weighting & !is.finite(diagnostics$ps.min)] <- default.ps.min
    diagnostics$ps.p01[mar.weighting & !is.finite(diagnostics$ps.p01)] <- default.ps.p01
    diagnostics$ps.p99[mar.weighting & !is.finite(diagnostics$ps.p99)] <- default.ps.p99
    diagnostics$ps.max[mar.weighting & !is.finite(diagnostics$ps.max)] <- default.ps.max
    diagnostics$smoking.prevalence <- mean(sim.dat$T == 1)
    list(estimates = estimates, diagnostics = diagnostics)
}
}, envir = envir)
}

.application_main <- function(args = commandArgs(TRUE)) {
steps <- c('propensity','descriptive','distribution','heterogeneity','ps_robustness','weight_stress','msm','lmz','lmz_gamma')
if (!length(args) || args[1] == '--list') { cat('Usage: Rscript application.R STEP CLEAN_CSV OUTPUT_DIR --run\nSteps:', paste(steps, collapse=', '), '\n'); return(invisible(NULL)) }
if(length(args)<3 || !args[1]%in%steps) stop('Specify STEP CLEAN_CSV OUTPUT_DIR [--run]')
step <- args[1]; input <- normalizePath(args[2], mustWork=TRUE); out <- normalizePath(args[3], mustWork=FALSE)
if(!'--run'%in%args) {cat('PLAN ONLY:',step,'\nInput:',input,'\nOutput:',out,'\nAdd --run to execute.\n');return(invisible(NULL))}
dir.create(out, recursive=TRUE, showWarnings=FALSE)
target <- file.path(out,'csl_smoking_birthweight_complete_case.csv')
if(!file.exists(target)) {if(!file.copy(input,target)) stop('Cannot stage input')} else if(!identical(read.csv(input),read.csv(target))) stop('Output directory contains different input data')
Sys.setenv(ELW_TASK=step,CSL_SCRIPT_DIR=out,CSL_SOURCE_DIR=out,CSL_COMPLETE_OUTPUT_DIR=out,CSL_COMPLETE_OUTPUT=out,CSL_DISTRIBUTIONAL_OUTPUT_DIR=out,CSL_HETEROGENEITY_OUTPUT_DIR=out,CSL_MANUSCRIPT_DIR=file.path(out,'FIGURES'),LMZ_APP_OUTPUT=file.path(out,'lmz'),CSL_INCLUDE_EXTRA_FLEXIBLE_PS='true')
.application_load_base(environment())
# ---- 01_prepare_csl.R ----
requested_steps <- step
run_step <- function(x) identical(x, step)
## ---- clean ----
if (run_step("clean")) local({
  ## ============================================================
  ## CSL smoking-birthweight analysis: data cleaning
  ## ============================================================
  
  
  input_file <- Sys.getenv(
    "CSL_RAW_FILE",
    unset = "CSLdata171003.csv"
  )
  output_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = "CSL-COMPLETE-PROCESS"
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  missing_threshold <- 0.60
  dat <- read.csv(input_file, check.names = FALSE)
  
  ## The first unnamed column is a sequential export row number, not a verified
  ## maternal or patient identifier.
  if (names(dat)[1] == "") names(dat)[1] <- "row_id"
  
  ## CSL special codes representing unknown, missing, or unavailable values.
  special_missing_codes <- list(
    Momrace = c(8, 9),
    Marital = c(8, 9),
    Education = c(8, 9),
    Hxpreterm = c(8, 9),
    Hxstillbirth = c(8, 9),
    Diabetes = c(8, 9),
    Hxasthma = c(8, 9),
    HIV = c(8, 9),
    Preeclampsia = c(8, 9),
    AnteGDM = c(8, 9),
    AnteSGA = c(8, 9),
    AnteLGA = c(8, 9),
    Sex = c(3, 8, 9)
  )
  
  for (variable in intersect(names(special_missing_codes), names(dat))) {
    dat[[variable]][dat[[variable]] %in% special_missing_codes[[variable]]] <- NA
  }
  
  missing_before <- data.frame(
    variable = names(dat),
    missing_n = vapply(dat, function(x) sum(is.na(x)), integer(1)),
    missing_rate = vapply(dat, function(x) mean(is.na(x)), numeric(1))
  )
  missing_before$drop_over_threshold <- (
    missing_before$missing_rate > missing_threshold
  )
  
  dropped_variables <- missing_before$variable[
    missing_before$drop_over_threshold
  ]
  
  ## Primary total-effect analysis: use only variables measured before (or fixed
  ## before) smoking could affect pregnancy outcomes. Pregnancy complications,
  ## admission weight, gestational age, and antenatal SGA/LGA assessments are not
  ## included as propensity-score covariates because they may be downstream of
  ## smoking. High-missingness baseline variables are removed by the rule above.
  candidate_baseline_covariates <- c(
    "Momrace", "Marital", "Education", "Momage", "Parity",
    "Hxpreterm", "Hxstillbirth", "Height", "PreBMI", "Diabetes",
    "Hxasthma", "HIV", "Alcohol", "Sex"
  )
  analysis_covariates <- setdiff(
    intersect(candidate_baseline_covariates, names(dat)), dropped_variables
  )
  analysis_variables <- c("row_id", "Smoke", "Birthweight", analysis_covariates)
  missing_analysis_variables <- setdiff(analysis_variables, names(dat))
  if (length(missing_analysis_variables)) {
    stop("Missing primary analysis variables: ",
         paste(missing_analysis_variables, collapse = ", "))
  }
  dat_retained <- dat[analysis_variables]
  
  complete <- complete.cases(dat_retained)
  dat_complete <- dat_retained[complete, , drop = FALSE]
  
  if (!all(c("Smoke", "Birthweight") %in% names(dat_complete))) {
    stop("Smoke and Birthweight must remain after variable filtering.")
  }
  if (!all(dat_complete$Smoke %in% c(1, 2))) {
    stop("Unexpected Smoke coding after cleaning; expected values 1 and 2.")
  }
  
  ## Original CSL coding: 1 = no smoking, 2 = smoking.
  dat_complete$treatment_smoking <- as.integer(dat_complete$Smoke == 2)
  dat_complete$outcome_birthweight <- as.numeric(dat_complete$Birthweight)
  
  front <- c("row_id", "treatment_smoking", "outcome_birthweight")
  dat_complete <- dat_complete[
    c(front, setdiff(names(dat_complete), front))
  ]
  
  sample_flow <- data.frame(
    step = c(
      "Raw records",
      "After identifying variables with >60% missingness",
      "Complete cases for primary baseline-covariate analysis"
    ),
    n_records = c(nrow(dat), nrow(dat), nrow(dat_complete)),
    n_variables = c(ncol(dat), ncol(dat) - length(dropped_variables),
                    ncol(dat_complete)),
    n_smokers = c(
      sum(dat$Smoke == 2, na.rm = TRUE),
      sum(dat_retained$Smoke == 2, na.rm = TRUE),
      sum(dat_complete$treatment_smoking == 1)
    ),
    n_non_smokers = c(
      sum(dat$Smoke == 1, na.rm = TRUE),
      sum(dat_retained$Smoke == 1, na.rm = TRUE),
      sum(dat_complete$treatment_smoking == 0)
    )
  )
  
  write.csv(
    dat_complete,
    file.path(output_dir, "csl_smoking_birthweight_complete_case.csv"),
    row.names = FALSE
  )
  write.csv(
    missing_before,
    file.path(output_dir, "csl_missingness_report.csv"),
    row.names = FALSE
  )
  write.csv(
    sample_flow,
    file.path(output_dir, "csl_cleaning_sample_flow.csv"),
    row.names = FALSE
  )
  writeLines(
    dropped_variables,
    file.path(output_dir, "csl_variables_dropped_high_missingness.txt")
  )
  writeLines(
    analysis_covariates,
    file.path(output_dir, "csl_primary_covariates.txt")
  )
  writeLines(
    intersect(
      c("Preeclampsia", "AnteGDM", "AnteSGA", "AnteLGA",
        "Admweight", "BESTGA"),
      names(dat)
    ),
    file.path(output_dir, "csl_post_exposure_variables_excluded.txt")
  )
  
  cat("Cleaning complete.\n")
  cat("Raw sample:", nrow(dat), "\n")
  cat("Variables with >60% missingness:",
      paste(dropped_variables, collapse = ", "), "\n")
  cat("Primary propensity-score covariates:",
      paste(analysis_covariates, collapse = ", "), "\n")
  cat("Complete-case sample:", nrow(dat_complete), "\n")
  cat(
    "Smoking groups: T=0", sum(dat_complete$treatment_smoking == 0),
    "; T=1", sum(dat_complete$treatment_smoking == 1), "\n"
  )
})

## ---- propensity ----
if (run_step("propensity")) local({
  ## ============================================================
  ## CSL smoking-birthweight analysis: tuned propensity scores
  ## ============================================================
  
  
  output_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = "CSL-COMPLETE-PROCESS"
  )
  input_file <- file.path(output_dir, "csl_smoking_birthweight_complete_case.csv")
  seed <- 20260615
  n_folds <- as.integer(Sys.getenv("CSL_PS_CV_FOLDS", "5"))
  rf_num_trees <- as.integer(Sys.getenv("CSL_RF_NUM_TREES", "300"))
  rf_num_threads <- as.integer(Sys.getenv(
    "CSL_RF_NUM_THREADS",
    as.character(max(1L, parallel::detectCores(logical = FALSE) - 1L))
  ))
  include_gam <- tolower(Sys.getenv("CSL_INCLUDE_GAM_PS", "false")) %in%
    c("true", "1", "yes")
  include_random_forest <- tolower(Sys.getenv("CSL_INCLUDE_RANDOM_FOREST_PS", "true")) %in%
    c("true", "1", "yes")
  include_extra_flexible <- tolower(Sys.getenv("CSL_INCLUDE_EXTRA_FLEXIBLE_PS", "false")) %in%
    c("true", "1", "yes")
  
  required_packages <- c("ggplot2", "gbm", "ranger", "kernlab", "nnet")
  if (include_gam) {
    required_packages <- c(required_packages, "mgcv")
  }
  if (include_random_forest) {
    required_packages <- c(required_packages, "ranger")
  }
  if (include_extra_flexible) {
    required_packages <- c(required_packages, "kernlab", "nnet")
  }
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  }
  suppressPackageStartupMessages(library(kernlab))
  
  dat <- read.csv(input_file, check.names = FALSE)
  invisible(NULL) # Kernel helpers are included above.
  A <- as.integer(dat$treatment_smoking)
  excluded <- c(
    "row_id", "Smoke", "Birthweight",
    "treatment_smoking", "outcome_birthweight"
  )
  covariates <- setdiff(names(dat), excluded)
  ps_data <- dat[covariates]
  
  categorical_variables <- intersect(
    c(
      "Momrace", "Marital", "Education", "Parity", "Hxpreterm",
      "Hxstillbirth", "Diabetes", "Hxasthma", "HIV", "Alcohol",
      "Preeclampsia", "AnteGDM", "AnteSGA", "AnteLGA", "Sex"
    ),
    names(ps_data)
  )
  ps_data[categorical_variables] <- lapply(ps_data[categorical_variables], factor)
  varying <- vapply(ps_data, function(x) length(unique(x)) > 1L, logical(1))
  dropped_constant <- names(ps_data)[!varying]
  ps_data <- ps_data[varying]
  
  standardize_matrix <- function(x) {
    means <- colMeans(x)
    sds <- apply(x, 2, stats::sd)
    keep <- is.finite(sds) & sds > 0
    out <- sweep(x[, keep, drop = FALSE], 2, means[keep], "-")
    sweep(out, 2, sds[keep], "/")
  }
  
  X_raw <- model.matrix(~ . - 1, data = ps_data)
  storage.mode(X_raw) <- "double"
  X_main <- standardize_matrix(X_raw)
  
  continuous <- intersect(c("Momage", "Height", "PreBMI"), colnames(X_main))
  
  make_folds <- function(y, k, seed) {
    set.seed(seed)
    fold <- integer(length(y))
    for (level in sort(unique(y))) {
      idx <- sample(which(y == level))
      fold[idx] <- rep(seq_len(k), length.out = length(idx))
    }
    fold
  }
  
  validate_probability <- function(p, name = "predicted probability") {
    p <- as.numeric(p)
    if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
      stop(name, " must contain finite values in [0, 1].")
    }
    p
  }
  
  kernel_index_predict <- function(u_train, y_train, u_test,
                                   bandwidth_constant = 1.06) {
    u_train <- as.numeric(u_train)
    u_test <- as.numeric(u_test)
    bandwidth <- bandwidth_constant * stats::sd(u_train) *
      length(u_train)^(-1 / 5)
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      return(rep(mean(y_train), length(u_test)))
    }
    nbin <- min(1000L, max(100L, ceiling(sqrt(length(u_train)) * 4)))
    breaks <- seq(min(u_train), max(u_train), length.out = nbin + 1L)
    bin <- findInterval(u_train, breaks, all.inside = TRUE)
    bin_count <- as.numeric(tabulate(bin, nbins = nbin))
    sum_by_bin <- rowsum(y_train, bin, reorder = FALSE)
    bin_sum <- numeric(nbin)
    matched <- match(rownames(sum_by_bin), as.character(seq_len(nbin)))
    bin_sum[matched] <- as.numeric(sum_by_bin[, 1])
    centers <- (breaks[-1L] + breaks[-length(breaks)]) / 2
    keep <- bin_count > 0
    centers <- centers[keep]
    bin_count <- bin_count[keep]
    bin_sum <- bin_sum[keep]
    z <- outer(centers, centers, "-") / bandwidth
    kern <- exp(-0.5 * z^2)
    denom <- as.vector(kern %*% bin_count)
    numer <- as.vector(kern %*% bin_sum)
    fitted_centers <- numer / denom
    out <- stats::approx(
      x = centers, y = fitted_centers, xout = u_test,
      method = "linear", rule = 2, ties = "ordered"
    )$y
    out[!is.finite(out)] <- mean(y_train)
    validate_probability(out, "kernel predicted probability")
  }
  
  fit_predict <- function(method, train, test, param, fit_seed) {
    X_train <- X_main[train, , drop = FALSE]
    X_test <- X_main[test, , drop = FALSE]
    d_train <- data.frame(A = A[train], X_train, check.names = TRUE)
    d_test <- data.frame(X_test, check.names = TRUE)
  
    if (method == "logistic_expanded") {
      ## Legacy method key retained for cache compatibility. This is a
      ## main-effects logistic regression with no interactions or nonlinear terms.
      fit <- suppressWarnings(glm(A ~ ., d_train, family = binomial()))
      p <- predict(fit, d_test, type = "response")
  
    } else if (method == "kernel_single_index") {
      p <- pca_kernel_predict(
        X_raw[train, , drop = FALSE], A[train],
        X_raw[test, , drop = FALSE], n_components = 3L,
        grid_bins = 45L
      )$probability
  
    } else if (method == "gam") {
      smooth <- intersect(make.names(continuous), names(d_train))
      linear <- setdiff(names(d_train), c("A", smooth))
      rhs <- c(paste0("s(", smooth, ", k=", param$k, ")"), linear)
      form <- as.formula(paste("A ~", paste(rhs, collapse = " + ")))
      fit <- suppressWarnings(mgcv::gam(
        form, d_train, family = binomial(), method = "REML", select = TRUE
      ))
      p <- predict(fit, d_test, type = "response")
  
    } else if (method == "random_forest") {
      d_train$A <- factor(d_train$A, levels = c(0, 1))
      fit <- ranger::ranger(
        A ~ ., d_train, probability = TRUE, num.trees = rf_num_trees,
        mtry = max(1L, round(param$mtry_fraction * ncol(X_train))),
        min.node.size = param$min_node_size, seed = fit_seed,
        num.threads = rf_num_threads
      )
      p <- predict(fit, data = d_test)$predictions[, "1"]
  
    } else if (method == "gbm") {
      set.seed(fit_seed)
      fit <- gbm::gbm.fit(
        x = X_train, y = A[train], distribution = "bernoulli",
        n.trees = param$n_trees,
        interaction.depth = param$interaction_depth,
        n.minobsinnode = param$min_node_size,
        shrinkage = param$shrinkage,
        bag.fraction = param$bag_fraction,
        nTrain = nrow(X_train), keep.data = FALSE, verbose = FALSE
      )
      p <- as.numeric(predict(
        fit, newdata = X_test, n.trees = param$n_trees,
        type = "response"
      ))
  
    } else if (method == "svm_radial") {
      ## Fixed-size training subset keeps the radial SVM computationally simple;
      ## sampling is stratified and uses a fixed seed, with no tuning search.
      set.seed(fit_seed)
      max_svm_train <- min(10000L, nrow(d_train))
      svm_rows <- unlist(lapply(sort(unique(A[train])), function(level) {
        candidates <- which(A[train] == level)
        take <- max(1L, round(as.numeric(max_svm_train) *
          as.numeric(length(candidates)) / nrow(d_train)))
        sample(candidates, min(take, length(candidates)))
      }), use.names = FALSE)
      d_train <- d_train[svm_rows, , drop = FALSE]
      svm_y <- as.integer(as.character(d_train$A))
      d_train$A <- factor(d_train$A, levels = c(0, 1))
      fit <- kernlab::ksvm(
        A ~ ., d_train, type = "C-svc",
        kernel = kernlab::rbfdot(sigma = param$sigma), C = param$C,
        prob.model = FALSE, scaled = FALSE
      )
      score_train <- as.numeric(predict(fit, d_train, type = "decision"))
      score_test <- as.numeric(predict(fit, d_test, type = "decision"))
      calibration <- suppressWarnings(glm(
        svm_y ~ score_train, family = binomial()
      ))
      p <- as.numeric(predict(
        calibration, newdata = data.frame(score_train = score_test),
        type = "response"
      ))
  
    } else if (method == "neural_network") {
      set.seed(fit_seed)
      fit <- nnet::nnet(
        A ~ ., d_train, size = param$size, decay = param$decay,
        maxit = 800, trace = FALSE, MaxNWts = 20000
      )
      p <- as.vector(predict(fit, d_test, type = "raw"))
  
    } else {
      stop("Unknown method: ", method)
    }
    validate_probability(p)
  }
  
  log_loss <- function(y, p) {
    mean(ifelse(y == 1, -log(p), -log1p(-p)))
  }
  brier_score <- function(y, p) mean((y - p)^2)
  auc_rank <- function(y, p) {
    n1 <- sum(y == 1)
    n0 <- sum(y == 0)
    (sum(rank(p, ties.method = "average")[y == 1]) - n1 * (n1 + 1) / 2) / (n1 * n0)
  }
  
  calibration_stats <- function(y, p) {
    lp <- qlogis(p)
    if (any(!is.finite(lp))) {
      return(c(intercept = NA_real_, slope = NA_real_))
    }
    intercept <- suppressWarnings(coef(glm(y ~ 1, family = binomial(), offset = lp))[1])
    slope <- suppressWarnings(coef(glm(y ~ lp, family = binomial()))[2])
    c(intercept = unname(intercept), slope = unname(slope))
  }
  
  parameter_grids <- list(
    logistic_expanded = data.frame(model = "all prespecified main effects"),
    kernel_single_index = data.frame(
      pca_components = 3L,
      bandwidth = "multivariate normal-reference"
    ),
    gam = data.frame(k = c(4L, 6L, 8L)),
    random_forest = expand.grid(
      mtry_fraction = 0.25,
      min_node_size = 50L
    ),
    gbm = data.frame(
      n_trees = 1000L, interaction_depth = 3L,
      min_node_size = 50L, shrinkage = 0.02, bag_fraction = 0.8
    ),
    svm_radial = data.frame(C = 1, sigma = 1 / ncol(X_main)),
    neural_network = data.frame(size = 5L, decay = 0.10)
  )
  
  method_labels <- c(
    logistic_expanded = "PLE",
    kernel_single_index = "PCA-Kernel (k = 3)",
    gam = "Generalized additive model",
    random_forest = "Random forest",
    gbm = "GBM",
    svm_radial = "Radial SVM",
    neural_network = "ANN"
  )
  methods <- c(
    "logistic_expanded", "kernel_single_index", "random_forest",
    "svm_radial", "neural_network", "gbm"
  )
  fold <- make_folds(A, n_folds, seed)
  
  tuning_rows <- list()
  best_parameters <- list()
  oof <- matrix(NA_real_, nrow(dat), length(methods), dimnames = list(NULL, methods))
  ps <- oof
  
  for (j in seq_along(methods)) {
    method <- methods[j]
    grid <- parameter_grids[[method]]
    candidate_predictions <- vector("list", nrow(grid))
    message("Cross-validating ", method_labels[[method]], "...")
  
    for (g in seq_len(nrow(grid))) {
      pred <- rep(NA_real_, nrow(dat))
      param <- as.list(grid[g, , drop = FALSE])
      for (v in seq_len(n_folds)) {
        test <- which(fold == v)
        train <- which(fold != v)
        pred[test] <- fit_predict(method, train, test, param, seed + 1000 * j + 10 * g + v)
      }
      candidate_predictions[[g]] <- pred
      parameter_text <- paste(
        paste(names(param), unlist(param), sep = "="), collapse = "; "
      )
      tuning_rows[[length(tuning_rows) + 1L]] <- data.frame(
        method = method, candidate = g,
        log_loss = log_loss(A, pred), brier = brier_score(A, pred),
        parameters = parameter_text
      )
    }
  
    method_tuning <- do.call(rbind, tuning_rows)
    method_tuning <- method_tuning[method_tuning$method == method, , drop = FALSE]
    best_candidate <- method_tuning$candidate[which.min(method_tuning$log_loss)]
    best_parameters[[method]] <- as.list(grid[best_candidate, , drop = FALSE])
    oof[, method] <- candidate_predictions[[best_candidate]]
    ps[, method] <- fit_predict(
      method, seq_len(nrow(dat)), seq_len(nrow(dat)),
      best_parameters[[method]], seed + 10000 + j
    )
  }
  
  tuning_results <- do.call(rbind, tuning_rows)
  diagnostics <- do.call(rbind, lapply(methods, function(method) {
    p <- oof[, method]
    p0 <- p[A == 0]
    p1 <- p[A == 1]
    calibration <- calibration_stats(A, p)
    data.frame(
      method = method,
      label = method_labels[[method]],
      oof_log_loss = log_loss(A, p),
      oof_brier = brier_score(A, p),
      oof_auc = auc_rank(A, p),
      calibration_intercept = unname(calibration["intercept"]),
      calibration_slope = unname(calibration["slope"]),
      T0_max = max(p0), T0_q99 = unname(quantile(p0, 0.99)),
      T0_ge_090 = sum(p0 >= 0.90), T0_ge_095 = sum(p0 >= 0.95),
      T1_min = min(p1), T1_q01 = unname(quantile(p1, 0.01)),
      T1_le_010 = sum(p1 <= 0.10), T1_le_005 = sum(p1 <= 0.05)
    )
  }))
  rownames(diagnostics) <- NULL
  
  individual_scores <- data.frame(
    row_id = dat$row_id,
    treatment_smoking = A,
    outcome_birthweight = dat$outcome_birthweight,
    ps,
    check.names = FALSE
  )
  oof_scores <- data.frame(row_id = dat$row_id, A = A, oof, check.names = FALSE)
  
  plot_data <- do.call(rbind, lapply(methods, function(method) {
    data.frame(
      method = method_labels[[method]], score = ps[, method],
      treatment = ifelse(A == 1, "T=1: smoking", "T=0: no smoking")
    )
  }))
  plot_data$method <- factor(plot_data$method, method_labels[methods])
  plot <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(score, ggplot2::after_stat(density), fill = treatment, color = treatment)
  ) +
    ggplot2::geom_histogram(bins = 50, position = "identity", alpha = 0.42, linewidth = 0.2) +
    ggplot2::geom_vline(xintercept = c(0.05, 0.95), linetype = "dotted", linewidth = 0.3) +
    ggplot2::geom_vline(xintercept = c(0.10, 0.90), linetype = "dashed", linewidth = 0.3) +
    ggplot2::facet_wrap(~method, ncol = 2) +
    ggplot2::coord_cartesian(xlim = c(0, 1)) +
    ggplot2::labs(
      x = "Estimated probability of smoking", y = "Density", fill = NULL, color = NULL
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(legend.position = "top")
  
  write.csv(individual_scores, file.path(output_dir, "csl_smoking_propensity_scores_full_sample.csv"), row.names = FALSE)
  write.csv(oof_scores, file.path(output_dir, "csl_smoking_propensity_scores_oof.csv"), row.names = FALSE)
  write.csv(tuning_results, file.path(output_dir, "csl_smoking_propensity_tuning_results.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(output_dir, "csl_smoking_propensity_model_diagnostics.csv"), row.names = FALSE)
  writeLines(
    c(
      paste("Covariates:", paste(covariates, collapse = ", ")),
      paste("Dropped constant variables:", paste(dropped_constant, collapse = ", ")),
      paste("Cross-validation folds:", n_folds),
      paste("Random forest trees:", rf_num_trees),
      paste("Random forest threads:", rf_num_threads),
      "PLE: all prespecified pretreatment covariate main effects",
      "PCA--Kernel: first 3 standardized principal components",
      "PCA--Kernel bandwidth: 3-dimensional normal-reference rule",
      "GBM: 1,000 trees, interaction depth 3, minimum node size 50, shrinkage 0.02, bag fraction 0.8",
      "Random forest: 300 trees, mtry fraction 0.25, minimum node size 50",
      "Radial SVM: fixed C=1 and sigma=1/p; stratified training subset capped at 10,000; training-fold Platt calibration",
      "ANN: one hidden layer with 5 units and weight decay 0.10",
      "Cross-validation is diagnostic only; it does not select covariates or tuning parameters",
      capture.output(print(best_parameters))
    ),
    file.path(output_dir, "csl_smoking_propensity_model_specification.txt")
  )
  ggplot2::ggsave(file.path(output_dir, "csl_smoking_propensity_histograms.pdf"), plot, width = 8.4, height = 7.2)
  
  cat("Tuned propensity-score analysis complete.\n")
  print(diagnostics)
  print(best_parameters)
})

## ---- sample_table ----
if (run_step("sample_table")) local({
  ## Build the supplementary included-versus-excluded comparison for the CSL analysis.
  
  input_file <- Sys.getenv(
    "CSL_RAW_FILE",
    unset = "CSLdata171003.csv"
  )
  output_dir <- Sys.getenv(
    "CSL_TABLE_OUTPUT_DIR",
    unset = "output"
  )
  
  dat <- read.csv(input_file, check.names = FALSE)
  
  special_missing_codes <- list(
    Momrace = c(8, 9), Marital = c(8, 9), Hxpreterm = c(8, 9),
    Diabetes = c(8, 9), Hxasthma = c(8, 9), HIV = c(8, 9),
    Sex = c(3, 8, 9)
  )
  for (v in intersect(names(special_missing_codes), names(dat))) {
    dat[[v]][dat[[v]] %in% special_missing_codes[[v]]] <- NA
  }
  
  analysis_vars <- c(
    "Smoke", "Birthweight", "Momrace", "Marital", "Momage", "Parity",
    "Hxpreterm", "Height", "PreBMI", "Diabetes", "Hxasthma", "HIV",
    "Alcohol", "Sex"
  )
  stopifnot(nrow(dat) == 213514L, all(analysis_vars %in% names(dat)))
  included <- complete.cases(dat[analysis_vars])
  stopifnot(sum(included) == 118523L, sum(!included) == 94991L)
  
  groups <- list(Included = dat[included, , drop = FALSE],
                 Excluded = dat[!included, , drop = FALSE])
  
  fmt_n <- function(x) format(as.integer(x), big.mark = ",", scientific = FALSE)
  fmt_missing <- function(x, n) sprintf("%s (%.1f\\%%)", fmt_n(x), 100 * x / n)
  fmt_cont <- function(x, digits = 1, bigmark = "") {
    q <- unname(quantile(x, c(.25, .5, .75), na.rm = TRUE))
    paste0(formatC(q[2], format = "f", digits = digits, big.mark = bigmark), " (",
           formatC(q[1], format = "f", digits = digits, big.mark = bigmark), ", ",
           formatC(q[3], format = "f", digits = digits, big.mark = bigmark), ")")
  }
  fmt_level <- function(x, value) {
    den <- sum(!is.na(x)); num <- sum(x == value, na.rm = TRUE)
    sprintf("%s (%.1f\\%%)", fmt_n(num), 100 * num / den)
  }
  
  rows <- list()
  add_row <- function(label, variable, summary_fun) {
    vals <- lapply(groups, function(g) summary_fun(g[[variable]]))
    miss <- lapply(groups, function(g) fmt_missing(sum(is.na(g[[variable]])), nrow(g)))
    rows[[length(rows) + 1L]] <<- data.frame(
      characteristic = label,
      included_summary = vals[[1]], included_missing = miss[[1]],
      excluded_summary = vals[[2]], excluded_missing = miss[[2]],
      stringsAsFactors = FALSE
    )
  }
  add_level <- function(label, variable, value) {
    vals <- lapply(groups, function(g) fmt_level(g[[variable]], value))
    rows[[length(rows) + 1L]] <<- data.frame(
      characteristic = paste0("\\quad ", label),
      included_summary = vals[[1]], included_missing = "",
      excluded_summary = vals[[2]], excluded_missing = "",
      stringsAsFactors = FALSE
    )
  }
  
  add_row("Recorded smoking", "Smoke", function(x) fmt_level(x, 2))
  add_row("Birthweight, g", "Birthweight", function(x) fmt_cont(x, 0, ","))
  add_row("Maternal age, years", "Momage", function(x) fmt_cont(x, 1))
  add_row("Maternal height, m", "Height", function(x) fmt_cont(x, 2))
  add_row("Prepregnancy BMI, kg/m$^2$", "PreBMI", function(x) fmt_cont(x, 1))
  
  add_row("\\textit{Maternal race}", "Momrace", function(x) "")
  add_level("White, non-Hispanic", "Momrace", 1)
  add_level("Black, non-Hispanic", "Momrace", 2)
  add_level("Hispanic", "Momrace", 3)
  add_level("Asian", "Momrace", 4)
  add_level("American Indian/Alaska Native", "Momrace", 5)
  
  add_row("\\textit{Marital status}", "Marital", function(x) "")
  add_level("Married", "Marital", 1)
  add_level("Widowed", "Marital", 2)
  add_level("Nonmarried (excluding widowed)", "Marital", 3)
  
  add_row("Parous", "Parity", function(x) fmt_level(x, 1))
  add_row("History of preterm birth", "Hxpreterm", function(x) fmt_level(x, 2))
  add_row("Preexisting diabetes", "Diabetes", function(x) fmt_level(x, 2))
  add_row("History of asthma", "Hxasthma", function(x) fmt_level(x, 2))
  add_row("HIV", "HIV", function(x) fmt_level(x, 2))
  add_row("Alcohol use", "Alcohol", function(x) fmt_level(x, 2))
  add_row("Male infant", "Sex", function(x) fmt_level(x, 2))
  
  out <- do.call(rbind, rows)
  write.csv(out, file.path(output_dir, "csl_included_excluded_comparison.csv"),
            row.names = FALSE)
  
  body <- apply(out, 1, function(z) paste(z, collapse = " & "))
  body <- paste0(body, " \\\\")
  tex <- c(
    "\\begin{table*}[htbp]",
    "\\TableFont\\centering",
    "\\caption{Comparison of Complete-Case and Excluded Records in the Available CSL Smoking-Analysis Extract}",
    "\\label{tab:app_included_excluded}",
    "\\begin{tabular*}{\\textwidth}{@{\\extracolsep{\\fill}}lcccc}",
    "\\toprule",
    "\\textbf{Characteristic} & \\multicolumn{2}{c}{\\textbf{Included ($N=118{,}523$)}} & \\multicolumn{2}{c}{\\textbf{Excluded ($N=94{,}991$)}} \\\\",
    "\\cmidrule(lr){2-3}\\cmidrule(lr){4-5}",
    " & \\textbf{Observed summary} & \\textbf{Missing, $n$ (\\%)} & \\textbf{Observed summary} & \\textbf{Missing, $n$ (\\%)} \\\\",
    "\\midrule", body, "\\bottomrule", "\\end{tabular*}",
    "\\TableNote{The comparison is restricted to the 213,514 records available in the smoking-analysis extract. Continuous variables are summarized as median (first quartile, third quartile); categorical variables as count (percentage among records with that variable observed). Missing percentages use the group total as denominator. Special unknown or unavailable codes were treated as missing. The 15,154-delivery difference between the published CSL cohort and the available extract cannot be characterized because individual-level records and exclusion reasons were unavailable.}",
    "\\end{table*}"
  )
  writeLines(tex, file.path(output_dir, "csl_included_excluded_comparison_table.tex"))
  
  cat("Included:", sum(included), "\nExcluded:", sum(!included), "\n")
})
# ---- 02_analyze_csl.R ----
requested_steps <- step
run_step <- function(x) identical(x, step)
## ---- descriptive ----
if (run_step("descriptive")) local({
  ## ============================================================
  ## CSL application: smoking and birthweight distribution
  ## ============================================================
  
  
  script_dir <- Sys.getenv(
    "CSL_SCRIPT_DIR",
    unset = "."
  )
  output_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = file.path(script_dir, "CSL-COMPLETE-PROCESS")
  )
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  manuscript_dir <- if (basename(manuscript_root) == "FIGURES") {
    manuscript_root
  } else {
    file.path(manuscript_root, "FIGURES")
  }
  dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
  .application_load_base(environment())
  
  required_packages <- c(
    "dplyr", "tidyr", "ggplot2", "patchwork", "MASS", "writexl"
  )
  missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
  ]
  if (length(missing_packages)) {
    stop("Missing packages: ", paste(missing_packages, collapse = ", "))
  }
  suppressPackageStartupMessages({
    library(dplyr)
    library(tidyr)
    library(ggplot2)
    library(patchwork)
  })
  
  set.seed(20260615)
  dat <- read.csv(
    file.path(output_dir, "csl_smoking_birthweight_complete_case.csv"),
    check.names = FALSE
  )
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  
  excluded <- c(
    "row_id", "Smoke", "Birthweight",
    "treatment_smoking", "outcome_birthweight"
  )
  covariates <- setdiff(names(dat), excluded)
  factor_vars <- intersect(
    c(
      "Momrace", "Marital", "Education", "Parity", "Hxpreterm",
      "Hxstillbirth", "Diabetes", "Hxasthma", "HIV", "Alcohol",
      "Preeclampsia", "AnteGDM", "AnteSGA", "AnteLGA", "Sex"
    ),
    covariates
  )
  ps_data <- dat[covariates]
  ps_data[factor_vars] <- lapply(ps_data[factor_vars], factor)
  X_balance <- model.matrix(~ . - 1, data = ps_data)
  storage.mode(X_balance) <- "double"
  
  ## Descriptive summary
  continuous <- c(
    Birthweight = "Birthweight, g",
    Momage = "Maternal age, years",
    Height = "Maternal height, m",
    PreBMI = "Prepregnancy BMI, kg/m$^2$",
    Admweight = "Admission weight, kg",
    BESTGA = "Gestational age, weeks"
  )
  continuous <- continuous[names(continuous) %in% names(dat)]
  binary <- c(
    Parity = "Parous",
    Hxpreterm = "History of preterm birth",
    Hxstillbirth = "History of stillbirth",
    Diabetes = "Preexisting diabetes",
    Hxasthma = "History of asthma",
    HIV = "HIV",
    Alcohol = "Alcohol use",
    Preeclampsia = "Preeclampsia",
    AnteGDM = "Gestational diabetes",
    AnteSGA = "Antepartum SGA",
    AnteLGA = "Antepartum LGA",
    Sex = "Male infant"
  )
  binary <- binary[names(binary) %in% names(dat)]
  
  desc_cont <- lapply(names(continuous), function(v) {
    x <- dat[[v]]
    data.frame(
      characteristic = continuous[[v]],
      overall = format_median_iqr(x, ifelse(v == "Birthweight", 0, 1)),
      T0 = format_median_iqr(x[A == 0], ifelse(v == "Birthweight", 0, 1)),
      T1 = format_median_iqr(x[A == 1], ifelse(v == "Birthweight", 0, 1)),
      p_value = format_p(wilcox.test(x ~ A, exact = FALSE)$p.value)
    )
  })
  desc_bin <- lapply(names(binary), function(v) {
    x <- dat[[v]]
    positive <- if (v == "Parity") 1 else 2
    data.frame(
      characteristic = binary[[v]],
      overall = format_binary(x, positive),
      T0 = format_binary(x[A == 0], positive),
      T1 = format_binary(x[A == 1], positive),
      p_value = format_p(
        suppressWarnings(chisq.test(table(x, A), correct = FALSE)$p.value)
      )
    )
  })
  descriptive <- bind_rows(desc_cont, desc_bin)
  write.csv(
    descriptive, file.path(output_dir, "csl_smoking_descriptive_summary.csv"),
    row.names = FALSE
  )
  write_descriptive_latex(
    descriptive, file.path(manuscript_dir, "csl_smoking_descriptive_table.tex"),
    nrow(dat), sum(A == 0), sum(A == 1)
  )
  
  ## Distributional shape of observed birthweight
  shape_statistics <- function(x) {
    n <- length(x)
    centered <- x - mean(x)
    second_moment <- mean(centered^2)
    skewness <- mean(centered^3) / second_moment^(3 / 2)
    pearson_kurtosis <- mean(centered^4) / second_moment^2
    skewness_p <- 2 * pnorm(-abs(skewness / sqrt(6 / n)))
    kurtosis_p <- 2 * pnorm(-abs((pearson_kurtosis - 3) / sqrt(24 / n)))
    c(
      mean = mean(x), sd = sd(x), median = median(x),
      skewness = skewness, pearson_kurtosis = pearson_kurtosis,
      skewness_p = skewness_p, kurtosis_p = kurtosis_p, n = n
    )
  }
  
  shape_matrix <- cbind(
    overall = shape_statistics(Y),
    no_smoking = shape_statistics(Y[A == 0]),
    smoking = shape_statistics(Y[A == 1])
  )
  shape_summary <- data.frame(
    statistic = rownames(shape_matrix), shape_matrix,
    row.names = NULL, check.names = FALSE
  )
  write.csv(
    shape_summary,
    file.path(output_dir, "csl_smoking_birthweight_shape_summary.csv"),
    row.names = FALSE
  )
  
  shape_labels <- c(
    mean = "Mean, g", sd = "Standard deviation, g", median = "Median, g",
    skewness = "Skewness", pearson_kurtosis = "Pearson kurtosis",
    skewness_p = "Zero-skewness test ($p$-value)",
    kurtosis_p = "Normal-kurtosis test ($p$-value)",
    n = "Number of observations"
  )
  format_shape <- function(value, statistic) {
    if (statistic == "n") return(format(round(value), big.mark = ","))
    if (statistic %in% c("skewness_p", "kurtosis_p")) {
      return(ifelse(value < 0.001, "$<0.001$", sprintf("%.3f", value)))
    }
    if (statistic %in% c("skewness", "pearson_kurtosis")) {
      return(sprintf("%.3f", value))
    }
    formatC(value, format = "f", digits = 1, big.mark = ",")
  }
  format_shape_plot <- function(value, statistic) {
    if (statistic == "n") return(format(round(value), big.mark = ","))
    if (statistic %in% c("skewness_p", "kurtosis_p")) {
      return(ifelse(value < 0.001, "<0.001", sprintf("%.3f", value)))
    }
    if (statistic %in% c("skewness", "pearson_kurtosis")) {
      return(sprintf("%.3f", value))
    }
    formatC(value, format = "f", digits = 1, big.mark = ",")
  }
  shape_display <- shape_summary[shape_summary$statistic != "n", ]
  shape_body <- vapply(seq_len(nrow(shape_display)), function(i) {
    statistic <- shape_display$statistic[i]
    paste0(
      shape_labels[[statistic]], " & ",
      format_shape(shape_display$overall[i], statistic), " & ",
      format_shape(shape_display$no_smoking[i], statistic), " & ",
      format_shape(shape_display$smoking[i], statistic), " \\\\"
    )
  }, character(1))
  writeLines(
    c(
      "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccc}",
      "\\toprule",
      "\\textbf{Statistic} & \\textbf{Overall} & \\textbf{Observed $Y(0)$} & \\textbf{Observed $Y(1)$} \\\\",
      sprintf(
        "& \\textbf{N = %s} & \\textbf{N = %s} & \\textbf{N = %s} \\\\",
        format_shape(shape_matrix["n", "overall"], "n"),
        format_shape(shape_matrix["n", "no_smoking"], "n"),
        format_shape(shape_matrix["n", "smoking"], "n")
      ),
      "\\midrule", shape_body, "\\bottomrule", "\\end{tabular*}"
    ),
    file.path(manuscript_dir, "csl_smoking_birthweight_shape_table.tex")
  )
  
  hist_data <- data.frame(
    birthweight = Y,
    group = ifelse(A == 1, "Smoking", "No smoking")
  )
  hist_data$group <- factor(
    hist_data$group, levels = c("No smoking", "Smoking")
  )
  info_colors <- c(
    "No smoking" = "#35B779", "Smoking" = "#F2A900"
  )
  bin_width <- 100
  bin_breaks <- seq(
    floor(min(hist_data$birthweight, na.rm = TRUE) / bin_width) * bin_width,
    ceiling(max(hist_data$birthweight, na.rm = TRUE) / bin_width) * bin_width,
    by = bin_width
  )
  hist_bins <- hist_data |>
    mutate(
      bin = cut(
        birthweight, breaks = bin_breaks, include.lowest = TRUE, right = FALSE
      )
    ) |>
    count(group, bin, name = "count") |>
    complete(group, bin, fill = list(count = 0)) |>
    group_by(group) |>
    mutate(
      frequency = count / sum(count),
      bin_mid = bin_breaks[as.integer(bin)] + bin_width / 2
    ) |>
    ungroup()
  outcome_plot <- ggplot(
    hist_bins,
    aes(bin_mid, frequency, fill = group, color = group)
  ) +
    geom_col(
      width = bin_width, position = "identity", alpha = 0.38, linewidth = 0.2
    ) +
    geom_vline(xintercept = 2500, linetype = "dashed", linewidth = 0.4) +
    annotate(
      "text", x = 2500, y = Inf, label = "2,500 g",
      hjust = -0.15, vjust = 1.25, fontface = "bold", size = 3.2
    ) +
    scale_fill_manual(values = info_colors) +
    scale_color_manual(values = info_colors, guide = "none") +
    scale_x_continuous(labels = function(x) formatC(x, format = "d", big.mark = ",")) +
    scale_y_continuous(labels = function(x) paste0(round(100 * x), "%")) +
    labs(x = NULL, y = "Within-group frequency", fill = NULL) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      plot.margin = margin(2, 8, 2, 8)
    )
  
  format_shape_strip <- function(value, statistic) {
    if (statistic %in% c("mean", "median")) {
      return(format(round(value), big.mark = ","))
    }
    sprintf("%.2f", value)
  }
  strip_groups <- data.frame(
    group = c("overall", "no_smoking", "smoking"),
    group_label = c("OVERALL", "NO SMOKING", "SMOKING"),
    x_center = c(0.5, 1.5, 2.5),
    group_color = c("grey45", info_colors[["No smoking"]], info_colors[["Smoking"]])
  )
  strip_stats <- data.frame(
    statistic = c("mean", "median", "skewness", "pearson_kurtosis"),
    stat_label = c("Mean", "Median", "Skewness", "Kurtosis"),
    x_offset = c(-0.35, -0.12, 0.12, 0.35)
  )
  stats_panel <- merge(strip_groups, strip_stats, all = TRUE) |>
    mutate(
      value = vapply(
        seq_len(n()),
        function(i) shape_summary[
          shape_summary$statistic == statistic[i], group[i]
        ],
        numeric(1)
      ),
      label = vapply(
        seq_len(n()),
        function(i) format_shape_strip(value[i], statistic[i]),
        character(1)
      ),
      x = x_center + x_offset
    )
  p_panel <- bind_rows(
    strip_groups |>
      mutate(
        test = "Skew",
        y = 0.14,
        p_value = c(
          shape_summary$overall[shape_summary$statistic == "skewness_p"],
          shape_summary$no_smoking[shape_summary$statistic == "skewness_p"],
          shape_summary$smoking[shape_summary$statistic == "skewness_p"]
        )
      ),
    strip_groups |>
      mutate(
        test = "Kurt",
        y = 0.00,
        p_value = c(
          shape_summary$overall[shape_summary$statistic == "kurtosis_p"],
          shape_summary$no_smoking[shape_summary$statistic == "kurtosis_p"],
          shape_summary$smoking[shape_summary$statistic == "kurtosis_p"]
        )
      )
  ) |>
    mutate(
      label = ifelse(
        p_value < 0.001,
        paste0(test, " p < 0.001"),
        paste0(test, " p = ", sprintf("%.3f", p_value))
      )
    )
  stats_plot <- ggplot(
    stats_panel,
    aes(x = x, y = 0)
  ) +
    annotate("segment", x = 0, xend = 3, y = 0.96, yend = 0.96,
             color = "grey85", linewidth = 0.35) +
    annotate("text", x = 1.5, y = 1.10, label = "Birthweight (g)",
             color = "grey45", size = 2.85) +
    geom_segment(
      data = data.frame(x = c(1, 2)),
      aes(x = x, xend = x, y = 0.02, yend = 0.96),
      inherit.aes = FALSE, color = "grey86", linewidth = 0.35
    ) +
    geom_text(
      data = strip_groups,
      aes(x = x_center, y = 0.78, label = group_label, color = group_color),
      inherit.aes = FALSE, size = 2.85, show.legend = FALSE
    ) +
    scale_color_identity() +
    geom_text(aes(label = label), y = 0.54, fontface = "bold",
              color = "grey10", size = 3.45) +
    geom_text(aes(label = stat_label), y = 0.34, color = "grey45", size = 2.55) +
    geom_label(
      data = p_panel,
      aes(x = x_center, y = y, label = label),
      inherit.aes = FALSE, fill = "#FDE4E2", color = "#CC2F2F",
      linewidth = 0, label.padding = unit(0.15, "lines"),
      fontface = "bold", size = 2.35
    ) +
    coord_cartesian(xlim = c(0, 3), ylim = c(-0.12, 1.14), clip = "off") +
    theme_void(base_size = 11) +
    theme(
      plot.margin = margin(2, 4, 4, 4)
    )
  
  birthweight_infographic <- outcome_plot / stats_plot +
    plot_layout(heights = c(1.65, 1.25))
  ggsave(
    file.path(manuscript_dir, "csl_smoking_birthweight_histogram.pdf"),
    outcome_plot, width = 8.4, height = 4.4
  )
  ggsave(
    file.path(manuscript_dir, "csl_smoking_birthweight_infographic.pdf"),
    birthweight_infographic, width = 9.4, height = 5.2
  )
  
  ## Covariate overlap motivation for the introduction.  Among continuous baseline
  ## covariates, plot the two variables with the largest raw treated-control SMDs.
  continuous_overlap_vars <- intersect(c("Momage", "Height", "PreBMI"), names(dat))
  continuous_overlap_labels <- c(
    Momage = "Maternal age, years",
    Height = "Maternal height, m",
    PreBMI = "Prepregnancy BMI, kg/m^2"
  )
  covariate_overlap_summary <- bind_rows(lapply(continuous_overlap_vars, function(v) {
    x1 <- as.numeric(dat[[v]][A == 1])
    x0 <- as.numeric(dat[[v]][A == 0])
    pooled_sd <- sqrt((var(x1, na.rm = TRUE) + var(x0, na.rm = TRUE)) / 2)
    data.frame(
      variable = v,
      label = continuous_overlap_labels[[v]],
      mean_smoking = mean(x1, na.rm = TRUE),
      mean_no_smoking = mean(x0, na.rm = TRUE),
      smd = (mean(x1, na.rm = TRUE) - mean(x0, na.rm = TRUE)) / pooled_sd,
      ks = suppressWarnings(as.numeric(stats::ks.test(x1, x0)$statistic))
    )
  }))
  write.csv(
    covariate_overlap_summary,
    file.path(output_dir, "csl_smoking_continuous_covariate_overlap.csv"),
    row.names = FALSE
  )
  overlap_focus <- covariate_overlap_summary |>
    slice_max(order_by = abs(smd), n = 2, with_ties = FALSE) |>
    arrange(desc(abs(smd)))
  overlap_plot_data <- bind_rows(lapply(seq_len(nrow(overlap_focus)), function(i) {
    v <- overlap_focus$variable[[i]]
    data.frame(
      variable = v,
      label = overlap_focus$label[[i]],
      value = as.numeric(dat[[v]]),
      exposure_group = factor(ifelse(A == 1, "Smoking", "No smoking"),
                              levels = c("No smoking", "Smoking"))
    )
  }))
  overlap_plot_data$label <- factor(overlap_plot_data$label, levels = overlap_focus$label)
  range_stats <- overlap_plot_data |>
    group_by(label, exposure_group) |>
    summarize(
      minimum = min(value, na.rm = TRUE),
      q05 = quantile(value, 0.05, na.rm = TRUE),
      mean = mean(value, na.rm = TRUE),
      q95 = quantile(value, 0.95, na.rm = TRUE),
      maximum = max(value, na.rm = TRUE),
      .groups = "drop"
    ) |>
    mutate(
      range_label = sprintf(
        "%s range: %.0f-%.0f", .data$exposure_group, .data$minimum, .data$maximum
      )
    )
  annotation_stats <- range_stats |>
    select(label, exposure_group, range_label) |>
    tidyr::pivot_wider(names_from = exposure_group, values_from = range_label) |>
    left_join(overlap_focus |> select(label, smd, ks), by = "label") |>
    mutate(
      value = Inf,
      y = 1.5,
      label_text = sprintf(
        "%s\n%s\nSMD = %.3f; KS = %.3f",
        Smoking, `No smoking`, smd, ks
      )
    )
  overlap_plot <- ggplot(
    overlap_plot_data,
    aes(x = value, y = exposure_group, fill = exposure_group, color = exposure_group)
  ) +
    geom_violin(trim = FALSE, scale = "width", alpha = 0.30, linewidth = 0.55) +
    geom_segment(
      data = range_stats,
      aes(x = minimum, xend = maximum, y = exposure_group, yend = exposure_group,
          color = exposure_group),
      inherit.aes = FALSE, linewidth = 1.2
    ) +
    geom_segment(
      data = range_stats,
      aes(x = q05, xend = q95, y = exposure_group, yend = exposure_group),
      inherit.aes = FALSE, color = "grey15", linewidth = 3.0, alpha = 0.75
    ) +
    geom_point(
      data = range_stats,
      aes(x = mean, y = exposure_group),
      inherit.aes = FALSE, shape = 21, fill = "white", color = "grey15", size = 2.4,
      stroke = 0.6
    ) +
    scale_fill_manual(values = info_colors) +
    scale_color_manual(values = info_colors) +
    labs(
      x = NULL,
      y = NULL,
      fill = NULL, color = NULL
    ) +
    geom_label(
      data = annotation_stats,
      aes(x = value, y = y, label = label_text),
      inherit.aes = FALSE, hjust = 1.05, vjust = 0.5,
      linewidth = 0, fill = "white", color = "grey20", size = 2.45,
      label.padding = unit(0.12, "lines")
    ) +
    facet_wrap(~label, ncol = 1, scales = "free_x") +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "none",
      axis.text.y = element_text(color = "grey15"),
      strip.background = element_rect(fill = "grey95", color = "grey80", linewidth = 0.3),
      strip.text = element_text(face = "bold", size = 9),
      panel.spacing.y = unit(0.55, "lines")
    )
  ggsave(
    file.path(manuscript_dir, "csl_smoking_covariate_overlap_ranges.pdf"),
    overlap_plot, width = 5.2, height = 5.0
  )
  
  ## Read tuned full-sample propensity scores generated by script 02.
  score_file <- file.path(output_dir, "csl_smoking_propensity_scores_full_sample.csv")
  diagnostic_file <- file.path(output_dir, "csl_smoking_propensity_model_diagnostics.csv")
  if (!file.exists(score_file) || !file.exists(diagnostic_file)) {
    stop("Run the propensity step of application.R before this analysis.")
  }
  score_data <- read.csv(score_file, check.names = FALSE)
  if (!identical(as.character(score_data$row_id), as.character(dat$row_id))) {
    stop("Propensity-score file does not match the analysis data order.")
  }
  
  candidate_methods <- c(
    "logistic_expanded", "gbm", "kernel_single_index", "gam",
    "random_forest", "svm_radial", "neural_network"
  )
  labels <- c(
    logistic_expanded = "PLE",
    gbm = "GBM",
    kernel_single_index = "PCA-Kernel (k = 3)",
    gam = "GAM",
    random_forest = "Random forest",
    svm_radial = "Radial SVM",
    neural_network = "ANN"
  )
  methods <- intersect(candidate_methods, names(score_data))
  if (!length(methods)) stop("No propensity-score columns found.")
  ps <- as.matrix(score_data[methods])
  ps_diagnostics <- read.csv(diagnostic_file, check.names = FALSE)
  
  ps_plot_data <- bind_rows(lapply(methods, function(m) {
    data.frame(
      method = labels[[m]], score = ps[, m],
      group = ifelse(A == 1, "Smoking", "No smoking")
    )
  }))
  ps_plot_data$method <- factor(ps_plot_data$method, labels[methods])
  ps_bin_width <- 0.02
  ps_breaks <- seq(0, 1, by = ps_bin_width)
  ps_hist <- ps_plot_data |>
    mutate(bin = cut(score, breaks = ps_breaks, include.lowest = TRUE, right = FALSE)) |>
    count(method, group, bin, name = "count") |>
    complete(method, group, bin, fill = list(count = 0)) |>
    group_by(method, group) |>
    mutate(
      density = count / (sum(count) * ps_bin_width),
      bin_mid = ps_breaks[as.integer(bin)] + ps_bin_width / 2,
      signed_density = ifelse(group == "Smoking", -density, density),
      signed_density_sqrt = sign(signed_density) * sqrt(abs(signed_density))
    ) |>
    ungroup()
  ps_thresholds <- data.frame(
    x = c(0.05, 0.10, 0.90, 0.95),
    threshold_type = factor(
      c("Outer", "Inner", "Inner", "Outer"),
      levels = c("Outer", "Inner")
    )
  )
  ps_density_breaks <- c(0, 1, 5, 10, 20, 30)
  ps_ybreaks <- c(-rev(sqrt(ps_density_breaks[-1])), sqrt(ps_density_breaks))
  ps_ylabels <- c(rev(ps_density_breaks[-1]), ps_density_breaks)
  ps_ymax <- max(
    sqrt(30),
    ceiling(max(abs(ps_hist$signed_density_sqrt), na.rm = TRUE) * 10) / 10
  )
  ps_plot <- ggplot(
    ps_hist, aes(bin_mid, signed_density_sqrt, fill = group, color = group)
  ) +
    geom_col(width = ps_bin_width, position = "identity", alpha = 0.45, linewidth = 0.2) +
    geom_vline(
      data = ps_thresholds,
      aes(xintercept = x, linetype = threshold_type),
      inherit.aes = FALSE, color = "grey10", linewidth = 0.35
    ) +
    facet_wrap(~method, ncol = 2) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    scale_y_continuous(
      limits = c(-ps_ymax, ps_ymax),
      breaks = ps_ybreaks[abs(ps_ybreaks) <= ps_ymax],
      labels = ps_ylabels[abs(ps_ybreaks) <= ps_ymax]
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1")
    ) +
    scale_fill_manual(values = info_colors) +
    scale_color_manual(values = c("No smoking" = "#218A5A", "Smoking" = "#C78300")) +
    scale_linetype_manual(values = c(Outer = "dotted", Inner = "longdash"), guide = "none") +
    labs(
      x = "Estimated probability of smoking",
      y = "Density (sqrt scale)",
      fill = NULL, color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      panel.spacing.x = unit(1.05, "lines"),
      axis.text = element_text(color = "grey25", size = 8.8),
      strip.background = element_rect(fill = "grey94", color = "grey80", linewidth = 0.3),
      strip.text = element_text(face = "bold")
    )
  ggsave(
    file.path(manuscript_dir, "csl_smoking_propensity_histograms.pdf"),
    ps_plot, width = 8.4, height = 7.2
  )
  
  gbm_overlap_plot <- ggplot(
    filter(ps_hist, method == "GBM"),
    aes(bin_mid, signed_density_sqrt, fill = group, color = group)
  ) +
    geom_col(width = ps_bin_width, position = "identity", alpha = 0.45, linewidth = 0.2) +
    geom_vline(
      data = ps_thresholds,
      aes(xintercept = x, linetype = threshold_type),
      inherit.aes = FALSE, color = "grey10", linewidth = 0.35
    ) +
    coord_cartesian(xlim = c(0, 1), clip = "off") +
    scale_y_continuous(
      limits = c(-ps_ymax, ps_ymax),
      breaks = ps_ybreaks[abs(ps_ybreaks) <= ps_ymax],
      labels = ps_ylabels[abs(ps_ybreaks) <= ps_ymax]
    ) +
    scale_x_continuous(
      breaks = seq(0, 1, by = 0.25),
      labels = c("0", "0.25", "0.50", "0.75", "1")
    ) +
    scale_fill_manual(values = info_colors) +
    scale_color_manual(values = c("No smoking" = "#218A5A", "Smoking" = "#C78300")) +
    scale_linetype_manual(values = c(Outer = "dotted", Inner = "longdash"), guide = "none") +
    labs(
      x = "Estimated probability of smoking",
      y = "Density (sqrt scale)",
      fill = NULL, color = NULL
    ) +
    theme_classic(base_size = 11) +
    theme(
      legend.position = "top",
      axis.text = element_text(color = "grey25", size = 8.8)
    )
  ggsave(
    file.path(manuscript_dir, "csl_smoking_propensity_gbm.pdf"),
    gbm_overlap_plot, width = 6.6, height = 4.6
  )
  
  ## MAR ELW distributional effects
  probs <- c(0.10, 0.25, 0.50, 0.75)
  quantile_labels <- c("10th", "25th", "50th", "75th")
  make_rows <- function(method, q1, q0) {
    data.frame(
      method = method, quantile = quantile_labels,
      tau = probs, q1 = q1, q0 = q0,
      difference = q1 - q0
    )
  }
  naive_q1 <- quantile(Y[A == 1], probs, type = 1, names = FALSE)
  naive_q0 <- quantile(Y[A == 0], probs, type = 1, names = FALSE)
  
  mar_results <- lapply(methods, function(m) {
    fit1 <- elw_group(A, ps[, m], "treated")
    fit0 <- elw_group(A, ps[, m], "control")
    list(
      fit1 = fit1, fit0 = fit0,
      q1 = weighted_quantile(Y[A == 1], fit1$weight, probs),
      q0 = weighted_quantile(Y[A == 0], fit0$weight, probs),
      lbw1 = weighted_probability(Y[A == 1], fit1$weight, function(x) x < 2500),
      lbw0 = weighted_probability(Y[A == 0], fit0$weight, function(x) x < 2500)
    )
  })
  names(mar_results) <- methods
  
  quantile_results <- bind_rows(
    make_rows("Naive", naive_q1, naive_q0),
    bind_rows(lapply(methods, function(m) {
      make_rows(paste0("ELW-", labels[[m]]),
                mar_results[[m]]$q1, mar_results[[m]]$q0)
    }))
  )
  weight_results <- bind_rows(lapply(methods, function(m) {
    r <- mar_results[[m]]
    n_group <- c(sum(A == 1), sum(A == 0))
    data.frame(
      method = labels[[m]], group = c("Smoking", "No smoking"),
      n = n_group,
      p99_weight = c(
        unname(quantile(r$fit1$weight, 0.99)),
        unname(quantile(r$fit0$weight, 0.99))
      ),
      max_weight = c(max(r$fit1$weight), max(r$fit0$weight)),
      max_relative_weight = n_group * c(max(r$fit1$weight), max(r$fit0$weight)),
      effective_sample_size = c(
        sum(r$fit1$weight)^2 / sum(r$fit1$weight^2),
        sum(r$fit0$weight)^2 / sum(r$fit0$weight^2)
      ),
      relative_effective_sample_size = c(
        sum(r$fit1$weight)^2 / sum(r$fit1$weight^2) / n_group[1],
        sum(r$fit0$weight)^2 / sum(r$fit0$weight^2) / n_group[2]
      )
    )
  }))
  
  weighted_mean <- function(x, w) sum(x * w) / sum(w)
  balance_results <- bind_rows(lapply(methods, function(m) {
    w1 <- mar_results[[m]]$fit1$weight
    w0 <- mar_results[[m]]$fit0$weight
    smd <- vapply(seq_len(ncol(X_balance)), function(j) {
      x1 <- X_balance[A == 1, j]
      x0 <- X_balance[A == 0, j]
      pooled_sd <- sqrt((var(x1) + var(x0)) / 2)
      if (!is.finite(pooled_sd) || pooled_sd == 0) return(0)
      (weighted_mean(x1, w1) - weighted_mean(x0, w0)) / pooled_sd
    }, numeric(1))
    data.frame(
      method = m,
      max_abs_smd = max(abs(smd)),
      mean_abs_smd = mean(abs(smd)),
      n_abs_smd_over_010 = sum(abs(smd) > 0.10)
    )
  }))
  
  weight_model_summary <- bind_rows(lapply(methods, function(m) {
    d <- weight_results[weight_results$method == labels[[m]], ]
    data.frame(
      method = m,
      maximum_elw_weight = max(d$max_weight),
      maximum_relative_elw_weight = max(d$max_relative_weight),
      minimum_relative_ess = min(d$relative_effective_sample_size)
    )
  }))
  model_evaluation <- ps_diagnostics %>%
    dplyr::select(
      method, label, oof_log_loss, oof_brier, oof_auc,
      calibration_intercept, calibration_slope,
      T0_max, T0_ge_090, T1_min, T1_le_010
    ) %>%
    left_join(balance_results, by = "method") %>%
    left_join(weight_model_summary, by = "method")
  
  eligible <- with(model_evaluation, max_abs_smd <= 0.10)
  candidate_rows <- which(eligible)
  if (length(candidate_rows)) {
    recommended_method <- model_evaluation$method[
      candidate_rows[which.min(model_evaluation$oof_log_loss[candidate_rows])]
    ]
  } else {
    ## When no model meets the conventional balance benchmark, prioritize the
    ## smallest residual imbalance and report the remaining diagnostics.
    recommended_method <- model_evaluation$method[
      which.min(model_evaluation$max_abs_smd)
    ]
  }
  ## The manuscript primary model follows the diagnostics unless explicitly
  ## overridden for a sensitivity run.
  primary_method <- Sys.getenv("CSL_PRIMARY_PS", recommended_method)
  if (!primary_method %in% methods) stop("Unknown CSL_PRIMARY_PS: ", primary_method)
  model_evaluation$recommended_primary <- model_evaluation$method == primary_method
  
  write.csv(
    model_evaluation,
    file.path(output_dir, "csl_smoking_propensity_model_evaluation.csv"),
    row.names = FALSE
  )
  ps_table_rows <- vapply(seq_len(nrow(model_evaluation)), function(i) {
    z <- model_evaluation[i, ]
    display_label <- labels[[z$method]]
    selected <- z$method == recommended_method
    label <- if (selected) {
      paste0("\\textbf{", display_label, "}")
    } else display_label
    sprintf(
      "%s & %.4f & %.4f & %.3f & %.3f & %.2f & %.3f & %s \\\\",
      label, z$oof_log_loss, z$oof_brier, z$calibration_slope,
      z$max_abs_smd, z$maximum_relative_elw_weight, z$minimum_relative_ess,
      ifelse(selected, "Yes", "")
    )
  }, character(1))
  writeLines(c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lccccccc}",
    "\\toprule",
    "\\textbf{Model} & \\textbf{Log loss} & \\textbf{Brier} & \\textbf{Cal. slope} & \\textbf{Max. $|\\mathrm{SMD}|$} & \\textbf{Max. $n_t w$} & \\textbf{Min. ESS/$n_t$} & \\textbf{Selected} \\\\",
    "\\midrule", ps_table_rows, "\\bottomrule", "\\end{tabular*}"
  ), file.path(manuscript_dir, "csl_smoking_ps_evaluation_table.tex"))
  lbw_results <- bind_rows(lapply(methods, function(m) {
    r <- mar_results[[m]]
    data.frame(
      method = paste0("ELW-", labels[[m]]), risk_smoking = r$lbw1,
      risk_no_smoking = r$lbw0,
      risk_difference = r$lbw1 - r$lbw0
    )
  }))
  
  write.csv(quantile_results, file.path(output_dir, "csl_smoking_elw_quantiles.csv"), row.names = FALSE)
  write.csv(weight_results, file.path(output_dir, "csl_smoking_elw_weights.csv"), row.names = FALSE)
  write.csv(lbw_results, file.path(output_dir, "csl_smoking_elw_low_birthweight.csv"), row.names = FALSE)
  write_quantile_latex(
    quantile_results,
    file.path(manuscript_dir, "csl_smoking_elw_quantile_point_table.tex")
  )
  write_dte_latex(
    lbw_results,
    file.path(manuscript_dir, "csl_smoking_elw_dte2500_table.tex")
  )
  
  ## Weight-stability comparison with normalized inverse probability weighting
  weight_comparison <- bind_rows(lapply(methods, function(m) {
    p <- ps[, m]
    ipw1 <- ((A == 1) / p)[A == 1]
    ipw0 <- ((A == 0) / (1 - p))[A == 0]
    ipw1 <- ipw1 / sum(ipw1)
    ipw0 <- ipw0 / sum(ipw0)
    data.frame(
      propensity_model = labels[[m]],
      estimator = rep(c("ELW", "Normalized IPW"), each = 2),
      group = rep(c("Smoking", "No smoking"), 2),
      effective_sample_size = c(
        effective_sample_size(mar_results[[m]]$fit1$weight),
        effective_sample_size(mar_results[[m]]$fit0$weight),
        effective_sample_size(ipw1),
        effective_sample_size(ipw0)
      ),
      maximum_normalized_weight = c(
        max(mar_results[[m]]$fit1$weight),
        max(mar_results[[m]]$fit0$weight),
        max(ipw1),
        max(ipw0)
      )
    )
  }))
  write.csv(
    weight_comparison,
    file.path(output_dir, "csl_smoking_elw_vs_ipw_weights.csv"),
    row.names = FALSE
  )
  
  primary <- mar_results[[primary_method]]
  ## The marginal-sensitivity and Shao--Wang analyses are run separately by
  ## The sensitivity analyses use these saved estimates as inputs.
  writexl::write_xlsx(
    list(
      descriptive = descriptive,
      outcome_shape = shape_summary,
      propensity_diagnostics = ps_diagnostics,
      propensity_model_evaluation = model_evaluation,
      quantiles = quantile_results,
      low_birthweight = lbw_results,
      weights = weight_results
    ),
    file.path(output_dir, "csl_smoking_birthweight_analysis_results.xlsx")
  )
  
  cat("Analysis complete. N =", nrow(dat), "\n")
  cat("Diagnostic-preferred propensity model:", recommended_method, "\n")
  cat("Primary propensity model:", primary_method, "\n")
  print(ps_diagnostics)
  print(model_evaluation)
  print(quantile_results)
  print(lbw_results)
  print(weight_results)
})

## ---- distribution ----
if (run_step("distribution")) local({
  ## Distributional inference for the primary ELW analysis.
  ## Implements the exchangeable bootstrap in Algorithm 2 of the ELW paper.
  
  
  source_dir <- Sys.getenv(
    "CSL_SCRIPT_DIR",
    unset = "."
  )
  input_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = file.path(source_dir, "CSL-COMPLETE-PROCESS")
  )
  output_dir <- Sys.getenv(
    "CSL_DISTRIBUTIONAL_OUTPUT_DIR",
    unset = input_dir
  )
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  manuscript_dir <- if (basename(manuscript_root) == "FIGURES") {
    manuscript_root
  } else {
    file.path(manuscript_root, "FIGURES")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
  primary_ps_method <- Sys.getenv("CSL_PRIMARY_PS_METHOD", "gbm")
  if (!requireNamespace("gbm", quietly = TRUE)) {
    stop("Package 'gbm' is required for the primary GBM analysis.")
  }
  gbm_spec <- list(
    n_trees = 1000L,
    interaction_depth = 3L,
    min_node_size = 50L,
    shrinkage = 0.02,
    bag_fraction = 0.8
  )
  ## Matches the full-sample GBM fit in 01_prepare_csl.R.
  gbm_point_seed <- 20270621L
  kernel_bandwidth_constant <- as.numeric(
    Sys.getenv("CSL_KERNEL_BW_CONSTANT", "1.06")
  )
  if (!is.finite(kernel_bandwidth_constant) || kernel_bandwidth_constant <= 0) {
    kernel_bandwidth_constant <- 1.06
  }
  kernel_implementation <- "pca3_grid_product_normal_reference_v2"
  invisible(NULL) # Kernel helpers are included above.
  
  dat <- read.csv(
    file.path(input_dir, "csl_smoking_birthweight_complete_case.csv"),
    check.names = FALSE
  )
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  N <- length(Y)
  
  excluded <- c(
    "row_id", "Smoke", "Birthweight",
    "treatment_smoking", "outcome_birthweight"
  )
  covariates <- setdiff(names(dat), excluded)
  factor_vars <- intersect(
    c(
      "Momrace", "Marital", "Education", "Parity", "Hxpreterm",
      "Hxstillbirth", "Diabetes", "Hxasthma", "HIV", "Alcohol",
      "Preeclampsia", "AnteGDM", "AnteSGA", "AnteLGA", "Sex"
    ),
    covariates
  )
  ps_data <- dat[covariates]
  ps_data[factor_vars] <- lapply(ps_data[factor_vars], factor)
  
  standardize_matrix <- function(x) {
    means <- colMeans(x)
    sds <- apply(x, 2, stats::sd)
    keep <- is.finite(sds) & sds > 0
    out <- sweep(x[, keep, drop = FALSE], 2, means[keep], "-")
    sweep(out, 2, sds[keep], "/")
  }
  
  X_raw <- model.matrix(~ . - 1, data = ps_data)
  storage.mode(X_raw) <- "double"
  X_main <- standardize_matrix(X_raw)
  design <- data.frame(A = A, X_main, check.names = TRUE)
  
  validate_probability <- function(p, name = "propensity score") {
    p <- as.numeric(p)
    if (any(!is.finite(p)) || any(p < 0 | p > 1)) {
      stop(name, " must contain finite values in [0, 1].")
    }
    p
  }
  
  kernel_index_predict <- function(u_train, y_train, u_test, weights,
                                   bandwidth_constant = 1.06) {
    u_train <- as.numeric(u_train)
    u_test <- as.numeric(u_test)
    weights <- as.numeric(weights)
    bandwidth <- bandwidth_constant * stats::sd(u_train) *
      length(u_train)^(-1 / 5)
    if (!is.finite(bandwidth) || bandwidth <= 0) {
      return(rep(stats::weighted.mean(y_train, weights), length(u_test)))
    }
    nbin <- min(1000L, max(100L, ceiling(sqrt(length(u_train)) * 4)))
    breaks <- seq(min(u_train), max(u_train), length.out = nbin + 1L)
    bin <- findInterval(u_train, breaks, all.inside = TRUE)
    bin_weight <- as.numeric(rowsum(weights, bin, reorder = FALSE)[, 1])
    bin_sum <- as.numeric(rowsum(weights * y_train, bin, reorder = FALSE)[, 1])
    present <- as.integer(rownames(rowsum(weights, bin, reorder = FALSE)))
    weight_by_bin <- numeric(nbin)
    sum_by_bin <- numeric(nbin)
    weight_by_bin[present] <- bin_weight
    sum_by_bin[present] <- bin_sum
    centers <- (breaks[-1L] + breaks[-length(breaks)]) / 2
    keep <- weight_by_bin > 0
    centers <- centers[keep]
    weight_by_bin <- weight_by_bin[keep]
    sum_by_bin <- sum_by_bin[keep]
    z <- outer(centers, centers, "-") / bandwidth
    kern <- exp(-0.5 * z^2)
    denom <- as.vector(kern %*% weight_by_bin)
    numer <- as.vector(kern %*% sum_by_bin)
    fitted_centers <- numer / denom
    out <- stats::approx(
      x = centers, y = fitted_centers, xout = u_test,
      method = "linear", rule = 2, ties = "ordered"
    )$y
    out[!is.finite(out)] <- stats::weighted.mean(y_train, weights)
    validate_probability(out, "kernel propensity score")
  }
  
  fit_propensity <- function(multiplier = rep(1, N),
                             method = primary_ps_method,
                             fit_seed = NULL) {
    if (method == "logistic_expanded") {
      ## Legacy method key retained for cache compatibility. This is a
      ## main-effects logistic regression with no interactions or nonlinear terms.
      fit <- suppressWarnings(glm(
        A ~ ., data = design, family = binomial(), weights = multiplier
      ))
      return(validate_probability(fitted(fit)))
    }
    if (method == "kernel_single_index") {
      return(pca_kernel_predict(
        X_raw, A, weights = multiplier, n_components = 3L,
        grid_bins = 45L
      )$probability)
    }
    if (method == "gbm") {
      if (!is.null(fit_seed)) set.seed(fit_seed)
      fit <- gbm::gbm.fit(
        x = X_main, y = A, distribution = "bernoulli", w = multiplier,
        n.trees = gbm_spec$n_trees,
        interaction.depth = gbm_spec$interaction_depth,
        n.minobsinnode = gbm_spec$min_node_size,
        shrinkage = gbm_spec$shrinkage,
        bag.fraction = gbm_spec$bag_fraction,
        nTrain = N, keep.data = FALSE, verbose = FALSE
      )
      return(validate_probability(as.numeric(predict(
        fit, newdata = X_main, n.trees = gbm_spec$n_trees,
        type = "response"
      )), "GBM propensity score"))
    }
    stop("Unknown CSL_PRIMARY_PS_METHOD: ", method)
  }
  
  elw_side <- function(pi_hat, treatment_level, multiplier = rep(1, N),
                       tol = 1e-10) {
    observed <- A == treatment_level
    e <- if (treatment_level == 1L) pi_hat else 1 - pi_hat
    total <- sum(multiplier)
    n1 <- sum(multiplier * A)
    r1 <- n1 / total
    r <- if (treatment_level == 1L) r1 else 1 - r1
    e_obs <- e[observed]
    m_obs <- multiplier[observed]
    lower <- min(e_obs)
    upper <- r + (1 - r) * lower
    root_fun <- function(alpha) {
      denom <- 1 + (1 - r) * (e_obs - alpha) / (r * (1 - alpha))
      sum(m_obs * (e_obs - alpha) / denom)
    }
    eps <- .Machine$double.eps
    lo <- lower
    hi <- upper - eps
    alpha <- uniroot(root_fun, c(lo, hi), tol = tol)$root
    denom <- 1 + (1 - r) * (e_obs - alpha) / (r * (1 - alpha))
    raw_weight <- m_obs / denom
    list(
      y = Y[observed],
      weight = raw_weight / sum(raw_weight),
      alpha = alpha
    )
  }
  
  weighted_cdf <- function(y, weight, grid) {
    ord <- order(y)
    approx(
      x = c(-Inf, y[ord], Inf),
      y = c(0, cumsum(weight[ord]), 1),
      xout = grid, method = "constant", f = 0, rule = 2,
      ties = "ordered"
    )$y
  }
  
  weighted_quantile <- function(y, weight, probs) {
    ord <- order(y)
    y <- y[ord]
    cw <- cumsum(weight[ord]) / sum(weight)
    vapply(probs, function(p) y[which(cw >= p)[1]], numeric(1))
  }
  
  grid <- sort(unique(Y))
  tau <- c(0.10, 0.25, 0.50, 0.75)
  intervals <- data.frame(
    label = c(
      "$(-\\infty,1{,}500)$",
      "$[1{,}500,2{,}500)$",
      "$[2{,}500,4{,}000)$",
      "$[4{,}000,\\infty)$"
    ),
    lower = c(-Inf, 1500, 2500, 4000),
    upper = c(1500, 2500, 4000, Inf),
    stringsAsFactors = FALSE
  )
  
  weighted_interval <- function(y, weight, lower, upper) {
    sum(weight[y >= lower & y < upper])
  }
  
  summarize_fit <- function(pi_hat, multiplier = rep(1, N)) {
    fit1 <- elw_side(pi_hat, 1L, multiplier)
    fit0 <- elw_side(pi_hat, 0L, multiplier)
    F1 <- weighted_cdf(fit1$y, fit1$weight, grid)
    F0 <- weighted_cdf(fit0$y, fit0$weight, grid)
    qte <- weighted_quantile(fit1$y, fit1$weight, tau) -
      weighted_quantile(fit0$y, fit0$weight, tau)
    pte <- vapply(seq_len(nrow(intervals)), function(j) {
      p1 <- weighted_interval(
        fit1$y, fit1$weight, intervals$lower[j], intervals$upper[j]
      )
      p0 <- weighted_interval(
        fit0$y, fit0$weight, intervals$lower[j], intervals$upper[j]
      )
      p1 - p0
    }, numeric(1))
    c(F1 = F1, F0 = F0, qte = qte, pte = pte)
  }
  
  set.seed(20260618)
  pi_hat <- fit_propensity(fit_seed = if (primary_ps_method == "gbm") {
    gbm_point_seed
  } else {
    NULL
  })
  point <- summarize_fit(pi_hat)
  n_grid <- length(grid)
  idx_F1 <- seq_len(n_grid)
  idx_F0 <- n_grid + seq_len(n_grid)
  idx_qte <- 2 * n_grid + seq_along(tau)
  idx_pte <- 2 * n_grid + length(tau) + seq_len(nrow(intervals))
  D_point <- point[idx_F1] - point[idx_F0]
  
  B <- as.integer(Sys.getenv("CSL_DISTRIBUTIONAL_BOOTSTRAP_B", "500"))
  detected_cores <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected_cores)) {
    detected_cores <- parallel::detectCores(logical = TRUE)
  }
  if (!is.finite(detected_cores)) {
    detected_cores <- suppressWarnings(as.integer(tryCatch(
      system("getconf _NPROCESSORS_ONLN", intern = TRUE)[1],
      error = function(e) NA_character_
    )))
  }
  if (!is.finite(detected_cores)) detected_cores <- 2L
  bootstrap_cores <- suppressWarnings(as.integer(Sys.getenv(
    "CSL_DISTRIBUTIONAL_BOOTSTRAP_CORES",
    as.character(max(1L, detected_cores - 1L))
  )))
  if (!is.finite(bootstrap_cores) || bootstrap_cores < 1L) bootstrap_cores <- 1L
  bootstrap_cores <- min(B, bootstrap_cores)
  bootstrap_file <- file.path(
    output_dir,
    paste0(
      "distributional_bootstrap_draws_", primary_ps_method, "_",
      if (primary_ps_method == "gbm") "fixed_spec_v1" else kernel_implementation,
      "_B", B, ".rds"
    )
  )
  reuse_bootstrap <- file.exists(bootstrap_file) &&
    !tolower(Sys.getenv("CSL_FORCE_BOOTSTRAP", "false")) %in% c("true", "1", "yes")
  if (reuse_bootstrap) {
    bootstrap_draws <- readRDS(bootstrap_file)
    reuse_bootstrap <- identical(bootstrap_draws$primary_ps_method, primary_ps_method)
    if (reuse_bootstrap && primary_ps_method == "gbm") {
      reuse_bootstrap <- isTRUE(all.equal(bootstrap_draws$gbm_spec, gbm_spec))
    } else if (reuse_bootstrap) {
      reuse_bootstrap <- isTRUE(all.equal(
        bootstrap_draws$kernel_bandwidth_constant,
        kernel_bandwidth_constant
      ))
    }
  }
  if (reuse_bootstrap) {
    point <- bootstrap_draws$point
    boot <- bootstrap_draws$boot
    D_point <- point[idx_F1] - point[idx_F0]
    if (!is.null(bootstrap_draws$primary_ps_method)) {
      primary_ps_method <- bootstrap_draws$primary_ps_method
    }
  } else {
    bootstrap_one <- function(b) {
      set.seed(20260618L + b)
      wb <- rexp(N)
      pi_b <- fit_propensity(wb)
      summarize_fit(pi_b, wb)
    }
    batch_size <- suppressWarnings(as.integer(Sys.getenv("CSL_BOOTSTRAP_BATCH_SIZE", "25")))
    if (!is.finite(batch_size) || batch_size < 1L) batch_size <- 25L
    partial_file <- paste0(bootstrap_file, ".partial.rds")
    boot <- matrix(NA_real_, B, length(point))
    completed <- 0L
    if (file.exists(partial_file)) {
      partial <- readRDS(partial_file)
      compatible <- identical(partial$B, B) &&
        isTRUE(all.equal(partial$gbm_spec, gbm_spec)) &&
        identical(partial$primary_ps_method, primary_ps_method)
      if (compatible) {
        boot <- partial$boot
        completed <- max(c(0L, which(rowSums(is.finite(boot)) == ncol(boot))))
        message("Resuming after ", completed, " completed bootstrap replicates.")
      }
    }
    message("Running ", B, " distributional bootstrap replicates in batches of ",
            batch_size, " on ", bootstrap_cores, " core(s).")
    while (completed < B) {
      ids <- seq.int(completed + 1L, min(B, completed + batch_size))
      values <- if (.Platform$OS.type != "windows" && bootstrap_cores > 1L) {
        parallel::mclapply(ids, bootstrap_one, mc.cores = min(bootstrap_cores, length(ids)),
                           mc.preschedule = TRUE, mc.set.seed = FALSE)
      } else lapply(ids, bootstrap_one)
      boot[ids, ] <- do.call(rbind, values)
      completed <- max(ids)
      saveRDS(list(B = B, gbm_spec = gbm_spec,
                   primary_ps_method = primary_ps_method, boot = boot), partial_file)
      message("Completed bootstrap replicate ", completed, "/", B)
    }
    message("Completed ", B, " distributional bootstrap replicates.")
    saveRDS(
      list(
        N = N, B = B, grid = grid, tau = tau, intervals = intervals,
        primary_ps_method = primary_ps_method,
        gbm_spec = gbm_spec,
        gbm_point_seed = gbm_point_seed,
        kernel_bandwidth_constant = kernel_bandwidth_constant,
        kernel_implementation = kernel_implementation,
        bootstrap_cores = bootstrap_cores,
        point = point, boot = boot,
        idx_F1 = idx_F1, idx_F0 = idx_F0, idx_qte = idx_qte, idx_pte = idx_pte
      ),
      bootstrap_file
    )
    unlink(partial_file)
  }
  
  ## Reconstruct interval probabilities from the cached CDF draws so that changes
  ## to clinical interval boundaries do not require refitting the bootstrap models.
  cdf_below <- function(F, cutoff) {
    if (is.infinite(cutoff) && cutoff < 0) return(0)
    if (is.infinite(cutoff) && cutoff > 0) return(1)
    idx <- which(grid < cutoff)
    if (length(idx)) F[max(idx)] else 0
  }
  pte_from_cdfs <- function(F1, F0) {
    vapply(seq_len(nrow(intervals)), function(j) {
      p1 <- cdf_below(F1, intervals$upper[j]) -
        cdf_below(F1, intervals$lower[j])
      p0 <- cdf_below(F0, intervals$upper[j]) -
        cdf_below(F0, intervals$lower[j])
      p1 - p0
    }, numeric(1))
  }
  point[idx_pte] <- pte_from_cdfs(point[idx_F1], point[idx_F0])
  boot[, idx_pte] <- t(vapply(seq_len(nrow(boot)), function(b) {
    pte_from_cdfs(boot[b, idx_F1], boot[b, idx_F0])
  }, numeric(nrow(intervals))))
  
  D_boot <- boot[, idx_F1, drop = FALSE] - boot[, idx_F0, drop = FALSE]
  centered_D <- sweep(D_boot, 2, D_point, "-")
  
  ## The ATE is recovered exactly from each discrete counterfactual CDF, so it
  ## uses the same bootstrap draws as the DTE, QTE, and PTE analyses.
  mean_from_cdf <- function(F) sum(grid * diff(c(0, F)))
  mean1_point <- mean_from_cdf(point[idx_F1])
  mean0_point <- mean_from_cdf(point[idx_F0])
  mean1_boot <- apply(boot[, idx_F1, drop = FALSE], 1, mean_from_cdf)
  mean0_boot <- apply(boot[, idx_F0, drop = FALSE], 1, mean_from_cdf)
  ate_boot <- mean1_boot - mean0_boot
  ate_ci <- quantile(ate_boot, c(0.025, 0.975), names = FALSE)
  ate_estimate <- mean1_point - mean0_point
  ate_se <- sd(ate_boot)
  ate_results <- data.frame(
    mean1 = mean1_point,
    mean0 = mean0_point,
    estimate = ate_estimate,
    se = ate_se,
    ci_lower = ate_ci[1],
    ci_upper = ate_ci[2],
    p_value = mean(abs(ate_boot - ate_estimate) >= abs(ate_estimate))
  )
  
  T_eq <- sqrt(N) * max(abs(D_point))
  T_eq_boot <- sqrt(N) * apply(abs(centered_D), 1, max)
  p_eq <- mean(T_eq_boot >= T_eq)
  
  ## One-sided harmful-effect test: H0 F1 <= F0 for all y.
  T_harm <- sqrt(N) * max(D_point)
  T_harm_boot <- sqrt(N) * apply(centered_D, 1, max)
  p_harm <- mean(T_harm_boot >= T_harm)
  
  tests <- data.frame(
    hypothesis = c(
      "F1(y) = F0(y) for all y",
      "F1(y) <= F0(y) for all y"
    ),
    statistic = c(T_eq, T_harm),
    p_value = c(p_eq, p_harm)
  )
  
  qte_boot <- boot[, idx_qte, drop = FALSE]
  qte_se <- apply(qte_boot, 2, sd)
  qte_ci <- t(apply(qte_boot, 2, quantile, c(0.025, 0.975), names = FALSE))
  qte_results <- data.frame(
    tau = tau, estimate = point[idx_qte], se = qte_se,
    ci_lower = qte_ci[, 1], ci_upper = qte_ci[, 2],
    p_value = vapply(seq_along(tau), function(j) {
      mean(abs(qte_boot[, j] - point[idx_qte][j]) >= abs(point[idx_qte][j]))
    }, numeric(1))
  )
  
  pte_boot <- boot[, idx_pte, drop = FALSE]
  pte_se <- apply(pte_boot, 2, sd)
  pte_ci <- t(apply(pte_boot, 2, quantile, c(0.025, 0.975), names = FALSE))
  interval_probability <- function(F, lower, upper) {
    cdf_below(F, upper) - cdf_below(F, lower)
  }
  pte_p1 <- vapply(seq_len(nrow(intervals)), function(j) {
    interval_probability(point[idx_F1], intervals$lower[j], intervals$upper[j])
  }, numeric(1))
  pte_p0 <- vapply(seq_len(nrow(intervals)), function(j) {
    interval_probability(point[idx_F0], intervals$lower[j], intervals$upper[j])
  }, numeric(1))
  pte_results <- data.frame(
    interval = intervals$label, lower = intervals$lower, upper = intervals$upper,
    P1 = pte_p1, P0 = pte_p0, estimate = point[idx_pte], se = pte_se,
    ci_lower = pte_ci[, 1], ci_upper = pte_ci[, 2],
    p_value = vapply(seq_len(nrow(intervals)), function(j) {
      mean(abs(pte_boot[, j] - point[idx_pte][j]) >= abs(point[idx_pte][j]))
    }, numeric(1))
  )
  
  ## Familywise inference for the four QTE and four PTE hypotheses. The absolute
  ## studentized maximum is already two-sided, so its Monte Carlo p-value is not
  ## multiplied by two.
  family_estimate <- c(qte_results$estimate, pte_results$estimate)
  family_boot <- cbind(qte_boot, pte_boot)
  family_se <- apply(family_boot, 2, sd)
  family_max_t <- apply(
    abs(sweep(sweep(family_boot, 2, family_estimate, "-"), 2, family_se, "/")),
    1, max
  )
  family_observed_t <- abs(family_estimate / family_se)
  family_adjusted_p <- vapply(family_observed_t, function(z) {
    sum(family_max_t >= z) / B
  }, numeric(1))
  family_results <- data.frame(
    estimand = c(paste0("QTE ", tau), intervals$label),
    estimate = family_estimate,
    se = family_se,
    absolute_t = family_observed_t,
    max_t_adjusted_p = family_adjusted_p
  )
  
  ## Pointwise and studentized uniform 95% bands use the IQR-scaled bootstrap
  ## empirical-process standard deviation. The raw pointwise scale can jump sharply
  ## between adjacent unique birthweight thresholds in sparse tails, so the DTE
  ## curve uses a light smoothing of log(scale). To avoid boundary smoothing
  ## artifacts in the studentized supremum, the stabilized scale is not allowed
  ## to fall below the raw IQR scale. Scalar QTE/PTE inference above remains
  ## based on the unsmoothed bootstrap draws.
  D_scale_raw <- apply(sqrt(N) * centered_D, 2, function(x) {
    diff(quantile(x, c(0.25, 0.75), names = FALSE)) /
      diff(qnorm(c(0.25, 0.75)))
  })
  D_scale_raw[D_scale_raw < 1e-8] <- NA_real_
  smooth_scale <- function(x, scale, df = 80) {
    keep <- is.finite(scale) & scale > 0
    if (sum(keep) < 10L) return(scale)
    fit <- stats::smooth.spline(
      x = x[keep], y = log(scale[keep]),
      df = min(df, sum(keep) - 1L)
    )
    out <- exp(stats::predict(fit, x = x)$y)
    out[!is.finite(out)] <- scale[!is.finite(out)]
    out
  }
  D_scale_smooth <- smooth_scale(grid, D_scale_raw)
  D_scale <- D_scale_smooth
  raw_keep <- is.finite(D_scale_raw) & D_scale_raw > 0
  D_scale[raw_keep] <- pmax(D_scale_smooth[raw_keep], D_scale_raw[raw_keep])
  t_boot <- apply(abs(sweep(sqrt(N) * centered_D, 2, D_scale, "/")), 1, max, na.rm = TRUE)
  critical <- quantile(t_boot, 0.95, names = FALSE)
  dte_results <- data.frame(
    y = grid, F1 = point[idx_F1], F0 = point[idx_F0], DTE = D_point,
    D_scale_raw = D_scale_raw, D_scale_smooth = D_scale_smooth, D_scale = D_scale,
    pointwise_lower = D_point - qnorm(0.975) * D_scale / sqrt(N),
    pointwise_upper = D_point + qnorm(0.975) * D_scale / sqrt(N),
    uniform_lower = D_point - critical * D_scale / sqrt(N),
    uniform_upper = D_point + critical * D_scale / sqrt(N)
  )
  
  write.csv(tests, file.path(output_dir, "distribution_tests.csv"), row.names = FALSE)
  write.csv(ate_results, file.path(output_dir, "ate_inference.csv"), row.names = FALSE)
  write.csv(qte_results, file.path(output_dir, "qte_inference.csv"), row.names = FALSE)
  write.csv(pte_results, file.path(output_dir, "pte_inference.csv"), row.names = FALSE)
  write.csv(
    family_results,
    file.path(output_dir, "qte_pte_max_t_familywise_inference.csv"),
    row.names = FALSE
  )
  write.csv(dte_results, file.path(output_dir, "dte_curve.csv"), row.names = FALSE)
  
  format_p <- function(p) ifelse(p == 0, "$<0.002$", sprintf("%.3f", p))
  format_tau <- function(x) sprintf("%.2f", x)
  format_integer <- function(x) format(
    round(x), big.mark = ",", scientific = FALSE, trim = TRUE
  )
  
  quantile_file <- file.path(input_dir, "csl_smoking_elw_quantiles.csv")
  if (file.exists(quantile_file)) {
    quantile_point <- read.csv(quantile_file, check.names = FALSE)
    target_method <- switch(
      primary_ps_method,
      logistic_expanded = "ELW-PLE",
      kernel_single_index = "ELW-Kernel",
      gbm = "ELW-GBM",
      stop("No quantile-table method label for ", primary_ps_method)
    )
    quantile_point <- quantile_point[quantile_point$method == target_method, ]
    qte_table <- merge(qte_results, quantile_point, by = "tau", all.x = TRUE)
    qte_table <- qte_table[order(qte_table$tau), ]
  } else {
    qte_table <- data.frame(
      tau = qte_results$tau, q1 = NA_real_, q0 = NA_real_,
      difference = qte_results$estimate
    )
    qte_table <- cbind(qte_table, qte_results[setdiff(names(qte_results), "tau")])
  }
  qte_lines <- vapply(seq_len(nrow(qte_results)), function(i) {
    z <- qte_table[i, ]
    sprintf(
      "%s & %s & %s & %s & %.1f & [%s, %s] & %s \\\\",
      format_tau(z$tau), format_integer(z$q1), format_integer(z$q0),
      format_integer(z$estimate), z$se, format_integer(z$ci_lower),
      format_integer(z$ci_upper), format_p(z$p_value)
    )
  }, character(1))
  writeLines(c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}ccccccc}",
    "\\toprule",
    "\\textbf{$\\tau$} & \\textbf{$\\hat q_1(\\tau)$} & \\textbf{$\\hat q_0(\\tau)$} & \\textbf{QTE} & \\textbf{SE} & \\textbf{95\\% CI} & \\textbf{$p$-value} \\\\",
    "\\midrule", qte_lines, "\\bottomrule", "\\end{tabular*}"
  ), file.path(manuscript_dir, "csl_smoking_elw_quantile_table.tex"))
  writeLines(c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}ccccccc}",
    "\\toprule",
    "\\textbf{$\\tau$} & \\textbf{$\\hat q_1(\\tau)$} & \\textbf{$\\hat q_0(\\tau)$} & \\textbf{QTE} & \\textbf{SE} & \\textbf{95\\% CI} & \\textbf{$p$-value} \\\\",
    "\\midrule", qte_lines, "\\bottomrule", "\\end{tabular*}"
  ), file.path(manuscript_dir, "csl_smoking_qte_inference_table.tex"))
  
  pte_lines <- vapply(seq_len(nrow(pte_results)), function(i) {
    z <- pte_results[i, ]
    sprintf(
      "%s & %.3f & %.3f & %.3f & %.3f & [%.3f, %.3f] & %s \\\\",
      z$interval, z$P1, z$P0, z$estimate, z$se, z$ci_lower, z$ci_upper,
      format_p(z$p_value)
    )
  }, character(1))
  writeLines(c(
    "\\begin{tabular*}{\\linewidth}{@{\\extracolsep{\\fill}}lcccccc}",
    "\\toprule",
    "\\textbf{Birthweight interval (g)} & \\textbf{$\\hat P_1$} & \\textbf{$\\hat P_0$} & \\textbf{PTE} & \\textbf{SE} & \\textbf{95\\% CI} & \\textbf{$p$-value} \\\\",
    "\\midrule", pte_lines, "\\bottomrule", "\\end{tabular*}"
  ), file.path(manuscript_dir, "csl_smoking_pte_inference_table.tex"))
  
  pdf(file.path(manuscript_dir, "csl_smoking_elw_dte_curve.pdf"), width = 7.4, height = 4.8)
  par(mar = c(4.2, 4.4, 0.8, 0.8), las = 1)
  dte_uniform_col <- "#D8EAF7"
  dte_pointwise_col <- "#7FB3D5"
  dte_line_col <- "#08306B"
  plot(
    grid, D_point, type = "n", xlab = "Birthweight threshold (g)",
    ylab = expression(hat(F)[1](y) - hat(F)[0](y)),
    xlim = quantile(Y, c(0.01, 0.99)),
    ylim = range(c(
      dte_results$uniform_lower, dte_results$uniform_upper,
      dte_results$pointwise_lower, dte_results$pointwise_upper
    ), na.rm = TRUE)
  )
  polygon(
    c(grid, rev(grid)),
    c(dte_results$uniform_lower, rev(dte_results$uniform_upper)),
    col = dte_uniform_col, border = NA
  )
  polygon(
    c(grid, rev(grid)),
    c(dte_results$pointwise_lower, rev(dte_results$pointwise_upper)),
    col = dte_pointwise_col, border = NA
  )
  lines(grid, D_point, col = dte_line_col, lwd = 1.8)
  abline(h = 0, lty = 2, lwd = 0.8, col = "gray35")
  abline(v = 2500, lty = 3, lwd = 0.8, col = "gray35")
  legend(
    "topleft", c("DTE", "95% pointwise CB", "95% uniform CB"),
    lty = c(1, NA, NA), lwd = c(1.6, NA, NA), pch = c(NA, 15, 15),
    col = c(dte_line_col, dte_pointwise_col, dte_uniform_col),
    pt.cex = c(NA, 2.2, 2.2), bty = "n"
  )
  dev.off()
  
  pdf(file.path(manuscript_dir, "csl_smoking_elw_distribution_curves.pdf"),
      width = 9.0, height = 4.1)
  par(mfrow = c(1, 2), mar = c(4.2, 4.6, 2.0, 1.0), las = 1)
  xlim_cdf <- c(500, 5000)
  xlim_dte <- xlim_cdf
  cdf_grid <- c(0, grid)
  F0_plot <- c(0, point[idx_F0])
  F1_plot <- c(0, point[idx_F1])
  plot(
    cdf_grid, F0_plot, type = "s", col = "#35B779", lwd = 1.8,
    xlab = "Birthweight (g)", ylab = "Estimated CDF",
    xlim = xlim_cdf, ylim = c(-0.04, 1.04), xaxs = "i", yaxs = "i",
    main = "Estimated counterfactual CDFs"
  )
  lines(cdf_grid, F1_plot, type = "s", col = "#F2A900", lwd = 1.8)
  abline(v = 2500, lty = 3, lwd = 0.8, col = "gray35")
  legend(
    "topleft", c(expression(hat(F)[0](y)), expression(hat(F)[1](y))),
    lty = 1, lwd = 1.8, col = c("#35B779", "#F2A900"), bty = "n",
    cex = 0.88
  )
  plot(
    grid, D_point, type = "n", xlab = "Birthweight threshold (g)",
    ylab = expression(hat(F)[1](y) - hat(F)[0](y)),
    xlim = xlim_dte, xaxs = "i",
    ylim = range(c(
      dte_results$uniform_lower, dte_results$uniform_upper,
      dte_results$pointwise_lower, dte_results$pointwise_upper
    ), na.rm = TRUE),
    main = "Distributional treatment effect"
  )
  polygon(
    c(grid, rev(grid)),
    c(dte_results$uniform_lower, rev(dte_results$uniform_upper)),
    col = dte_uniform_col, border = NA
  )
  polygon(
    c(grid, rev(grid)),
    c(dte_results$pointwise_lower, rev(dte_results$pointwise_upper)),
    col = dte_pointwise_col, border = NA
  )
  lines(grid, D_point, col = dte_line_col, lwd = 1.8)
  abline(h = 0, lty = 2, lwd = 0.8, col = "gray35")
  abline(v = 2500, lty = 3, lwd = 0.8, col = "gray35")
  legend(
    "topleft", c("DTE", "95% pointwise CB", "95% uniform CB"),
    lty = c(1, NA, NA), lwd = c(1.6, NA, NA), pch = c(NA, 15, 15),
    col = c(dte_line_col, dte_pointwise_col, dte_uniform_col),
    pt.cex = c(NA, 2.0, 2.0), bty = "n", cex = 0.82
  )
  dev.off()
  
  saveRDS(
    list(
      primary_ps_method = primary_ps_method,
      tests = tests, qte = qte_results, pte = pte_results, dte = dte_results
    ),
    file.path(output_dir, "distributional_inference.rds")
  )
  
  print(tests)
  print(ate_results)
  print(qte_results)
  print(pte_results)
})

## ---- heterogeneity ----
if (run_step("heterogeneity")) local({
  ## Prespecified heterogeneity analyses by infant sex and maternal age.
  ## One global GBM fit is shared by all subgroups within each bootstrap draw.
  
  
  source_dir <- Sys.getenv("CSL_SCRIPT_DIR", unset = ".")
  input_dir <- Sys.getenv("CSL_COMPLETE_OUTPUT_DIR", unset = file.path(source_dir, "CSL-COMPLETE-PROCESS"))
  output_dir <- Sys.getenv("CSL_HETEROGENEITY_OUTPUT_DIR", unset = input_dir)
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  figure_dir <- if (basename(manuscript_root) == "FIGURES") manuscript_root else file.path(manuscript_root, "FIGURES")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  
  if (!requireNamespace("gbm", quietly = TRUE)) stop("Package 'gbm' is required.")
  gbm_spec <- list(
    n_trees = 1000L, interaction_depth = 3L, min_node_size = 50L,
    shrinkage = 0.02, bag_fraction = 0.8
  )
  
  dat <- read.csv(file.path(input_dir, "csl_smoking_birthweight_complete_case.csv"), check.names = FALSE)
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  N <- length(Y)
  
  excluded <- c("row_id", "Smoke", "Birthweight", "treatment_smoking", "outcome_birthweight")
  covariates <- setdiff(names(dat), excluded)
  factor_vars <- intersect(
    c("Momrace", "Marital", "Education", "Parity", "Hxpreterm", "Hxstillbirth",
      "Diabetes", "Hxasthma", "HIV", "Alcohol", "Preeclampsia", "AnteGDM",
      "AnteSGA", "AnteLGA", "Sex"),
    covariates
  )
  ps_data <- dat[covariates]
  ps_data[factor_vars] <- lapply(ps_data[factor_vars], factor)
  X_raw <- model.matrix(~ . - 1, data = ps_data)
  storage.mode(X_raw) <- "double"
  means <- colMeans(X_raw)
  sds <- apply(X_raw, 2, stats::sd)
  keep <- is.finite(sds) & sds > 0
  X_main <- sweep(X_raw[, keep, drop = FALSE], 2, means[keep], "-")
  X_main <- sweep(X_main, 2, sds[keep], "/")
  
  fit_propensity <- function(multiplier = rep(1, N), seed = 20270621L) {
    set.seed(seed)
    fit <- gbm::gbm.fit(
      x = X_main, y = A, distribution = "bernoulli", w = multiplier,
      n.trees = gbm_spec$n_trees,
      interaction.depth = gbm_spec$interaction_depth,
      n.minobsinnode = gbm_spec$min_node_size,
      shrinkage = gbm_spec$shrinkage,
      bag.fraction = gbm_spec$bag_fraction,
      nTrain = N, keep.data = FALSE, verbose = FALSE
    )
    p <- as.numeric(predict(fit, newdata = X_main, n.trees = gbm_spec$n_trees, type = "response"))
    if (any(!is.finite(p)) || any(p < 0 | p > 1)) stop("Invalid GBM propensity score.")
    p
  }
  
  grid <- sort(unique(Y))
  plot_grid <- seq(500, 5000, by = 10)
  intervals <- data.frame(
    lower = c(-Inf, 1500, 2500, 4000), upper = c(1500, 2500, 4000, Inf),
    label = c("(-Inf, 1,500)", "[1,500, 2,500)", "[2,500, 4,000)", "[4,000, Inf)"),
    stringsAsFactors = FALSE
  )
  
  group_sets <- list(
    sex = list(
      labels = c("Female", "Male"),
      index = list(dat$Sex == 1, dat$Sex == 2)
    ),
    age = list(
      labels = c("<20", "20--34", ">=35"),
      index = list(dat$Momage < 20, dat$Momage >= 20 & dat$Momage < 35, dat$Momage >= 35)
    )
  )
  
  group_signature <- lapply(group_sets, function(gs) {
    list(labels = gs$labels, sizes = vapply(gs$index, sum, numeric(1)))
  })
  
  weighted_cdf <- function(y, weight, at) {
    ord <- order(y)
    approx(
      x = c(-Inf, y[ord], Inf), y = c(0, cumsum(weight[ord]), 1),
      xout = at, method = "constant", f = 0, rule = 2, ties = "ordered"
    )$y
  }
  
  elw_side <- function(pi_hat, treatment_level, subgroup, multiplier, tol = 1e-10) {
    eligible <- subgroup
    observed <- eligible & A == treatment_level
    e <- if (treatment_level == 1L) pi_hat else 1 - pi_hat
    total <- sum(multiplier[eligible])
    r1 <- sum(multiplier[eligible] * A[eligible]) / total
    r <- if (treatment_level == 1L) r1 else 1 - r1
    e_obs <- e[observed]
    m_obs <- multiplier[observed]
    lower <- min(e_obs)
    upper <- r + (1 - r) * lower
    root_fun <- function(alpha) {
      denom <- 1 + (1 - r) * (e_obs - alpha) / (r * (1 - alpha))
      sum(m_obs * (e_obs - alpha) / denom)
    }
    alpha <- uniroot(root_fun, c(lower, upper - .Machine$double.eps), tol = tol)$root
    denom <- 1 + (1 - r) * (e_obs - alpha) / (r * (1 - alpha))
    raw <- m_obs / denom
    list(y = Y[observed], weight = raw / sum(raw))
  }
  
  summarize_group <- function(pi_hat, subgroup, multiplier) {
    f1 <- elw_side(pi_hat, 1L, subgroup, multiplier)
    f0 <- elw_side(pi_hat, 0L, subgroup, multiplier)
    F1 <- weighted_cdf(f1$y, f1$weight, plot_grid)
    F0 <- weighted_cdf(f0$y, f0$weight, plot_grid)
    pte <- vapply(seq_len(nrow(intervals)), function(j) {
      p1 <- sum(f1$weight[f1$y >= intervals$lower[j] & f1$y < intervals$upper[j]])
      p0 <- sum(f0$weight[f0$y >= intervals$lower[j] & f0$y < intervals$upper[j]])
      p1 - p0
    }, numeric(1))
    c(dte = F1 - F0, pte = pte)
  }
  
  summarize_all <- function(pi_hat, multiplier) {
    unlist(lapply(group_sets, function(gs) {
      unlist(lapply(gs$index, function(idx) summarize_group(pi_hat, idx, multiplier)), use.names = FALSE)
    }), use.names = FALSE)
  }
  
  component_index <- function(dimension, group, component = c("dte", "pte")) {
    component <- match.arg(component)
    dims <- names(group_sets)
    offset <- 0L
    for (d in dims) {
      ng <- length(group_sets[[d]]$labels)
      if (d == dimension) break
      offset <- offset + ng * (length(plot_grid) + nrow(intervals))
    }
    block <- length(plot_grid) + nrow(intervals)
    start <- offset + (group - 1L) * block
    if (component == "dte") start + seq_along(plot_grid) else start + length(plot_grid) + seq_len(nrow(intervals))
  }
  
  point_pi <- fit_propensity(seed = 20270621L)
  point <- summarize_all(point_pi, rep(1, N))
  
  B <- as.integer(Sys.getenv("CSL_HETEROGENEITY_BOOTSTRAP_B", "500"))
  cores <- suppressWarnings(as.integer(Sys.getenv(
    "CSL_HETEROGENEITY_BOOTSTRAP_CORES",
    as.character(max(1L, parallel::detectCores(logical = FALSE) - 1L))
  )))
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  cores <- min(cores, B)
  cache_file <- file.path(output_dir, paste0("csl_heterogeneity_bootstrap_gbm_B", B, ".rds"))
  reuse <- file.exists(cache_file) && !tolower(Sys.getenv("CSL_FORCE_BOOTSTRAP", "false")) %in% c("true", "1", "yes")
  
  if (reuse) {
    cache <- readRDS(cache_file)
    reuse <- identical(cache$gbm_spec, gbm_spec) && identical(cache$B, B) &&
      identical(cache$plot_grid, plot_grid) && identical(cache$group_signature, group_signature)
  }
  if (reuse) {
    point <- cache$point
    boot <- cache$boot
  } else {
    bootstrap_one <- function(b) {
      set.seed(20260802L + b)
      mult <- rexp(N)
      pi_b <- fit_propensity(mult, seed = 20260802L + b)
      summarize_all(pi_b, mult)
    }
    batch_size <- suppressWarnings(as.integer(Sys.getenv("CSL_BOOTSTRAP_BATCH_SIZE", "25")))
    if (!is.finite(batch_size) || batch_size < 1L) batch_size <- 25L
    partial_file <- paste0(cache_file, ".partial.rds")
    boot <- matrix(NA_real_, B, length(point))
    completed <- 0L
    if (file.exists(partial_file)) {
      partial <- readRDS(partial_file)
      compatible <- identical(partial$B, B) &&
        isTRUE(all.equal(partial$gbm_spec, gbm_spec)) &&
        identical(partial$plot_grid, plot_grid) &&
        identical(partial$group_signature, group_signature)
      if (compatible) {
        boot <- partial$boot
        completed <- max(c(0L, which(rowSums(is.finite(boot)) == ncol(boot))))
        message("Resuming after ", completed, " completed heterogeneity replicates.")
      }
    }
    message("Running ", B, " joint heterogeneity bootstrap draws in batches of ",
            batch_size, " on ", cores, " cores.")
    while (completed < B) {
      ids <- seq.int(completed + 1L, min(B, completed + batch_size))
      draws <- if (.Platform$OS.type != "windows" && cores > 1L) {
        parallel::mclapply(ids, bootstrap_one, mc.cores = min(cores, length(ids)),
                           mc.preschedule = TRUE, mc.set.seed = FALSE)
      } else lapply(ids, bootstrap_one)
      boot[ids, ] <- do.call(rbind, draws)
      completed <- max(ids)
      saveRDS(list(B = B, boot = boot, plot_grid = plot_grid, gbm_spec = gbm_spec,
                   group_signature = group_signature), partial_file)
      message("Completed heterogeneity replicate ", completed, "/", B)
    }
    saveRDS(list(B = B, point = point, boot = boot, plot_grid = plot_grid, gbm_spec = gbm_spec,
                 group_signature = group_signature), cache_file)
    unlink(partial_file)
  }
  
  # Centered bootstrap test, as in Section 2.3 and the overall analysis.
  tail_p <- function(draws, estimate) mean(abs(draws - estimate) >= abs(estimate))
  all_dte <- list()
  all_pte <- list()
  tests <- list()
  pte_tests <- list()
  
  for (dimension in names(group_sets)) {
    gs <- group_sets[[dimension]]
    ng <- length(gs$labels)
    dte_point <- sapply(seq_len(ng), function(g) point[component_index(dimension, g, "dte")])
    dte_boot <- lapply(seq_len(ng), function(g) boot[, component_index(dimension, g, "dte"), drop = FALSE])
    pte_point <- sapply(seq_len(ng), function(g) point[component_index(dimension, g, "pte")])
    pte_boot <- lapply(seq_len(ng), function(g) boot[, component_index(dimension, g, "pte"), drop = FALSE])
  
    for (g in seq_len(ng)) {
      lo <- apply(dte_boot[[g]], 2, quantile, probs = 0.025, names = FALSE)
      hi <- apply(dte_boot[[g]], 2, quantile, probs = 0.975, names = FALSE)
      all_dte[[length(all_dte) + 1L]] <- data.frame(
        dimension = dimension, subgroup = gs$labels[g], y = plot_grid,
        estimate = dte_point[, g], ci_lower = lo, ci_upper = hi
      )
      for (j in seq_len(nrow(intervals))) {
        draws <- pte_boot[[g]][, j]
        ci <- quantile(draws, c(0.025, 0.975), names = FALSE)
        all_pte[[length(all_pte) + 1L]] <- data.frame(
          dimension = dimension, subgroup = gs$labels[g], interval = intervals$label[j],
          estimate = pte_point[j, g], se = sd(draws), ci_lower = ci[1], ci_upper = ci[2],
          p_value = tail_p(draws, pte_point[j, g])
        )
      }
    }
  
    pairs <- combn(seq_len(ng), 2L)
    point_diffs <- lapply(seq_len(ncol(pairs)), function(k) dte_point[, pairs[1, k]] - dte_point[, pairs[2, k]])
    boot_diffs <- lapply(seq_len(ncol(pairs)), function(k) dte_boot[[pairs[1, k]]] - dte_boot[[pairs[2, k]]])
    stat <- sqrt(N) * max(vapply(point_diffs, function(z) max(abs(z)), numeric(1)))
    stat_b <- vapply(seq_len(B), function(b) {
      max(vapply(seq_along(boot_diffs), function(k) {
        max(abs(boot_diffs[[k]][b, ] - point_diffs[[k]]))
      }, numeric(1))) * sqrt(N)
    }, numeric(1))
    tests[[length(tests) + 1L]] <- data.frame(
      dimension = dimension,
      hypothesis = "All subgroup DTE curves are equal",
      statistic = stat,
      p_value = mean(stat_b >= stat)
    )
  
    family_size <- ncol(pairs) * nrow(intervals)
    for (k in seq_len(ncol(pairs))) {
      g1 <- pairs[1, k]
      g2 <- pairs[2, k]
      for (j in seq_len(nrow(intervals))) {
        estimate <- pte_point[j, g2] - pte_point[j, g1]
        draws <- pte_boot[[g2]][, j] - pte_boot[[g1]][, j]
        ci <- quantile(draws, c(0.025, 0.975), names = FALSE)
        p_raw <- tail_p(draws, estimate)
        pte_tests[[length(pte_tests) + 1L]] <- data.frame(
          dimension = dimension,
          contrast = paste0(gs$labels[g2], " minus ", gs$labels[g1]),
          interval = intervals$label[j], estimate = estimate, se = sd(draws),
          ci_lower = ci[1], ci_upper = ci[2], p_value = p_raw,
          bonferroni_p = min(1, family_size * p_raw), family_size = family_size
        )
      }
    }
  
    for (k in seq_len(ncol(pairs))) {
      g1 <- pairs[1, k]
      g2 <- pairs[2, k]
      estimate <- sum(pte_point[1:2, g2]) - sum(pte_point[1:2, g1])
      draws <- rowSums(pte_boot[[g2]][, 1:2, drop = FALSE]) -
        rowSums(pte_boot[[g1]][, 1:2, drop = FALSE])
      ci <- quantile(draws, c(0.025, 0.975), names = FALSE)
      pte_tests[[length(pte_tests) + 1L]] <- data.frame(
        dimension = dimension,
        contrast = paste0(gs$labels[g2], " minus ", gs$labels[g1]),
        interval = "(-Inf, 2,500)", estimate = estimate, se = sd(draws),
        ci_lower = ci[1], ci_upper = ci[2], p_value = tail_p(draws, estimate),
        bonferroni_p = NA_real_, family_size = NA_integer_
      )
    }
  }
  
  dte_results <- do.call(rbind, all_dte)
  pte_results <- do.call(rbind, all_pte)
  test_results <- do.call(rbind, tests)
  pte_test_results <- do.call(rbind, pte_tests)
  write.csv(dte_results, file.path(output_dir, "csl_heterogeneity_dte.csv"), row.names = FALSE)
  write.csv(pte_results, file.path(output_dir, "csl_heterogeneity_pte.csv"), row.names = FALSE)
  write.csv(test_results, file.path(output_dir, "csl_heterogeneity_tests.csv"), row.names = FALSE)
  write.csv(pte_test_results, file.path(output_dir, "csl_heterogeneity_pte_interaction_tests.csv"), row.names = FALSE)
  
  plot_dimension <- function(dimension, filename, colors) {
    gs <- group_sets[[dimension]]
    dd <- dte_results[dte_results$dimension == dimension, ]
    pp <- pte_results[pte_results$dimension == dimension, ]
    pdf(file.path(figure_dir, filename), width = 9.0, height = 4.3)
    op <- par(mfrow = c(1, 2), mar = c(4.5, 4.6, 2.1, 0.7), mgp = c(2.6, 0.8, 0), las = 1, family = "serif")
    ylim <- range(dd$ci_lower, dd$ci_upper)
    plot(range(plot_grid), ylim, type = "n", xlab = "Birthweight (g)", ylab = expression(F[1](y) - F[0](y)))
    abline(h = 0, col = "grey45", lty = 2)
    abline(v = 2500, col = "grey65", lty = 3)
    for (g in seq_along(gs$labels)) {
      z <- dd[dd$subgroup == gs$labels[g], ]
      polygon(c(z$y, rev(z$y)), c(z$ci_lower, rev(z$ci_upper)), col = adjustcolor(colors[g], alpha.f = 0.15), border = NA)
      lines(z$y, z$estimate, col = colors[g], lwd = 1.8)
    }
    legend("topright", legend = gs$labels, col = colors, lwd = 2, bty = "n", cex = 0.88)
    mtext("(a) DTE", side = 3, line = 0.5, cex = 1.02)
  
    mat <- sapply(gs$labels, function(lbl) pp$estimate[pp$subgroup == lbl])
    lo <- sapply(gs$labels, function(lbl) pp$ci_lower[pp$subgroup == lbl])
    hi <- sapply(gs$labels, function(lbl) pp$ci_upper[pp$subgroup == lbl])
    ylim_p <- range(lo, hi, 0)
    ylim_p <- ylim_p + c(-0.08, 0.08) * diff(ylim_p)
    bp <- barplot(t(mat), beside = TRUE, col = colors, border = NA, ylim = ylim_p,
                  names.arg = c("<1,500", "1,500-2,499", "2,500-3,999", ">=4,000"),
                  las = 2, ylab = "Probability treatment effect", cex.names = 0.78)
    abline(h = 0, col = "grey35")
    arrows(bp, t(lo), bp, t(hi), angle = 90, code = 3, length = 0.035, lwd = 1.0)
    legend("topright", legend = gs$labels, fill = colors, border = NA, bty = "n", cex = 0.88)
    mtext("(b) PTE", side = 3, line = 0.5, cex = 1.02)
    par(op)
    dev.off()
  }
  
  plot_dimension("sex", "csl_smoking_heterogeneity_sex.pdf", c("#C85A3A", "#3C78A8"))
  plot_dimension("age", "csl_smoking_heterogeneity_age.pdf", c("#3C78A8", "#4F9B67", "#C85A3A"))
  
  message("Heterogeneity analysis complete.")
})

## ---- ps_robustness ----
if (run_step("ps_robustness")) local({
  ## Section 5.6.1: descriptive robustness to propensity-score specification.
  ## No new bootstrap inference is performed here.
  
  script_dir <- Sys.getenv("CSL_SCRIPT_DIR", unset = ".")
  output_dir <- Sys.getenv("CSL_COMPLETE_OUTPUT_DIR", unset = file.path(script_dir, "CSL-COMPLETE-PROCESS"))
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  figure_dir <- if (basename(manuscript_root) == "FIGURES") manuscript_root else file.path(manuscript_root, "FIGURES")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  .application_load_base(environment())
  
  dat <- read.csv(file.path(output_dir, "csl_smoking_birthweight_complete_case.csv"), check.names = FALSE)
  scores <- read.csv(file.path(output_dir, "csl_smoking_propensity_scores_full_sample.csv"), check.names = FALSE)
  stopifnot(identical(as.character(dat$row_id), as.character(scores$row_id)))
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  methods <- c("gbm", "logistic_expanded", "neural_network")
  labels <- c(gbm = "GBM", logistic_expanded = "PLE", neural_network = "ANN")
  colors <- c(GBM = "#222222", PLE = "#3C78A8", ANN = "#C85A3A")
  stopifnot(all(methods %in% names(scores)))
  
  grid <- seq(500, 5000, by = 10)
  taus <- c(0.10, 0.25, 0.50, 0.75)
  intervals <- data.frame(
    interval = c("(-Inf, 1,500)", "[1,500, 2,500)", "[2,500, 4,000)", "[4,000, Inf)"),
    lower = c(-Inf, 1500, 2500, 4000), upper = c(1500, 2500, 4000, Inf),
    stringsAsFactors = FALSE
  )
  
  weighted_cdf <- function(y, w, at) {
    ord <- order(y)
    approx(c(-Inf, y[ord], Inf), c(0, cumsum(w[ord]), 1), xout = at,
           method = "constant", f = 0, rule = 2, ties = "ordered")$y
  }
  
  estimate_one <- function(method, subgroup = rep(TRUE, length(A))) {
    aa <- A[subgroup]
    yy <- Y[subgroup]
    pp <- as.numeric(scores[[method]][subgroup])
    f1 <- elw_group(aa, pp, "treated")
    f0 <- elw_group(aa, pp, "control")
    y1 <- yy[aa == 1]
    y0 <- yy[aa == 0]
    F1 <- weighted_cdf(y1, f1$weight, grid)
    F0 <- weighted_cdf(y0, f0$weight, grid)
    q1 <- weighted_quantile(y1, f1$weight, taus)
    q0 <- weighted_quantile(y0, f0$weight, taus)
    pte <- vapply(seq_len(nrow(intervals)), function(j) {
      sum(f1$weight[y1 >= intervals$lower[j] & y1 < intervals$upper[j]]) -
        sum(f0$weight[y0 >= intervals$lower[j] & y0 < intervals$upper[j]])
    }, numeric(1))
    list(F1 = F1, F0 = F0, dte = F1 - F0, q1 = q1, q0 = q0,
         qte = q1 - q0, pte = pte)
  }
  
  fits <- setNames(lapply(methods, estimate_one), labels[methods])
  dte <- do.call(rbind, lapply(names(fits), function(m) data.frame(
    method = m, y = grid, F1 = fits[[m]]$F1, F0 = fits[[m]]$F0,
    estimate = fits[[m]]$dte
  )))
  qte <- do.call(rbind, lapply(names(fits), function(m) data.frame(
    method = m, tau = taus, q1 = fits[[m]]$q1, q0 = fits[[m]]$q0,
    estimate = fits[[m]]$qte
  )))
  pte <- do.call(rbind, lapply(names(fits), function(m) data.frame(
    method = m, interval = intervals$interval, estimate = fits[[m]]$pte
  )))
  
  subgroups <- list(
    Female = dat$Sex == 1, Male = dat$Sex == 2,
    `Age <20` = dat$Momage < 20,
    `Age 20--34` = dat$Momage >= 20 & dat$Momage < 35,
    `Age >=35` = dat$Momage >= 35
  )
  heterogeneity <- do.call(rbind, lapply(names(subgroups), function(g) {
    do.call(rbind, lapply(methods, function(m) {
      z <- estimate_one(m, subgroups[[g]])
      data.frame(
        subgroup = g, method = labels[[m]], n = sum(subgroups[[g]]),
        n_smokers = sum(A[subgroups[[g]]] == 1),
        dte_2500_minus = sum(z$pte[1:2]),
        max_dte = max(z$dte), y_at_max = grid[which.max(z$dte)]
      )
    }))
  }))
  
  write.csv(dte, file.path(output_dir, "csl_ps_robustness_dte.csv"), row.names = FALSE)
  write.csv(qte, file.path(output_dir, "csl_ps_robustness_qte.csv"), row.names = FALSE)
  write.csv(pte, file.path(output_dir, "csl_ps_robustness_pte.csv"), row.names = FALSE)
  write.csv(heterogeneity, file.path(output_dir, "csl_ps_robustness_heterogeneity.csv"), row.names = FALSE)
  
  pdf(file.path(figure_dir, "csl_ps_robustness_dte_pte.pdf"), width = 9.0, height = 4.3, family = "serif")
  op <- par(mfrow = c(1, 2), mar = c(4.6, 4.7, 2.1, 0.8), mgp = c(2.7, 0.8, 0), las = 1, family = "serif")
  plot(range(grid), range(dte$estimate), type = "n", xlab = "Birthweight (g)",
       ylab = expression(F[1](y) - F[0](y)))
  abline(h = 0, col = "grey50", lty = 2)
  abline(v = 2500, col = "grey65", lty = 3)
  for (m in names(fits)) lines(grid, fits[[m]]$dte, col = colors[m], lwd = ifelse(m == "GBM", 2.3, 1.7), lty = ifelse(m == "GBM", 1, 2))
  legend("topright", legend = names(fits), col = colors[names(fits)], lwd = c(2.3, 1.7, 1.7), lty = c(1, 2, 2), bty = "n", cex = 0.86)
  mtext("(a) DTE", side = 3, line = 0.5, cex = 1.02)
  
  mat <- sapply(names(fits), function(m) fits[[m]]$pte)
  bp <- barplot(t(mat), beside = TRUE, col = colors[names(fits)], border = NA,
                names.arg = c("<1,500", "1,500-2,499", "2,500-3,999", ">=4,000"),
                las = 2, ylim = c(-0.055, 0.055), ylab = "PTE", cex.names = 0.78)
  abline(h = 0, col = "grey35")
  legend("topright", legend = names(fits), fill = colors[names(fits)], border = NA, bty = "n", cex = 0.86)
  mtext("(b) PTE", side = 3, line = 0.5, cex = 1.02)
  par(op)
  dev.off()
  
  message("Propensity-score robustness sensitivity complete.")
})

## ---- weight_stress ----
if (run_step("weight_stress")) local({
  ## CSL diagnostic: ELW versus SIPW and illustrative propensity-score trimming
  
  
  script_dir <- Sys.getenv(
    "CSL_SCRIPT_DIR",
    unset = "/Volumes/POST/elw-application/r code"
  )
  output_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = file.path(script_dir, "CSL-COMPLETE-PROCESS")
  )
  
  .application_load_base(environment())
  
  dat <- read.csv(
    file.path(output_dir, "csl_smoking_birthweight_complete_case.csv"),
    check.names = FALSE
  )
  scores <- read.csv(
    file.path(output_dir, "csl_smoking_propensity_scores_full_sample.csv"),
    check.names = FALSE
  )
  stopifnot(identical(as.character(dat$row_id), as.character(scores$row_id)))
  
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  p <- as.numeric(scores$svm_radial)
  
  ess <- function(w) sum(w)^2 / sum(w^2)
  weighted_q <- function(y, w, tau) {
    ord <- order(y)
    y <- y[ord]
    w <- w[ord] / sum(w)
    y[which(cumsum(w) >= tau)[1L]]
  }
  
  summarize_two_group <- function(label, keep, w1, w0) {
    a <- A[keep]
    y <- Y[keep]
    w1 <- w1 / sum(w1)
    w0 <- w0 / sum(w0)
    data.frame(
      estimator = label,
      retained_n = sum(keep),
      excluded_n = sum(!keep),
      minimum_group_ess = min(ess(w1), ess(w0)),
      maximum_normalized_weight = max(w1, w0),
      low_birthweight_risk_difference =
        weighted.mean(y[a == 1] < 2500, w1) -
        weighted.mean(y[a == 0] < 2500, w0),
      qte_0_10 =
        weighted_q(y[a == 1], w1, 0.10) -
        weighted_q(y[a == 0], w0, 0.10)
    )
  }
  
  keep_full <- rep(TRUE, length(A))
  fit1 <- elw_group(A, p, "treated")
  fit0 <- elw_group(A, p, "control")
  elw_row <- summarize_two_group(
    "ELW, no trimming", keep_full, fit1$weight, fit0$weight
  )
  
  sipw1 <- 1 / p[A == 1]
  sipw0 <- 1 / (1 - p[A == 0])
  sipw_row <- summarize_two_group(
    "SIPW, no trimming", keep_full, sipw1, sipw0
  )
  
  trim_cutoff <- 0.005
  keep_trim <- p >= trim_cutoff & p <= 1 - trim_cutoff
  a_trim <- A[keep_trim]
  p_trim <- p[keep_trim]
  trim_row <- summarize_two_group(
    "SIPW, 0.005 <= ehat <= 0.995",
    keep_trim,
    1 / p_trim[a_trim == 1],
    1 / (1 - p_trim[a_trim == 0])
  )
  
  out <- rbind(elw_row, sipw_row, trim_row)
  write.csv(
    out,
    file.path(output_dir, "csl_elw_sipw_trimming_stress.csv"),
    row.names = FALSE
  )
  print(out)
})
# ---- 03_csl_sensitivity.R ----
requested_steps <- step
run_step <- function(x) identical(x, step)
## ---- msm ----
if (run_step("msm")) local({
  ## ============================================================
  ## CSL marginal odds-ratio sensitivity analysis
  ## ============================================================
  
  
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_path <- if (length(file_arg)) {
    gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
  } else {
    "03_csl_sensitivity.R"
  }
  script_path <- gsub("~\\+~", " ", script_path)
  script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
  
  output_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = file.path(script_dir, "CSL-COMPLETE-PROCESS")
  )
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  manuscript_dir <- if (basename(manuscript_root) == "FIGURES") {
    manuscript_root
  } else {
    file.path(manuscript_root, "FIGURES")
  }
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(manuscript_dir, recursive = TRUE, showWarnings = FALSE)
  
  .application_load_base(environment())
  
  data_file <- file.path(output_dir, "csl_smoking_birthweight_complete_case.csv")
  score_file <- file.path(output_dir, "csl_smoking_propensity_scores_full_sample.csv")
  evaluation_file <- file.path(output_dir, "csl_smoking_propensity_model_evaluation.csv")
  qte_file <- file.path(output_dir, "csl_smoking_elw_quantiles.csv")
  lbw_file <- file.path(output_dir, "csl_smoking_elw_low_birthweight.csv")
  required_files <- c(data_file, score_file, evaluation_file, qte_file, lbw_file)
  if (any(!file.exists(required_files))) {
    stop("Run the current CSL propensity-score and application scripts first.")
  }
  
  dat <- read.csv(data_file, check.names = FALSE)
  score_data <- read.csv(score_file, check.names = FALSE)
  evaluation <- read.csv(evaluation_file, check.names = FALSE)
  if (!identical(as.character(score_data$row_id), as.character(dat$row_id))) {
    stop("Propensity-score file does not match the analysis data order.")
  }
  
  primary_method <- Sys.getenv("CSL_PRIMARY_PS", "auto")
  if (primary_method == "auto") {
    selected <- evaluation$method[as.logical(evaluation$recommended_primary)]
    if (length(selected) != 1L) {
      stop("The cached propensity evaluation must identify one primary model.")
    }
    primary_method <- selected
  }
  if (!primary_method %in% names(score_data)) {
    stop("Primary propensity score is unavailable: ", primary_method)
  }
  
  method_labels <- c(
    gbm = "GBM",
    logistic_expanded = "PLE",
    kernel_single_index = "Kernel",
    gam = "GAM",
    random_forest = "Random forest",
    svm_radial = "Support vector machine",
    neural_network = "Neural network"
  )
  if (!primary_method %in% names(method_labels)) {
    stop("No display label is defined for: ", primary_method)
  }
  primary_label <- paste0("ELW-", method_labels[[primary_method]])
  
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  propensity <- as.numeric(score_data[[primary_method]])
  probs <- c(0.10, 0.25, 0.50, 0.75)
  quantile_labels <- c("10th", "25th", "50th", "75th")
  Lambda_grid <- c(seq(1, 1.5, by = 0.05), 2, 3)
  
  msm_sensitivity <- elw_msm_sensitivity(
    A, Y, propensity, Lambda_grid, probs, threshold = 2500
  )
  msm_sensitivity$quantile <- quantile_labels[
    match(msm_sensitivity$tau, probs)
  ]
  msm_sensitivity$contains_zero <- with(
    msm_sensitivity, effect_lower <= 0 & effect_upper >= 0
  )
  
  ## Guard against mixing sensitivity results with a different application fit.
  qte_reference <- read.csv(qte_file, check.names = FALSE)
  qte_reference <- qte_reference[qte_reference$method == primary_label, ]
  qte_reference <- qte_reference[order(qte_reference$tau), ]
  lbw_reference <- read.csv(lbw_file, check.names = FALSE)
  lbw_reference <- lbw_reference[lbw_reference$method == primary_label, ]
  baseline <- msm_sensitivity[msm_sensitivity$Lambda == 1, ]
  baseline <- baseline[order(baseline$tau), ]
  baseline_ok <- nrow(qte_reference) == length(probs) &&
    nrow(lbw_reference) == 1L &&
    isTRUE(all.equal(baseline$q1_lower, qte_reference$q1, tolerance = 0)) &&
    isTRUE(all.equal(baseline$q1_upper, qte_reference$q1, tolerance = 0)) &&
    isTRUE(all.equal(baseline$q0_lower, qte_reference$q0, tolerance = 0)) &&
    isTRUE(all.equal(baseline$q0_upper, qte_reference$q0, tolerance = 0)) &&
    isTRUE(all.equal(
      baseline$lbw_rd_lower,
      rep(lbw_reference$risk_difference, length(probs)),
      tolerance = 1e-12
    )) &&
    isTRUE(all.equal(
      baseline$lbw_rd_upper,
      rep(lbw_reference$risk_difference, length(probs)),
      tolerance = 1e-12
    ))
  if (!baseline_ok) {
    stop("Lambda = 1 does not reproduce the current ", primary_label, " analysis.")
  }
  
  critical_rows <- lapply(probs, function(tau) {
    x <- msm_sensitivity[msm_sensitivity$tau == tau, ]
    hit <- which(x$contains_zero)
    critical <- if (length(hit)) min(x$Lambda[hit]) else NA_real_
    data.frame(
      tau = tau,
      quantile = quantile_labels[match(tau, probs)],
      critical_Lambda = critical
    )
  })
  msm_critical <- do.call(rbind, critical_rows)
  
  interpolate_crossing <- function(x, value, side = c("upper", "lower")) {
    side <- match.arg(side)
    x <- x[order(x$Lambda), ]
    hit <- if (side == "upper") which(value(x) >= 0) else which(value(x) <= 0)
    if (!length(hit)) return(NA_real_)
    j <- hit[1]
    if (j == 1L) return(x$Lambda[j])
    x0 <- x$Lambda[j - 1L]
    x1 <- x$Lambda[j]
    y0 <- value(x)[j - 1L]
    y1 <- value(x)[j]
    if (!is.finite(y0) || !is.finite(y1) || y1 == y0) return(x1)
    x0 + (0 - y0) * (x1 - x0) / (y1 - y0)
  }
  critical_interpolated <- rbind(
    do.call(rbind, lapply(probs, function(tau) {
      x <- msm_sensitivity[msm_sensitivity$tau == tau, ]
      data.frame(
        estimand = paste0("QTE tau=", sprintf("%.2f", tau)),
        critical_Lambda = interpolate_crossing(x, function(z) z$effect_upper, "upper")
      )
    })),
    data.frame(
      estimand = "Low-birthweight risk difference",
      critical_Lambda = interpolate_crossing(
        msm_sensitivity[msm_sensitivity$tau == probs[1], ],
        function(z) z$lbw_rd_lower, "lower"
      )
    )
  )
  
  write.csv(
    msm_sensitivity,
    file.path(output_dir, "csl_smoking_msm_sensitivity.csv"),
    row.names = FALSE
  )
  write.csv(
    msm_critical,
    file.path(output_dir, "csl_smoking_msm_critical_lambda.csv"),
    row.names = FALSE
  )
  write.csv(
    critical_interpolated,
    file.path(output_dir, "csl_smoking_msm_critical_lambda_interpolated.csv"),
    row.names = FALSE
  )
  selected_rows <- msm_sensitivity[
    msm_sensitivity$Lambda %in% c(1, 1.25, 1.5, 2, 3),
  ]
  write_msm_latex(
    selected_rows,
    file.path(output_dir, "csl_smoking_msm_sensitivity_table.tex")
  )
  write_msm_latex(
    selected_rows,
    file.path(manuscript_dir, "csl_smoking_msm_sensitivity_table.tex")
  )
  
  
  cat("Marginal sensitivity analysis complete.\n")
  cat("Primary propensity model:", primary_method, "\n")
  print(baseline[, c(
    "Lambda", "tau", "q1_lower", "q0_lower", "effect_lower",
    "lbw_rd_lower"
  )])
})

## ---- msm_plot ----
if (run_step("msm_plot")) local({
  ## Redraw marginal sensitivity figure from saved results.
  ## This script does not rerun the sensitivity analysis.
  
  if (!requireNamespace("ggplot2", quietly = TRUE) ||
      !requireNamespace("patchwork", quietly = TRUE)) {
    stop("Packages 'ggplot2' and 'patchwork' are required.")
  }
  
  cmd_args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", cmd_args, value = TRUE)
  script_path <- if (length(file_arg)) {
    gsub("~+~", " ", sub("^--file=", "", file_arg[[1]]), fixed = TRUE)
  } else {
    "03_csl_sensitivity.R"
  }
  script_path <- gsub("~\\+~", " ", script_path)
  script_dir <- dirname(normalizePath(script_path, mustWork = FALSE))
  
  main_dir <- script_dir
  output_dir <- Sys.getenv(
    "CSL_COMPLETE_OUTPUT_DIR",
    unset = file.path(main_dir, "CSL-COMPLETE-PROCESS")
  )
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  manuscript_dir <- if (basename(manuscript_root) == "FIGURES") {
    manuscript_root
  } else {
    file.path(manuscript_root, "FIGURES")
  }
  
  msm_path <- file.path(output_dir, "csl_smoking_msm_sensitivity.csv")
  if (!file.exists(msm_path)) {
    stop("Cannot find saved sensitivity results: ", msm_path)
  }
  
  msm_sensitivity <- read.csv(msm_path)
  
  ## Do not redraw a sensitivity cache whose MAR baseline comes from another fit.
  evaluation_path <- file.path(output_dir, "csl_smoking_propensity_model_evaluation.csv")
  qte_path <- file.path(output_dir, "csl_smoking_elw_quantiles.csv")
  if (file.exists(evaluation_path) && file.exists(qte_path)) {
    evaluation <- read.csv(evaluation_path, check.names = FALSE)
    selected <- Sys.getenv("CSL_PRIMARY_PS", "auto")
    if (selected == "auto") {
      selected <- evaluation$method[as.logical(evaluation$recommended_primary)]
    }
    method_labels <- c(
      gbm = "GBM",
      logistic_expanded = "PLE", kernel_single_index = "Kernel",
      gam = "GAM", random_forest = "Random forest",
      svm_radial = "Support vector machine", neural_network = "Neural network"
    )
    if (length(selected) != 1L || !selected %in% names(method_labels)) {
      stop("Cannot identify the current primary propensity model.")
    }
    qte_reference <- read.csv(qte_path, check.names = FALSE)
    qte_reference <- qte_reference[
      qte_reference$method == paste0("ELW-", method_labels[[selected]]),
    ]
    qte_reference <- qte_reference[order(qte_reference$tau), ]
    baseline <- msm_sensitivity[msm_sensitivity$Lambda == 1, ]
    baseline <- baseline[order(baseline$tau), ]
    if (nrow(qte_reference) != nrow(baseline) ||
        !isTRUE(all.equal(
          baseline$effect_lower, qte_reference$difference, tolerance = 0
        )) ||
        !isTRUE(all.equal(
          baseline$effect_upper, qte_reference$difference, tolerance = 0
        ))) {
      stop("The cached Lambda = 1 sensitivity results do not match the current primary analysis.")
    }
  }
  
  msm_sensitivity <- msm_sensitivity[msm_sensitivity$Lambda >= 1 & msm_sensitivity$Lambda <= 1.50, ]
  msm_sensitivity$quantile <- factor(
    msm_sensitivity$quantile,
    levels = c("10th", "25th", "50th", "75th")
  )
  
  interpolate_grid <- seq(1, 1.50, by = 0.0025)
  msm_plot_data <- do.call(rbind, lapply(
    split(msm_sensitivity, msm_sensitivity$quantile),
    function(x) {
      x <- x[order(x$Lambda), ]
      data.frame(
        Lambda = interpolate_grid,
        tau = x$tau[1],
        quantile = x$quantile[1],
        effect_lower = approx(x$Lambda, x$effect_lower, xout = interpolate_grid)$y,
        effect_upper = approx(x$Lambda, x$effect_upper, xout = interpolate_grid)$y
      )
    }
  ))
  msm_plot_data$quantile <- factor(
    msm_plot_data$quantile,
    levels = levels(msm_sensitivity$quantile)
  )
  
  dte_data <- msm_sensitivity[
    !duplicated(msm_sensitivity$Lambda),
    c("Lambda", "lbw_rd_lower", "lbw_rd_upper")
  ]
  dte_data <- dte_data[order(dte_data$Lambda), ]
  dte_plot_data <- data.frame(
    Lambda = interpolate_grid,
    effect_lower = approx(
      dte_data$Lambda, dte_data$lbw_rd_lower, xout = interpolate_grid
    )$y,
    effect_upper = approx(
      dte_data$Lambda, dte_data$lbw_rd_upper, xout = interpolate_grid
    )$y
  )
  dte_hit <- which(dte_data$lbw_rd_lower <= 0 & dte_data$lbw_rd_upper >= 0)
  if (!length(dte_hit)) stop("The low-birthweight DTE interval does not cross zero.")
  dte_j <- dte_hit[1]
  if (dte_j == 1L) {
    dte_critical <- dte_data$Lambda[dte_j]
  } else {
    dte_x0 <- dte_data$Lambda[dte_j - 1L]
    dte_x1 <- dte_data$Lambda[dte_j]
    dte_y0 <- dte_data$lbw_rd_lower[dte_j - 1L]
    dte_y1 <- dte_data$lbw_rd_lower[dte_j]
    dte_critical <- dte_x0 + (0 - dte_y0) *
      (dte_x1 - dte_x0) / (dte_y1 - dte_y0)
  }
  
  first_zero <- do.call(rbind, lapply(
    split(msm_sensitivity, msm_sensitivity$quantile),
    function(x) {
      x <- x[order(x$Lambda), ]
      hit <- which(x$effect_lower <= 0 & x$effect_upper >= 0)
      if (!length(hit)) return(NULL)
      j <- hit[1]
      if (j == 1L) {
        critical <- x$Lambda[j]
      } else {
        x0 <- x$Lambda[j - 1L]
        x1 <- x$Lambda[j]
        y0 <- x$effect_upper[j - 1L]
        y1 <- x$effect_upper[j]
        critical <- if (is.finite(y0) && is.finite(y1) && y1 != y0) {
          x0 + (0 - y0) * (x1 - x0) / (y1 - y0)
        } else {
          x1
        }
      }
      data.frame(Lambda = critical, quantile = x$quantile[j])
    }
  ))
  first_zero$quantile <- factor(
    first_zero$quantile,
    levels = levels(msm_sensitivity$quantile)
  )
  first_zero$label <- sprintf("%.3f", first_zero$Lambda)
  plot_range <- range(c(
    msm_plot_data$effect_lower, msm_plot_data$effect_upper, 0
  ), na.rm = TRUE)
  plot_span <- diff(plot_range)
  ## Place crossing-value labels just above the x-axis, away from the curves and zero line.
  first_zero$label_y <- rep(plot_range[1] + 0.045 * plot_span, nrow(first_zero))
  
  quantile_colors <- c(
    "10th" = "#0072B2",
    "25th" = "#009E73",
    "50th" = "#D55E00",
    "75th" = "#CC79A7"
  )
  quantile_linetypes <- c(
    "10th" = "solid",
    "25th" = "dashed",
    "50th" = "dotdash",
    "75th" = "longdash"
  )
  
  msm_plot <- ggplot2::ggplot(
    msm_plot_data,
    ggplot2::aes(
      x = Lambda,
      color = quantile,
      linetype = quantile,
      group = quantile
    )
  ) +
    ggplot2::geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    ggplot2::geom_line(ggplot2::aes(y = effect_lower), linewidth = 0.65) +
    ggplot2::geom_line(ggplot2::aes(y = effect_upper), linewidth = 0.65) +
    ggplot2::geom_vline(
      data = first_zero,
      ggplot2::aes(xintercept = Lambda, color = quantile),
      linetype = "dotted", linewidth = 0.45, show.legend = FALSE
    ) +
    ggplot2::geom_point(
      data = first_zero,
      ggplot2::aes(x = Lambda, y = 0, color = quantile),
      size = 1.8, show.legend = FALSE
    ) +
    ggplot2::geom_text(
      data = first_zero,
      ggplot2::aes(x = Lambda, y = label_y, label = label, color = quantile),
      hjust = 0.5, vjust = 0, size = 3.1, fontface = "bold", show.legend = FALSE
    ) +
    ggplot2::scale_color_manual(
      values = quantile_colors,
      labels = c("10th" = "0.10", "25th" = "0.25", "50th" = "0.50", "75th" = "0.75")
    ) +
    ggplot2::scale_linetype_manual(
      values = quantile_linetypes,
      labels = c("10th" = "0.10", "25th" = "0.25", "50th" = "0.50", "75th" = "0.75")
    ) +
    ggplot2::scale_x_continuous(
      limits = c(1, 1.50),
      breaks = seq(1, 1.50, by = 0.10),
      expand = ggplot2::expansion(mult = c(0.01, 0.08))
    ) +
    ggplot2::scale_y_continuous(
      labels = function(x) formatC(x, format = "d", big.mark = ",")
    ) +
    ggplot2::labs(
      title = "(b) QTE",
      x = expression(Lambda),
      y = "ELW quantile-effect sensitivity interval",
      color = expression(tau),
      linetype = expression(tau)
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      legend.position = "top",
      legend.title = ggplot2::element_text(size = 10),
      legend.key.width = grid::unit(1.7, "lines"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
  
  dte_range <- range(c(dte_plot_data$effect_lower, dte_plot_data$effect_upper, 0))
  dte_label_y <- dte_range[1] + 0.06 * diff(dte_range)
  dte_plot <- ggplot2::ggplot(dte_plot_data, ggplot2::aes(x = Lambda)) +
    ggplot2::geom_hline(yintercept = 0, color = "grey45", linewidth = 0.35) +
    ggplot2::geom_line(
      ggplot2::aes(y = effect_lower), color = "#0072B2", linewidth = 0.75
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = effect_upper), color = "#0072B2", linewidth = 0.75
    ) +
    ggplot2::geom_vline(
      xintercept = dte_critical, color = "#0072B2",
      linetype = "dotted", linewidth = 0.45
    ) +
    ggplot2::annotate(
      "point", x = dte_critical, y = 0,
      color = "#0072B2", size = 1.8
    ) +
    ggplot2::annotate(
      "text", x = dte_critical, y = dte_label_y,
      label = sprintf("%.3f", dte_critical), color = "#0072B2",
      vjust = 0, size = 3.1, fontface = "bold"
    ) +
    ggplot2::scale_x_continuous(
      limits = c(1, 1.50), breaks = seq(1, 1.50, by = 0.10),
      expand = ggplot2::expansion(mult = c(0.01, 0.08))
    ) +
    ggplot2::labs(
      title = "(a) DTE",
      x = expression(Lambda),
      y = "ELW low-birthweight DTE sensitivity interval"
    ) +
    ggplot2::theme_classic(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 11),
      plot.margin = ggplot2::margin(6, 10, 6, 6)
    )
  
  combined_plot <- dte_plot + msm_plot +
    patchwork::plot_layout(widths = c(0.42, 0.56))
  
  for (d in c(output_dir, manuscript_dir)) {
    if (!dir.exists(d)) dir.create(d, recursive = TRUE)
    ggplot2::ggsave(
      file.path(d, "csl_smoking_msm_sensitivity.pdf"),
      msm_plot, width = 6.8, height = 4.8
    )
    ggplot2::ggsave(
      file.path(d, "csl_smoking_msm_sensitivity_dte.pdf"),
      dte_plot, width = 5.2, height = 4.8
    )
    ggplot2::ggsave(
      file.path(d, "csl_smoking_msm_sensitivity_combined.pdf"),
      combined_plot, width = 9.0, height = 4.3
    )
  }
})

## ---- lmz ----
if (run_step("lmz")) local({
  ## CSL application: LMZ shadow-variable MNAR sensitivity analysis.
  
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(file_arg)) {
    dirname(normalizePath(gsub("~+~", " ", sub("^--file=", "", file_arg[[1L]]), fixed = TRUE)))
  } else getwd()
  
  .application_load_lmz(environment())
  input_file <- file.path(
    Sys.getenv("CSL_COMPLETE_OUTPUT", file.path(script_dir, "CSL-COMPLETE-PROCESS")),
    "csl_smoking_birthweight_complete_case.csv"
  )
  output_dir <- Sys.getenv(
    "LMZ_APP_OUTPUT",
    file.path(script_dir, "lmz-application")
  )
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  dat <- read.csv(input_file, check.names = FALSE)
  A <- as.integer(dat$treatment_smoking)
  Y <- as.numeric(dat$outcome_birthweight)
  
  selection_vars <- c(
    "Momrace", "Marital", "Momage", "Parity", "Hxpreterm", "Height",
    "PreBMI", "Diabetes", "Hxasthma", "HIV", "Alcohol"
  )
  factor_vars <- c(
    "Momrace", "Marital", "Parity", "Hxpreterm", "Diabetes", "Hxasthma",
    "HIV", "Alcohol"
  )
  design_data <- dat[selection_vars]
  design_data[factor_vars] <- lapply(design_data[factor_vars], factor)
  X <- model.matrix(~ ., data = design_data)
  ps_glm <- suppressWarnings(glm.fit(X, A, family = binomial()))
  ell <- as.numeric(ps_glm$linear.predictors)
  
  lmz_data <- data.frame(T = A, Y = Y, Z = factor(dat$Sex), ell = ell)
  ps_fit <- list(linear_predictor = ell)
  lmz <- estimate_lmz_probabilities(lmz_data, ps_fit)
  
  i1 <- A == 1L
  i0 <- A == 0L
  rho1 <- lmz$fit1$rho[i1]
  rho0 <- lmz$fit0$rho[i0]
  simulation_validate_probability(c(rho1, rho0), "observed LMZ probability")
  
  p1 <- rep(0.5, length(A)); p1[i1] <- rho1
  p0 <- rep(0.5, length(A)); p0[i0] <- 1 - rho0
  elw1 <- elw_group(A, p1, "treated")
  elw0 <- elw_group(A, p0, "control")
  
  alpha1 <- mean(i1); alpha0 <- mean(i0)
  sipw1 <- alpha1 / rho1; sipw1 <- sipw1 / sum(sipw1)
  sipw0 <- alpha0 / rho0; sipw0 <- sipw0 / sum(sipw0)
  probs <- c(0.10, 0.25, 0.50, 0.75)
  
  summarize_method <- function(method, w1, w0) {
    q1 <- weighted_quantile(Y[i1], w1, probs)
    q0 <- weighted_quantile(Y[i0], w0, probs)
    data.frame(
      method = method,
      F1_2500 = sum(w1 * (Y[i1] < 2500)) / sum(w1),
      F0_2500 = sum(w0 * (Y[i0] < 2500)) / sum(w0),
      RD_2500 = sum(w1 * (Y[i1] < 2500)) / sum(w1) -
        sum(w0 * (Y[i0] < 2500)) / sum(w0),
      qte_10 = q1[1] - q0[1], qte_25 = q1[2] - q0[2],
      qte_50 = q1[3] - q0[3], qte_75 = q1[4] - q0[4]
    )
  }
  
  estimates <- rbind(
    summarize_method("ELW--LMZ (MNAR)", elw1$weight, elw0$weight),
    summarize_method("SIPW--LMZ (MNAR)", sipw1, sipw0)
  )
  
  side_diagnostic <- function(side, fit, rho, weight, elw = NULL) {
    data.frame(
      side = side, convergence = fit$convergence,
      beta_intercept = fit$beta[1], beta_index = fit$beta[2],
      gamma = paste(signif(fit$gamma, 6), collapse = ","),
      rho_min = min(rho), rho_q01 = quantile(rho, 0.01),
      rho_median = median(rho), rho_max = max(rho),
      ess = sum(weight)^2 / sum(weight^2), max_weight = max(weight / sum(weight)),
      normalization_error = if (is.null(elw)) NA_real_ else elw$normalization_error,
      calibration_error = if (is.null(elw)) NA_real_ else elw$calibration_error
    )
  }
  diagnostics <- rbind(
    side_diagnostic("Y(1)", lmz$fit1, rho1, elw1$weight, elw1),
    side_diagnostic("Y(0)", lmz$fit0, rho0, elw0$weight, elw0)
  )
  
  write.csv(estimates, file.path(output_dir, "csl_lmz_mnar_estimates.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(output_dir, "csl_lmz_mnar_diagnostics.csv"), row.names = FALSE)
  saveRDS(
    list(estimates = estimates, diagnostics = diagnostics, lmz = lmz,
         elw = list(side1 = elw1, side0 = elw0),
         settings = list(shadow = "Sex", selection_vars = selection_vars)),
    file.path(output_dir, "csl_lmz_mnar_results.rds")
  )
  print(estimates)
  print(diagnostics)
})

## ---- missing_data ----
if (run_step("missing_data")) local({
  ## Section 5.6.3: multiply imputed expanded-sample sensitivity analysis.
  ## The primary expanded sample keeps all 211,059 records with observed birthweight.
  ## No bootstrap inference or conformal analysis is performed here.
  
  if (!requireNamespace("mice", quietly = TRUE)) stop("Package 'mice' is required.")
  if (!requireNamespace("gbm", quietly = TRUE)) stop("Package 'gbm' is required.")
  
  raw_file <- Sys.getenv(
    "CSL_RAW_FILE",
    unset = "CSLdata171003.csv"
  )
  script_dir <- Sys.getenv("CSL_SCRIPT_DIR", unset = ".")
  output_dir <- Sys.getenv("CSL_COMPLETE_OUTPUT_DIR", unset = file.path(script_dir, "CSL-COMPLETE-PROCESS"))
  manuscript_root <- Sys.getenv(
    "CSL_MANUSCRIPT_DIR",
    unset = "output"
  )
  figure_dir <- if (basename(manuscript_root) == "FIGURES") manuscript_root else file.path(manuscript_root, "FIGURES")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
  .application_load_base(environment())
  
  M <- as.integer(Sys.getenv("CSL_MI_M", "20"))
  maxit <- as.integer(Sys.getenv("CSL_MI_MAXIT", "5"))
  cores <- as.integer(Sys.getenv("CSL_MI_ANALYSIS_CORES", "5"))
  if (!is.finite(M) || M < 2L) M <- 20L
  if (!is.finite(maxit) || maxit < 1L) maxit <- 5L
  if (!is.finite(cores) || cores < 1L) cores <- 1L
  
  raw <- read.csv(raw_file, check.names = FALSE)
  if (names(raw)[1] == "") names(raw)[1] <- "row_id"
  special_missing_codes <- list(
    Momrace = c(8, 9), Marital = c(8, 9), Hxpreterm = c(8, 9),
    Diabetes = c(8, 9), Hxasthma = c(8, 9), HIV = c(8, 9), Sex = c(3, 8, 9)
  )
  for (v in names(special_missing_codes)) {
    raw[[v]][raw[[v]] %in% special_missing_codes[[v]]] <- NA
  }
  vars <- c(
    "row_id", "Smoke", "Birthweight", "Momrace", "Marital", "Momage",
    "Parity", "Hxpreterm", "Height", "PreBMI", "Diabetes", "Hxasthma",
    "HIV", "Alcohol", "Sex"
  )
  dat <- raw[!is.na(raw$Birthweight), vars]
  stopifnot(nrow(dat) == 211059L, all(dat$Smoke %in% c(1, 2)))
  missing_mask <- is.na(dat)
  
  factor_vars <- c(
    "Smoke", "Momrace", "Marital", "Parity", "Hxpreterm", "Diabetes",
    "Hxasthma", "HIV", "Alcohol", "Sex"
  )
  dat[factor_vars] <- lapply(dat[factor_vars], factor)
  
  method <- mice::make.method(dat)
  method[c("row_id", "Smoke", "Birthweight", "Parity", "Alcohol")] <- ""
  method[c("Momage", "Height", "PreBMI")] <- "pmm"
  method[c("Momrace", "Marital")] <- "polyreg"
  method[c("Hxpreterm", "Diabetes", "Hxasthma", "HIV", "Sex")] <- "logreg"
  predictor_matrix <- mice::make.predictorMatrix(dat)
  predictor_matrix[, "row_id"] <- 0
  predictor_matrix["row_id", ] <- 0
  diag(predictor_matrix) <- 0
  
  cache_file <- file.path(output_dir, paste0("csl_expanded_mice_M", M, "_maxit", maxit, ".rds"))
  if (file.exists(cache_file) && !tolower(Sys.getenv("CSL_FORCE_MI", "false")) %in% c("true", "1", "yes")) {
    imp <- readRDS(cache_file)
  } else {
    set.seed(20260803)
    imp <- mice::mice(
      dat, m = M, maxit = maxit, method = method,
      predictorMatrix = predictor_matrix, seed = 20260803,
      printFlag = TRUE
    )
    saveRDS(imp, cache_file, compress = FALSE)
  }
  
  grid <- seq(500, 5000, by = 10)
  taus <- c(0.10, 0.25, 0.50, 0.75)
  intervals <- data.frame(
    interval = c("(-Inf, 1,500)", "[1,500, 2,500)", "[2,500, 4,000)", "[4,000, Inf)"),
    lower = c(-Inf, 1500, 2500, 4000), upper = c(1500, 2500, 4000, Inf),
    stringsAsFactors = FALSE
  )
  gbm_spec <- list(n_trees = 1000L, interaction_depth = 3L, min_node_size = 50L,
                   shrinkage = 0.02, bag_fraction = 0.8)
  
  weighted_cdf <- function(y, w, at) {
    ord <- order(y)
    approx(c(-Inf, y[ord], Inf), c(0, cumsum(w[ord]), 1), xout = at,
           method = "constant", f = 0, rule = 2, ties = "ordered")$y
  }
  
  analyze_completed <- function(completed, imputation, scenario, bmi_delta = 0, height_delta = 0) {
    ## Delta shifts apply only to values that were originally missing.
    completed$PreBMI[missing_mask[, "PreBMI"]] <- completed$PreBMI[missing_mask[, "PreBMI"]] + bmi_delta
    completed$Height[missing_mask[, "Height"]] <- completed$Height[missing_mask[, "Height"]] + height_delta
    A <- as.integer(as.character(completed$Smoke)) - 1L
    Y <- as.numeric(completed$Birthweight)
    X <- completed[setdiff(vars, c("row_id", "Smoke", "Birthweight"))]
    X[factor_vars[factor_vars %in% names(X)]] <- lapply(X[factor_vars[factor_vars %in% names(X)]], factor)
    Xmat <- model.matrix(~ . - 1, data = X)
    storage.mode(Xmat) <- "double"
    set.seed(20260803L + imputation + match(scenario, c("MAR MI", "Lower delta", "Higher delta")) * 1000L)
    fit <- gbm::gbm.fit(
      x = Xmat, y = A, distribution = "bernoulli",
      n.trees = gbm_spec$n_trees, interaction.depth = gbm_spec$interaction_depth,
      n.minobsinnode = gbm_spec$min_node_size, shrinkage = gbm_spec$shrinkage,
      bag.fraction = gbm_spec$bag_fraction, nTrain = nrow(completed),
      keep.data = FALSE, verbose = FALSE
    )
    p <- as.numeric(predict(fit, Xmat, n.trees = gbm_spec$n_trees, type = "response"))
    f1 <- elw_group(A, p, "treated")
    f0 <- elw_group(A, p, "control")
    y1 <- Y[A == 1]
    y0 <- Y[A == 0]
    F1 <- weighted_cdf(y1, f1$weight, grid)
    F0 <- weighted_cdf(y0, f0$weight, grid)
    q1 <- weighted_quantile(y1, f1$weight, taus)
    q0 <- weighted_quantile(y0, f0$weight, taus)
    pte <- vapply(seq_len(nrow(intervals)), function(j) {
      sum(f1$weight[y1 >= intervals$lower[j] & y1 < intervals$upper[j]]) -
        sum(f0$weight[y0 >= intervals$lower[j] & y0 < intervals$upper[j]])
    }, numeric(1))
    group_index <- list(
      Female = as.character(completed$Sex) == "1",
      Male = as.character(completed$Sex) == "2",
      `Age <25` = completed$Momage < 25,
      `Age 25--34` = completed$Momage >= 25 & completed$Momage < 35,
      `Age >=35` = completed$Momage >= 35
    )
    subgroup <- do.call(rbind, lapply(names(group_index), function(g) {
      idx <- group_index[[g]]
      aa <- A[idx]
      yy <- Y[idx]
      pp <- p[idx]
      z1 <- elw_group(aa, pp, "treated")
      z0 <- elw_group(aa, pp, "control")
      d2500 <- sum(z1$weight[yy[aa == 1] < 2500]) - sum(z0$weight[yy[aa == 0] < 2500])
      data.frame(subgroup = g, dte_2500 = d2500)
    }))
    list(
      scalar = data.frame(
        imputation = imputation, scenario = scenario,
        mean1 = sum(y1 * f1$weight), mean0 = sum(y0 * f0$weight),
        ate = sum(y1 * f1$weight) - sum(y0 * f0$weight),
        risk_difference = pte[1] + pte[2]
      ),
      dte = data.frame(imputation = imputation, scenario = scenario, y = grid,
                       F1 = F1, F0 = F0, estimate = F1 - F0),
      qte = data.frame(imputation = imputation, scenario = scenario, tau = taus,
                       q1 = q1, q0 = q0, estimate = q1 - q0),
      pte = data.frame(imputation = imputation, scenario = scenario,
                       interval = intervals$interval, estimate = pte),
      subgroup = transform(subgroup, imputation = imputation, scenario = scenario)
    )
  }
  
  ## A quarter-SD shift is used as a moderate stress test in both directions.
  bmi_shift <- 0.25 * sd(as.numeric(dat$PreBMI), na.rm = TRUE)
  height_shift <- 0.25 * sd(as.numeric(dat$Height), na.rm = TRUE)
  scenarios <- data.frame(
    scenario = c("MAR MI", "Lower delta", "Higher delta"),
    bmi_delta = c(0, -bmi_shift, bmi_shift),
    height_delta = c(0, -height_shift, height_shift),
    stringsAsFactors = FALSE
  )
  
  one_imputation <- function(m) {
    completed <- mice::complete(imp, m)
    lapply(seq_len(nrow(scenarios)), function(s) {
      analyze_completed(completed, m, scenarios$scenario[s],
                        scenarios$bmi_delta[s], scenarios$height_delta[s])
    })
  }
  message("Analyzing ", M, " imputations under three scenarios on ", cores, " cores.")
  results <- if (.Platform$OS.type != "windows" && cores > 1L) {
    parallel::mclapply(seq_len(M), one_imputation, mc.cores = min(cores, M),
                       mc.preschedule = TRUE, mc.set.seed = FALSE)
  } else lapply(seq_len(M), one_imputation)
  results <- unlist(results, recursive = FALSE)
  
  bind_component <- function(name) do.call(rbind, lapply(results, `[[`, name))
  scalar_raw <- bind_component("scalar")
  dte_raw <- bind_component("dte")
  qte_raw <- bind_component("qte")
  pte_raw <- bind_component("pte")
  subgroup_raw <- bind_component("subgroup")
  
  pool_numeric <- function(data, keys, values) {
    means <- aggregate(data[values], data[keys], mean)
    sds <- aggregate(data[values], data[keys], sd)
    names(means)[names(means) %in% values] <- paste0(values, "_mean")
    names(sds)[names(sds) %in% values] <- paste0(values, "_between_sd")
    merge(means, sds, by = keys, sort = FALSE)
  }
  scalar_pooled <- pool_numeric(scalar_raw, "scenario", c("mean1", "mean0", "ate", "risk_difference"))
  dte_pooled <- pool_numeric(dte_raw, c("scenario", "y"), c("F1", "F0", "estimate"))
  qte_pooled <- pool_numeric(qte_raw, c("scenario", "tau"), c("q1", "q0", "estimate"))
  pte_pooled <- pool_numeric(pte_raw, c("scenario", "interval"), "estimate")
  subgroup_pooled <- pool_numeric(subgroup_raw, c("scenario", "subgroup"), "dte_2500")
  
  flatten <- function(x) data.frame(x, check.names = FALSE)
  scalar_pooled <- flatten(scalar_pooled)
  dte_pooled <- flatten(dte_pooled)
  qte_pooled <- flatten(qte_pooled)
  pte_pooled <- flatten(pte_pooled)
  subgroup_pooled <- flatten(subgroup_pooled)
  
  write.csv(scalar_raw, file.path(output_dir, "csl_mi_scalar_raw.csv"), row.names = FALSE)
  write.csv(dte_raw, file.path(output_dir, "csl_mi_dte_raw.csv"), row.names = FALSE)
  write.csv(qte_raw, file.path(output_dir, "csl_mi_qte_raw.csv"), row.names = FALSE)
  write.csv(pte_raw, file.path(output_dir, "csl_mi_pte_raw.csv"), row.names = FALSE)
  write.csv(subgroup_raw, file.path(output_dir, "csl_mi_subgroup_raw.csv"), row.names = FALSE)
  write.csv(scalar_pooled, file.path(output_dir, "csl_mi_scalar_pooled.csv"), row.names = FALSE)
  write.csv(dte_pooled, file.path(output_dir, "csl_mi_dte_pooled.csv"), row.names = FALSE)
  write.csv(qte_pooled, file.path(output_dir, "csl_mi_qte_pooled.csv"), row.names = FALSE)
  write.csv(pte_pooled, file.path(output_dir, "csl_mi_pte_pooled.csv"), row.names = FALSE)
  write.csv(subgroup_pooled, file.path(output_dir, "csl_mi_subgroup_pooled.csv"), row.names = FALSE)
  write.csv(scenarios, file.path(output_dir, "csl_mi_delta_scenarios.csv"), row.names = FALSE)
  
  ## Compare with the complete-case GBM point estimates created in 5.6.1.
  cc_dte <- read.csv(file.path(output_dir, "csl_ps_robustness_dte.csv"))
  cc_dte <- cc_dte[cc_dte$method == "GBM", ]
  cc_pte <- read.csv(file.path(output_dir, "csl_ps_robustness_pte.csv"))
  cc_pte <- cc_pte[cc_pte$method == "GBM", ]
  scenario_colors <- c("MAR MI" = "#3C78A8", "Lower delta" = "#4F9B67", "Higher delta" = "#C85A3A")
  
  pdf(file.path(figure_dir, "csl_missing_data_sensitivity.pdf"), width = 9.0, height = 4.3, family = "serif")
  op <- par(mfrow = c(1, 2), mar = c(4.6, 4.7, 2.1, 0.8), mgp = c(2.7, 0.8, 0), las = 1, family = "serif")
  plot(range(grid), range(c(cc_dte$estimate, dte_pooled$estimate_mean)), type = "n",
       xlab = "Birthweight (g)", ylab = expression(F[1](y) - F[0](y)))
  abline(h = 0, col = "grey50", lty = 2)
  abline(v = 2500, col = "grey65", lty = 3)
  lines(cc_dte$y, cc_dte$estimate, col = "#222222", lwd = 2.3)
  for (s in scenarios$scenario) {
    z <- dte_pooled[dte_pooled$scenario == s, ]
    lines(z$y, z$estimate_mean, col = scenario_colors[s], lwd = 1.7, lty = 2)
  }
  legend("topright", legend = c("Complete case", scenarios$scenario),
         col = c("#222222", scenario_colors), lwd = c(2.3, rep(1.7, 3)),
         lty = c(1, rep(2, 3)), bty = "n", cex = 0.78)
  mtext("(a) DTE", side = 3, line = 0.5, cex = 1.02)
  
  mat <- cbind(
    `Complete case` = cc_pte$estimate,
    sapply(scenarios$scenario, function(s) pte_pooled$estimate_mean[pte_pooled$scenario == s])
  )
  barplot(t(mat), beside = TRUE, col = c("#222222", scenario_colors), border = NA,
          names.arg = c("<1,500", "1,500-2,499", "2,500-3,999", ">=4,000"),
          las = 2, ylab = "Probability treatment effect", cex.names = 0.78)
  abline(h = 0, col = "grey35")
  legend("topright", legend = colnames(mat), fill = c("#222222", scenario_colors),
         border = NA, bty = "n", cex = 0.72)
  mtext("(b) PTE", side = 3, line = 0.5, cex = 1.02)
  par(op)
  dev.off()
  
  message("Multiple-imputation sensitivity complete.")
})
if(step == 'lmz_gamma') local({
commandArgs <- function(trailingOnly=FALSE) c(input, file.path(out,'lmz-gamma'), Sys.getenv('CSL_LMZ_GAMMA_B','500'), Sys.getenv('CSL_LMZ_GAMMA_CORES','1'))
# Exchangeably weighted bootstrap for the standardized LMZ selection coefficients gamma_t.
# For the linear g_t and normal observed-outcome model, gamma_t has the closed form
#   gamma_t = (2t - 1) * (sex coefficient of the marginal response logit) * SD_t / delta_t,
# where SD_t is the observed-side birthweight SD and delta_t the observed male-female mean
# difference. Each draw refits the treatment index, the marginal response logit, and the
# observed-side moments with Exp(1) multipliers.
# Usage:
#   Rscript lmz_gamma_bootstrap.R COMPLETE_CASE_CSV OUTPUT_DIR [B] [CORES]
args <- commandArgs(TRUE)
if (length(args) < 2L) stop("Usage: Rscript lmz_gamma_bootstrap.R COMPLETE_CASE_CSV OUTPUT_DIR [B] [CORES]")
d <- read.csv(args[1])
out <- args[2]; dir.create(out, recursive = TRUE, showWarnings = FALSE)
B <- if (length(args) >= 3L) as.integer(args[3]) else 500L
cores <- if (length(args) >= 4L) as.integer(args[4]) else max(1L, parallel::detectCores(logical = FALSE) - 1L)

a <- d$treatment_smoking
y <- d$outcome_birthweight
n <- length(a)
vars <- c("Momrace", "Marital", "Momage", "Parity", "Hxpreterm", "Height", "PreBMI",
          "Diabetes", "Hxasthma", "HIV", "Alcohol")
fac <- c("Momrace", "Marital", "Parity", "Hxpreterm", "Diabetes", "Hxasthma", "HIV", "Alcohol")
x <- d[vars]; x[fac] <- lapply(x[fac], factor)
X <- model.matrix(~ ., x)
Z <- model.matrix(~ Z, data = data.frame(Z = factor(d$Sex)))

fit_gamma <- function(w) {
  w <- w / mean(w)
  first <- suppressWarnings(glm.fit(X, a, weights = w, family = binomial()))
  stopifnot(first$converged, all(is.finite(first$coefficients)))
  selection <- suppressWarnings(glm.fit(cbind(1, first$linear.predictors, Z[, 2]), a,
                                        weights = w, family = binomial()))
  stopifnot(selection$converged, all(is.finite(selection$coefficients)))
  vapply(1:0, function(t) {
    r <- a == t
    delta <- lm.wfit(Z[r, , drop = FALSE], y[r], w[r])$coefficients[2]
    mu <- weighted.mean(y[r], w[r])
    s <- sqrt(sum(w[r] * (y[r] - mu)^2) / (sum(w[r]) - 1))
    unname((2 * t - 1) * selection$coefficients[3] * s / delta)
  }, numeric(1))
}

point <- setNames(fit_gamma(rep(1, n)), c("T1", "T0"))
draws <- do.call(rbind, parallel::mclapply(seq_len(B), function(b) {
  set.seed(20260618L + b)
  fit_gamma(rexp(n))
}, mc.cores = min(cores, B), mc.set.seed = FALSE))
colnames(draws) <- c("T1", "T0")
ci <- apply(draws, 2, quantile, c(0.025, 0.975))
result <- data.frame(t = 1:0, estimate = point, bootstrap_SE = apply(draws, 2, sd),
                     lower = ci[1, ], upper = ci[2, ], B = B)
write.csv(data.frame(b = seq_len(B), draws), file.path(out, "gamma-bootstrap-draws.csv"), row.names = FALSE)
write.csv(result, file.path(out, "gamma-bootstrap-ci.csv"), row.names = FALSE)
print(result, digits = 4, row.names = FALSE)
})
invisible(out)
}

if (sys.nframe() == 0L) .application_main()
