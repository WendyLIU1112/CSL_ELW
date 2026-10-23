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
