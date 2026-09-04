# Empirical likelihood weighting for counterfactual distributions
#
# This file is self-contained. The main entry point is elw().

.elw_arm <- function(outcome, observed, rho, level, tol) {
  observed <- as.logical(observed)
  y <- as.numeric(outcome[observed])
  rho_observed <- as.numeric(rho[observed])
  n <- length(outcome)
  n_t <- sum(observed)
  r_t <- n_t / n

  if (!n_t || n_t == n) {
    stop("Each exposure level must contain at least one observed and one unobserved outcome.")
  }
  if (anyNA(y) || any(!is.finite(y))) stop("Observed outcomes must be finite.")
  if (length(rho) != n || anyNA(rho_observed) ||
      any(!is.finite(rho_observed)) || any(rho_observed < 0 | rho_observed > 1)) {
    stop("Selection probabilities must be finite values in [0, 1].")
  }

  lower <- min(rho_observed)
  upper <- r_t + (1 - r_t) * lower
  denominator <- function(alpha) {
    r_t + (1 - r_t) * rho_observed - alpha
  }
  estimating_function <- function(alpha) {
    sum((rho_observed - alpha) / denominator(alpha))
  }

  f_lower <- estimating_function(lower)
  if (abs(f_lower) <= tol) {
    alpha <- lower
  } else {
    upper_search <- upper - max(.Machine$double.eps, abs(upper) * .Machine$double.eps)
    f_upper <- estimating_function(upper_search)
    if (!is.finite(f_lower) || !is.finite(f_upper) || f_lower * f_upper > 0) {
      stop("The ELW calibration root could not be bracketed for exposure level ", level, ".")
    }
    alpha <- stats::uniroot(
      estimating_function, c(lower, upper_search),
      f.lower = f_lower, f.upper = f_upper, tol = tol
    )$root
  }

  denom <- 1 + ((1 - r_t) / r_t) *
    (rho_observed - alpha) / (1 - alpha)
  tilting_weight <- 1 / denom
  probability_mass <- tilting_weight / n_t
  probability_mass <- probability_mass / sum(probability_mass)
  tilting_weight <- n_t * probability_mass

  structure(list(
    level = level,
    observed = which(observed),
    outcome = y,
    selection_probability = rho_observed,
    alpha = alpha,
    tilting_weight = tilting_weight,
    probability_mass = probability_mass,
    effective_sample_size = 1 / sum(probability_mass^2),
    normalization_error = abs(mean(tilting_weight) - 1),
    calibration_error = abs(mean(tilting_weight * (rho_observed - alpha)))
  ), class = "elw_arm")
}

#' Fit empirical likelihood weights for two counterfactual distributions
#'
#' @param outcome Numeric outcome vector.
#' @param treatment Binary exposure vector coded 0 and 1.
#' @param propensity_score Optional P(T = 1 | X) vector. Use this under
#'   conditional exchangeability.
#' @param selection_probability Optional two-column matrix or list containing
#'   P(T = 0 | Y(0), X) and P(T = 1 | Y(1), X). Supply either this argument or
#'   propensity_score.
#' @param tol Numerical tolerance for the scalar calibration root.
#' @return An object of class elw_fit.
elw <- function(outcome, treatment, propensity_score = NULL,
                selection_probability = NULL, tol = 1e-10) {
  outcome <- as.numeric(outcome)
  treatment <- as.integer(treatment)
  n <- length(outcome)
  if (length(treatment) != n || anyNA(treatment) ||
      !all(treatment %in% c(0L, 1L))) {
    stop("treatment must have the same length as outcome and be coded 0 and 1.")
  }
  if (is.null(propensity_score) == is.null(selection_probability)) {
    stop("Supply exactly one of propensity_score or selection_probability.")
  }

  if (!is.null(propensity_score)) {
    propensity_score <- as.numeric(propensity_score)
    if (length(propensity_score) != n || anyNA(propensity_score) ||
        any(!is.finite(propensity_score)) ||
        any(propensity_score < 0 | propensity_score > 1)) {
      stop("propensity_score must contain n finite values in [0, 1].")
    }
    rho <- list(`0` = 1 - propensity_score, `1` = propensity_score)
    probability_source <- "propensity score"
  } else {
    if (is.data.frame(selection_probability) || is.matrix(selection_probability)) {
      if (nrow(selection_probability) != n || ncol(selection_probability) != 2L) {
        stop("selection_probability must have n rows and two columns.")
      }
      rho <- list(`0` = selection_probability[, 1L],
                  `1` = selection_probability[, 2L])
    } else if (is.list(selection_probability) && length(selection_probability) == 2L) {
      rho <- list(`0` = selection_probability[[1L]],
                  `1` = selection_probability[[2L]])
    } else {
      stop("selection_probability must be a two-column object or a list of two vectors.")
    }
    if (any(vapply(rho, length, integer(1L)) != n)) {
      stop("Each selection-probability vector must have length n.")
    }
    probability_source <- "supplied selection probability"
  }

  arms <- list(
    `0` = .elw_arm(outcome, treatment == 0L, rho[["0"]], 0L, tol),
    `1` = .elw_arm(outcome, treatment == 1L, rho[["1"]], 1L, tol)
  )
  structure(list(
    call = match.call(), n = n, treatment = treatment,
    probability_source = probability_source, arms = arms
  ), class = "elw_fit")
}

print.elw_fit <- function(x, ...) {
  cat("Empirical likelihood weighting fit\n")
  cat("Observations:", x$n, "\n")
  cat("Selection probabilities:", x$probability_source, "\n")
  for (level in c("0", "1")) {
    arm <- x$arms[[level]]
    cat("Level", level, ": n =", length(arm$outcome),
        ", alpha =", format(arm$alpha, digits = 5),
        ", ESS =", format(arm$effective_sample_size, digits = 6), "\n")
  }
  invisible(x)
}

elw_cdf <- function(object, y, level) {
  stopifnot(inherits(object, "elw_fit"), level %in% c(0, 1))
  arm <- object$arms[[as.character(level)]]
  vapply(y, function(cutoff) {
    sum(arm$probability_mass[arm$outcome <= cutoff])
  }, numeric(1L))
}

elw_quantile <- function(object, probs, level) {
  stopifnot(inherits(object, "elw_fit"), level %in% c(0, 1))
  if (anyNA(probs) || any(probs < 0 | probs > 1)) stop("probs must lie in [0, 1].")
  arm <- object$arms[[as.character(level)]]
  ord <- order(arm$outcome)
  y <- arm$outcome[ord]
  cumulative_mass <- cumsum(arm$probability_mass[ord])
  vapply(probs, function(probability) {
    y[which(cumulative_mass >= probability)[1L]]
  }, numeric(1L))
}

elw_probability <- function(object, interval, level) {
  stopifnot(inherits(object, "elw_fit"), level %in% c(0, 1))
  if (length(interval) != 2L || interval[1L] > interval[2L]) {
    stop("interval must be c(lower, upper).")
  }
  arm <- object$arms[[as.character(level)]]
  sum(arm$probability_mass[
    arm$outcome >= interval[1L] & arm$outcome < interval[2L]
  ])
}

elw_dte <- function(object, y) {
  elw_cdf(object, y, 1L) - elw_cdf(object, y, 0L)
}

elw_qte <- function(object, probs) {
  elw_quantile(object, probs, 1L) - elw_quantile(object, probs, 0L)
}

elw_pte <- function(object, interval) {
  elw_probability(object, interval, 1L) -
    elw_probability(object, interval, 0L)
}
