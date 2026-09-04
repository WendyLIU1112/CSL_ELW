# CSL-calibrated data-generating process used in the manuscript simulations
#
# This file is self-contained and uses only base R and stats. Prepare a
# calibration object with make_csl_simulation_base(), then generate a sample
# with generate_csl_simulation_data(). The restricted CSL data are not included.

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
  X.main <- as.data.frame(X.std, check.names = TRUE)
  ps.full.formula <- stats::as.formula(
    paste("A ~", paste(names(X.main), collapse = " + "))
  )
  ps.full.fit <- suppressWarnings(stats::glm(
    ps.full.formula, data = data.frame(A = A, X.main), family = stats::binomial()
  ))
  ## The simulation calibration retains the entire prespecified basis.  There
  ## is deliberately no AIC or other data-adaptive variable selection.
  ps.selected.fit <- ps.full.fit
  ps.selected.terms <- attr(stats::terms(ps.selected.fit), "term.labels")
  ell <- as.vector(stats::predict(ps.selected.fit, type = "link"))
  shadow.cols <- grepl("^Sex", colnames(X.std))
  X.no.shadow <- X.std[, !shadow.cols, drop = FALSE]
  ps.no.shadow.fit <- suppressWarnings(stats::glm.fit(
    cbind("(Intercept)" = 1, X.no.shadow), A, family = stats::binomial()
  ))
  ell.no.shadow <- as.vector(ps.no.shadow.fit$linear.predictors)

  outcome.data <- data.frame(Y = Y[A == 0], X.main[A == 0, , drop = FALSE])
  outcome.full.formula <- stats::as.formula(
    paste("Y ~", paste(names(X.main), collapse = " + "))
  )
  outcome.full.fit <- suppressWarnings(stats::lm(
    outcome.full.formula, data = outcome.data
  ))
  outcome.selected.fit <- outcome.full.fit
  outcome.selected.terms <- attr(stats::terms(outcome.selected.fit), "term.labels")
  m0 <- as.numeric(stats::predict(
    outcome.selected.fit,
    newdata = data.frame(X.main, check.names = TRUE)
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
    X.no.shadow = X.no.shadow,
    X.expanded = X.expanded,
    ell = ell,
    ell.no.shadow = ell.no.shadow,
    m0 = m0,
    sigma0 = sigma0,
    eps.pool = eps.pool,
    lbw.risk = lbw.risk,
    risk.percentile = risk.percentile,
    smoking.prevalence = mean(A == 1),
    ps.selected.terms = ps.selected.terms,
    ps.selected.formula = stats::formula(ps.selected.fit),
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
                                         treatment.nonlinear = FALSE,
                                         overlap.step = FALSE,
                                         overlap.low.prob = 1e-4,
                                         overlap.high.quantile = 0.90,
                                         seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  is.mnar <- gamma0 != 0 || gamma1 != 0
  id <- sample.int(nrow(base$dat), n, replace = TRUE)
  if (is.mnar) {
    ## LMZ-compatible CSL design. Infant sex is resampled jointly with the
    ## remaining CSL covariates and acts as the shadow variable: it shifts the
    ## outcome distribution but is excluded from h(U).
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

    ## If q_0(y|Z) is normal and g(Y(1)) is linear, exponential tilting
    ## leaves both complete-case outcome distributions normal. This matches
    ## the working LMZ likelihood while retaining the CSL outcome scale.
    log.m <- gamma * (mu0 + tau - y.center) / sigma0 + 0.5 * gamma^2
    h <- b * u
    a <- calibrate_logit_intercept(h + log.m, base$smoking.prevalence)
    p.treated <- stats::plogis(a + h + log.m)
    T <- stats::rbinom(n, 1L, p.treated)
    Y0 <- stats::rnorm(n, mu0 + T * gamma * sigma0, sigma0)
    Y1 <- Y0 + tau
    Y <- ifelse(T == 1L, Y1, Y0)
    pi.true <- stats::plogis(a + h + gamma * (Y1 - y.center) / sigma0)
    r <- rank(mu0, ties.method = "average") / n

    return(data.frame(
      id = id, Z = Z, T = T, Y = Y, Y0 = Y0, Y1 = Y1,
      pi.true = pi.true, rho1.true = pi.true, rho0.true = 1 - pi.true,
      ell = ell, overlap_stratum = 0, risk.percentile = r,
      delta = tau, check.names = FALSE
    ) |>
      cbind(as.data.frame(X, check.names = FALSE)) |>
      structure(X.expanded = X.expanded, X.main = X.no.shadow,
                ps.selected.terms = colnames(X.no.shadow), intercept.a = a,
                b = b, gamma0 = gamma0, gamma1 = gamma1,
                outcome.nonlinearity = 0, effect.nonlinearity = 0,
                lmz.compatible = TRUE, shadow.variable = "Sex"))
  }
  ell <- base$ell[id]
  X <- base$X[id, , drop = FALSE]
  X.expanded <- base$X.expanded[id, , drop = FALSE]
  m0 <- base$m0[id]
  sampled.lbw.risk <- base$lbw.risk[id]
  r <- rank(sampled.lbw.risk, ties.method = "average") / n
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
  overlap.stratum <- as.numeric(ell > stats::quantile(
    base$ell, overlap.high.quantile, names = FALSE
  ))
  if (overlap.step && !is.mnar) {
    share.high <- mean(overlap.stratum)
    high.prob <- (base$smoking.prevalence -
                    (1 - share.high) * overlap.low.prob) / share.high
    high.prob <- min(max(high.prob, overlap.low.prob), 1 - 1e-8)
    a <- stats::qlogis(overlap.low.prob)
    treatment.lp <- (stats::qlogis(high.prob) - a) * overlap.stratum
  } else {
    treatment.lp <- b * ell
    if (treatment.nonlinear) {
      ## Fixed nonlinear MAR mechanism.  These terms are not included in the
      ## main-effects PLE fitted below, but can be learned by the tree methods.
      age <- if ("Momage" %in% colnames(X)) X[, "Momage"] else 0
      bmi <- if ("PreBMI" %in% colnames(X)) X[, "PreBMI"] else 0
      height <- if ("Height" %in% colnames(X)) X[, "Height"] else 0
      treatment.lp <- treatment.lp +
        0.60 * (age^2 - 1) + 0.60 * sin(bmi) + 0.50 * age * bmi -
        0.40 * as.numeric(height > 0.75)
    }
  }
  if (!overlap.step && isTRUE(all.equal(b, 1)) && gamma0 == 0 && gamma1 == 0) {
    a <- 0
  } else if (!overlap.step || is.mnar) {
    a <- calibrate_logit_intercept(
      treatment.lp + gamma1 * z1 + gamma0 * z0,
      base$smoking.prevalence
    )
  }
  if (overlap.step && !is.mnar) {
    ## The intercept and high-stratum coefficient above already preserve the
    ## target prevalence; no replication-specific performance tuning is used.
    a <- a
  } else {
    a <- a
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
    overlap_stratum = overlap.stratum,
    risk.percentile = r,
    delta = delta,
    check.names = FALSE
  ) |>
    cbind(as.data.frame(X, check.names = FALSE)) |>
    structure(X.expanded = X.expanded,
              X.main = X,
              ps.selected.terms = base$ps.selected.terms,
              intercept.a = a, b = b, gamma0 = gamma0, gamma1 = gamma1,
              outcome.nonlinearity = outcome.nonlinearity,
              effect.nonlinearity = effect.nonlinearity)
}
