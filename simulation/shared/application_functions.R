validate_probability <- function(p, name="probability") { stopifnot(all(is.finite(p)), all(p >= 0 & p <= 1)); as.numeric(p) }
gbm_spec <- list(n_trees=1000L, interaction_depth=3L, min_node_size=50L, shrinkage=0.02, bag_fraction=0.8)
primary_ps_method <- "gbm"
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
  
