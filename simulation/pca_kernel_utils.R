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
                               grid_bins = 45L,
                               kernel = c("gaussian", "epanechnikov"),
                               bandwidth_multiplier = 1) {
  kernel <- match.arg(kernel)
  if (is.null(weights)) weights <- rep(1, nrow(X_train))
  pca <- weighted_pca_transform(
    X_train, X_test, weights = weights, n_components = n_components
  )
  z_train <- pca$train
  z_test <- pca$test
  dimension <- ncol(z_train)
  ## The bandwidth rule is fixed by the nominal sample size; multiplier sums
  ## fluctuate across draws but do not redefine the bootstrap sample size.
  bandwidth <- bandwidth_multiplier *
    normal_reference_bandwidth(nrow(X_train), dimension)

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
  standardized_distance <- outer(centers, centers, "-") / bandwidth
  kernel_matrix <- if (kernel == "gaussian") {
    exp(-0.5 * standardized_distance^2)
  } else {
    matrix(
      0.75 * pmax(0, as.numeric(1 - standardized_distance^2)),
      nrow = nrow(standardized_distance), ncol = ncol(standardized_distance)
    )
  }
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
    grid_bins = grid_bins,
    kernel = kernel
  )
}
