requested_steps <- "estimation"
run_step <- function(step) step %in% requested_steps

if (run_step("estimation")) local({
  args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed = TRUE)
  file_arg <- grep("^--file=", args, value = TRUE)
  script_dir <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[[1L]]))) else getwd()
  source(file.path(script_dir, "csl_smoking_all_functions.R"))
  source(file.path(script_dir, "pca_kernel_utils.R"))
  
  dat <- read.csv(
    Sys.getenv("CSL_SIM_DATA", file.path(script_dir, "data", "csl_smoking_birthweight_complete_case.csv")),
    check.names = FALSE
  )
  base <- make_csl_simulation_base(dat)
  R <- as.integer(Sys.getenv("SIM_R", "1000"))
  iter.start <- as.integer(Sys.getenv("SIM_ITER_START", "1"))
  n <- as.integer(Sys.getenv("SIM_N", "10000"))
  detected.cores <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected.cores)) {
    detected.cores <- suppressWarnings(as.integer(system("getconf _NPROCESSORS_ONLN", intern = TRUE)))
  }
  cores <- as.integer(Sys.getenv("SIM_CORES", as.character(max(1L, detected.cores - 1L))))
  out.dir <- Sys.getenv("SIM_OUT", file.path(script_dir, "results", "estimation-mar"))
  gbm.n.trees <- as.integer(Sys.getenv("GBM_TREES", "500"))
  gbm.depth <- as.integer(Sys.getenv("GBM_DEPTH", "2"))
  gbm.min.node <- as.integer(Sys.getenv("GBM_MIN_NODE", "100"))
  gbm.shrinkage <- as.numeric(Sys.getenv("GBM_SHRINKAGE", "0.02"))
  gbm.bag <- as.numeric(Sys.getenv("GBM_BAG", "0.5"))
  mnar.gamma <- as.numeric(Sys.getenv("SIM_MNAR_GAMMA", "-0.20"))
  instrument.effect <- as.numeric(Sys.getenv("SIM_INSTRUMENT_EFFECT", "300"))
  heteroskedastic <- tolower(Sys.getenv("SIM_HETEROSKEDASTIC", "false")) %in%
    c("true", "t", "1", "yes", "y")
  dir.create(out.dir, recursive = TRUE, showWarnings = FALSE)
  
  scenarios <- data.frame(
    scenario = c("MAR", "MAR", "MAR", "MAR--NL", "MAR--NL", "MAR--NL"),
    b = rep(1, 6),
    gamma1 = rep(0, 6),
    ps.method = c("PLE", "GBM", "honest_tree", "PLE", "GBM", "honest_tree"),
    stringsAsFactors = FALSE
  )
  job.filter <- Sys.getenv("SIM_FILTER", "")
  if (nzchar(job.filter)) {
    wanted <- strsplit(job.filter, ",", fixed = TRUE)[[1]]
    scenario.key <- paste(scenarios$scenario, scenarios$ps.method, sep = ":")
    scenarios <- scenarios[scenario.key %in% wanted, , drop = FALSE]
  }
  jobs <- merge(scenarios, data.frame(iter = seq.int(iter.start, length.out = R)), all = TRUE)
  
  one <- function(j) {
    sc <- jobs[j, ]
    z <- run_one_csl_simulation(
      base = base, iter = sc$iter, n = n, b = sc$b,
      gamma1 = sc$gamma1, probs = c(0.10, 0.25, 0.50, 0.75),
      seed.offset = 2026081100L,
      heteroskedastic = heteroskedastic,
      outcome.nonlinearity = 0,
      effect.nonlinearity = 0,
      treatment.nonlinear = sc$scenario == "MAR--NL",
      instrument.effect = instrument.effect,
      overlap.step = FALSE,
      overlap.low.prob = 0.001,
      overlap.high.quantile = 0.95,
      include_mnar_sipw = sc$scenario == "MNAR",
      include_mnar_elw = sc$scenario == "MNAR",
      gbm.n.trees = gbm.n.trees,
      gbm.interaction.depth = gbm.depth,
      gbm.min.node.size = gbm.min.node,
      gbm.shrinkage = gbm.shrinkage,
      gbm.bag.fraction = gbm.bag,
      ps.method = sc$ps.method
    )
    z$estimates$scenario <- sc$scenario
    z$diagnostics$scenario <- sc$scenario
    z
  }
  
  if (.Platform$OS.type != "windows" && cores > 1L) {
    res <- parallel::mclapply(
      seq_len(nrow(jobs)), one, mc.cores = cores,
      mc.preschedule = FALSE, mc.set.seed = FALSE
    )
  } else if (cores > 1L) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Install package 'future.apply' for parallel execution on Windows.")
    }
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)
    res <- future.apply::future_lapply(seq_len(nrow(jobs)), one, future.seed = FALSE)
  } else {
    res <- lapply(seq_len(nrow(jobs)), one)
  }
  est <- simulation_bind_rows(lapply(res, `[[`, "estimates"))
  dia <- simulation_bind_rows(lapply(res, `[[`, "diagnostics"))
  write.csv(est, file.path(out.dir, "estimates_all.csv"), row.names = FALSE)
  write.csv(dia, file.path(out.dir, "diagnostics_all.csv"), row.names = FALSE)
  
  keep.est <- est[est$estimand %in% c("qte", "lbw_rd"), ]
  keep.est$target <- ifelse(
    keep.est$estimand == "qte",
    sprintf("QTE(%.2f)", keep.est$tau),
    "DTE(2500-)"
  )
  keep.est$error.value <- keep.est$estimate - keep.est$truth
  group <- interaction(keep.est$scenario, keep.est$ps.method,
                       keep.est$method, keep.est$target, drop = TRUE)
  scalar <- do.call(rbind, lapply(split(keep.est, group), function(x) {
    ok <- is.finite(x$error.value)
    e <- x$error.value[ok]
    bias <- mean(e)
    sd.r <- sqrt(mean((e - bias)^2))  # Empirical SD, divisor m in every scenario.
    data.frame(
      scenario = x$scenario[1], ps.method = x$ps.method[1],
      method = x$method[1], target = x$target[1],
      truth = mean(x$truth), bias = bias, sd = sd.r,
      rmse = sqrt(mean(e^2)), mae = stats::median(abs(e)),
      unavailable = mean(!ok), n.available = sum(ok),
      identity.error = abs(mean(e^2) - bias^2 - sd.r^2)
    )
  }))
  
  dgroup <- interaction(dia$scenario, dia$ps.method, dia$method, drop = TRUE)
  cdf <- do.call(rbind, lapply(split(dia, dgroup), function(x) data.frame(
    scenario = x$scenario[1], ps.method = x$ps.method[1], method = x$method[1],
    mean.ise = mean(x$ise, na.rm = TRUE),
    unavailable = mean(!is.finite(x$ise)),
    min.observed.prob = min(x$ps.min, na.rm = TRUE)
  )))
  write.csv(scalar, file.path(out.dir, "scalar_summary.csv"), row.names = FALSE)
  write.csv(cdf, file.path(out.dir, "cdf_summary.csv"), row.names = FALSE)
  saveRDS(list(settings = list(R = R, n = n, scenarios = scenarios,
                               outcome.nonlinearity = 0,
                               effect.nonlinearity = 0,
                               seed.offset = 2026081100L),
               estimates = est, diagnostics = dia, scalar = scalar, cdf = cdf),
          file.path(out.dir, "simulation_results.rds"))
  print(scalar, row.names = FALSE)
  print(cdf, row.names = FALSE)
})

