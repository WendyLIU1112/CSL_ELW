requested_steps <- "lmz"
run_step <- function(step) step %in% requested_steps

if (run_step("lmz")) local({
  
  args <- gsub("~+~", " ", commandArgs(trailingOnly = FALSE), fixed = TRUE)
  file_arg <- grep("^--file=", args, value = TRUE)
  main_dir <- if (length(file_arg)) {
    dirname(normalizePath(sub("^--file=", "", file_arg[1L])))
  } else getwd()
  function_file <- file.path(main_dir, "csl_smoking_all_functions_lmz.R")
  if (!file.exists(function_file)) function_file <- file.path(main_dir, "csl_smoking_all_functions.R")
  source(function_file)
  
  input_file <- Sys.getenv(
    "CSL_SIM_DATA",
    file.path(main_dir, "data", "csl_smoking_birthweight_complete_case.csv")
  )
  dat <- read.csv(input_file, check.names = FALSE)
  base <- make_csl_simulation_base(dat)
  
  reps <- as.integer(Sys.getenv("LMZ_SIM_REPS", "1000"))
  n_sim <- as.integer(Sys.getenv("LMZ_SIM_N", "10000"))
  detected.cores <- parallel::detectCores(logical = FALSE)
  if (!is.finite(detected.cores)) {
    detected.cores <- suppressWarnings(as.integer(system("getconf _NPROCESSORS_ONLN", intern = TRUE)))
  }
  cores <- as.integer(Sys.getenv("LMZ_SIM_CORES", as.character(max(1L, detected.cores - 1L))))
  shadow_effect <- as.numeric(Sys.getenv("LMZ_SHADOW_EFFECT", "300"))
  mnar_gamma <- as.numeric(Sys.getenv("LMZ_MNAR_GAMMA", "-0.30"))
  output_dir <- Sys.getenv("LMZ_SIM_OUTPUT", file.path(main_dir, "results", "estimation-mnar-lmz"))
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  run_rep <- function(iter) {
    run_one_csl_simulation(
      base = base, iter = iter, n = n_sim, b = 1,
      gamma0 = 0, gamma1 = mnar_gamma,
      instrument.effect = shadow_effect,
      heteroskedastic = TRUE, ps.method = "PLE",
      include_mnar_sipw = TRUE, include_mnar_elw = TRUE,
      include_mnar_lql = FALSE
    )
  }
  
  if (.Platform$OS.type != "windows" && cores > 1L) {
    result <- parallel::mclapply(seq_len(reps), run_rep, mc.cores = cores,
                                 mc.preschedule = TRUE, mc.set.seed = FALSE)
  } else if (cores > 1L) {
    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop("Install package 'future.apply' for parallel execution on Windows.")
    }
    future::plan(future::multisession, workers = cores)
    on.exit(future::plan(future::sequential), add = TRUE)
    result <- future.apply::future_lapply(seq_len(reps), run_rep, future.seed = FALSE)
  } else {
    result <- lapply(seq_len(reps), run_rep)
  }
  
  estimates <- simulation_bind_rows(lapply(result, `[[`, "estimates"))
  diagnostics <- simulation_bind_rows(lapply(result, `[[`, "diagnostics"))
  summary <- summarize_csl_simulation(estimates, diagnostics)
  
  write.csv(estimates, file.path(output_dir, "simulation_estimates.csv"), row.names = FALSE)
  write.csv(diagnostics, file.path(output_dir, "simulation_diagnostics.csv"), row.names = FALSE)
  write.csv(summary$estimates, file.path(output_dir, "simulation_summary_estimates.csv"), row.names = FALSE)
  write.csv(summary$diagnostics, file.path(output_dir, "simulation_summary_diagnostics.csv"), row.names = FALSE)
  saveRDS(list(estimates = estimates, diagnostics = diagnostics, summary = summary),
          file.path(output_dir, "simulation_results.rds"))
  
  print(subset(summary$estimates, method %in% c("ELW-LMZ", "SIPW-LMZ")))
})

