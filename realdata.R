#' Reproduce the CBASP data application
#'
#' Run this file with Rscript realdata.R, or source it and call rd_reproduce().
#' Defaults match the supplied manuscript: 453 fixed observed baseline profiles,
#' noise SD 2, 1000 independent simulation seeds, quadratic age splines with a
#' fixed interior knot zero, and four age strata for comparison allocations.
#' Table SD is the empirical SD of treatment estimates across replications.
#' This is a simulation calibrated to real data; observed-outcome analyses are
#' exported separately. No published table numbers are inserted into outputs.
#'
#' @examples
#' source("realdata.R")
#' ans <- rd_reproduce(data_file = "CBASP.csv", repetitions = 1000)
#' ans$table
#' rd_initialize("CBASP.csv")
#' one_seed <- delta(1) # optional full inference, including bootstrap B = 500

# Resolve the source location before leaving source()'s evaluation frame. This
# lets the script find its neighboring CSV when launched from another folder.
.rd_source_file <- tryCatch(sys.frame(1)$ofile, error = function(e) NULL)
if (is.null(.rd_source_file)) {
  .rd_file_arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(.rd_file_arg)) .rd_source_file <- sub("^--file=", "", .rd_file_arg[1])
}
.rd_source_dir <- if (is.null(.rd_source_file)) getwd() else
  dirname(normalizePath(.rd_source_file, mustWork = FALSE))

#' Prepare the observed baseline profiles for the trial procedures
#'
#' Read the original CBASP CSV, retain NEF and combination participants, and
#' standardize age and HAMD24 using all eligible baseline profiles. Install the
#' legacy data objects in this script's environment for delta() and test methods.
#' @param data_file Path to CBASP.csv with A2, AGE, HAMD24, and FinalHAMD columns.
#' @return List containing profiles, original treatment codes, eligibility and
#'   complete-outcome counts, centering/scaling constants, and input checksum.
#' @details A2 = -1 denotes NEF, 0 denotes combination, and 1 denotes CBT alone.
#'   Missing outcomes are retained for simulations but excluded from observed
#'   outcome fits. Missing baseline covariates cause an explicit error.
rd_initialize <- function(data_file = file.path(.rd_source_dir, "CBASP.csv")) {
  if (!file.exists(data_file)) stop("Cannot find CBASP.csv: ", data_file)
  raw <- read.csv(data_file, check.names = TRUE)
  needed <- c("A2", "AGE", "HAMD24", "FinalHAMD")
  if (!all(needed %in% names(raw)))
    stop("Missing CSV columns: ", paste(setdiff(needed, names(raw)), collapse = ", "))
  if (!all(vapply(raw[needed], is.numeric, logical(1))))
    stop("A2, AGE, HAMD24, and FinalHAMD must be numeric.")
  if (anyNA(raw$A2) || any(!raw$A2 %in% c(-1, 0, 1)))
    stop("Expected A2 codes -1 (NEF), 0 (combination), and 1 (CBT).")
  d <- raw[raw$A2 != 1, , drop = FALSE]
  if (!nrow(d) || any(!is.finite(as.matrix(d[c("AGE", "HAMD24")]))))
    stop("Eligible participants must have finite AGE and HAMD24 values.")
  if (any(!is.na(d$FinalHAMD) & !is.finite(d$FinalHAMD)))
    stop("FinalHAMD must be finite or missing.")
  original_A2 <- d$A2
  d$A2 <- d$A2 + 1 # 0 = NEF, 1 = combination; same coding as the original fit.
  centers <- vapply(d[c("AGE", "HAMD24")], mean, numeric(1))
  scales <- vapply(d[c("AGE", "HAMD24")], sd, numeric(1))
  if (any(scales <= 0)) stop("Both baseline covariates must have positive SD.")
  d$AGE2 <- d$AGE^2
  d$HAMD242 <- d$HAMD24^2
  d$AGE_STD <- (d$AGE - centers["AGE"]) / scales["AGE"]
  d$HAMD24_STD <- (d$HAMD24 - centers["HAMD24"]) / scales["HAMD24"]
  d$AGE_STD2 <- d$AGE_STD^2
  env <- environment(rd_initialize)
  # mgcv recognizes the bare s() formula special; bind it without attaching packages.
  assign("s", mgcv::s, envir = env)
  assign("df2", d, envir = env)
  assign("sample_size", nrow(d), envir = env)
  assign("BLCOV", d[c("AGE", "HAMD24")], envir = env)
  assign("BLCOVSTD", d[c("AGE_STD", "HAMD24_STD")], envir = env)
  list(profiles = d, original_A2 = original_A2, eligible_n = nrow(d),
       observed_n = sum(complete.cases(d[needed])), centers = centers,
       scales = scales, data_file = normalizePath(data_file),
       data_md5 = unname(tools::md5sum(data_file)))
}

#' Analyze observed outcomes and reproduce the smooth-effect figure
#'
#' Fit the quadratic calibration regression and an additive GAM to the observed
#' outcomes. Export coefficients, model summaries, and SmoothEffect.pdf/png.
#' @param profile List returned by rd_initialize().
#' @param output_dir Existing directory for generated model and figure files.
#' @return List containing quadratic and additive model objects and the number
#'   of complete observations used. These are actual observed-data analyses.
rd_fit_observed <- function(profile, output_dir) {
  d <- profile$profiles
  d <- d[complete.cases(d[c("A2", "AGE", "HAMD24", "FinalHAMD")]), ]
  fit <- lm(FinalHAMD ~ A2 + AGE + I(AGE^2) + HAMD24, data = d)
  additive <- mgcv::gam(FinalHAMD ~ A2 + s(AGE) + s(HAMD24), data = d)
  cmat <- summary(fit)$coefficients
  write.csv(data.frame(term = rownames(cmat), cmat, row.names = NULL,
                       check.names = FALSE),
            file.path(output_dir, "observed_quadratic_coefficients.csv"),
            row.names = FALSE)
  capture.output(summary(fit), file = file.path(output_dir, "observed_quadratic_fit.txt"))
  capture.output(summary(additive), file = file.path(output_dir, "observed_additive_fit.txt"))
  for (ext in c("pdf", "png")) {
    path <- file.path(output_dir, paste0("SmoothEffect.", ext))
    if (ext == "pdf") pdf(path, width = 10, height = 5) else
      png(path, width = 1800, height = 900, res = 180, type = "cairo")
    tryCatch({
      par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))
      plot(additive, select = 1, shade = FALSE, xlab = "Age", ylab = "s(AGE)")
      plot(additive, select = 2, shade = FALSE, xlab = "Baseline HAMD24",
           ylab = "s(HAMD24)")
    }, finally = dev.off())
  }
  list(quadratic = fit, additive = additive, n = nrow(d))
}

#' Obtain the six reported treatment estimates for one simulated trial seed
#'
#' Call the existing working-model estimators before dataset-specific variance
#' calculations. Every method sets the same seed independently, as in delta().
#' @param seed Positive integer replication seed.
#' @param sigma Positive normal noise SD; manuscript default 2.
#' @param strata Number of standard-normal age bins for comparator allocations.
#' @return Named length-six numeric vector of treatment estimates; four table
#'   methods followed by the proposed two-sample and partially linear methods.
#' @details Requires rd_initialize() first. Repeated trials independently
#'   randomize treatment and generate noise on the same fixed baseline profiles.
#'   Empirical SD therefore measures repeated-trial variation conditional on
#'   those profiles. No bootstrap is needed to calculate this empirical SD.
rd_estimates_one <- function(seed, sigma = 2, strata = 4) {
  c("Proposed-lm" = test2.ind(seed, sample_size, 500, sigma, estimate_only = TRUE),
    "Ma2015lm" = fftest.bnz2(seed, sample_size, strata, sigma, estimate_only = TRUE),
    "Ma2022mr" = ftest2.p.c2(seed, sample_size, strata, sigma, estimate_only = TRUE),
    "Ye2018plm" = testye.bnz2.c2(seed, sample_size, strata, sigma, estimate_only = TRUE),
    "Proposed-t-test" = test1.ind(seed, sample_size, 500, sigma, estimate_only = TRUE),
    "Proposed-plm" = unname(test3.ind(seed, sample_size, sigma)[1]))
}

#' Simulate independent trials for the real-data-based table
#'
#' Collect treatment estimates using consecutive deterministic seeds. Save one
#' row per method and seed, without running unneeded bootstrap variance fits.
#' @param repetitions Integer number of trials, at least 2; default 1000.
#' @param seed_start Positive first seed; default 1.
#' @param sigma Positive normal noise SD; default 2.
#' @param strata Integer number of comparator age bins; default 4.
#' @param verbose Logical; print progress every 100 trials.
#' @return Data frame with columns seed, sigma, method, and estimate.
rd_run_simulations <- function(repetitions = 1000, seed_start = 1,
                               sigma = 2, strata = 4, verbose = TRUE) {
  seeds <- seq.int(seed_start, length.out = repetitions)
  draws <- vector("list", repetitions)
  for (j in seq_along(seeds)) {
    values <- rd_estimates_one(seeds[j], sigma, strata)
    if (any(!is.finite(values))) stop("Nonfinite estimate at seed ", seeds[j])
    draws[[j]] <- data.frame(seed = seeds[j], sigma = sigma,
                             method = names(values), estimate = unname(values))
    if (verbose && (j %% 100L == 0L || j == repetitions))
      message("Completed ", j, " of ", repetitions, " trial seeds.")
  }
  do.call(rbind, draws)
}

#' Read treatment estimates from existing legacy seed.RData files
#'
#' Load each file in a private environment and extract the six reported methods
#' from the length-60 delta() result at the requested noise setting.
#' @param input_dir Directory containing seed1.RData, seed2.RData, etc.
#' @param expected Positive seed numbers to use; extra files are ignored.
#' @param sigma One of sqrt(3), 2, or sqrt(5), matching legacy delta() settings.
#' @return Data frame with seed, sigma, method, and estimate columns.
#' @details Legacy unnamed vectors do not record knots, strata, sample size, or
#'   package versions. The caller must verify that the batch used the desired
#'   settings; changing this script cannot change previously simulated results.
rd_read_legacy_seeds <- function(input_dir, expected = 1:1000, sigma = 2) {
  h <- which(abs(c(sqrt(3), 2, sqrt(5)) - sigma) < 1e-10)
  if (length(h) != 1L) stop("Legacy files support only noise SD sqrt(3), 2, or sqrt(5).")
  if (!dir.exists(input_dir)) stop("Seed directory does not exist: ", input_dir)
  files <- list.files(input_dir, pattern = "^seed[0-9]+[.]RData$", full.names = TRUE)
  ids <- as.integer(sub("^seed([0-9]+)[.]RData$", "\\1", basename(files)))
  if (anyDuplicated(ids)) stop("Duplicate seed IDs in ", input_dir)
  missing <- setdiff(expected, ids)
  if (length(missing)) stop("Missing seed files: ", paste(head(missing, 20), collapse = ", "))
  files <- files[match(expected, ids)]
  # Offsets follow delta()'s nine method blocks, with three noise settings each.
  starts <- c("Proposed-lm" = 10, "Ma2015lm" = 31, "Ma2022mr" = 43,
              "Ye2018plm" = 55, "Proposed-t-test" = 1, "Proposed-plm" = 19)
  widths <- c(3, 2, 2, 2, 3, 2)
  index <- starts + (h - 1L) * widths
  rows <- lapply(seq_along(files), function(j) {
    e <- new.env(parent = emptyenv())
    names_loaded <- load(files[j], envir = e)
    if ("res1" %in% names_loaded) value <- e$res1 else {
      candidates <- names_loaded[vapply(names_loaded, function(nm)
        is.numeric(e[[nm]]) && length(e[[nm]]) == 60L, logical(1))]
      if (length(candidates) != 1L)
        stop("Need one numeric length-60 delta() result in ", files[j])
      value <- e[[candidates]]
    }
    if (!is.numeric(value) || length(value) != 60L || any(!is.finite(value[index])))
      stop("Invalid treatment estimates in ", files[j])
    data.frame(seed = expected[j], sigma = sigma, method = names(starts),
               estimate = unname(value[index]))
  })
  do.call(rbind, rows)
}

#' Summarize repeated-trial estimates and build the paper table
#'
#' Compute each method's mean treatment estimate and empirical SD. The SE of
#' the mean is SD divided by sqrt(R); SD itself is not an average fitted SE.
#' @param estimates Data frame returned by rd_run_simulations() or the legacy reader.
#' @return List with an unrounded per-method summary and a two-row four-method
#'   table in the manuscript order. The first column identifies the statistic.
rd_summarize <- function(estimates) {
  order <- c("Proposed-lm", "Ma2015lm", "Ma2022mr", "Ye2018plm",
             "Proposed-t-test", "Proposed-plm")
  if (!all(order %in% estimates$method) || length(unique(estimates$sigma)) != 1L)
    stop("Need all six methods at one noise setting.")
  if (anyDuplicated(estimates[c("seed", "method")]) ||
      any(!is.finite(estimates$estimate))) stop("Duplicate or invalid estimates.")
  seed_sets <- lapply(order, function(m) sort(estimates$seed[estimates$method == m]))
  if (length(seed_sets[[1]]) < 2L ||
      !all(vapply(seed_sets, identical, logical(1), seed_sets[[1]])))
    stop("Each method must contain the same set of at least two seeds.")
  s <- do.call(rbind, lapply(order, function(m) {
    v <- estimates$estimate[estimates$method == m]
    data.frame(method = m, repetitions = length(v), sigma = unique(estimates$sigma),
               mean_estimate = mean(v), empirical_sd = sd(v),
               mc_se_mean = sd(v) / sqrt(length(v)))
  }))
  tab <- data.frame(Statistic = c("Treatment estimate", "SD"),
                    check.names = FALSE)
  for (j in 1:4) tab[[s$method[j]]] <- c(s$mean_estimate[j], s$empirical_sd[j])
  list(summary = s, table = tab)
}

#' Export the table and its supplementary results
#'
#' Write full-precision CSV/RDS summaries, a booktabs LaTeX table fragment, an
#' HTML preview, and PDF/PNG table images. Values are computed from estimates.
#' @param summaries List returned by rd_summarize().
#' @param output_dir Existing directory for generated files.
#' @return Invisibly, the unrounded four-method table.
rd_write_table <- function(summaries, output_dir) {
  tab <- summaries$table
  s <- summaries$summary
  write.csv(tab, file.path(output_dir, "Table_D1.csv"), row.names = FALSE)
  write.csv(s, file.path(output_dir, "all_methods_summary.csv"), row.names = FALSE)
  write.csv(s[5:6, ], file.path(output_dir, "supplementary_methods.csv"), row.names = FALSE)
  caption <- sprintf(paste0("Treatment effect estimates under linear adjustment. ",
    "Mean estimates and empirical standard deviations across %d independent ",
    "trials on the fixed observed baseline profiles; noise SD = %g."),
    s$repetitions[1], s$sigma[1])
  cells <- rbind(sprintf("%.3f", as.numeric(tab[1, -1])),
                 sprintf("%.3f", as.numeric(tab[2, -1])))
  tex <- c("\\begin{table}[!h]", paste0("\\caption{", caption, "}"),
    "\\label{tb:realdata}", "\\centering", "\\begin{tabular}{lcccc}",
    "\\toprule", paste0("Methods & ", paste(sprintf("\\texttt{%s}", names(tab)[-1]),
                                           collapse = " & "), " \\\\"),
    "\\midrule", paste0("$\\hat u_1-\\hat u_2$ & ", paste(cells[1, ], collapse = " & "), " \\\\"),
    paste0("SD & ", paste(cells[2, ], collapse = " & "), " \\\\"),
    "\\bottomrule", "\\end{tabular}", "\\end{table}")
  writeLines(tex, file.path(output_dir, "Table_D1.tex"))
  html <- c("<!doctype html><meta charset='utf-8'><title>Table D.1</title>",
    "<style>body{font-family:serif;max-width:900px;margin:40px auto}table{border-collapse:collapse;width:100%}th,td{padding:12px;text-align:center}thead{border-top:2px solid;border-bottom:1px solid}tbody{border-bottom:2px solid}td:first-child{text-align:left}</style>",
    paste0("<p>", caption, "</p><table><thead><tr><th>Methods</th>",
            paste0("<th>", names(tab)[-1], "</th>", collapse = ""), "</tr></thead><tbody>"))
  for (j in 1:2) html <- c(html, paste0("<tr><td>", tab$Statistic[j], "</td>",
                                  paste0("<td>", cells[j, ], "</td>", collapse = ""), "</tr>"))
  writeLines(c(html, "</tbody></table>"), file.path(output_dir, "Table_D1.html"))
  for (ext in c("pdf", "png")) {
    path <- file.path(output_dir, paste0("Table_D1.", ext))
    if (ext == "pdf") pdf(path, width = 10, height = 2.6) else
      png(path, width = 1800, height = 470, res = 180, type = "cairo")
    tryCatch({
      par(mar = c(1, 0.5, 1, 0.5), family = "serif")
      plot.new(); plot.window(xlim = c(0, 5), ylim = c(0, 3))
      x <- seq(0.5, 4.5, by = 1)
      segments(0.1, c(2.75, 1.85, 0.15), 4.9, c(2.75, 1.85, 0.15))
      text(x, 2.3, c("Methods", names(tab)[-1]), cex = 1.05)
      text(x[1], c(1.3, 0.5), c(expression(hat(u)[1] - hat(u)[2]), "SD"))
      text(x[-1], 1.3, cells[1, ], cex = 1.1)
      text(x[-1], 0.5, cells[2, ], cex = 1.1)
    }, finally = dev.off())
  }
  invisible(tab)
}

#' Run the complete real-data example and table reproduction
#'
#' Prepare the CSV, analyze observed outcomes, generate repeated simulated
#' trials, and export the table, supplementary summaries, figures, and audit
#' settings. Called automatically when realdata.R is executed with Rscript.
#' @param data_file Original CBASP.csv path; default beside this script.
#' @param repetitions Integer trial count, at least 2; manuscript default 1000.
#' @param seed_start Positive first simulation seed; default 1.
#' @param sigma Positive noise SD; manuscript default 2.
#' @param strata Integer comparator age-bin count; manuscript default 4.
#' @param output_dir Output directory; default realdata_results beside this script.
#' @param results_dir Optional directory of existing length-60 seed.RData files.
#'   NULL generates fresh trials; old files must have compatible settings.
#' @param verbose Logical; report progress during fresh simulation.
#' @return Invisibly, a list containing estimates, summary, table, observed
#'   models, settings, and output_dir. Also saves reproduction.rds and CSV files.
#' @details Published numbers are Monte Carlo realizations. A fresh run need
#'   not match their final decimals if original seeds or package versions differ.
#'   Simulation outcomes use raw AGE/HAMD24 and the manuscript coefficients;
#'   standardized covariates are used in randomization and working regressions.
rd_reproduce <- function(data_file = file.path(.rd_source_dir, "CBASP.csv"),
                         repetitions = 1000, seed_start = 1, sigma = 2, strata = 4,
                         output_dir = file.path(.rd_source_dir, "realdata_results"),
                         results_dir = NULL, verbose = TRUE) {
  for (nm in c("repetitions", "seed_start", "strata")) {
    v <- get(nm)
    if (!is.numeric(v) || length(v) != 1L || !is.finite(v) || v != floor(v) || v < 1)
      stop(nm, " must be a positive integer.")
  }
  if (repetitions < 2 || strata < 2) stop("Need at least two repetitions and strata.")
  if (seed_start + repetitions - 1 > .Machine$integer.max)
    stop("The requested seeds exceed the R integer range.")
  if (!is.numeric(sigma) || length(sigma) != 1L || !is.finite(sigma) || sigma <= 0)
    stop("sigma must be a positive finite noise SD.")
  pkgs <- c("mgcv", "splines2", if (is.null(results_dir)) "carat")
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required R packages first: ", paste(missing, collapse = ", "))
  profile <- rd_initialize(data_file)
  if (profile$eligible_n != 453L)
    warning("This CSV has ", profile$eligible_n, " eligible profiles; the manuscript uses 453.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  if (!dir.exists(output_dir)) stop("Cannot create output directory: ", output_dir)
  observed <- rd_fit_observed(profile, output_dir)
  seeds <- seq.int(seed_start, length.out = repetitions)
  if (is.null(results_dir)) estimates <- rd_run_simulations(repetitions, seed_start, sigma, strata, verbose) else {
    warning("Reading legacy vectors: verify the original knots, strata, sample size and package versions; these were not saved in the files.")
    estimates <- rd_read_legacy_seeds(results_dir, seeds, sigma)
  }
  summaries <- rd_summarize(estimates)
  rd_write_table(summaries, output_dir)
  settings <- list(
    mode = if (is.null(results_dir)) "fresh manuscript-settings simulation" else "legacy files; original settings unverified",
    data_file = profile$data_file, data_md5 = profile$data_md5,
    eligible_profiles = profile$eligible_n, observed_complete_outcomes = observed$n,
    repetitions = repetitions, seeds = seeds, noise_sd = sigma,
    u1 = 20.32581, u2 = 25.196228, true_contrast = 20.32581 - 25.196228,
    age_linear = -0.972179, age_quadratic = 0.011102, hamd_linear = 0.391615,
    proposed = list(order = 3, internal_knots = 0, rho = 0.9,
                    boundary_knots = range(profile$profiles$AGE_STD)),
    comparator_strata = strata, baseline_centers = profile$centers,
    baseline_scales = profile$scales,
    sd_definition = "sd(treatment estimates across independent simulation seeds)",
    results_dir = results_dir,
    legacy_bootstrap_B = 500,
    table_bootstrap = "not required for point estimates or empirical SD",
    package_versions = setNames(vapply(pkgs, function(p) as.character(packageVersion(p)), character(1)), pkgs),
    R_version = R.version.string, RNG_kind = RNGkind())
  ans <- c(summaries, list(estimates = estimates, observed = observed,
                           settings = settings, output_dir = normalizePath(output_dir)))
  saveRDS(ans, file.path(output_dir, "reproduction.rds"))
  write.csv(estimates, file.path(output_dir, "trial_estimates.csv"), row.names = FALSE)
  write.csv(data.frame(seed = seeds), file.path(output_dir, "simulation_seeds.csv"), row.names = FALSE)
  capture.output(str(settings, max.level = 3), file = file.path(output_dir, "settings.txt"))
  capture.output(sessionInfo(), file = file.path(output_dir, "sessionInfo.txt"))
  if (verbose) message("Results saved in ", normalizePath(output_dir))
  invisible(ans)
}

#' Parse command-line arguments for the reproducible example
#'
#' Translate --name=value options into arguments for rd_reproduce(). Unknown or
#' repeated options cause an error rather than silently changing the analysis.
#' @param args Character command-line arguments; defaults to trailing arguments.
#' @return Named argument list for do.call(rd_reproduce, ...).
rd_parse_args <- function(args = commandArgs(trailingOnly = TRUE)) {
  allowed <- c("data-file", "repetitions", "seed-start", "sigma", "strata",
               "output-dir", "results-dir")
  out <- list()
  for (a in args) {
    if (!grepl("^--[a-z-]+=.+$", a)) stop("Use command-line options of the form --name=value.")
    key <- sub("^--([^=]+)=.*$", "\\1", a)
    value <- sub("^--[^=]+=", "", a)
    if (!key %in% allowed) stop("Unknown option: ", key)
    key <- gsub("-", "_", key)
    if (key %in% names(out)) stop("Repeated option: ", key)
    out[[key]] <- if (key %in% c("repetitions", "seed_start", "sigma", "strata"))
      suppressWarnings(as.numeric(value)) else value
  }
  out
}

#' Assign treatments with local univariate spline balancing
#'
#' Evaluate fixed-knot splines once, then update a separate signed imbalance vector in
#' each interval as participants arrive.
#'
#' @param x Finite numeric vector of standardized age values in arrival order.
#' @param J Number of local intervals; must equal length(internal_knots) + 1. Default 2.
#' @param k B-spline degree; default 2 means order 3.
#' @param rho Probability favoring the treatment that reduces local imbalance; default
#'   0.9.
#' @param internal_knots Interior spline knots on the standardized age scale; default 0,
#'   as in the manuscript.
#' @return Numeric vector of length length(x), containing 0/1 treatment assignments.
#' @details Features are an intercept plus the marginal splines returned with intercept =
#'   FALSE. Boundary knots use the observed range. Ties use probability 0.5; draws
#'   consume the current RNG state.
covRcpp_uni_b <- function(x, J = 2, k = 2, rho = 0.9, internal_knots = 0) {
  # x: Univariate vector of covariates
  # J: Number of intervals
  # k: Spline degree; default 2 means order 3 (quadratic).

  n <- length(x)

  # The manuscript uses a fixed knot at standardized age zero, giving two
  # local intervals. Replacing one marginal spline by the intercept preserves
  # the complete univariate spline span; there is no tensor basis in this example.
  if (!is.numeric(x) || !length(x) || any(!is.finite(x)))
    stop("x must be a nonempty finite numeric vector.")
  if (!is.numeric(internal_knots) || any(!is.finite(internal_knots)))
    stop("internal_knots must be finite numeric values.")
  internal_knots <- sort(unique(internal_knots))
  if (length(J) != 1L || J != length(internal_knots) + 1L)
    stop("J must equal the number of internal knots plus one.")
  if (any(internal_knots <= min(x) | internal_knots >= max(x)))
    stop("Every internal knot must lie strictly inside the observed range.")
  if (length(rho) != 1L || !is.finite(rho) || rho < 0.5 || rho > 1)
    stop("rho must be between 0.5 and 1.")
  boundary_knots <- c(min(x), internal_knots, max(x))

  # 2. Generate B-spline Basis
  # This creates the matrix of basis functions B_{j,k}(X)
  bs_basis <- splines2::bSpline(x, knots = internal_knots, degree = k, intercept = FALSE,
                      Boundary.knots = range(boundary_knots))

  # Add the constant term phi0 = 1
  # phi(X) = (1, B_{1,k}, ..., B_{d,k})
  phi_full <- cbind(1, bs_basis)
  d_dim <- ncol(phi_full)

  # 3. Initialize Storage
  T_ls <- numeric(n)

  # We use a matrix to store the history of imbalances.
  # Rows correspond to intervals (1 to J).
  # Columns correspond to the basis functions (1 to d_dim).
  Lambda_history <- matrix(0, nrow = J, ncol = d_dim)

  for (i in 1:n) {
    # A. Identify the interval j where X_i falls
    # returns index 1 to J
    interval_idx <- findInterval(x[i], boundary_knots, rightmost.closed = TRUE, all.inside = TRUE)

    # B. Retrieve Lambda^j (The imbalance specific to this interval)
    Lambda_j <- Lambda_history[interval_idx, ]

    # C. Get phi^j(X)
    phi_i <- phi_full[i, ]

    # D. Calculate Inner Product
    test_val <- sum(Lambda_j * phi_i)

    # E. Assignment Decision
    if (i == 1) {
      T_ls[i] <- rbinom(1, 1, 0.5)
    } else {
      if (test_val < 0) {
        T_ls[i] <- rbinom(1, 1, rho)      # Favor Treatment 1
      } else if (test_val > 0) {
        T_ls[i] <- rbinom(1, 1, 1 - rho)  # Favor Treatment 0
      } else {
        T_ls[i] <- rbinom(1, 1, 0.5)
      }
    }

    # F. Update ONLY the Lambda for interval j
    Lambda_new <- Lambda_j + (2 * T_ls[i] - 1) * phi_i

    # Save back to history
    Lambda_history[interval_idx, ] <- Lambda_new
  }

  return(T_ls)
}


#' Discretize observations at standard-normal quantiles
#'
#' Convert a continuous covariate into bins for comparator randomization and regression.
#'
#' @param x Numeric vector of observed covariate values.
#' @param b Integer number of bins, at least 2; this is a bin count, not a bootstrap
#'   count.
#' @return Numeric vector of bin codes 1 through b in the same order as x; missing
#'   observations remain missing.
#' @details The breaks are -Inf, qnorm(1/b), ..., qnorm((b-1)/b), Inf. cut() uses its
#'   default right-closed intervals. Although the factor labels are 0 through b - 1,
#'   as.numeric() returns their one-based factor codes.
disc2 <- function(x, b){
  bins = c(-Inf,qnorm(p = seq(1, b - 1, by = 1) / b), Inf)
  a = as.numeric(cut(x, breaks = bins, labels = c(seq(0, b - 1, by = 1))))
  return(a)
}
########
# inference

#' Simulate the proposed two-sample inference procedure
#'
#' Balance standardized age, generate a calibrated outcome, estimate the mean contrast,
#' and calculate adjusted and naive normal tests.
#'
#' @param i Integer trial seed passed to set.seed().
#' @param n Trial sample size; must equal sample_size established by rd_initialize().
#' @param b Integer bootstrap count; delta() uses 500.
#' @param sigma Positive normal error standard deviation in the calibrated simulation.
#' @param estimate_only Logical; TRUE returns the treatment contrast immediately after
#'   its working regression, skipping later inference. Default FALSE retains the original
#'   inferential path.
#' @return Length-3 numeric vector: treatment contrast, adjusted rejection indicator,
#'   naive rejection indicator. With estimate_only = TRUE, a single unname()d numeric
#'   treatment contrast.
#' @details Uses df2 and BLCOV installed in this source environment by rd_initialize().
#'   The outcome uses the manuscript coefficients, u1 = 20.32581 and u2 = 25.196228,
#'   giving true contrast -4.870418. This simulates outcomes on observed baselines;
#'   rd_fit_observed() analyzes observed FinalHAMD. Normal-test indicators use the 0.975
#'   cutoff and, with this nonzero contrast, describe power rather than Type I error. The
#'   legacy inference path resamples outcomes/covariates and refits GAMs. Allocation
#'   retains the fixed internal knot zero and uses each sample range as its spline
#'   boundaries. Table D.1 uses only the pre-bootstrap point estimates and their
#'   empirical SD across independent trials.
test1.ind <- function(i, n, b, sigma, estimate_only = FALSE){
  
  set.seed(i)
  
  x = df2$AGE_STD; z = df2$HAMD24_STD
  treat <- covRcpp_uni_b(df2$AGE_STD)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)

  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y=ifelse(treat==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  
  md.l <- lm(Y ~ 0 + treat.mx)
  hat.t = md.l$coefficients
   # Table D.1 uses empirical SD across trials, so its estimates need no
   # dataset-specific variance fit or bootstrap. The default retains inference.
   if (isTRUE(estimate_only)) return(unname(hat.t[1] - hat.t[2]))
  sig.epsi2 <- sum(md.l$residuals^2) / (n-2)
  
  bootf = rep(0, b)
  vzx = rep(0, b)
  for(j in 1: b){
    resample = sample(1: n, replace = T)
    ynew = Y[resample]; xnew = x[resample]; znew = z[resample]
   
    treatnew = covRcpp_uni_b(xnew)
    treatnew2 = 1 - treatnew
    mdnew = mgcv::gam(ynew ~ 0 + treatnew + treatnew2 + s(xnew) + s(znew))
    xp = predict(mdnew, type = "lpmatrix")
    fz = xp[, 12: 20] %*% matrix(mdnew$coefficients[12 : 20], ncol = 1)
    
    fzmd = mgcv::gam(fz ~ s(xnew))
    efz = as.vector(predict(fzmd))
    
    bootf[j] = sum((2 * treatnew - 1) * efz) / sqrt(n)
    vzx[j] = var((fz - efz))
  }
  sigma.m2 = var(bootf)
  evzx = mean(vzx)
  
  md = mgcv::gam(Y ~ 0 + treat + treat2 + s(x) + s(z))
  xp = predict(md, type = "lpmatrix")
  hats = xp[, 3: 20] %*% matrix(md$coefficients[3 : 20], ncol = 1)
  
   # Fitted prognostic variance is removed from the working residual component.
  vgf = var(hats)
  
  
  sige2 <- sig.epsi2 - vgf + evzx + sigma.m2
 
  s <- (hat.t[1] - hat.t[2]) / sqrt(sige2 * 4 / n)
  s.n <- (hat.t[1] - hat.t[2]) / sqrt(sig.epsi2 * 4 / n)
  
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  return(c((hat.t[1] - hat.t[2]), prob, prob.n))
} 

#' Simulate the proposed linear working-model procedure
#'
#' Balance standardized age and estimate the contrast with a linear regression on
#' treatment, standardized age, and baseline HAMD24.
#'
#' @return Length-3 numeric vector: treatment contrast, adjusted rejection indicator,
#'   naive rejection indicator. With estimate_only = TRUE, a single unname()d numeric
#'   treatment contrast.
test2.ind <- function(i, n, b, sigma, estimate_only = FALSE){
  
  set.seed(i)

  x = df2$AGE_STD; z = df2$HAMD24_STD
  treat <- covRcpp_uni_b(df2$AGE_STD)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)

  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y=ifelse(treat==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  
  md.l <- lm(Y ~ 0 + treat.mx + x + z)
  hat.t = md.l$coefficients
   # Table D.1 uses empirical SD across trials, so its estimates need no
   # dataset-specific variance fit or bootstrap. The default retains inference.
   if (isTRUE(estimate_only)) return(unname(hat.t[1] - hat.t[2]))
  sig.epsi2 <- sum(md.l$residuals^2) / (n-2)
  
  bootf = rep(0, b)
  vzx = rep(0, b)
  
  md = mgcv::gam(Y ~ 0 + treat + treat2 + s(x) + s(z))
  xp = predict(md, type = "lpmatrix")
  hats = as.vector(xp[, 3: 20] %*% matrix(md$coefficients[3 : 20], ncol = 1))
  
  mnew = cbind(x, z); cm = cov(mnew)
  v = matrix(c(mean(x * hats), mean(z * hats)), ncol = 1)
   # Retain the original linear projection correction used by this script.
  e1 = t(v) %*% solve(cm) %*% v
  
   # Fitted prognostic variance is removed from the working residual component.
  vgf = var(hats)
  
  
  for(j in 1: b){
    resample = sample(1: n, replace = T)
    ynew = Y[resample]; xnew = x[resample]; znew = z[resample]

    treatnew = covRcpp_uni_b(xnew)
    treatnew2 = 1 - treatnew
    
    mdnew = mgcv::gam(ynew ~ 0 + treatnew + treatnew2 + s(xnew) + s(znew))
    xp = predict(mdnew, type = "lpmatrix")
    fz = xp[, 12: 20] %*% matrix(mdnew$coefficients[12 : 20], ncol = 1)
    
    md.ln = lm(ynew ~ 0 + treatnew + treatnew2 + xnew + znew)
    beta0 = md.ln$coefficients[4]
    
    fzminus = (fz - znew %*% matrix(beta0, ncol = 1))
    fzmd = mgcv::gam(fzminus ~ s(xnew))
    efz = as.vector(predict(fzmd))
    
    bootf[j] = sum((2 * treatnew - 1) * efz) / sqrt(n)
    vzx[j] = var(fzminus - efz)
  }
  sigma.m2 = var(bootf)
  evzx = mean(vzx)
  
  sige2 <- sig.epsi2 - vgf + e1 + evzx + sigma.m2

  s <- (hat.t[1] - hat.t[2]) / sqrt(sige2 * 4 / n)
  s.n <- (hat.t[1] - hat.t[2]) / sqrt(sig.epsi2 * 4 / n)
  
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  return(c(hat.t[1] - hat.t[2], prob, prob.n))
} 

#' Simulate the proposed partially linear procedure
#'
#' Balance standardized age and estimate the treatment contrast from an outcome GAM with
#' smooth age and HAMD24 effects.
#'
#' @return Length-2 numeric vector: treatment contrast and rejection indicator.
#' @details Uses df2 and BLCOV installed in this source environment by rd_initialize().
#'   The outcome uses the manuscript coefficients, u1 = 20.32581 and u2 = 25.196228,
#'   giving true contrast -4.870418. This simulates outcomes on observed baselines;
#'   rd_fit_observed() analyzes observed FinalHAMD. Normal-test indicators use the 0.975
#'   cutoff and, with this nonzero contrast, describe power rather than Type I error.
test3.ind <- function(i, n, sigma){
  
  set.seed(i)
   
  x = df2$AGE_STD; z = df2$HAMD24_STD
  treat <- covRcpp_uni_b(df2$AGE_STD)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  
  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y=ifelse(treat==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  
  md = mgcv::gam(Y ~ 0 + treat + treat2 + s(x) + s(z))
  xp = predict(md, type = "lpmatrix")
  
  hat.t = md$coefficients[1:2]
  sig.epsi2 <- sum(md$residuals^2) / (n-2)
  
  s <- (hat.t[1] - hat.t[2]) / sqrt(sig.epsi2 * 4 / n)
  
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  return(c((hat.t[1] - hat.t[2]), prob))
} 

#' Simulate the Ma2015 comparison variant 1
#'
#' Use carat::StrBCD() on age bins and apply the original tau2.f variance correction to
#' the fitted treatment contrast.
#'
#' @return Length-2 numeric vector: treatment contrast and rejection indicator.
fftest.bnz1 <- function(i, n, b, sigma){
  set.seed(i)
  # data generation
  x = df2$AGE_STD; z = df2$HAMD24_STD
  
  xd = disc2(x, b); zd = disc2(z, b)
  X.d = data.frame(xd = xd, extra = rep(1, n))
  
  trt <- ifelse(carat::StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)

  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y.f=ifelse(trt==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  

  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / (n-2)
  
  strx = model.matrix( ~ factor(xd) - 1); strz = model.matrix( ~ factor(zd) - 1)
  varx = 0; varz = 0
  for (k in 1: ncol(strx)) {
    varx = varx + var(x[strx[,k] == 1]) * sum(strx[,k] == 1) / n
  }
  for (k in 1: ncol(strz)) {
    varz = varz + var(z[strz[,k] == 1]) * sum(strz[,k] == 1) / n
  }
  
  e2.f = beta2[3]^2 * varx + beta2[4]^2 * varz
  e1.f = beta2[3]^2 * var(x) + beta2[4]^2 * var(z)
  
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  M = matrix(c(trt, 1 - trt), ncol = 2)
  L = matrix(c(1, -1), ncol = 2)
  temp = L %*% solve(t(M) %*% M) %*% t(L)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * temp * tau2.f)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(c((hat.t.f[1] - hat.t.f[2]), prob.f))
}

#' Simulate the Ma2015 comparison variant 2
#'
#' Use carat::StrBCD() on age bins and apply the original tau2.f variance correction to
#' the fitted treatment contrast.
#'
#' @param estimate_only Logical; TRUE returns the treatment contrast immediately after
#'   its working regression, skipping later inference. Default FALSE retains the original
#'   inferential path.
#' @return Length-2 numeric vector: treatment contrast and rejection indicator. With
#'   estimate_only = TRUE, a single unname()d numeric treatment contrast.
fftest.bnz2 <- function(i, n, b, sigma, estimate_only = FALSE){
  
  set.seed(i)
  x = df2$AGE_STD; z = df2$HAMD24_STD
  
  xd = disc2(x, b); zd = disc2(z, b)
  X.d = data.frame(xd = xd, extra = rep(1, n))
  
  trt <- ifelse(carat::StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y.f=ifelse(trt==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx + z)
  hat.t.f = md.l.f$coefficients
   # Table D.1 uses empirical SD across trials, so its estimates need no
   # dataset-specific variance fit or bootstrap. The default retains inference.
   if (isTRUE(estimate_only)) return(unname(hat.t.f[1] - hat.t.f[2]))
  sig.epsi2 <- sum(md.l.f$residuals^2) / (n-2)
  
  strx = model.matrix( ~ factor(xd) - 1)
  varx = 0
  for (k in 1: ncol(strx)) {
    varx = varx + var(x[strx[,k] == 1]) * sum(strx[,k] == 1) / n
  }
  
  e2.f = beta2[3]^2 * varx 
  e1.f = beta2[3]^2 * var(x)
  
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * tau2.f * 4 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(c((hat.t.f[1] - hat.t.f[2]), prob.f))
}

#' Simulate the Ma2022 comparison variant 1
#'
#' Use carat::StrPBR() and the within-/between-stratum variance decomposition.
#'
#' @return Length-2 numeric vector: treatment contrast and rejection indicator.
ftest1.p.c2 <- function(i, n, b, sigma){
  set.seed(i)
  # data generation
  x = df2$AGE_STD; z = df2$HAMD24_STD
  
  dis1 <- disc2(x, b)
  bs =  model.matrix( ~ factor(dis1) - 1)

  X.d = data.frame(xd = dis1, zd = rep(1, n))        
  
  trt <- ifelse(carat::StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y.f=ifelse(trt==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  
  yk1 = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    yk1[i] = mean(Y.f[trt == 1 & bs[,i] == 1])
  }
  
  yk0 = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    yk0[i] = mean(Y.f[trt == 0 & bs[,i] == 1])
  }
  
  yk1.square = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    yk1.square[i] = sum((Y.f[trt == 1 & bs[,i] == 1] - yk1[i])^2) * (sum( bs[,i] == 1) / n) / sum(trt == 1 & bs[,i] == 1)
  }
  
  yk0.square = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    yk0.square[i] = sum((Y.f[trt == 0 & bs[,i] == 1] - yk0[i])^2) * (sum( bs[,i] == 1) / n) / sum(trt == 0 & bs[,i] == 1)
  }
  
  sig1.square = 2 * sum(yk1.square, na.rm = T) + 2 * sum(yk0.square, na.rm = T)
  
  ay = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    ay[i] = (yk1[i] - mean(Y.f[trt == 1]) - (yk0[i] - mean(Y.f[trt == 0])))^2 * (sum( bs[,i] == 1) / n)
  }
  
  sig2.square = sum(ay, na.rm = T)
  md.l.f <- lm(Y.f ~ trt)
  
  hat.t.f = md.l.f$coefficients[2]
  
  
  s.f <- hat.t.f / sqrt((sig1.square + sig2.square) / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(c(hat.t.f, prob.f))
} 

#' Simulate the Ma2022 comparison variant 2
#'
#' Use carat::StrBCD() and regress on treatment, age-stratum indicators, and standardized
#' baseline HAMD24.
#'
#' @param i Integer trial seed passed to set.seed().
#' @param n Trial sample size; must equal sample_size established by rd_initialize().
#' @param b Integer number of standard-normal quantile bins for comparator allocation;
#'   delta() uses 4. This is not a bootstrap count.
#' @param sigma Positive normal error standard deviation in the calibrated simulation.
#' @param estimate_only Logical; TRUE returns the treatment contrast immediately after
#'   its working regression, skipping later inference. Default FALSE retains the original
#'   inferential path.
#' @return Length-2 numeric vector: treatment contrast and rejection indicator. With
#'   estimate_only = TRUE, a single unname()d numeric treatment contrast.
ftest2.p.c2 <- function(i, n, b, sigma, estimate_only = FALSE){
  set.seed(i)
  # data generation
  x = df2$AGE_STD; z = df2$HAMD24_STD
  
  dis1 <- disc2(x, b)
  bs =  model.matrix(~ factor(dis1) - 1)

  xd = dis1
  X.d = data.frame(xd = xd, e = rep(1, n))
  
  trt <- ifelse(carat::StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y.f=ifelse(trt==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  
  md.l.f <- lm(Y.f ~ trt + bs + z)
  hat.t.f = md.l.f$coefficients[2]
   # Table D.1 uses empirical SD across trials, so its estimates need no
   # dataset-specific variance fit or bootstrap. The default retains inference.
   if (isTRUE(estimate_only)) return(unname(hat.t.f))
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(c(hat.t.f, prob.f))
} 

#' Simulate the Ye2018 comparison variant 1
#'
#' Use carat::AdjBCD() on age strata and the original size-weighted within-stratum
#' outcome variance.
#'
#' @return Length-2 numeric vector: treatment contrast and rejection indicator.
testye.bnz1.c2 <- function(i, n, b, sigma){
  set.seed(i)
  # data generation
  x = df2$AGE_STD; z = df2$HAMD24_STD
  
  dis1 <- disc2(x, b)
  bs =  model.matrix(~ factor(dis1) - 1)
  
  xd = dis1
  X.d = data.frame(xd = xd, zd = rep(1, n))
  
  trt = ifelse(carat::AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y.f=ifelse(trt==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(dis1))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(c((beta.hat[1] - beta.hat[2]), prob.f))
} 

#' Simulate the Ye2018 comparison variant 2
#'
#' Use carat::AdjBCD() on age strata and the original size-weighted within-stratum
#' outcome variance.
#'
#' @return Length-2 numeric vector: treatment contrast and rejection indicator. With
#'   estimate_only = TRUE, a single unname()d numeric treatment contrast.
testye.bnz2.c2 <- function(i, n, b, sigma, estimate_only = FALSE){
  set.seed(i)
  # data generation
  x = df2$AGE_STD; z = df2$HAMD24_STD
  
  dis1 <- disc2(x, b)
  bs =  model.matrix(~ factor(dis1) - 1)
  
  xd = dis1
  X.d = data.frame(xd = xd, zd = rep(1, n))
  
  trt = ifelse(carat::AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  # Quadratic outcome model calibrated to the observed CBASP data.
  u0=25.196228
  u1=20.32581
  Y.f=ifelse(trt==0,u0,u1)-0.972179*BLCOV[,1]+ 0.011102*BLCOV[,1]^2+0.391615*BLCOV[,2]+rnorm(n,mean=0,sd=sigma)
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + xd)
  beta.hat <- md.l2$coefficients
   # Table D.1 uses empirical SD across trials, so its estimates need no
   # dataset-specific variance fit or bootstrap. The default retains inference.
   if (isTRUE(estimate_only)) return(unname(beta.hat[1] - beta.hat[2]))
  
  stratas = split(Y.f, list(dis1))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(c((beta.hat[1] - beta.hat[2]), prob.f))
} 

##########




#' Run the legacy real-data-based simulation for one seed
#'
#' Evaluate nine inference variants at error SDs sqrt(3), 2, and sqrt(5), using the
#' initialized observed baseline profiles.
#'
#' @param i Positive integer outer trial seed.
#' @return Length-60 numeric vector in the legacy order: test1 (1:9), test2 (10:18),
#'   test3 (19:24), fftest1 (25:30), fftest2 (31:36), ftest1 (37:42), ftest2 (43:48),
#'   testye1 (49:54), testye2 (55:60).
#' @details Within each method block, noise SDs follow sqrt(3), 2, sqrt(5). test1/test2
#'   store contrast/adjusted/naive per setting; other methods store contrast/rejection.
#'   Requires rd_initialize() first. Comparator allocations use four age strata.
#'   rd_reproduce() obtains table and supplementary estimates directly, without
#'   unnecessary bootstrap fits.
delta <- function(i){
  
  # Each callback takes noise SD h and returns the named trial procedure's contrast/test vector.
  # case 2
  power.bnz1.ind = sapply(c(sqrt(3), 2, sqrt(5)), function(h) test1.ind(i, sample_size, 500, h))
  power.bnz2.ind = sapply(c(sqrt(3), 2, sqrt(5)), function(h) test2.ind(i, sample_size, 500, h))
  power.bnz3.ind = sapply(c(sqrt(3), 2, sqrt(5)), function(h) test3.ind(i, sample_size, h))
  f1.p.c2 = sapply(c(sqrt(3), 2, sqrt(5)), function(h) fftest.bnz1(i, sample_size, 4, h))
  
  f2.p.c2 = sapply(c(sqrt(3), 2, sqrt(5)), function(h) fftest.bnz2(i, sample_size, 4, h))
  
  powerf1.p.c2 = sapply(c(sqrt(3), 2, sqrt(5)), function(h) ftest1.p.c2(i, sample_size, 4, h))
  
  powerf2.p.c2 = sapply(c(sqrt(3), 2, sqrt(5)), function(h) ftest2.p.c2(i, sample_size, 4, h))
  
  powery1.p.c2 = sapply(c(sqrt(3), 2, sqrt(5)), function(h) testye.bnz1.c2(i, sample_size, 4, h))
  
  powery2.p.c2 = sapply(c(sqrt(3), 2, sqrt(5)), function(h) testye.bnz2.c2(i, sample_size, 4, h))

  
  return(c(power.bnz1.ind, power.bnz2.ind, power.bnz3.ind, f1.p.c2, f2.p.c2, powerf1.p.c2, powerf2.p.c2, powery1.p.c2, powery2.p.c2))
  
}



# Executing the script runs the complete example; source() defines its functions.
if (sys.nframe() == 0L) {
  result <- do.call(rd_reproduce, rd_parse_args())
  print(result$table, row.names = FALSE, digits = 6)
}
