# Figure S4: recompute the simulation called in the original yebin.Rmd.
# Source this file, then run ans <- s4_reproduce(); ans$g4.
# All plotted rates are calculated from saved trials. Historical values appear
# only in a separate comparison CSV, never as inputs to the plotted rates.
# To reproduce Figure S4, input:
# source("reproduce_figure_S4_original.R")
# ans <- s4_reproduce(repetitions = 1000)
# g4 <- ans$g4
# g4
# ans$historical_comparison

#' Allocate treatment using the original global spline-balancing rule
#' @param x Numeric n-by-p feature matrix, including all marginal spline columns.
#' @param w0 Overall treatment-count weight before normalization; default 1.
#' @param w1 Feature-imbalance weight before normalization; default 1.
#' @param rho Biased-coin probability; default 0.9.
#' @return Length-n treatment vector coded 0/1.
#' @details Preserve squared-norm comparisons and runif() draws, since a
#'   dot-product shortcut can change numerical ties and seed-specific results.
s4_allocate <- function(x, w0 = 1, w1 = 1, rho = .9) {
  w <- w0 + w1; w0 <- w0 / w; w1 <- w1 / w
  n <- nrow(x); p <- ncol(x)
  T_ls <- numeric(n); Lambda <- rep(0, 1 + p)
  for (i in seq_len(n)) {
    phi <- c(sqrt(w0), sqrt(w1) * x[i, ])
    Test1 <- sum((Lambda + phi)^2); Test0 <- sum((Lambda - phi)^2)
    if (Test1 < Test0) T_ls[i] <- sum(runif(1) < rho) else
      if (Test1 > Test0) T_ls[i] <- sum(runif(1) < (1 - rho)) else
        T_ls[i] <- sum(runif(1) < .5)
    Lambda <- Lambda + (2 * T_ls[i] - 1) * phi
  }
  T_ls
}

#' Discretize one covariate by empirical quantiles
#' @param x Numeric covariate vector.
#' @param b Number of bins per covariate.
#' @return Numeric bin codes 1,...,b, using quantile()'s default type 7.
#' @details Cuts are sample quantiles, not qnorm() values.
s4_discretize <- function(x, b) {
  cuts <- c(-Inf, quantile(x, probs = seq_len(b - 1L) / b), Inf)
  as.numeric(cut(x, breaks = cuts, labels = seq.int(0L, b - 1L)))
}

#' Calculate the original comparator's joint-stratum variance
#' @param Y Outcome vector or n-by-1 matrix.
#' @param xbin First covariate's empirical bin codes.
#' @param zbin Second covariate's empirical bin codes.
#' @return List containing contrast variance and stratum-size diagnostics.
#' @details Strata with <=1 subject contribute zero, as in yebin.Rmd. No
#'   small-stratum correction is added. Preserve the original summation order.
s4_ye_variance <- function(Y, xbin, zbin) {
  strata <- split(Y, list(xbin, zbin)); counts <- lengths(strata); sig2 <- 0
  for (j in seq_along(strata)) {
    if (counts[j] > 1L) sig2 <- sig2 + var(strata[[j]]) * counts[j]
  }
  sig2 <- sig2 / length(Y)
  list(variance = 4 * sig2 / length(Y), empty_strata = sum(counts == 0L),
       singleton_strata = sum(counts == 1L), small_strata = sum(counts < 5L),
       occupied_strata = sum(counts > 0L))
}

#' Compute one test using the original simulation settings
#' @param i Positive replication number, also used directly as the random seed.
#' @param n Sample size, originally 600.
#' @param u True contrast u1-u2, originally 0 or 0.25.
#' @param method Proposed-t-test or Ye2018plm (the original comparator label).
#' @param b Bins per covariate for Ye; NA for the proposed method.
#' @return One row containing the estimate, variance, rejection and diagnostics.
#' @details Reset set.seed(i) separately for every method/bin/effect. Generate
#'   latent normals with MASS::mvrnorm(), clip Z to [-3,3], allocate treatment,
#'   then draw epsilon. Y=1+u*T+exp(X)/2+log(Z+4)+epsilon. Proposed features are
#'   cbind(B_X,B_Z), degree 2 (order 3), no interior knots, intercept=TRUE and
#'   default sample-range boundaries. The two-arm lm supplies the contrast;
#'   the additive GAM is used only for the variance correction.
s4_trial <- function(i, n, u, method, b = NA_integer_) {
  set.seed(i)
  covariates <- MASS::mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x <- covariates[, 1]; z <- covariates[, 2]
  z[z < -3] <- -3; z[z > 3] <- 3
  if (method == "Proposed-t-test") {
    bs <- splines2::bSpline(x, degree = 2, knots = NULL, intercept = TRUE)
    bs <- cbind(bs, splines2::bSpline(z, degree = 2, knots = NULL, intercept = TRUE))
    T <- s4_allocate(bs)
  } else {
    xbin <- s4_discretize(x, b); zbin <- s4_discretize(z, b)
    X.d <- data.frame(dis1 = xbin, dis2 = zbin)
    # The original AdjBCD() default is a=3; make that setting explicit.
    T <- ifelse(carat::AdjBCD(X.d, a = 3)$assignments == "A", 1, 0)
  }
  trt.mx <- matrix(c(T, 1 - T), ncol = 2)
  beta <- matrix(c(1 + u, 1), ncol = 1)
  if (method == "Proposed-t-test") {
    # The original proposed function first forms g, then adds treatment/noise.
    # Keep that arithmetic order as well as the original random draw order.
    g <- exp(x) / 2 + log(z + 4)
    Y <- trt.mx %*% beta + g + rnorm(n)
  } else Y <- trt.mx %*% beta + exp(x) / 2 + log(z + 4) + rnorm(n)
  fit <- lm(Y ~ 0 + trt.mx)
  estimate <- coef(fit)[1] - coef(fit)[2]
  sig.t <- sum(residuals(fit)^2) / (n - 2); sig.g <- NA_real_
  if (method == "Proposed-t-test") {
    gamfit <- mgcv::gam(Y ~ 0 + trt.mx + s(x) + s(z))
    if (!isTRUE(gamfit$converged)) stop("GAM did not converge at seed ", i)
    hhat <- fitted(gamfit) - trt.mx %*% matrix(coef(gamfit)[1:2], ncol = 1)
    sig.g <- var(hhat)
    diag <- list(variance = as.numeric(4 * (sig.t - sig.g) / n),
      empty_strata = NA_integer_, singleton_strata = NA_integer_,
      small_strata = NA_integer_, occupied_strata = NA_integer_)
  } else diag <- s4_ye_variance(Y, xbin, zbin)
  valid <- is.finite(diag$variance) && diag$variance > 0 && is.finite(estimate)
  data.frame(replicate = i, n = n, method = method, bins = b,
    effect_u1_minus_u2 = u, estimate_u1_minus_u2 = unname(estimate),
    treatment_variance = diag$variance,
    reject = if (valid) as.numeric(abs(estimate) / sqrt(diag$variance) >= qnorm(.975)) else NA_real_,
    status = if (valid) "ok" else "invalid_variance", seed = i,
    pooled_residual_variance = as.numeric(sig.t), fitted_covariate_variance = as.numeric(sig.g),
    empty_strata = diag$empty_strata, singleton_strata = diag$singleton_strata,
    small_strata = diag$small_strata, occupied_strata = diag$occupied_strata,
    singleton_policy = "legacy_zero", stringsAsFactors = FALSE)
}

#' Run one Figure S4 replication over all methods, bins and effects
#' @param i Positive replication identifier and seed.
#' @param n Sample size; default 600.
#' @param bins Bins per covariate; plotted defaults 2,4,...,18.
#' @param effects True contrasts u1-u2; default c(0,0.25).
#' @return Data frame containing 20 trial rows at the defaults.
#' @details The old simulation also computed b=20, but its plot selected only
#'   the first nine bin counts. Supply seq(2,20,2) to compute that extra bin.
s4_delta <- function(i, n = 600L, bins = seq(2L, 18L, 2L), effects = c(0, .25)) {
  if (length(i) != 1L || !is.finite(i) || i < 1 || i != floor(i) ||
      i > .Machine$integer.max || length(n) != 1L || n < 30 || n != floor(n) ||
      !length(bins) || any(!is.finite(bins) | bins < 2 | bins != floor(bins)) ||
      anyDuplicated(bins) || !length(effects) || any(!is.finite(effects)) || anyDuplicated(effects))
    stop("Invalid S4 replication, sample size, bins or effects.")
  rows <- list(); k <- 0L
  for (method in c("Proposed-t-test", "Ye2018plm")) {
    for (b in if (method == "Proposed-t-test") NA_integer_ else bins) {
      for (u in effects) {
        k <- k + 1L; rows[[k]] <- s4_trial(i, n, u, method, b)
      }
    }
  }
  do.call(rbind, rows)
}

#' Return the historical values explicitly inserted in the old plotting chunk
#' @return Data frame used only to compare fresh calculations with the old plot.
#' @details These are archived plotted rates, not saved original trial outcomes.
s4_historical_rates <- function() {
  data.frame(method = c(rep("Proposed-t-test", 2), rep("Ye2018plm", 18)),
    bins = c(NA, NA, rep(seq(2, 18, 2), 2)), effect = c(0, .25, rep(c(0, .25), each = 9)),
    historical_rate = c(.067, .837,
      .047, .063, .066, .073, .072, .078, .110, .133, .166,
      .626, .685, .713, .714, .725, .720, .711, .746, .754))
}

#' Summarize and plot completed Figure S4 trials
#' @param trials Combined data frames returned by s4_delta().
#' @param output_dir Destination directory for summaries, figure and ggplot object.
#' @return Invisibly, list with summary, expanded plot_data and g4.
#' @details Proposed estimates are evaluated once per seed/effect and repeated
#'   across bin positions only for drawing horizontal reference lines.
s4_plot <- function(trials, output_dir) {
  if (anyDuplicated(trials[c("replicate", "method", "bins", "effect_u1_minus_u2")]))
    stop("Duplicate S4 trial rows.")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  trials$group_bins <- ifelse(is.na(trials$bins), 0, trials$bins)
  groups <- split(trials, interaction(trials$method, trials$effect_u1_minus_u2,
                                     trials$group_bins, drop = TRUE))
  summary <- do.call(rbind, lapply(groups, function(z) {
    valid <- is.finite(z$reject); total <- sum(valid); rejected <- sum(z$reject[valid])
    ci <- if (total) binom.test(rejected, total)$conf.int else c(NA, NA)
    data.frame(method = z$method[1], bins = z$bins[1], effect = z$effect_u1_minus_u2[1],
      loaded = nrow(z), evaluable = total, invalid = sum(!valid), rejected = rejected,
      rejection_rate = if (total) rejected / total else NA_real_,
      lower = ci[1], upper = ci[2], mean_estimate = mean(z$estimate_u1_minus_u2), mean_variance = mean(z$treatment_variance),
      mean_singleton_strata = if (all(is.na(z$singleton_strata))) NA_real_ else mean(z$singleton_strata),
      mean_empty_strata = if (all(is.na(z$empty_strata))) NA_real_ else mean(z$empty_strata))
  }))
  rownames(summary) <- NULL
  if (any(summary$invalid)) warning("Some S4 tests are invalid; inspect evaluable/invalid counts.")
  bins <- sort(unique(trials$bins[!is.na(trials$bins)]))
  plot_data <- do.call(rbind, lapply(seq_len(nrow(summary)), function(j) {
    z <- summary[j, ]; if (is.na(z$bins)) {z <- z[rep(1, length(bins)), ]; z$bins <- bins}; z
  }))
  plot_data$curve <- paste0(plot_data$method, ", u1-u2 = ", format(plot_data$effect, trim = TRUE))
  lev <- unique(plot_data$curve[order(plot_data$method, plot_data$effect, decreasing = FALSE)])
  lev <- c(lev[grepl("^Proposed", lev)], lev[!grepl("^Proposed", lev)])
  plot_data$curve <- factor(plot_data$curve, levels = lev)
  g4 <- ggplot2::ggplot(plot_data, ggplot2::aes(bins, rejection_rate, colour = curve,
                                              linetype = curve, group = curve)) +
    ggplot2::geom_line(linewidth = .65) +
    ggplot2::geom_point(data = plot_data[plot_data$method == "Ye2018plm", ], size = 1.7) +
    ggplot2::scale_colour_manual(values = rep(c("#5DADE2", "#4ECDC4", "#8C9CA9", "#264653"), length.out = length(lev))) +
    ggplot2::scale_linetype_manual(values = ifelse(grepl("^Proposed", lev), "solid", "dotdash")) +
    ggplot2::scale_x_continuous(breaks = c(5, 10, 15)) +
    ggplot2::scale_y_continuous(breaks = seq(0, 1, .25), limits = c(0, 1)) +
    ggplot2::labs(x = "Bins per covariate", y = "Rejection probability", colour = "Method", linetype = "Method") +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(legend.position = "right", panel.grid.minor = ggplot2::element_line(colour = "#F5F5F5"))
  comparison <- merge(summary, s4_historical_rates(), by = c("method", "bins", "effect"), all.x = TRUE)
  comparison$difference <- comparison$rejection_rate - comparison$historical_rate
  write.csv(comparison, file.path(output_dir, "Figure_S4_historical_comparison.csv"), row.names = FALSE)
  write.csv(summary, file.path(output_dir, "Figure_S4_summary.csv"), row.names = FALSE)
  write.csv(plot_data, file.path(output_dir, "Figure_S4_plot_data.csv"), row.names = FALSE)
  write.csv(trials, file.path(output_dir, "Figure_S4_trial_results.csv"), row.names = FALSE)
  for (ext in c("png", "pdf")) ggplot2::ggsave(file.path(output_dir, paste0("Figure_S4.", ext)),
    g4, width = 8.5, height = 4, dpi = 300, bg = "white")
  saveRDS(g4, file.path(output_dir, "g4.rds"))
  invisible(list(summary = summary, plot_data = plot_data, historical_comparison = comparison, g4 = g4))
}

#' Reproduce Figure S4 directly without the package master script
#' @param repetitions Number of replications; default 1000.
#' @param output_dir Relative destination for seeds, CSVs, figures and g4.rds.
#' @param first_seed First seed/replication identifier; default 1.
#' @param n Sample size; default 600.
#' @param bins Bin counts per covariate; default 2,4,...,18.
#' @param effects True treatment contrasts u1-u2; default c(0,0.25).
#' @return Invisibly, a list with summary, plot_data, historical_comparison, g4.
#' @details Compatible saved seeds are reused. R/package versions, generator
#'   and simulation-function fingerprints are checked before resuming.
s4_reproduce <- function(repetitions = 1000L, output_dir = "output/figure_S4_original",
                         first_seed = 1L, n = 600L, bins = seq(2L, 18L, 2L), effects = c(0, .25)) {
  if (length(repetitions) != 1L || !is.finite(repetitions) || repetitions < 1 || repetitions != floor(repetitions) ||
      length(first_seed) != 1L || !is.finite(first_seed) || first_seed < 1 || first_seed != floor(first_seed))
    stop("repetitions and first_seed must be positive integers.")
  packages <- c("MASS", "mgcv", "splines2", "carat", "ggplot2")
  for (pkg in packages) if (!requireNamespace(pkg, quietly = TRUE)) stop("Install package ", pkg)
  old_rng <- RNGkind(); on.exit(do.call(RNGkind, as.list(old_rng)), add = TRUE)
  RNGkind("Mersenne-Twister", "Inversion", "Rejection")
  seed_dir <- file.path(output_dir, "seeds"); dir.create(seed_dir, recursive = TRUE, showWarnings = FALSE)
  fingerprint <- tempfile(); on.exit(unlink(fingerprint), add = TRUE)
  saveRDS(lapply(list(s4_delta, s4_trial, s4_allocate, s4_discretize, s4_ye_variance), body), fingerprint)
  recipe <- list(setting = "original_yebin_Rmd", n = n, bins = bins, effects_u1_minus_u2 = effects,
    R = as.character(getRversion()), RNG = RNGkind(), implementation_md5 = unname(tools::md5sum(fingerprint)),
    packages = setNames(vapply(packages, function(p) as.character(utils::packageVersion(p)), character(1)), packages))
  settings_file <- file.path(output_dir, "settings.rds")
  if (file.exists(settings_file) && !identical(readRDS(settings_file), recipe))
    stop("Saved S4 settings/software differ; choose a new output_dir instead of mixing runs.")
  saveRDS(recipe, settings_file)
  ids <- seq.int(first_seed, length.out = repetitions)
  for (i in ids) {
    path <- file.path(seed_dir, paste0("seed", i, ".rds"))
    if (!file.exists(path)) saveRDS(s4_delta(i, n, bins, effects), path)
    if (i %% 100 == 0 || i == tail(ids, 1)) message("S4: seed ", i)
  }
  trials <- do.call(rbind, lapply(file.path(seed_dir, paste0("seed", ids, ".rds")), readRDS))
  if (nrow(trials) != repetitions * length(effects) * (1 + length(bins))) stop("Incomplete S4 results.")
  writeLines(capture.output(sessionInfo()), file.path(output_dir, "sessionInfo.txt"))
  s4_plot(trials, output_dir)
}
