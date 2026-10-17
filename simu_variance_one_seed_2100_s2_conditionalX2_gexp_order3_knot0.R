# File guide: preparation -> frozen-fit bootstrap -> independent reference -> seed output.
# Component quantities base, sigma_E, and sigma_m are unscaled; treatment_variance includes 4/n.
# Functions and local closures are documented below with roxygen2-style input/output headers.

# Scenario (ii): X ~ N(0,1) conditional on -2 < X < 2; g(X)=exp(X).
# This is continuous truncation, NOT clipping/winsorization: no endpoint atoms.
# X is not standardized or rescaled after sampling. Z is still clipped to [-3,3],
# f(Z)=4*log(Z+4), and the error variance in Scenario (ii) remains one.
# For EVERY n, use one fixed interior knot at 0 and boundaries at -2 and 2.
# Balancing splines retain ORDER 3 = DEGREE 2 (quadratic) and the existing
# constant-plus-all-but-the-first-spline parameterization (4 columns).
# Both original-data preparation functions and the independent reference
# generator use vc_scenario2_rX(), vc_scenario2_design(), and vc_fun("ii").
# Oracle X moments, support checks and Scenario (ii) metadata match the new law.
# Variance formulas, the X-only empirical bootstrap, nuisance fits, sample-size
# and seed grids, and Scenarios (i)/(iii) are unchanged from the supplied file.
# In particular, Scenario (iii)'s X coordinates retain their ORIGINAL clipping.
# Main output columns are unchanged; new metadata identify the bounds and g.
# Save results in a fresh directory, separate from earlier simulation settings.
# An Rmd that validates the old [-1,1] study name or support string must be
# updated for this version. Do not bypass checks to pool different settings.

# Empirical-bootstrap component AND treatment-effect variance study.
# Sample sizes: 300, 600, 900, 1200, 1500, 1800, 2100. B = 500 by default.
# The t-test and X-and-Z linear-model formulas match the corrected
# simu_updated_empirical_bootstrap.R; only the imbalance component is bootstrapped.
# In each scenario/n/method, treatment_variance estimates Var(tau_hat).
# Across seeds, compare mean(treatment_variance) with var(mc_treatment_estimate).
# One seed provides ONE independent MC treatment estimate, not an MC variance.
# Keep this version's results in a new folder: the 7-size RNG grid differs from v1.
# Source this file, then call delta(i) for ONE replicate across all sample sizes.
# Sourcing defines the workers; delta(i) runs one replicate and the writers save it.
# Workers need mgcv, MASS, splines2 and optionally Rcpp; plots need ggplot2.


library(mgcv)
library(MASS)
library(splines2)


# Draw independent X ~ N(0,1) conditional on -2 < X < 2.
# qnorm(U), U ~ Uniform(Phi(-2), Phi(2)), has this conditional distribution.
# Never clip the result or map rejected observations onto the endpoints.
#' Sample the Scenario (ii) conditional-normal covariate
#'
#' Draw X directly by inverse-CDF sampling from N(0,1) conditional on -2 < X < 2.
#'
#' @param n Positive integer number of observations to draw.
#' @return Numeric length-n vector strictly inside (-2,2); invalid counts or numerical
#'   endpoint draws cause an error.
#' @details Uses the current RNG state and adds no endpoint atoms or standardization.
vc_scenario2_rX <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != floor(n)) {
    stop("n must be a positive integer.")
  }
  p_lower <- stats::pnorm(-2)
  p_upper <- stats::pnorm(2)
  x <- stats::qnorm(stats::runif(n, min = p_lower, max = p_upper))
  if (any(!is.finite(x)) || any(x <= -2 | x >= 2)) {
    stop("Conditional-normal sampling produced a value outside (-2,2); no clipping was applied.")
  }
  x
}

# Scenario (ii): one fixed interior knot at zero for every sample size.
# Order 3 = polynomial degree 2; preserve the existing feature parameterization.
#' Record the fixed Scenario (ii) knot count
#'
#' Validate the trial size and return the number of interior knots used in its balancing
#' design.
#'
#' @param n Positive integer trial sample size.
#' @return Integer 1L for every valid n; invalid counts cause an error.
#' @details The knot is always at zero; n is validated for metadata consistency and does
#'   not change the knot rule.
vc_scenario2_nknots <- function(n) {
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) ||
      n < 1 || n != floor(n)) {
    stop("n must be a positive integer.")
  }
  1L
}

#' Evaluate the fixed quadratic Scenario (ii) CAR design
#'
#' Construct order-3 splines with boundaries [-2,2], one interior knot at zero, and two
#' local allocation cells.
#'
#' @param x Nonempty finite numeric vector, or one-column matrix, of observations in
#'   [-2,2].
#' @return List with phi (n-by-4 intercept-plus-splines matrix), cell (length-n labels
#'   1/2), J (2L), interior_knots (0), boundary_knots (-2,2), degree (2L), and order
#'   (3L).
#' @details This evaluates observed raw X without drawing or transforming it. The full
#'   marginal spline span is parameterized by a constant plus the three splines retained
#'   by intercept = FALSE. Zero belongs to the upper cell.
vc_scenario2_design <- function(x) {
  # X is drawn from the conditional normal law in the original/reference generator.
  # This function only evaluates the basis; it does not transform or resample X.
  if (!is.numeric(x) || length(x) == 0L ||
      (!is.null(dim(x)) && !(is.matrix(x) && ncol(x) == 1L))) {
    stop("Scenario (ii) requires a numeric vector or a one-column numeric matrix.")
  }
  x <- as.numeric(x)
  if (any(!is.finite(x)) || any(x < -2 | x > 2)) {
    stop("Scenario (ii) X must be finite and within [-2, 2].")
  }

  # Fixed boundaries and interior knot; no sample-quantile calculation.
  breaks <- c(-2, 0, 2)
  bs <- splines2::bSpline(
    x,
    knots = 0,
    Boundary.knots = c(-2, 2),
    degree = 2L,
    intercept = FALSE
  )
  # Preserve the existing convention: constant + all but the first spline.
  # Bins are [-2, 0) and [0, 2]; the sampling law has no endpoint atoms.
  cell <- findInterval(x, breaks, rightmost.closed = TRUE, all.inside = TRUE)
  list(
    phi = cbind(1, bs),
    cell = as.integer(cell),
    J = 2L,
    interior_knots = 0,
    boundary_knots = c(-2, 2),
    degree = 2L,
    order = 3L
  )
}

#' Construct fixed balancing features and cell labels
#'
#' Build the design once from observed covariates so original and bootstrap CAR paths use
#' identical feature definitions.
#'
#' @param x Observed X: a vector for dis/ind or an n-by-p matrix for ts.
#' @param type Design identifier: dis for stratum indicators, ind for quantile-knot
#'   marginal splines, or ts for zero-split tensor splines.
#' @param J Number of empirical quantile intervals for type = ind; default 2.
#' @param k B-spline degree for ind and ts; order is k + 1. Default 3; Scenario (iii)
#'   calls this with k = 2.
#' @return List with phi (n-by-d feature matrix), cell (length-n integer labels), and J
#'   (number of imbalance vectors to allocate).
#' @details The dis design uses one global cell. The ind design uses observed quantile
#'   boundaries. The ts design uses boundary knots [-1,1] and an interior knot at zero in
#'   each coordinate, forms the complete tensor, then replaces its first column with an
#'   intercept, preserving the span. The code checks distinct quantile knots and
#'   observations on both sides of zero.
eb_car_design <- function(x, type, J = 2, k = 3) {
  if (type == "dis") {
    strata = interaction(as.data.frame(x), drop = TRUE, sep = ".")
    phi = cbind(1, model.matrix(~ strata - 1))
    # The discrete rule balances a global intercept and stratum indicators.
    return(list(phi = phi, cell = rep(1L, nrow(phi)), J = 1L))
  }
  if (type == "ind") {
    knots = as.numeric(quantile(x, seq(0, 1, length.out = J + 1)))
    if (any(diff(knots) <= 0)) stop("CAR quantile knots must be distinct.")
    bs = bSpline(x, knots = knots[-c(1, length(knots))], degree = k,
                 intercept = FALSE, Boundary.knots = range(knots))
    cell = findInterval(x, knots, rightmost.closed = TRUE, all.inside = TRUE)
    return(list(phi = cbind(1, bs), cell = cell, J = J))
  }
  x = as.matrix(x)
  p = ncol(x)
  # Knot callback: d is a coordinate index; returns ordered observed endpoints plus zero.
  knots = lapply(seq_len(p), function(d) unique(sort(c(range(x[, d]), 0))))
  # Basis callback: d is a coordinate index; returns its full n-row marginal spline
  # matrix.
  bs_list = lapply(seq_len(p), function(d) {
    rng = range(x[, d])
    if (rng[1] >= 0 || rng[2] <= 0)
      stop("The multivariate CAR basis requires observations on both sides of zero.")
    bSpline(x[, d], knots = 0, Boundary.knots = c(-1, 1), degree = k, intercept = TRUE)
  })
  if (p == 1) {
    phi = bs_list[[1]]
  } else {
    terms = paste0("bs_list[[", seq_len(p), "]]")
    phi = model.matrix(as.formula(paste("~", paste(terms, collapse = ":"), "- 1")))
  }
  # Cell callbacks: i selects a row and returns its character cell key;
  # the inner callback takes coordinate d and returns its integer interval label.
  ids = vapply(seq_len(nrow(x)), function(i) {
    paste(vapply(seq_len(p), function(d)
      findInterval(x[i, d], knots[[d]], rightmost.closed = TRUE,
                   all.inside = TRUE), integer(1)), collapse = "-")
  }, character(1))
  cell = match(ids, unique(ids))
  # Form the COMPLETE tensor first, then replace ONE tensor column by 1.
  # Because the full tensor columns sum to one, this preserves its full span.
  # Dropping a marginal column BEFORE taking products would lose other terms.
  # degree 2 (= order 3) gives 1 + 15 = 16 features for two covariates.
  list(phi = cbind(1, phi[, -1, drop = FALSE]), cell = cell, J = max(cell))
}

#' Restart CAR on an original or resampled row sequence
#'
#' Start all cell imbalances at zero and sequentially allocate treatments using the
#' supplied fixed design.
#'
#' @param design List returned by eb_car_design(), containing phi, cell, and J.
#' @param index Integer vector of donor row indices in arrival order; duplicates are
#'   allowed. Defaults to all original rows.
#' @param rho Probability of choosing the treatment that reduces local imbalance; default
#'   0.9.
#' @return Numeric 0/1 treatment vector of length length(index), in the supplied sequence
#'   order.
#' @details Only rows are selected: knots and balancing features are never re-estimated.
#'   Ties use probability 0.5. Each call consumes random draws from the current RNG
#'   state.
eb_assign_car <- function(design, index = seq_len(nrow(design$phi)), rho = 0.9) {
  phi = design$phi[index, , drop = FALSE]
  cell = design$cell[index]
  # One imbalance row per local cell; each allocation call starts a fresh path.
  Lambda = matrix(0, design$J, ncol(phi))
  treat = numeric(length(index))
  for (i in seq_along(index)) {
    j = cell[i]
    test_val = sum(Lambda[j, ] * phi[i, ])
    prob = if (test_val < 0) rho else if (test_val > 0) 1 - rho else 0.5
    treat[i] = rbinom(1, 1, prob)
    Lambda[j, ] = Lambda[j, ] + (2 * treat[i] - 1) * phi[i, ]
  }
  treat
}

#' Fit the conditional mean on the mixed support of clipped X
#'
#' Estimate E(response | X) separately on the interior, edges, and corner atoms of the
#' Scenario (iii) square.
#'
#' @param X Numeric n-by-2 matrix with coordinates clipped to [-1,1].
#' @param response Numeric length-n response whose conditional mean is to be estimated.
#' @param k_multi Upper limit on each marginal tensor-smooth basis dimension in the
#'   interior; default 6.
#' @return List with fitted.values (length-n conditional-mean estimates in original row
#'   order) and fits (a named list of fitted GAM objects or scalar group means).
#' @details Corner atoms and groups with fewer than 12 rows use their observed mean.
#'   Interior groups use a two-dimensional smooth; edges use a univariate smooth. GAM
#'   fits use REML and must converge.
eb_fit_mixed_mean <- function(X, response, k_multi = 6L) {
  X = as.matrix(X)
  status = ifelse(X == -1, -1L, ifelse(X == 1, 1L, 0L))
  key = paste(status[, 1], status[, 2], sep = ":")
  efz = numeric(nrow(X))
  fits = list()
  for (group in unique(key)) {
    rows = which(key == group)
    free = which(status[rows[1], ] == 0L)
    nd = data.frame(response = response[rows], x1 = X[rows, 1], x2 = X[rows, 2])
    # Sparse components and corner atoms use their observed group mean.
    if (length(free) == 0L || length(rows) < 12L) {
      fits[[group]] = mean(response[rows])
      efz[rows] = fits[[group]]
    } else {
      if (length(free) == 2L) {
        kk = max(3L, min(k_multi, floor(sqrt(length(rows) / 3))))
        fit = gam(response ~ te(x1, x2, k = c(kk, kk), bs = c("cr", "cr")),
                  data = nd, method = "REML")
      } else {
        nd$t = X[rows, free]
        kk = min(8L, length(unique(nd$t)), max(3L, floor(length(rows) / 4)))
        fit = gam(response ~ s(t, k = kk, bs = "cr"), data = nd, method = "REML")
      }
      if (!isTRUE(fit$converged)) stop("Conditional-mean GAM did not converge.")
      fits[[group]] = fit
      efz[rows] = as.vector(fitted(fit))
    }
  }
  list(fitted.values = efz, fits = fits)
}

#' Prepare one two-sample outer trial (dis)
#'
#' Generate one trial, fit the outcome and conditional-mean models once, and retain the
#' quantities needed for the frozen-fit empirical bootstrap.
#'
#' @param i Integer original-trial seed passed to set.seed().
#' @param n Positive integer trial sample size.
#' @param u Scalar true contrast between treatment means 1+u and 1.
#' @param fun Vectorized g(X)+f(Z) function accepting observed x,z.
#' @return List with X (observed n-by-p covariates), z (length-n vector), design (fixed
#'   CAR features/cells), m (length-n fitted conditional means), sigma_E2 (conditional
#'   residual variance), base (working residual variance minus fitted Vgf), estimate
#'   (treatment contrast), Y (outcomes), treat (0/1 assignments), and beta_Z (NA for this
#'   treatment-only working model).
#' @details Scenario (i) uses discretized correlated X and clipped Z. No bootstrap is run
#'   here. The conditional response is fitted f(Z). Nuisance means are sample-centered;
#'   the population oracle is centered separately.
vc_prepare1.dis <- function(i, n, u, fun){
  
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0.7, 0.7, 1), ncol = 2))
  x = data[,1]
  z = data[,2]
  x[x < -2] = -2
  x[x > 2] = 2
  x[(x > -2) & (x <0)] = -1
  x[x == 0] = 0
  x[(x > 0) & (x < 2)] = 1
  z[z < -3] = -3
  z[z > 3] = 3
  
  
  g = fun(x, z)
  
  design = eb_car_design(x, "dis")
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  epsi <- rnorm(n)
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx)
  hat.t = coef(md.l)
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + x + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # Remove fitted treatment means to isolate the prognostic g(X)+f(Z) signal.
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  vgf = var(hats)
  fzmd = gam(fz ~ factor(x))
  if (!isTRUE(fzmd$converged)) stop("Conditional-mean GAM did not converge.")
  efz = as.vector(fitted(fzmd))

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fz - efz)

  # Keep fitted means exactly as in the empirical-bootstrap simulation file.
  # Smooth terms are sample-centered; oracle moments below are population-centered.
  return(list(X = as.matrix(x), z = z, design = design, m = efz,
              sigma_E2 = evzx, base = sig.epsi2 - vgf,
              estimate = unname(hat.t[1] - hat.t[2]), Y = as.vector(Y),
              treat = treat, beta_Z = NA_real_))
}

#' Prepare one linear working-model outer trial (dis)
#'
#' Generate one trial, fit the outcome and conditional-mean models once, and retain the
#' quantities needed for the frozen-fit empirical bootstrap.
#'
#' @return List with X (observed n-by-p covariates), z (length-n vector), design (fixed
#'   CAR features/cells), m (length-n fitted conditional means), sigma_E2 (conditional
#'   residual variance), base (working residual variance minus fitted Vgf plus fitted Q),
#'   estimate (treatment contrast), Y (outcomes), treat (0/1 assignments), and beta_Z
#'   (fitted working Z slope).
vc_prepare2.dis <- function(i, n, u, fun){
  
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0.7, 0.7, 1), ncol = 2))
  x = data[,1]
  z = data[,2]
  x[x < -2] = -2
  x[x > 2] = 2
  x[(x > -2) & (x <0)] = -1
  x[x == 0] = 0
  x[(x > 0) & (x < 2)] = 1
  z[z < -3] = -3
  z[z > 3] = 3
  
  
  g = fun(x, z)
  
  design = eb_car_design(x, "dis")
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  epsi <- rnorm(n)
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx + x + z)
  if (anyNA(coef(md.l)) || df.residual(md.l) <= 0)
    stop("Proposed-lm working model is rank deficient or has no residual degrees of freedom.")
  hat.t = coef(md.l)[c("treat.mx1", "treat.mx2")]
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + x + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # hhat estimates g(X) + f(Z); use U = (X, Z) in the TABLE's correction.
  # Remove fitted treatment means to isolate the prognostic g(X)+f(Z) signal.
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  vgf = var(hats)
  mnew = cbind(x, z)
  cm = cov(mnew)
  v = cov(mnew, hats)
  # Q = v^T Cov(U)^(-1) v adds back the fitted linear projection for raw U=(X,Z).
  e1 = as.numeric(crossprod(v, solve(cm, v)))
  # beta_Z is the Z coefficient after adjusting for ALL X coordinates.
  beta.work = coef(md.l)[colnames(mnew)]
  if (any(!is.finite(beta.work))) stop("Nonfinite Proposed-lm covariate coefficients.")
  beta0 = unname(beta.work["z"])
  fzminus = fz - beta0 * z
  fzmd = gam(fzminus ~ factor(x))
  if (!isTRUE(fzmd$converged)) stop("Conditional-mean GAM did not converge.")
  efz = as.vector(fitted(fzmd))

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fzminus - efz)

  # Keep fitted means exactly as in the empirical-bootstrap simulation file.
  # Smooth terms are sample-centered; oracle moments below are population-centered.
  return(list(X = as.matrix(x), z = z, design = design, m = efz,
              sigma_E2 = evzx, base = sig.epsi2 - vgf + e1,
              estimate = unname(hat.t[1] - hat.t[2]), Y = as.vector(Y),
              treat = treat, beta_Z = beta0))
}

#' Prepare one two-sample outer trial (ind)
#'
#' Generate one trial, fit the outcome and conditional-mean models once, and retain the
#' quantities needed for the frozen-fit empirical bootstrap.
#'
#' @return List with X (observed n-by-p covariates), z (length-n vector), design (fixed
#'   CAR features/cells), m (length-n fitted conditional means), sigma_E2 (conditional
#'   residual variance), base (working residual variance minus fitted Vgf), estimate
#'   (treatment contrast), Y (outcomes), treat (0/1 assignments), and beta_Z (NA for this
#'   treatment-only working model).
#' @details Scenario (ii) samples conditional-normal X independently of clipped Z and
#'   uses the fixed order-3 zero-knot design. No bootstrap is run here. The conditional
#'   response is fitted f(Z). Nuisance means are sample-centered; the population oracle
#'   is centered separately.
vc_prepare1.ind <- function(i, n, u, fun){
  
  set.seed(i)
  # Scenario (ii): X is conditional-normal, independent of the latent Z*.
  # Draw n observations directly; no endpoint clipping or random sample-size loss.
  x = vc_scenario2_rX(n)
  z = stats::rnorm(n)
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  # The same fixed knot at zero is used for every trial size.
  design = vc_scenario2_design(x)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  epsi <- rnorm(n)
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx)
  hat.t = coef(md.l)
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + s(x) + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # Remove fitted treatment means to isolate the prognostic g(X)+f(Z) signal.
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  vgf = var(hats)
  fzmd = gam(fz ~ s(x))
  if (!isTRUE(fzmd$converged)) stop("Conditional-mean GAM did not converge.")
  efz = as.vector(fitted(fzmd))

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fz - efz)

  # Keep fitted means exactly as in the empirical-bootstrap simulation file.
  # Smooth terms are sample-centered; oracle moments below are population-centered.
  return(list(X = as.matrix(x), z = z, design = design, m = efz,
              sigma_E2 = evzx, base = sig.epsi2 - vgf,
              estimate = unname(hat.t[1] - hat.t[2]), Y = as.vector(Y),
              treat = treat, beta_Z = NA_real_))
}

#' Prepare one linear working-model outer trial (ind)
#'
#' Generate one trial, fit the outcome and conditional-mean models once, and retain the
#' quantities needed for the frozen-fit empirical bootstrap.
#'
#' @return List with X (observed n-by-p covariates), z (length-n vector), design (fixed
#'   CAR features/cells), m (length-n fitted conditional means), sigma_E2 (conditional
#'   residual variance), base (working residual variance minus fitted Vgf plus fitted Q),
#'   estimate (treatment contrast), Y (outcomes), treat (0/1 assignments), and beta_Z
#'   (fitted working Z slope).
#' @details Scenario (ii) samples conditional-normal X independently of clipped Z and
#'   uses the fixed order-3 zero-knot design. No bootstrap is run here. The conditional
#'   response is fitted f(Z) minus beta_Z*Z, and Q uses all raw X coordinates and Z.
#'   Nuisance means are sample-centered; the population oracle is centered separately.
vc_prepare2.ind <- function(i, n, u, fun){
  
  set.seed(i)
  # Scenario (ii): X is conditional-normal, independent of the latent Z*.
  # Draw n observations directly; no endpoint clipping or random sample-size loss.
  x = vc_scenario2_rX(n)
  z = stats::rnorm(n)
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  # The same fixed knot at zero is used for every trial size.
  design = vc_scenario2_design(x)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  epsi <- rnorm(n)
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx + x + z)
  if (anyNA(coef(md.l)) || df.residual(md.l) <= 0)
    stop("Proposed-lm working model is rank deficient or has no residual degrees of freedom.")
  hat.t = coef(md.l)[c("treat.mx1", "treat.mx2")]
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + s(x) + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # hhat estimates g(X) + f(Z); use U = (X, Z) in the TABLE's correction.
  # Remove fitted treatment means to isolate the prognostic g(X)+f(Z) signal.
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  vgf = var(hats)
  mnew = cbind(x, z)
  cm = cov(mnew)
  v = cov(mnew, hats)
  # Q = v^T Cov(U)^(-1) v adds back the fitted linear projection for raw U=(X,Z).
  e1 = as.numeric(crossprod(v, solve(cm, v)))
  # beta_Z is the Z coefficient after adjusting for ALL X coordinates.
  beta.work = coef(md.l)[colnames(mnew)]
  if (any(!is.finite(beta.work))) stop("Nonfinite Proposed-lm covariate coefficients.")
  beta0 = unname(beta.work["z"])
  fzminus = fz - beta0 * z
  fzmd = gam(fzminus ~ s(x))
  if (!isTRUE(fzmd$converged)) stop("Conditional-mean GAM did not converge.")
  efz = as.vector(fitted(fzmd))

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fzminus - efz)

  # Keep fitted means exactly as in the empirical-bootstrap simulation file.
  # Smooth terms are sample-centered; oracle moments below are population-centered.
  return(list(X = as.matrix(x), z = z, design = design, m = efz,
              sigma_E2 = evzx, base = sig.epsi2 - vgf + e1,
              estimate = unname(hat.t[1] - hat.t[2]), Y = as.vector(Y),
              treat = treat, beta_Z = beta0))
}

#' Prepare one two-sample outer trial (ts)
#'
#' Generate one trial, fit the outcome and conditional-mean models once, and retain the
#' quantities needed for the frozen-fit empirical bootstrap.
#'
#' @param i Integer original-trial seed passed to set.seed().
#' @param n Positive integer trial sample size.
#' @param u Scalar true contrast between treatment means 1+u and 1.
#' @param fun Vectorized g(X)+f(Z) function accepting observed x1,x2,z.
#' @return List with X (observed n-by-p covariates), z (length-n vector), design (fixed
#'   CAR features/cells), m (length-n fitted conditional means), sigma_E2 (conditional
#'   residual variance), base (working residual variance minus fitted Vgf), estimate
#'   (treatment contrast), Y (outcomes), treat (0/1 assignments), and beta_Z (NA for this
#'   treatment-only working model).
#' @details Scenario (iii) clips both X coordinates to [-1,1], uses the complete order-3
#'   tensor span for allocation, and has heteroscedastic normal noise. No bootstrap is
#'   run here. The conditional response is fitted f(Z). Nuisance means are sample-
#'   centered; the population oracle is centered separately.
vc_prepare1.ts <- function(i, n, u, fun){
  
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  x = cbind(x1, x2)
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  xs1 = bSpline(x1, degree = 2, knots = NULL, intercept = T)
  xs2 = bSpline(x2, degree = 2, knots = NULL, intercept = T)
  
  # Tensor callback: i is a row index; returns the full product of the two marginal spline
  # rows.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  # Remove the redundant constant direction in the OUTCOME regression.
  xs = xs[, -1, drop = FALSE]
  g = fun(x1, x2, z)
  
  design = eb_car_design(x, "ts", k = 2)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  # Noise callback: i is a row index; returns one zero-mean draw with that row's
  # conditional variance.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx)
  hat.t = coef(md.l)
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + xs + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # Remove fitted treatment means to isolate the prognostic g(X)+f(Z) signal.
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  vgf = var(hats)
  fzmd = eb_fit_mixed_mean(x, fz)
  efz = fzmd$fitted.values

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fz - efz)

  # Keep fitted means exactly as in the empirical-bootstrap simulation file.
  # Smooth terms are sample-centered; oracle moments below are population-centered.
  return(list(X = as.matrix(x), z = z, design = design, m = efz,
              sigma_E2 = evzx, base = sig.epsi2 - vgf,
              estimate = unname(hat.t[1] - hat.t[2]), Y = as.vector(Y),
              treat = treat, beta_Z = NA_real_))
}

#' Prepare one linear working-model outer trial (ts)
#'
#' Generate one trial, fit the outcome and conditional-mean models once, and retain the
#' quantities needed for the frozen-fit empirical bootstrap.
#'
#' @return List with X (observed n-by-p covariates), z (length-n vector), design (fixed
#'   CAR features/cells), m (length-n fitted conditional means), sigma_E2 (conditional
#'   residual variance), base (working residual variance minus fitted Vgf plus fitted Q),
#'   estimate (treatment contrast), Y (outcomes), treat (0/1 assignments), and beta_Z
#'   (fitted working Z slope).
#' @details Scenario (iii) clips both X coordinates to [-1,1], uses the complete order-3
#'   tensor span for allocation, and has heteroscedastic normal noise. No bootstrap is
#'   run here. The conditional response is fitted f(Z) minus beta_Z*Z, and Q uses all raw
#'   X coordinates and Z. Nuisance means are sample-centered; the population oracle is
#'   centered separately.
vc_prepare2.ts <- function(i, n, u, fun){
  
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  x = cbind(x1, x2)
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  xs1 = bSpline(x1, degree = 2, knots = NULL, intercept = T)
  xs2 = bSpline(x2, degree = 2, knots = NULL, intercept = T)
  
  # Tensor callback: i is a row index; returns the full product of the two marginal spline
  # rows.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  # Remove the redundant constant direction in the OUTCOME regression.
  xs = xs[, -1, drop = FALSE]
  
  g = fun(x1, x2, z)
  
  design = eb_car_design(x, "ts", k = 2)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  # Noise callback: i is a row index; returns one zero-mean draw with that row's
  # conditional variance.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx + x1 + x2 + z)
  if (anyNA(coef(md.l)) || df.residual(md.l) <= 0)
    stop("Proposed-lm working model is rank deficient or has no residual degrees of freedom.")
  hat.t = coef(md.l)[c("treat.mx1", "treat.mx2")]
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + xs + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # hhat estimates g(X) + f(Z); use U = (X, Z) in the TABLE's correction.
  # Remove fitted treatment means to isolate the prognostic g(X)+f(Z) signal.
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  vgf = var(hats)
  mnew = cbind(x1, x2, z)
  cm = cov(mnew)
  v = cov(mnew, hats)
  # Q = v^T Cov(U)^(-1) v adds back the fitted linear projection for raw U=(X,Z).
  e1 = as.numeric(crossprod(v, solve(cm, v)))
  # beta_Z is the Z coefficient after adjusting for ALL X coordinates.
  beta.work = coef(md.l)[colnames(mnew)]
  if (any(!is.finite(beta.work))) stop("Nonfinite Proposed-lm covariate coefficients.")
  beta0 = unname(beta.work["z"])
  fzminus = fz - beta0 * z
  fzmd = eb_fit_mixed_mean(x, fzminus)
  efz = fzmd$fitted.values

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fzminus - efz)

  # Keep fitted means exactly as in the empirical-bootstrap simulation file.
  # Smooth terms are sample-centered; oracle moments below are population-centered.
  return(list(X = as.matrix(x), z = z, design = design, m = efz,
              sigma_E2 = evzx, base = sig.epsi2 - vgf + e1,
              estimate = unname(hat.t[1] - hat.t[2]), Y = as.vector(Y),
              treat = treat, beta_Z = beta0))
}
.vc_engine = new.env(parent = globalenv())
#' Initialize and cache the optional C++ simulation engine
#'
#' Compile CAR allocation and joint-row bootstrap routines into .vc_engine, or retain the
#' R fallback when Rcpp is unavailable.
#'
#' @return Logical TRUE when the cached CAR routine exists or compilation succeeds; FALSE
#'   when Rcpp is absent or compilation fails.
#' @details Mutates .vc_engine. vc_car_cpp takes phi, one-based cell labels, J, donor
#'   index, and rho, returning a 0/1 assignment vector. vc_boot_cpp takes phi, cell, J, a
#'   reward matrix, B, and rho, returning a B-by-reward-column matrix of scaled
#'   imbalances. Compilation errors emit a message; the embedded C++ source is preserved.
vc_enable_fast = function() {
  if (exists('vc_car_cpp', .vc_engine, inherits=FALSE)) return(TRUE)
  if (!requireNamespace('Rcpp',quietly=TRUE)) return(FALSE)
  tryCatch({
    Rcpp::cppFunction(env=.vc_engine, code='
    Rcpp::NumericVector vc_car_cpp(Rcpp::NumericMatrix phi, Rcpp::IntegerVector cell,
      int J, Rcpp::IntegerVector index, double rho) {
      int n=index.size(), p=phi.ncol();
      std::vector<double> L(J*p,0.0); Rcpp::NumericVector T(n);
      for(int i=0;i<n;i++) {
        int k=index[i]-1, j=cell[k]-1; long double dot=0;
        for(int q=0;q<p;q++) dot += (double)(L[j*p+q]*phi(k,q));
        double pr=dot<0?rho:(dot>0?1-rho:0.5);
        T[i]=R::rbinom(1,pr); double a=2*T[i]-1;
        for(int q=0;q<p;q++) L[j*p+q]+=a*phi(k,q);
      }
      return T;
    }')
    Rcpp::cppFunction(env=.vc_engine, code='
    Rcpp::NumericMatrix vc_boot_cpp(Rcpp::NumericMatrix phi, Rcpp::IntegerVector cell,
      int J, Rcpp::NumericMatrix reward, int B, double rho) {
      int n=phi.nrow(), p=phi.ncol(), M=reward.ncol();
      Rcpp::NumericMatrix out(B,M); double rootn=std::sqrt((double)n);
      for(int b=0;b<B;b++) {
        std::vector<double> L(J*p,0.0), totals(M,0.0);
        for(int i=0;i<n;i++) {
          int k=(int)std::floor(R::runif(0,1)*n), j=cell[k]-1;
          long double dot=0;
          for(int q=0;q<p;q++) dot+=(double)(L[j*p+q]*phi(k,q));
          double pr=dot<0?rho:(dot>0?1-rho:0.5), a=2*R::rbinom(1,pr)-1;
          for(int q=0;q<p;q++) L[j*p+q]+=a*phi(k,q);
          for(int q=0;q<M;q++) totals[q]+=a*reward(k,q);
        }
        for(int q=0;q<M;q++) out(b,q)=totals[q]/rootn;
      }
      return out;
    }')
    TRUE
  },error=function(e) {message('Rcpp unavailable; using slower R engine: ',conditionMessage(e)); FALSE})
}
# Retain the earlier pure-R allocator before installing the dispatch wrapper.
vc_assign_R = eb_assign_car
#' Dispatch CAR allocation to the selected engine
#'
#' Use the cached C++ allocator when present, otherwise the original R implementation
#' saved as vc_assign_R.
#'
#' @param design List returned by eb_car_design(), containing phi, cell, and J.
#' @param index Integer vector of donor row indices in arrival order; duplicates are
#'   allowed. Defaults to all original rows.
#' @param rho Probability of choosing the treatment that reduces local imbalance; default
#'   0.9.
#' @return Numeric 0/1 treatment vector of length length(index), in the supplied sequence
#'   order.
#' @details Preserves fixed donor-row features and resets all local imbalances on each
#'   call. vc_assign_R is the earlier R eb_assign_car definition; this later definition
#'   provides backend dispatch.
eb_assign_car = function(design,index=seq_len(nrow(design$phi)),rho=.9) {
  if(exists('vc_car_cpp',.vc_engine,inherits=FALSE))
    .vc_engine$vc_car_cpp(design$phi,as.integer(design$cell),as.integer(design$J),as.integer(index),rho)
  else vc_assign_R(design,index,rho)
}
#' Bootstrap imbalance statistics with fixed features and rewards
#'
#' Resample observed rows jointly, restart local CAR from zero on each path, and compute
#' signed reward totals divided by sqrt(n).
#'
#' @param design List containing fixed phi (feature matrix), cell (integer row labels),
#'   and J (number of imbalance vectors), returned by a CAR design helper.
#' @param reward Finite numeric n-by-M matrix of original fitted or oracle rewards; rows
#'   correspond to design$phi.
#' @param B Integer number of bootstrap paths; delta() requires at least 2.
#' @param rho Biased-coin probability favoring reduced cell imbalance; default 0.9.
#' @return Numeric B-by-M matrix; column j contains n^(-1/2) times the signed imbalance
#'   for reward column j, one value per bootstrap path.
#' @details Uses cached vc_boot_cpp when available, otherwise the R engine. Reward fits
#'   and feature definitions remain fixed, and all reward columns share the same donor
#'   sequence and allocation within a path.
vc_boot = function(design,reward,B,rho=.9) {
  stopifnot(nrow(reward)==nrow(design$phi),all(is.finite(reward)))
  if(exists('vc_boot_cpp',.vc_engine,inherits=FALSE))
    return(.vc_engine$vc_boot_cpp(design$phi,as.integer(design$cell),as.integer(design$J),reward,as.integer(B),rho))
  n=nrow(reward); ans=matrix(0,B,ncol(reward))
  for(b in seq_len(B)) {
    index=sample.int(n,n,replace=TRUE)
    a=2*eb_assign_car(design,index,rho)-1
    ans[b,]=as.vector(crossprod(a,reward[index,,drop=FALSE]))/sqrt(n)
  }
  ans
}
#' Select the scenario-specific prognostic mean function
#'
#' Return the vectorized g(X)+f(Z) function used by both original and reference trials.
#'
#' @param scenario Scenario identifier 'i', 'ii', or 'iii'.
#' @return A function taking observed x,z for i/ii or x1,x2,z for iii and returning
#'   numeric per-row prognostic means; unsupported identifiers return NULL.
#' @details Defaults are 2*x+exp(z)/2.5, exp(x)+4*log(z+4), and 2*x1*x2^2+exp(z)/2.5 for
#'   i, ii, and iii, respectively.
vc_fun = function(scenario) switch(scenario,
  #' Evaluate the Scenario (i) prognostic mean
  #'
  #' Supply the vectorized outcome mean selected by vc_fun().
  #'
  #' @param x Observed X vector after scenario preprocessing.
  #' @param z Observed clipped Z vector.
  #' @return Numeric vector 2*x+exp(z)/2.5, one value per participant.
  i=function(x,z) 2*x+exp(z)/2.5,
  #' Evaluate the Scenario (ii) prognostic mean
  #'
  #' Supply the vectorized outcome mean selected by vc_fun().
  #'
  #' @param x Observed X vector after scenario preprocessing.
  #' @param z Observed clipped Z vector.
  #' @return Numeric vector exp(x)+4*log(z+4), one value per participant.
  ii=function(x,z) exp(x)+4*log(z+4),
  #' Evaluate the Scenario (iii) prognostic mean
  #'
  #' Supply the vectorized outcome mean selected by vc_fun().
  #'
  #' @param x1 Observed first X-coordinate vector.
  #' @param x2 Observed second X-coordinate vector.
  #' @param z Observed clipped Z vector.
  #' @return Numeric vector 2*x1*x2^2+exp(z)/2.5, one value per participant.
  iii=function(x1,x2,z) 2*x1*x2^2+exp(z)/2.5)
#' Translate scenario identifiers to preparation-function suffixes
#'
#' Provide the suffix used to look up vc_prepare1 and vc_prepare2 for a selected
#' scenario.
#'
#' @param scenario Scenario identifier 'i', 'ii', or 'iii'.
#' @return Character 'dis', 'ind', or 'ts' for i, ii, or iii; unsupported identifiers
#'   return NULL.
vc_suffix = function(scenario) switch(scenario,i='dis',ii='ind',iii='ts')
#' Generate one independent reference trial
#'
#' Draw observed covariates, assign treatments with the scenario design, and generate
#' outcomes with the requested effect.
#'
#' @param seed Integer RNG seed for this independent trial; the function calls
#'   set.seed(seed).
#' @param n Positive integer trial size.
#' @param scenario Scenario identifier 'i', 'ii', or 'iii'.
#' @param effect Scalar true treatment contrast; default zero.
#' @return List with X (n-by-1 or n-by-2 observed covariates), z (length-n clipped Z), T
#'   (length-n 0/1 assignments), Y (length-n outcomes), and design (fixed CAR features
#'   and cells).
#' @details Scenario (i) discretizes X; (ii) conditions normal X on (-2,2), independently
#'   of Z; (iii) clips both correlated X coordinates to [-1,1]. Z is clipped to [-3,3].
#'   Noise variance is one for i/ii and log(1+2*X1^2+1.5*Z^2) for iii. This function fits
#'   no nuisance models.
vc_dataset = function(seed,n,scenario,effect=0) {
  set.seed(seed)
  if(scenario=='iii') {
    latent=MASS::mvrnorm(n,c(0,0,0),matrix(c(1,.5,.7,.5,1,0,.7,0,1),3))
    X=latent[,1:2]; X[]=pmin(1,pmax(-1,X)); z=pmin(3,pmax(-3,latent[,3]))
    design=eb_car_design(X,'ts',k=2)
    g=vc_fun(scenario)(X[,1],X[,2],z)
  } else {
    if(scenario=='ii') {
      # Same conditional X distribution and Z distribution as in both preparation functions.
      x=vc_scenario2_rX(n)
      z=pmin(3,pmax(-3,stats::rnorm(n)))
    } else {
      # Scenario (i) is unchanged.
      r=if(scenario=='i') .7 else 0
      latent=MASS::mvrnorm(n,c(0,0),matrix(c(1,r,r,1),2))
      x=latent[,1]; z=pmin(3,pmax(-3,latent[,2]))
      if(scenario=='i') x=ifelse(x<= -2,-2,ifelse(x<0,-1,ifelse(x<2,1,2)))
    }
    X=matrix(x,ncol=1)
    # Use the same fixed-zero-knot rule for the independent reference.
    if (scenario == "ii") {
      design = vc_scenario2_design(x)
    } else {
      design = eb_car_design(x, vc_suffix(scenario))
    }
    g=vc_fun(scenario)(x,z)
  }
  T=eb_assign_car(design)
  sd=if(scenario=='iii') sqrt(log(1+2*X[,1]^2+1.5*z^2)) else rep(1,n)
  Y=1+effect*T+g+rnorm(n,sd=sd)
  list(X=X,z=z,T=T,Y=Y,design=design)
}
#' Compute independent reference contrasts and oracle imbalances
#'
#' Analyze one reference trial with the two-sample and raw-covariate linear working
#' models, and evaluate the population conditional-mean imbalance.
#'
#' @param seed Integer RNG seed for this independent trial; the function calls
#'   set.seed(seed).
#' @param n Positive integer reference trial size.
#' @param scenario Scenario identifier 'i', 'ii', or 'iii'.
#' @param oracle Population-reference list returned by vc_make_oracle(), especially its
#'   two-column m(X) function.
#' @param effect Scalar true treatment contrast; default zero.
#' @return Two-row data frame with method, estimate (treatment contrast), and imbalance
#'   (scaled oracle m_a(X) imbalance), in two-sample then linear-model order.
#' @details Reference observations are independent of the outer trial when a separate
#'   seed is supplied. A sampling variance is obtained only after collecting these rows
#'   across repetitions; rank-deficient reference regression causes an error.
vc_reference_one = function(seed,n,scenario,oracle,effect=0) {
  d=vc_dataset(seed,n,scenario,effect)
  U=cbind(d$X,d$z); A=2*d$T-1
  fit1=lm.fit(cbind(d$T,1-d$T),d$Y)
  fit2=lm.fit(cbind(d$T,1-d$T,U),d$Y)
  if(anyNA(fit2$coefficients)) stop('Rank-deficient Monte Carlo working regression.')
  mstat=as.vector(crossprod(A,oracle$m(d$X)))/sqrt(n)
  data.frame(method=c('Proposed-t-test','Proposed-lm'),
    estimate=c(diff(fit1$coefficients[2:1]),diff(fit2$coefficients[2:1])),
    imbalance=mstat)
}
#' Estimate feasible variance components for one paired outer trial
#'
#' Prepare both working models on the same trial and bootstrap their fixed fitted
#' conditional means together with the population-reference means.
#'
#' @param seed Integer RNG seed for this independent trial; the function calls
#'   set.seed(seed).
#' @param boot_seed Separate integer seed for joint donor-row bootstrap paths.
#' @param n Positive integer outer trial size.
#' @param scenario Scenario identifier 'i', 'ii', or 'iii'.
#' @param oracle Population-reference list returned by vc_make_oracle(), especially its
#'   two-column m(X) function.
#' @param B Integer bootstrap count, at least 2 in delta().
#' @param effect Scalar true treatment contrast; default zero.
#' @return Two-row data frame: method, estimate (contrast), sigma_m (fitted-mean
#'   bootstrap variance), sigma_m_oracle (oracle-mean bootstrap variance), sigma_E
#'   (original conditional residual variance), base (remaining variance before sigma_E
#'   and sigma_m), variance (4/n times their sum), and bootstrap_mcse (SE of sigma_m from
#'   bootstrap draws).
#' @details Checks that outcomes and assignments match across the two preparations. Only
#'   sigma_m is bootstrapped; nuisance fits and sigma_E stay fixed. Oracle bootstrap
#'   components are diagnostics and do not enter the feasible variance.
vc_outer_one = function(seed,boot_seed,n,scenario,oracle,B,effect=0) {
  suf=vc_suffix(scenario); fun=vc_fun(scenario)
  # Preparation callback: a is method index 1/2; returns the corresponding original-data
  # preparation list.
  p=lapply(1:2,function(a) get(paste0('vc_prepare',a,'.',suf))(seed,n,effect,fun))
  stopifnot(identical(p[[1]]$Y,p[[2]]$Y),identical(p[[1]]$treat,p[[2]]$treat))
  rewards=cbind(p[[1]]$m,p[[2]]$m,oracle$m(p[[1]]$X))
  set.seed(boot_seed)
  draws=vc_boot(p[[1]]$design,rewards,B)
  vm=apply(draws,2,var)
  data.frame(method=c('Proposed-t-test','Proposed-lm'),
    # Extraction callbacks take one preparation list x and return its scalar
    # contrast, sigma_E2, base, or base+sigma_E2 in the matching named column.
    estimate=vapply(p,function(x)x$estimate,numeric(1)),
    sigma_m=vm[1:2], sigma_m_oracle=vm[3:4],
    sigma_E=vapply(p,function(x)x$sigma_E2,numeric(1)),
    base=vapply(p,function(x)x$base,numeric(1)),
    variance=4/n*(vapply(p,function(x)x$base+x$sigma_E2,numeric(1))+vm[1:2]),
    # MCSE callback: j selects a fitted-reward column; returns the SE of its bootstrap
    # sample variance.
    bootstrap_mcse=vapply(1:2,function(j)vc_variance_se(draws[,j]),numeric(1)))
}
#' Approximate the Monte Carlo standard error of a sample variance
#'
#' Use the empirical fourth central moment and the finite-sample variance-of-variance
#' formula.
#'
#' @param x Numeric vector of at least two independent finite draws of the target
#'   statistic.
#' @return Scalar plug-in standard error of var(x), with the estimated squared SE bounded
#'   below by zero.
#' @details This is the SE of a sample variance, not the standard deviation of x and not
#'   the SE of its mean. No missing-value filtering is performed here.
vc_variance_se = function(x) {
  n=length(x); v=var(x)
  sqrt(max(0,(mean((x-mean(x))^4)-(n-3)/(n-1)*v^2)/n))
}
# Population oracle used only for the Monte Carlo reference / diagnostic.
# It is never used to fit the feasible bootstrap estimator.
# Z = clip(Z*, -3, 3); Scenario (iii) Xj are still clipped to [-1,1].
# Scenario (ii) X instead has the continuous conditional N(0,1) law on [-2,2].
# f is identified by population centering: E[f(Z)] is subtracted.

#' Construct Gauss-Legendre quadrature nodes on the unit interval
#'
#' Use the symmetric Jacobi-matrix eigenvalues and eigenvectors to produce an integration
#' rule on [0,1].
#'
#' @param k Integer quadrature order, at least 8; default 64.
#' @return List with u (length-k nodes in ascending order) and w (length-k integration
#'   weights).
vc_gauss01 <- function(k = 64L) {
  stopifnot(length(k) == 1L, is.finite(k), k >= 8L, k == as.integer(k))
  j = seq_len(k - 1L)
  off = j / sqrt(4 * j^2 - 1)
  mat = matrix(0, k, k)
  mat[cbind(j, j + 1L)] = off
  eig = eigen(mat + t(mat), symmetric = TRUE)
  ord = order(eig$values)
  list(u = (eig$values[ord] + 1) / 2, w = eig$vectors[1, ord]^2)
}

#' Clip observations to specified bounds
#'
#' Apply lower and upper caps while preserving values inside the interval.
#'
#' @param x Numeric observations, possibly a vector or matrix.
#' @param lo Scalar lower clipping bound.
#' @param hi Scalar upper clipping bound.
#' @return Numeric values pmin(hi, pmax(lo, x)); with scalar bounds, matrix inputs are
#'   returned as a vector.
vc_clip <- function(x, lo, hi) pmin(hi, pmax(lo, x))

# Weighted nodes over a (possibly truncated) normal distribution.
# Tail intervals are integrated over finite latent bounds +/- 12 standard
# deviations. Their omitted probability is < 4e-33; integrating on the latent
# scale avoids the singular derivatives of inverse-normal tail quadrature.
#' Create weighted nodes for a normal interval integral
#'
#' Map the unit-interval quadrature rule to a latent normal interval and multiply weights
#' by its density.
#'
#' @param mu Scalar mean of the latent normal distribution.
#' @param sd Positive scalar normal standard deviation.
#' @param lo Lower integration endpoint; may be -Inf.
#' @param hi Upper integration endpoint; may be Inf.
#' @param gl List with u and w from vc_gauss01().
#' @return List with x (integration nodes) and w (unnormalized normal-probability
#'   weights); an empty interval returns x = mu and w = 0.
#' @details Infinite tails are approximated by mu +/- 12*sd. Weights integrate interval
#'   mass and are not normalized to sum to one.
vc_normal_nodes <- function(mu, sd, lo, hi, gl) {
  lower = max(lo, mu - 12 * sd)
  upper = min(hi, mu + 12 * sd)
  if (upper <= lower) return(list(x = mu, w = 0))
  x = lower + (upper - lower) * gl$u
  list(x = x, w = (upper - lower) * gl$w * dnorm(x, mu, sd))
}

# Exact conditional moments of clipped normal Z and f(Z)=exp(Z)/2.5.
#' Compute moments of clipped normal Z and its exponential mean term
#'
#' Include both tail atoms and the interior of Z clipped to [-3,3], with f(Z)=exp(Z)/2.5.
#'
#' @param mu Numeric latent normal mean(s), scalar or vector.
#' @param sd Positive latent normal standard deviation(s), scalar or conformable vector.
#' @return Numeric matrix with columns ef = E[f(Z)], ez = E[Z], ez2 = E[Z^2], and ezf =
#'   E[Z*f(Z)], one row per evaluated normal law.
#' @details These are uncentered moments. Scenario (ii), whose f differs, integrates its
#'   log mean term separately in vc_make_oracle().
vc_normal_clip_moments <- function(mu, sd) {
  lo = (-3 - mu) / sd
  hi = (3 - mu) / sd
  pl = pnorm(lo)
  ph = pnorm(hi, lower.tail = FALSE)
  pint = pnorm(hi) - pnorm(lo)
  ez = -3 * pl + 3 * ph + mu * pint + sd * (dnorm(lo) - dnorm(hi))
  ez2 = 9 * (pl + ph) + (mu^2 + sd^2) * pint +
    2 * mu * sd * (dnorm(lo) - dnorm(hi)) +
    sd^2 * (lo * dnorm(lo) - hi * dnorm(hi))
  ls = lo - sd
  hs = hi - sd
  ee = exp(mu + sd^2 / 2) * (pnorm(hs) - pnorm(ls))
  ef = (exp(-3) * pl + exp(3) * ph + ee) / 2.5
  ezf = (-3 * exp(-3) * pl + 3 * exp(3) * ph +
    (mu + sd^2) * ee +
    sd * exp(mu + sd^2 / 2) * (dnorm(ls) - dnorm(hs))) / 2.5
  cbind(ef = ef, ez = ez, ez2 = ez2, ezf = ezf)
}

#' Build population conditional means for Monte Carlo diagnostics
#'
#' Compute the population linear projection and conditional moments for the observed,
#' discretized or clipped support, using quadrature and boundary interpolation as needed.
#'
#' @param scenario Scenario 'i', 'ii', or 'iii'; numeric/character 1, 2, 3 are also
#'   mapped to these identifiers.
#' @param quadrature Integer Gauss-Legendre order, at least 8; default 64.
#' @param edge_grid Integer interpolation-grid size per clipped-X edge, at least 3;
#'   default 1001. delta() requires at least 11.
#' @return List with m (function returning n-by-2 centered oracle means a1/a2), beta
#'   (population Z projection slope), mean_f (E[f(Z)]), coefficients (all U slopes),
#'   cov_u (Cov(U)), cov_uh (Cov(U,g+f)), conditional_moments (function returning ef/ez),
#'   accuracy (integration diagnostics), scenario, quadrature, and edge_grid.
#' @details m1(X)=E[f(Z)|X]-E[f(Z)] and m2(X)=m1(X)-beta*E[Z|X]. Conditional independence
#'   makes both exactly zero in Scenario (ii). Other scenarios integrate the latent
#'   events represented by observed strata, edges, and corners. These references do not
#'   fit or replace the feasible nuisance models.
vc_make_oracle <- function(scenario, quadrature = 64L, edge_grid = 1001L) {
  scenario = as.character(scenario)
  scenario = switch(scenario, "1" = "i", "2" = "ii", "3" = "iii", scenario)
  if (!scenario %in% c("i", "ii", "iii")) stop("Unknown scenario.")
  if (length(edge_grid) != 1L || !is.finite(edge_grid) || edge_grid < 3L ||
      edge_grid != as.integer(edge_grid)) stop("edge_grid must be an integer >= 3.")
  gl = vc_gauss01(quadrature)
  marginal = vc_normal_clip_moments(0, 1)
  accuracy = list(quadrature = quadrature, edge_grid = edge_grid,
                  latent_tail_sd = 12, quadrature_method = "Gauss-Legendre on latent normal intervals")

  if (scenario == "ii") {
    z = c(-3, -3 + 6 * gl$u, 3)
    w = c(pnorm(-3), 6 * gl$w * dnorm(-3 + 6 * gl$u), pnorm(-3))
    f = 4 * log(z + 4)
    mean_f = sum(w * f)
    # X ~ N(0,1) | (-2 < X < 2), independent of the latent Z*.
    # Integrate against phi(x)/(Phi(2)-Phi(-2)) on (-2,2); NO endpoint atoms.
    x_mass = stats::pnorm(2) - stats::pnorm(-2)
    x_nodes = -2 + 4 * gl$u
    x_weights = 4 * gl$w * stats::dnorm(x_nodes) / x_mass
    ex2 = sum(x_weights * x_nodes^2)  # E[X]=0 exactly by symmetry.
    exg = sum(x_weights * x_nodes * exp(x_nodes))
    cov_u = diag(c(ex2, sum(w * z^2)))
    dimnames(cov_u) = list(c("x", "z"), c("x", "z"))
    cov_uh = c(x = exg, z = sum(w * z * f))
    coefficients = solve(cov_u, cov_uh)
    #' Evaluate population conditional moments in Scenario (ii)
    #'
    #' Use independence to return constant E[f(Z)] and zero E[Z|X].
    #'
    #' @param X Observed numeric matrix with one conditional-normal X column within
    #'   [-2,2].
    #' @return Numeric n-by-2 matrix with ef = E[f(Z)|X] and ez = E[Z|X].
    #' @details Local closure captures the population calculations from vc_make_oracle();
    #'   it is returned as conditional_moments.
    pred = function(X) {
      X = as.matrix(X)
      if (ncol(X) != 1L || any(!is.finite(X)) || any(abs(X) > 2))
        stop("Scenario (ii) oracle requires one conditional-normal X column in [-2,2].")
      # Conditioning only X preserves independence from Z: centered m_a(X) are zero.
      cbind(ef = rep(mean_f, nrow(X)), ez = rep(0, nrow(X)))
    }
    accuracy$probability_error = sum(w) - 1
    accuracy$x_probability_error = sum(x_weights) - 1
    accuracy$x_mean_error = sum(x_weights * x_nodes)
    accuracy$x_variance = ex2
    # Analytic variance of a symmetric conditional standard normal on [-2,2].
    accuracy$x_variance_error = ex2 - (1 - 4 * stats::dnorm(2) / x_mass)
    accuracy$x_normalizing_probability = x_mass
    accuracy$x_boundary_mass = 0
    accuracy$x_g_covariance = exg
    accuracy$x_transform = "conditional_normal[-2,2]"
    accuracy$g_function = "exp(X)"
    accuracy$centered_m_exactly_zero = TRUE
  } else if (scenario == "i") {
    levels = c(-2, -1, 1, 2)
    cuts = c(-Inf, -2, 0, 2, Inf)
    prob = diff(pnorm(cuts))
    # Stratum callback: j selects one discrete latent interval; returns normalized
    # conditional moment integrals.
    rows = lapply(seq_len(4L), function(j) {
      nd = vc_normal_nodes(0, 1, cuts[j], cuts[j + 1L], gl)
      colSums(vc_normal_clip_moments(0.7 * nd$x, sqrt(0.51)) * nd$w) / prob[j]
    })
    tab = do.call(rbind, rows)
    #' Evaluate population conditional moments in Scenario (i)
    #'
    #' Look up moments integrated over the latent normal event for each observed X
    #' stratum.
    #'
    #' @param X Observed numeric matrix with one discrete X column in {-2,-1,1,2}.
    #' @return Numeric n-by-2 matrix with ef = E[f(Z)|X] and ez = E[Z|X].
    #' @details Local closure captures the population calculations from vc_make_oracle();
    #'   it is returned as conditional_moments.
    pred = function(X) {
      idx = match(as.matrix(X)[, 1], levels)
      if (anyNA(idx)) stop("Scenario 1 oracle requires X in {-2,-1,1,2}.")
      tab[idx, c("ef", "ez"), drop = FALSE]
    }
    ex2 = sum(prob * levels^2)
    xz = sum(prob * levels * tab[, "ez"])
    xf = sum(prob * levels * tab[, "ef"])
    cov_u = matrix(c(ex2, xz, xz, marginal[1, "ez2"]), 2L,
                   dimnames = list(c("x", "z"), c("x", "z")))
    cov_uh = c(x = 2 * ex2 + xf, z = 2 * xz + marginal[1, "ezf"])
    coefficients = solve(cov_u, cov_uh)
    mean_f = unname(marginal[1, "ef"])
    accuracy$mean_f_error = sum(prob * tab[, "ef"]) - mean_f
    accuracy$mean_z_error = sum(prob * tab[, "ez"])
    accuracy$conditional_z_symmetry_error = max(abs(tab[, "ez"] + rev(tab[, "ez"])))
  } else {
    # Latent covariance: cor(X1*,X2*)=.5, cor(X1*,Z*)=.7,
    # cor(X2*,Z*)=0. Boundary X values represent tail events, not X*=+/-1.
    b = solve(matrix(c(1, 0.5, 0.5, 1), 2L), c(0.7, 0))
    sd_z = sqrt(1 - sum(b * c(0.7, 0)))
    #' Evaluate conditional normal moments at latent X coordinates
    #'
    #' Apply the population conditional law of Z given both latent X coordinates.
    #'
    #' @param x1 First latent coordinate(s), scalar or vector.
    #' @param x2 Conformable second latent coordinate(s).
    #' @return Numeric matrix with ef and ez columns from vc_normal_clip_moments().
    #' @details Captures the conditional regression coefficients b and residual sd_z;
    #'   also used as an integrand over clipped-X tail events.
    interior = function(x1, x2) {
      vc_normal_clip_moments(b[1] * x1 + b[2] * x2, sd_z)[, c("ef", "ez"), drop = FALSE]
    }
    #' Integrate population moments on one clipped-X edge
    #'
    #' Average conditional moments over the latent tail represented by the clipped
    #' coordinate.
    #'
    #' @param q Clipped coordinate index, 1 or 2.
    #' @param sign Boundary sign, -1 or 1.
    #' @param t Observed value of the other, free coordinate.
    #' @return Length-2 numeric vector of ef and ez conditional moments for that edge
    #'   point.
    #' @details Captures gl and interior(); interval probability weights are normalized
    #'   within the tail event.
    edge = function(q, sign, t) {
      lo = if (sign < 0) -Inf else 1
      hi = if (sign < 0) -1 else Inf
      nd = vc_normal_nodes(0.5 * t, sqrt(0.75), lo, hi, gl)
      vals = if (q == 1L) interior(nd$x, t) else interior(t, nd$x)
      colSums(vals * nd$w) / sum(nd$w)
    }
    grid = seq(-1, 1, length.out = edge_grid)
    edges = list()
    for (q in 1:2) for (sg in c(-1, 1)) {
      # Edge callback: t is a free-coordinate grid point; returns its two edge moments
      # before interpolation.
      val = t(vapply(grid, function(t) edge(q, sg, t), numeric(2)))
      edges[[paste(q, sg)]] = list(ef = approxfun(grid, val[, 1], rule = 2),
                                   ez = approxfun(grid, val[, 2], rule = 2))
    }

    # Integrate all nine latent rectangles. Moments use CLIPPED X in U and g.
    #' Integrate joint moments over a latent X rectangle
    #'
    #' Compute probability and clipped-observation moments needed for population
    #' projections and corner conditioning.
    #'
    #' @param lo1 First-coordinate lower latent bound.
    #' @param hi1 First-coordinate upper latent bound.
    #' @param lo2 Second-coordinate lower latent bound.
    #' @param hi2 Second-coordinate upper latent bound.
    #' @return Named numeric vector of unnormalized rectangle integrals: prob, ef, ez, X
    #'   means and second/cross moments, x1z/x2z, x1h/x2h for h=g+f, and gz.
    #' @details Bounds may be infinite. Uses captured gl and interior() with the
    #'   correlated normal X law; the integrands use X clipped to [-1,1].
    rectangle = function(lo1, hi1, lo2, hi2) {
      nx = vc_normal_nodes(0, 1, lo1, hi1, gl)
      ans = c(prob = 0, ef = 0, ez = 0, x1 = 0, x2 = 0,
              x1sq = 0, x2sq = 0, x1x2 = 0, x1z = 0, x2z = 0,
              x1h = 0, x2h = 0, gz = 0)
      for (j in seq_along(nx$x)) {
        ny = vc_normal_nodes(0.5 * nx$x[j], sqrt(0.75), lo2, hi2, gl)
        mm = interior(nx$x[j], ny$x)
        w = nx$w[j] * ny$w
        x1 = vc_clip(nx$x[j], -1, 1)
        x2 = vc_clip(ny$x, -1, 1)
        g = 2 * x1 * x2^2
        ans = ans + c(sum(w), sum(w * mm[, "ef"]), sum(w * mm[, "ez"]),
          sum(w * x1), sum(w * x2), sum(w * x1^2), sum(w * x2^2),
          sum(w * x1 * x2), sum(w * x1 * mm[, "ez"]),
          sum(w * x2 * mm[, "ez"]), sum(w * x1 * (g + mm[, "ef"])),
          sum(w * x2 * (g + mm[, "ef"])), sum(w * g * mm[, "ez"]))
      }
      ans
    }
    cuts = c(-Inf, -1, 1, Inf)
    corners = list()
    all_rect = list()
    for (j in 1:3) for (k in 1:3) {
      rr = rectangle(cuts[j], cuts[j + 1L], cuts[k], cuts[k + 1L])
      all_rect[[paste(j, k)]] = rr
      if (j != 2L && k != 2L)
        corners[[paste(if (j == 1L) -1 else 1, if (k == 1L) -1 else 1)]] =
          rr[c("ef", "ez")] / rr["prob"]
    }
    mm = Reduce(`+`, all_rect)
    cov_u = matrix(c(mm["x1sq"], mm["x1x2"], mm["x1z"],
                     mm["x1x2"], mm["x2sq"], mm["x2z"],
                     mm["x1z"], mm["x2z"], marginal[1, "ez2"]), 3L,
                   dimnames = list(c("x1", "x2", "z"), c("x1", "x2", "z")))
    cov_uh = setNames(c(mm["x1h"], mm["x2h"], mm["gz"] + marginal[1, "ezf"]),
                      c("x1", "x2", "z"))
    coefficients = solve(cov_u, cov_uh)
    mean_f = unname(marginal[1, "ef"])
    accuracy$probability_error = unname(mm["prob"] - 1)
    accuracy$mean_f_error = unname(mm["ef"] - mean_f)
    accuracy$zero_mean_errors = mm[c("x1", "x2", "ez")]
    accuracy$independent_x2_z_cov_error = unname(mm["x2z"])
    # Symmetry callback: key identifies a corner; returns the absolute sum of its Z mean
    # and the opposite corner's mean.
    accuracy$corner_z_symmetry_error = max(vapply(names(corners), function(key) {
      signs = as.numeric(strsplit(key, " ", fixed = TRUE)[[1]])
      abs(corners[[key]]["ez"] + corners[[paste(-signs[1], -signs[2])]]["ez"])
    }, numeric(1)))
    #' Evaluate population conditional moments in Scenario (iii)
    #'
    #' Combine interior conditional-normal moments with interpolated edge and integrated
    #' corner moments.
    #'
    #' @param X Observed numeric matrix with two clipped X columns within [-1,1].
    #' @return Numeric n-by-2 matrix with ef = E[f(Z)|X] and ez = E[Z|X].
    #' @details Local closure captures the population calculations from vc_make_oracle();
    #'   it is returned as conditional_moments.
    pred = function(X) {
      X = as.matrix(X)
      if (ncol(X) != 2L || any(!is.finite(X)) || any(abs(X) > 1))
        stop("Scenario 3 oracle requires two clipped X columns in [-1,1].")
      x1 = X[, 1]
      x2 = X[, 2]
      out = interior(x1, x2)
      for (q in 1:2) for (sg in c(-1, 1)) {
        idx = which(X[, q] == sg & abs(X[, 3L - q]) < 1)
        e = edges[[paste(q, sg)]]
        out[idx, ] = cbind(e$ef(X[idx, 3L - q]), e$ez(X[idx, 3L - q]))
      }
      for (s1 in c(-1, 1)) for (s2 in c(-1, 1)) {
        idx = which(x1 == s1 & x2 == s2)
        out[idx, ] = matrix(rep(corners[[paste(s1, s2)]], each = length(idx)), ncol = 2L)
      }
      out
    }
  }
  beta = unname(coefficients["z"])
  #' Evaluate the two centered population bootstrap rewards
  #'
  #' Convert population conditional moments into the two working-model imbalance means.
  #'
  #' @param X Observed n-by-p covariate matrix on the selected scenario support.
  #' @return Numeric n-by-2 matrix: a1 = E[f(Z)|X]-E[f(Z)] and a2 = a1-beta*E[Z|X].
  #' @details Captures pred(), mean_f, and the population Z projection slope beta from
  #'   vc_make_oracle().
  m = function(X) {
    mm = pred(X)
    m1 = mm[, "ef"] - mean_f
    cbind(a1 = m1, a2 = m1 - beta * mm[, "ez"])
  }
  list(m = m, beta = beta, mean_f = mean_f, coefficients = coefficients,
       cov_u = cov_u, cov_uh = cov_uh, conditional_moments = pred,
       accuracy = accuracy, scenario = scenario, quadrature = quadrature, edge_grid = edge_grid)
}
#' Compare outer variance estimates with independent Monte Carlo references
#'
#' Summarize fitted bootstrap, oracle bootstrap, and total adjusted variances separately
#' for each working model.
#'
#' @param outer Outer-trial data frame for one scenario/n cell, with method, sigma_m,
#'   sigma_m_oracle, and variance.
#' @param reference Independent reference data frame for the same cell, with method,
#'   estimate, and imbalance.
#' @param scenario Scenario identifier attached to summary rows.
#' @param n Trial size attached to summary rows.
#' @return Data frame with means, mean SEs, reference sample variances and their SEs,
#'   differences and difference SEs, ratios, relative bias, RMSE, repetition counts, and
#'   counts of nonpositive/negative estimates for each method and target.
#' @details Total variance uses reference treatment contrasts; imbalance components use
#'   reference oracle imbalances. Ratios are NA for references <= 1e-14. No failures or
#'   missing values are removed here; use validated inputs.
vc_summarize = function(outer,reference,scenario,n) {
  ans=list(); k=0
  for(method in unique(outer$method)) {
    d=outer[outer$method==method,]; r=reference[reference$method==method,]
    for(target in c('Bootstrap component','Oracle bootstrap component','Total adjusted variance')) {
      est=switch(target,'Bootstrap component'=d$sigma_m,
        'Oracle bootstrap component'=d$sigma_m_oracle,'Total adjusted variance'=d$variance)
      values=if(target=='Total adjusted variance') r$estimate else r$imbalance
      truth=var(values); mse=vc_variance_se(values)
      k=k+1
      ans[[k]]=data.frame(scenario=scenario,n=n,method=method,target=target,
        outer_repetitions=nrow(d),mc_repetitions=nrow(r),mean_estimate=mean(est),
        mean_estimate_se=sd(est)/sqrt(length(est)),mc_variance=truth,mc_variance_se=mse,
        difference=mean(est)-truth, difference_se=sqrt(var(est)/length(est)+mse^2),
        variance_ratio=if(truth>1e-14) mean(est)/truth else NA_real_,
        relative_bias=if(truth>1e-14) (mean(est)-truth)/truth else NA_real_,
        rmse_vs_mc=sqrt(mean((est-truth)^2)),
        nonpositive_estimates=sum(est<=0),
        negative_estimates=sum(est<0))
    }
  }
  do.call(rbind,ans)
}
#' Save component and treatment-variance level plots
#'
#' Compare estimated mean variances with Monte Carlo references, with a panel for each
#' working-model/scenario combination.
#'
#' @param summary Data frame returned by vc_summarize(), possibly combined over trial
#'   sizes/scenarios.
#' @param output_dir Character destination directory, created recursively if necessary.
#' @return The supplied summary invisibly; writes PNG/PDF figures for every target.
#' @details Total contrast variances are multiplied by n/4 for display only. Error bars
#'   use approximate 95% Monte Carlo intervals. These level plots are separate from the
#'   companion Rmd variance-ratio and Type I error figure. Requires ggplot2.
vc_plot = function(summary,output_dir) {
  if(!requireNamespace('ggplot2',quietly=TRUE)) stop('Install ggplot2 to plot results.')
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  for(target in unique(summary$target)) {
    d=summary[summary$target==target,]
    # n/4 scaling makes total variance visible across sample sizes; tables stay raw.
    scale=if(target=='Total adjusted variance') d$n/4 else rep(1,nrow(d))
    a=data.frame(d[,c('scenario','n','method')],series='Mean estimated variance',
      value=d$mean_estimate*scale,se=d$mean_estimate_se*scale)
    b=data.frame(d[,c('scenario','n','method')],series='Monte Carlo variance',
      value=d$mc_variance*scale,se=d$mc_variance_se*scale)
    z=rbind(a,b); z$scenario=factor(z$scenario,levels=c('i','ii','iii'),
      labels=c('Scenario (i)','Scenario (ii)','Scenario (iii)'))
    z$method=factor(z$method,levels=c('Proposed-t-test','Proposed-lm'))
    z$lower=pmax(0,z$value-1.96*z$se); z$upper=z$value+1.96*z$se
    # Show negative full-variance means honestly if one occurs.
    z$lower[z$value<0]=z$value[z$value<0]-1.96*z$se[z$value<0]
    p=ggplot2::ggplot(z,ggplot2::aes(n,value,colour=series,linetype=series,shape=series))+
      ggplot2::geom_errorbar(ggplot2::aes(ymin=lower,ymax=upper),width=25,linewidth=.35)+
      ggplot2::geom_line(linewidth=.65)+ggplot2::geom_point(size=2)+
      ggplot2::facet_wrap(ggplot2::vars(method,scenario),ncol=3,scales='free_y')+
      ggplot2::scale_x_continuous(breaks=sort(unique(d$n)))+
      ggplot2::scale_colour_manual(values=c('Mean estimated variance'='#0072B2','Monte Carlo variance'='#D55E00'))+
      ggplot2::scale_linetype_manual(values=c('Mean estimated variance'='solid','Monte Carlo variance'='dashed'))+
      ggplot2::scale_shape_manual(values=c('Mean estimated variance'=16,'Monte Carlo variance'=17))+
      ggplot2::labs(x='Sample size',y=if(target=='Total adjusted variance') 'Variance multiplied by n / 4' else expression(sigma[m]^2),
        title=target,
        subtitle=paste0(unique(d$outer_repetitions), ' outer datasets; ',unique(d$mc_repetitions),' independent Monte Carlo datasets'),
        colour=NULL,linetype=NULL,shape=NULL,
        caption='Bars: approximate 95% Monte Carlo uncertainty intervals. Panel y-scales differ.')+
      ggplot2::theme_bw(base_size=11)+ggplot2::theme(legend.position='bottom',strip.background=ggplot2::element_rect(fill='grey93'))
    name=switch(target,'Bootstrap component'='bootstrap_component',
      'Oracle bootstrap component'='oracle_bootstrap_component','Total adjusted variance'='total_adjusted_variance')
    ggplot2::ggsave(file.path(output_dir,paste0(name,'.png')),p,width=11,height=6.4,dpi=180)
    ggplot2::ggsave(file.path(output_dir,paste0(name,'.pdf')),p,width=11,height=6.4)
  }
  invisible(summary)
}
# One outer dataset AND one independent reference dataset per scenario/n cell.
# Across i=1:1000 these become 1000 outer estimates and 1000 MC observations.
#' Return the fixed trial-size grid
#'
#' Define the supported sample sizes shared by simulation, validation, and seed encoding.
#'
#' @return Integer vector c(300,600,900,1200,1500,1800,2100).
vc_sample_sizes = function() seq(300L,2100L,by=300L)

#' Construct the method-by-size-by-scenario study grid
#'
#' Validate requested subsets against the fixed study grid, without changing seed slots.
#'
#' @param n_values Distinct supported trial sizes from vc_sample_sizes(); defaults to all
#'   seven.
#' @param scenarios Distinct scenario identifiers from c('i','ii','iii'); defaults to all
#'   three.
#' @return Data frame with method, n, and scenario for every combination; two methods per
#'   scenario/n cell.
vc_grid = function(n_values=vc_sample_sizes(),scenarios=c('i','ii','iii')) {
  stopifnot(length(n_values)>0,length(scenarios)>0,!anyDuplicated(n_values),
    !anyDuplicated(scenarios),all(n_values %in% vc_sample_sizes()),
    all(scenarios %in% c('i','ii','iii')))
  expand.grid(method=c('Proposed-t-test','Proposed-lm'),n=n_values,
              scenario=scenarios,stringsAsFactors=FALSE)
}
#' Fingerprint the executable study functions
#'
#' Hash deparsed definitions from the fixed .vc_function_names list to identify
#' incompatible simulation implementations.
#'
#' @return Character scalar MD5 code identifier.
#' @details Reads functions in the enclosing source environment, includes the saved R
#'   allocation fallback, writes a temporary file, and removes it on exit. Source
#'   comments are not part of the deparsed definitions.
vc_code_id = function() {
  nm=.vc_function_names
  # Fingerprint callbacks take function name k; the predicate returns a logical
  # and the following callback returns one deparsed definition string.
  nm=nm[vapply(nm,function(k)is.function(get(k)),logical(1))]
  code=setNames(lapply(nm,function(k)paste(deparse(get(k),width.cutoff=500L),collapse='\n')),nm)
  f=tempfile(fileext='.rds'); on.exit(unlink(f))
  text=paste(names(code),unlist(code,use.names=FALSE),sep='\n',collapse='\n')
  writeBin(charToRaw(enc2utf8(text)),f)
  unname(tools::md5sum(f))
}
#' Record simulation settings and package versions
#'
#' Attach a study identifier and code fingerprint so plotting code can check that saved
#' batches are compatible.
#'
#' @param b Integer bootstrap count.
#' @param effect True treatment contrast.
#' @param base_seed Integer base for seed encoding.
#' @param quadrature Population quadrature order.
#' @param edge_grid Population edge interpolation-grid size.
#' @param backend Selected engine label 'R' or 'Rcpp'.
#' @return One-row data frame containing study, code_id, B, effect, base_seed,
#'   quadrature, edge_grid, backend, and R/mgcv/MASS/splines2/Rcpp versions.
#' @details Rcpp_version is recorded as unused for the R backend. Version lookups require
#'   the relevant installed packages.
vc_metadata = function(b,effect,base_seed,quadrature,edge_grid,backend) {
  data.frame(study='empirical_X_only_variance_v8_2100_s2_conditionalX2_gexp_order3_knot0',code_id=vc_code_id(),
    B=as.integer(b),effect=effect,base_seed=as.integer(base_seed),
    quadrature=as.integer(quadrature),edge_grid=as.integer(edge_grid),backend=backend,
    R_version=as.character(getRversion()),
    mgcv_version=as.character(packageVersion('mgcv')),
    MASS_version=as.character(packageVersion('MASS')),
    splines2_version=as.character(packageVersion('splines2')),
    Rcpp_version=if(backend=='Rcpp')as.character(packageVersion('Rcpp')) else 'unused',
    stringsAsFactors=FALSE)
}
#' Select the R or optional Rcpp engine
#'
#' Initialize the requested fast engine when possible, or remove cached compiled routines
#' to force the R implementation.
#'
#' @param fast Logical flag requesting the Rcpp engine.
#' @return Character backend label 'Rcpp' or 'R'.
#' @details Mutates .vc_engine. Requesting fast does not guarantee compilation succeeds;
#'   the selected backend is recorded in metadata.
vc_set_backend = function(fast) {
  backend=if(fast && vc_enable_fast()) 'Rcpp' else 'R'
  if(backend=='R') {
    found=intersect(c('vc_car_cpp','vc_boot_cpp'),ls(.vc_engine))
    if(length(found)) rm(list=found,envir=.vc_engine)
  }
  backend
}
#' Capture warnings and failures without losing a study cell
#'
#' Evaluate a zero-argument callback, collect unique warning messages, and turn an error
#' into an explicit result object.
#'
#' @param FUN Zero-argument function evaluating one oracle build or one trial-cell
#'   calculation.
#' @return List with value (callback result, or a vc_error list with message) and
#'   warnings (unique warning text joined by a separator).
#' @details Warning handlers muffle warnings after recording their text. The error
#'   handler preserves the condition message for delta() failure rows.
vc_capture = function(FUN) {
  ww=character()
  # Condition callbacks: w is a warning condition, recorded and muffled;
  # e is an error condition, returned as a vc_error object with its message.
  value=tryCatch(withCallingHandlers(FUN(),warning=function(w) {
    ww<<-c(ww,conditionMessage(w));invokeRestart('muffleWarning')
  }),error=function(e)structure(list(message=conditionMessage(e)),class='vc_error'))
  list(value=value,warnings=paste(unique(ww),collapse=' | '))
}

#' Run one outer replicate over the requested study grid
#'
#' Generate paired feasible working-model results and an independent reference trial for
#' every scenario/sample-size cell, with fixed seed slots and explicit failure records.
#'
#' @param i Positive integer outer repetition identifier.
#' @param b Integer bootstrap count at least 2; default 500.
#' @param n_values Distinct supported trial sizes; defaults to 300 through 2100 in steps
#'   of 300.
#' @param scenarios Distinct scenario identifiers; defaults to c('i','ii','iii').
#' @param effect Finite scalar true treatment contrast; default zero for null
#'   simulations.
#' @param base_seed Positive integer seed-encoding base; default 20260921.
#' @param quadrature Integer population quadrature order at least 8; default 64.
#' @param edge_grid Integer population edge-grid size at least 11; default 1001.
#' @param fast Logical request for the Rcpp engine; default TRUE.
#' @param verbose Logical flag to print per-cell progress; default TRUE.
#' @return Data frame with two method rows per scenario/n cell (42 rows at defaults).
#'   Includes replicate, scenario, n, treatment_estimate, treatment_variance,
#'   adjusted_sigma2, mc_treatment_estimate, base, sigma_E, sigma_m, sigma_m_oracle,
#'   bootstrap_mcse, MC imbalance, legacy aliases, three seed roles, settings/package
#'   metadata, status, error, and warnings.
#' @details treatment_variance = 4*(base+sigma_E+sigma_m)/n, while adjusted_sigma2 is
#'   their unscaled sum. One reference contrast per method/cell is saved, so a reference
#'   variance requires multiple repetitions. A fixed 63-seed block per i separates
#'   original/bootstrap/reference roles and preserves seeds when the grid is subsetted or
#'   reordered. Resets RNGkind; failures yield two explicit NA rows. Returns data only;
#'   vc_write_seed()/vc_write_rdata() save it. The companion report computes null
#'   rejection rates.
delta = function(i, b=500, n_values=vc_sample_sizes(),
                 scenarios=c('i','ii','iii'), effect=0, base_seed=20260921L,
                 quadrature=64L, edge_grid=1001L, fast=TRUE, verbose=TRUE) {
  grid=vc_grid(n_values,scenarios)
  ints=c(i,b,base_seed,quadrature,edge_grid)
  stopifnot(length(i)==1,length(b)==1,length(base_seed)==1,
    length(quadrature)==1,length(edge_grid)==1,
    all(is.finite(ints)),all(ints==floor(ints)),i>=1,b>=2,base_seed>=1,
    quadrature>=8,edge_grid>=11,length(effect)==1,is.finite(effect))
  # Fixed 63-seed block for each i: 3 scenarios x 7 n values x 3 roles.
  # Never use length(n_values) here: subsetting/reordering must preserve each seed.
  full_grid=vc_sample_sizes()
  seed_block=3L*length(full_grid)*3L
  if(as.double(base_seed)+seed_block*as.double(i)>.Machine$integer.max)
    stop('Replicate identifier too large for the seed encoding.')
  RNGkind('Mersenne-Twister','Inversion','Rejection')
  backend=vc_set_backend(fast)
  meta=vc_metadata(b,effect,base_seed,quadrature,edge_grid,backend)
  rows=list(); k=0L
  for(sc in scenarios) {
    # Zero-argument oracle callback returns the population-reference list for sc.
    oracle_result=vc_capture(function()vc_make_oracle(sc,quadrature,edge_grid))
    for(n in n_values) {
      if(verbose) message('Replicate ',i,': Scenario ',sc,', n=',n,', B=',b)
      cell=(match(sc,c('i','ii','iii'))-1L)*length(full_grid)+match(n,full_grid)-1L
      start=as.double(base_seed)+seed_block*(as.double(i)-1)+3*cell
      if(inherits(oracle_result$value,'vc_error')) {
        z=list(value=oracle_result$value,warnings=oracle_result$warnings)
      # Zero-argument cell callback returns the paired outer/reference data frame;
      # vc_capture preserves warnings or a failure message for this study cell.
      } else z=vc_capture(function() {
        oracle=oracle_result$value
        outer=vc_outer_one(start,start+1,n,sc,oracle,b,effect)
        ref=vc_reference_one(start+2,n,sc,oracle,effect)
        stopifnot(identical(outer$method,ref$method))
        value=data.frame(outer,mc_estimate=ref$estimate,mc_imbalance=ref$imbalance)
        numbers=value[vapply(value,is.numeric,logical(1))]
        # Finiteness callback takes one numeric result column x and returns whether every
        # entry is finite.
        if(!all(vapply(numbers,function(x)all(is.finite(x)),logical(1))))
          stop('Nonfinite result; no variance values may be silently omitted.')
        value
      })
      if(inherits(z$value,'vc_error')) {
        value=data.frame(method=c('Proposed-t-test','Proposed-lm'),estimate=NA_real_,
          sigma_m=NA_real_,sigma_m_oracle=NA_real_,sigma_E=NA_real_,base=NA_real_,
          variance=NA_real_,bootstrap_mcse=NA_real_,mc_estimate=NA_real_,mc_imbalance=NA_real_)
        status='error'; error=z$value$message
      } else {value=z$value;status='ok';error=''}
      # Explicit treatment-effect columns, with legacy names kept for clarity.
      # variance already includes 4/n. Do not scale it by 4/n again.
      value$treatment_estimate=value$estimate
      value$treatment_variance=value$variance
      value$adjusted_sigma2=value$base+value$sigma_E+value$sigma_m
      value$mc_treatment_estimate=value$mc_estimate
      k=k+1L
      rows[[k]]=data.frame(replicate=as.integer(i),scenario=sc,n=n,value,
        scenario2_interior_knots=if(sc=="ii") vc_scenario2_nknots(n) else NA_integer_,
        scenario2_knot_location=if(sc=="ii") 0 else NA_real_,
        scenario2_knot_rule=if(sc=="ii") "fixed_zero" else NA_character_,
        scenario2_X_transform=if(sc=="ii") "conditional_normal[-2,2]" else NA_character_,
        scenario2_X_lower=if(sc=="ii") -2 else NA_real_,
        scenario2_X_upper=if(sc=="ii") 2 else NA_real_,
        scenario2_g_function=if(sc=="ii") "exp(X)" else NA_character_,
        scenario2_spline_order=if(sc=="ii") 3L else NA_integer_,
        scenario2_spline_degree=if(sc=="ii") 2L else NA_integer_,
        scenario2_feature_dimension=if(sc=="ii") 4L else NA_integer_,
        outer_seed=start,bootstrap_seed=start+1,reference_seed=start+2,
        status=status,error=error,
        warnings=paste(unique(c(oracle_result$warnings,z$warnings)[nzchar(c(oracle_result$warnings,z$warnings))]),collapse=' | '),
        meta,stringsAsFactors=FALSE)
    }
  }
  ans=do.call(rbind,rows);rownames(ans)=NULL
  ans
}

#' Write one completed replicate as an atomic CSV file
#'
#' Save all rows for one repetition using a temporary file and final rename.
#'
#' @param result Data frame from delta() containing exactly one replicate identifier.
#' @param output_dir Character destination directory; defaults to the named seed-results
#'   folder.
#' @param overwrite Logical flag permitting deliberate replacement of an existing result;
#'   default FALSE.
#' @return Normalized path to variance_seed_XXXX.csv; creates the output directory and
#'   writes the CSV.
#' @details Missing values are written as NA. Refuses an existing file unless overwrite
#'   is TRUE and removes a leftover temporary file on exit.
vc_write_seed = function(result,output_dir='seed_results_2100_s2_conditionalX2_gexp_order3_knot0',overwrite=FALSE) {
  stopifnot(is.data.frame(result),length(unique(result$replicate))==1L)
  i=unique(result$replicate)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  path=file.path(output_dir,sprintf('variance_seed_%04d.csv',i))
  if(file.exists(path)&&!overwrite) stop('Result already exists: ',path,
    '. Use overwrite=TRUE only when deliberately replacing this seed.')
  # Atomic final rename prevents the combiner from reading a partial file.
  tmp=tempfile(pattern='.writing_',tmpdir=output_dir,fileext='.csv')
  on.exit(if(file.exists(tmp))unlink(tmp))
  write.csv(result,tmp,row.names=FALSE,na='NA')
  if(!file.rename(tmp,path)) stop('Could not rename completed result to ',path)
  normalizePath(path)
}

# Save in the same seed1.RData ... seed1000.RData format used previously.
#' Write one completed replicate as an atomic RData file
#'
#' Save the result data frame as object res1 in the seedN.RData format read by the
#' companion report.
#'
#' @param result Nonempty data frame from delta(), with exactly one finite replicate
#'   identifier.
#' @param output_dir Character destination directory; defaults to the named seed-results
#'   folder.
#' @param overwrite Logical flag permitting deliberate replacement of an existing result;
#'   default FALSE.
#' @return Normalized path to seedN.RData; creates the output directory and saves res1.
#' @details A temporary file and final rename prevent readers from loading partially
#'   written results. Refuses existing files by default and cleans temporary files on
#'   exit.
vc_write_rdata = function(result,output_dir='seed_results_2100_s2_conditionalX2_gexp_order3_knot0',overwrite=FALSE) {
  stopifnot(is.data.frame(result),nrow(result)>0L,
    length(unique(result$replicate))==1L,all(is.finite(result$replicate)))
  i=unique(result$replicate)
  dir.create(output_dir,recursive=TRUE,showWarnings=FALSE)
  path=file.path(output_dir,paste0('seed',i,'.RData'))
  if(file.exists(path)&&!overwrite) stop('Result already exists: ',path,
    '. Use overwrite=TRUE only when deliberately replacing this seed.')
  res1=result
  tmp=tempfile(pattern='.writing_',tmpdir=output_dir,fileext='.RData')
  on.exit(if(file.exists(tmp))unlink(tmp))
  save(res1,file=tmp)
  if(!file.rename(tmp,path)) stop('Could not rename completed result to ',path)
  normalizePath(path)
}

# Fixed function list prevents unrelated caller definitions from changing the study code
# fingerprint.
.vc_function_names = c("vc_scenario2_rX","vc_scenario2_nknots","vc_scenario2_design","vc_sample_sizes","vc_write_rdata","delta","eb_assign_car","eb_car_design","eb_fit_mixed_mean","vc_assign_R","vc_boot","vc_capture","vc_clip","vc_code_id","vc_dataset","vc_enable_fast","vc_fun","vc_gauss01","vc_grid","vc_make_oracle","vc_metadata","vc_normal_clip_moments","vc_normal_nodes","vc_outer_one","vc_plot","vc_prepare1.dis","vc_prepare1.ind","vc_prepare1.ts","vc_prepare2.dis","vc_prepare2.ind","vc_prepare2.ts","vc_reference_one","vc_set_backend","vc_suffix","vc_summarize","vc_variance_se","vc_write_seed")
