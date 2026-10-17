# File guide: helpers construct distributions, allocation features, and tests.
# Test suffixes dis, ind/c2, and ts/c3 identify Scenarios (i), (ii), and (iii).
# Functions return replicate-level decisions; the Rmd averages saved decisions.
# All function interfaces are documented below in roxygen2-style comments.
# These are standalone scripts; no documentation-generation package is needed.
# Source simulation scripts separately because they define overlapping names.
# The current delta() defaults use g(X)=2*X in (i), exp(X) in (ii), and
# 2*X1*X2^2 with X clipped to [-1,1] in (iii). Z is clipped to [-3,3].
# f(Z)=exp(Z)/2.5 in (i)/(iii) and 4*log(Z+4) in (ii).
# Treatment labels are 1 and 0; the treatment-mean contrast is u.
# The randomization, fit, and variance calculations are in the test bodies.

# Scenario (ii): X ~ N(0,1) conditional on -2 < X < 2; g(X)=exp(X).
# This is continuous truncation, NOT clipping/winsorization: no endpoint atoms.
# All nine Scenario (ii) methods use the same conditional-X helper. X is not
# standardized after sampling. The default delta() sets g(X)+f(Z) to
# exp(X)+4*log(Z+4) for EVERY Scenario (ii) proposed/comparison method.
# The original latent (X*,Z*) draw is retained; only X* is mapped to the new
# conditional distribution. No random draws are added; Z remains clipped to
# [-3,3]. Scenario (ii) noise and the conditional independence of X and Z remain.
# The original spline orders/degrees, knot rules, comparator bins, nuisance
# fitting, variance formulas, bootstrap logic, function arguments, output
# ordering, and Scenarios (i)/(iii) are unchanged from the supplied source.
# Low-level test functions still accept an explicit fun argument; for this
# setting, pass function(x,z) exp(x)+4*log(z+4).
# Save results in a fresh directory and source the three scripts separately,
# since they define overlapping function names, including delta().

# If X* ~ N(0,1), Phi(X*) ~ Uniform(0,1). The returned X has conditional CDF
# (Phi(x)-Phi(-2))/(Phi(2)-Phi(-2)), for -2 < x < 2.
# This helper does not consume additional random numbers.


#' Map latent normal draws to the Scenario (ii) conditional distribution
#'
#' Apply the inverse-CDF transform for a standard normal conditioned on -2 < X < 2, using
#' the supplied latent draws and no additional randomness.
#'
#' @param x_star Nonempty numeric vector of finite latent standard-normal draws.
#' @return Numeric vector of the same length, with every value strictly between -2 and 2.
#'   Invalid inputs or numerical endpoint values cause an error.
#' @details Phi(x_star) supplies uniform quantiles; rescaling those quantiles between
#'   Phi(-2) and Phi(2) gives the conditional distribution.
s2_conditional_x <- function(x_star) {
  if (!is.numeric(x_star) || is.complex(x_star) ||
      length(x_star) == 0L || any(!is.finite(x_star))) {
    stop("x_star must contain finite standard-normal draws.")
  }
  # Rescale the latent normal quantile to the conditional probability interval.
  p_lower <- stats::pnorm(-2)
  p_upper <- stats::pnorm(2)
  u <- stats::pnorm(as.numeric(x_star))
  x <- stats::qnorm(p_lower + (p_upper - p_lower) * u)
  if (any(!is.finite(x)) || any(x <= -2 | x >= 2)) {
    stop("Conditional-normal inverse CDF reached an endpoint numerically; no clipping was applied.")
  }
  x
}

# Proposed-lm X-and-Z revision (2026-09-21).
# All test2 functions adjust for ALL X coordinates and Z.
# Their variance follows the manuscript table with U=(X,Z) and the full
# covariance projection term. The bootstrap response is fhat(Z)-beta_Z_hat*Z.
# All slopes, fitted means and conditional residual components are frozen.
# Proposed-plm uses SSE/n as specified in the same table.
# All nine Setting 3 methods use error variance log(1+2*x1^2+1.5*z^2).
# Source this complete file instead of the previous empirical-bootstrap file.
# Original public test names and delta output ordering are retained.
# Source this file; delta(i) retains n=1000, b=500 and four effect sizes.
# Quadratic = degree 2 = order 3. Scenario 3 uses the COMPLETE tensor basis.
# Only sigma_m^2 is bootstrapped; sigma_E^2 and all nuisance fits use original data.
library(mgcv)
library(MASS)
library(splines2)
library(carat)

# Empirical covariate bootstrap: no histogram, jitter, or outcome resampling.
# Every bootstrap path starts at zero and uses the original balancing features.
#' Check the dependency for the comparison methods
#'
#' Check whether carat can be loaded before delta() runs its comparison methods.
#'
#' @return NULL when the package is available; otherwise an error explaining the missing
#'   dependency.
eb_require_carat <- function() {
  if (!requireNamespace("carat", quietly = TRUE))
    stop("The comparison methods require package 'carat'. Install it before calling delta().")
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
  # Knot callback: input d is a coordinate index; output is its ordered
  # observed endpoints and zero, used to identify the local allocation cell.
  knots = lapply(seq_len(p), function(d) unique(sort(c(range(x[, d]), 0))))
  # Basis callback: input d is a coordinate index; output is the complete
  # n-row marginal B-spline matrix for that coordinate.
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
  # Cell callbacks: outer input i is a row index and output is its character
  # cell key. Inner input d is a coordinate index and output is the integer
  # knot-interval number for x[i,d]. Joining these numbers defines one cell.
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
  # Each row stores one cell imbalance; every call starts a fresh CAR path.
  Lambda = matrix(0, design$J, ncol(phi))
  treat = numeric(length(index))
  for (i in seq_along(index)) {
    j = cell[i]
    test_val = sum(Lambda[j, ] * phi[i, ])
  # A negative inner product favors adding phi (treatment 1); ties are fair.
    prob = if (test_val < 0) rho else if (test_val > 0) 1 - rho else 0.5
    treat[i] = rbinom(1, 1, prob)
    Lambda[j, ] = Lambda[j, ] + (2 * treat[i] - 1) * phi[i, ]
  }
  treat
}

# X in scenario 3 has an interior, four edges and four corner atoms.
# Estimate E(response | X) separately on these observed support components.
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

#' Evaluate a two-sided normal test for a treatment contrast
#'
#' Compare the absolute standardized contrast with the 0.975 normal quantile.
#'
#' @param estimate Scalar estimated treatment contrast under the null of zero contrast.
#' @param variance Scalar variance of estimate, already including the 4/n scaling.
#' @return Numeric rejection indicator: 1 for rejection, 0 otherwise; NA_real_ with a
#'   warning for nonpositive or nonfinite variance.
eb_reject <- function(estimate, variance) {
  if (!is.finite(variance) || variance <= 0) {
    warning("Nonpositive or nonfinite variance estimate: returning NA; do not silently omit it.")
    return(NA_real_)
  }
  as.numeric(abs(estimate) / sqrt(variance) >= qnorm(0.975))
}

#' Validate the number of empirical bootstrap replicates
#'
#' Require enough bootstrap paths to compute a sample variance.
#'
#' @param b Scalar finite integer bootstrap count, at least 2.
#' @return NULL for a valid count; otherwise an error.
eb_check_bootstrap <- function(b) {
  if (length(b) != 1L || !is.finite(b) || b < 2 || b != as.integer(b))
    stop("b must be an integer of at least 2.")
}


#' Assign treatments by global weighted covariate balancing
#'
#' Sequentially favor the treatment that reduces the squared imbalance of an intercept
#' and the supplied features; maintain one global imbalance vector.
#'
#' @param x Numeric n-by-p feature matrix, in participant arrival order.
#' @param w0 Nonnegative weight for the overall treatment-count imbalance.
#' @param w1 Nonnegative weight for the feature imbalance; w0 + w1 must be positive.
#' @param rho Biased-coin probability of choosing the better-balanced treatment; default
#'   0.9.
#' @return Numeric vector of length n with treatment labels 0 or 1, in row order.
#' @details Weights are normalized by their sum, features are weighted by their square
#'   roots, and ties use probability 0.5. Despite the historical name, the implementation
#'   is R code. Random draws use the current RNG state.
covRcpp11 <- function(x, w0 = 1, w1 = 1, rho = 0.9) {
  w <- w0 + w1
  w0 <- w0 / w
  w1 <- w1 / w
  
  n <- nrow(x)
  p <- ncol(x)
  T_ls <- numeric(n)
  # Lambda accumulates signed treatment imbalance in the weighted features.
  Lambda <- rep(0, 1 + p)
  
  Test <- 0
  
  for (i in 1:n) {
    phi0 <- sqrt(w0)
    phi1 <- sqrt(w1) * x[i, ]
    phi <- c(phi0, phi1)
    
    Test1 <- sum((Lambda + phi)^2)
    Test0 <- sum((Lambda - phi)^2)
    
    if (Test1 < Test0) {
      T_ls[i] <- sum(runif(1) < rho)
    } else if (Test1 > Test0) {
      T_ls[i] <- sum(runif(1) < (1 - rho))
    } else if (Test1 == Test0) {
      T_ls[i] <- sum(runif(1) < 0.5)
    }
    
    Lambda <- Lambda + (2 * T_ls[i] - 1) * phi
  }
  
  return(T_ls)
}

#' Allocate using the multivariate tensor-spline CAR design
#'
#' Construct a ts design and allocate participants through eb_assign_car().
#'
#' @param x Numeric n-by-p matrix on [-1,1], with values on both sides of zero in each
#'   coordinate.
#' @param degree Marginal B-spline degree; default 2 means order 3.
#' @param rho Probability favoring the treatment that reduces local imbalance; default
#'   0.9.
#' @return Numeric length-n vector of treatment labels 0 or 1.
#' @details Uses zero-split cells and the complete tensor span in the intercept-plus-all-
#'   but-first-column parameterization of eb_car_design().
covRcpp_multi <- function(x, degree = 2, rho = 0.9) {
  # Complete tensor-product basis on [-1,1]^p; knot 0 in each coordinate.
  # degree = 2 means quadratic (order 3), not order 2 (linear).
  # Separate imbalance vectors are maintained for the cells split at zero.
  eb_assign_car(eb_car_design(x, "ts", k = degree), rho = rho)
}

#' Allocate using the univariate quantile-spline CAR design
#'
#' Construct an ind design and allocate participants through eb_assign_car().
#'
#' @param x Numeric vector of observed covariates in arrival order.
#' @param J Number of empirical quantile intervals; default 2.
#' @param k Marginal B-spline degree; default 3 means order 4.
#' @param rho Probability favoring the treatment that reduces local imbalance; default
#'   0.9.
#' @return Numeric vector of length length(x), with treatment labels 0 or 1.
covRcpp_uni_b <- function(x, J = 2, k = 3, rho = 0.9) {
  eb_assign_car(eb_car_design(x, "ind", J = J, k = k), rho = rho)
}

#' Allocate using the discrete-stratum CAR design
#'
#' Construct a dis design and balance the global intercept and stratum indicators.
#'
#' @param x Observed discrete covariate vector or matrix in arrival order.
#' @param rho Probability favoring the treatment that reduces imbalance; default 0.9.
#' @return Numeric vector of treatment labels 0 or 1, one per observation.
covRcpp_uni_s <- function(x, rho = 0.9) {
  eb_assign_car(eb_car_design(x, "dis"), rho = rho)
}

#' Choose a sign-dependent allocation probability
#'
#' Map a scalar imbalance sign to the supplied positive- and negative-sign probabilities.
#'
#' @param x Scalar imbalance statistic.
#' @param p Probability returned when x < 0.
#' @param q Probability returned when x > 0.
#' @return Scalar q for positive x, p for negative x, or 0.5 for a tie.
#' @details This is a deterministic helper; despite its historical name, the
#'   implementation is R code.
g_cpp <- function(x, p, q) {
  if (x > 0) {
    y <- q
  } else if (x < 0) {
    y <- p
  } else {
    y <- 0.5
  }
  return(y)
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

#' Simulate the proposed two-sample test in Scenario (i)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Fit original-data nuisance models once, resample observed donor rows, and restart CAR
#' to estimate the variance of the fitted conditional-mean imbalance.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param b Bootstrap count, integer at least 2, for estimating sigma_m^2.
#' @param fun Vectorized prognostic mean function fun(x, z), returning one g(X) + f(Z)
#'   value per participant.
#' @param details Logical; FALSE returns rejection decisions only, TRUE returns the
#'   diagnostic list described below.
#' @return By default, a length-2 numeric vector c(adjusted, naive), with 1 for
#'   rejection, 0 otherwise, and possible NA for invalid variances. With details = TRUE,
#'   a list containing rejection, estimate (treatment contrast), variance and
#'   variance_naive (contrast variances), sigma_E2 (conditional residual variance),
#'   sigma_m2 (bootstrap imbalance variance), noise_component, fitted_m (length-n fixed
#'   fitted means), and bootf (length-b imbalance statistics).
test1.dis <- function(i, n, u, b, fun, details = FALSE){
  
  eb_check_bootstrap(b)
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
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  # Vgf is the fitted prognostic variance to subtract from the working residual variance.
  vgf = var(hats)
  fzmd = gam(fz ~ factor(x))
  if (!isTRUE(fzmd$converged)) stop("Conditional-mean GAM did not converge.")
  efz = as.vector(fitted(fzmd))

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fz - efz)

  bootf = rep(0, b)
  for (j in seq_len(b)) {
    resample = sample.int(n, size = n, replace = TRUE)
    # The same donor row selects X, its FIXED features, and fitted m(X).
    # With two covariates this is joint row resampling, not independent margins.
    treatnew = eb_assign_car(design, index = resample)
    efznew = efz[resample]
    bootf[j] = sum((2 * treatnew - 1) * efznew) / sqrt(n)
  }
  # Sample variance over b independent CAR bootstrap paths estimates sigma_m^2.
  sigma.m2 = var(bootf)
  sige2 = sig.epsi2 - vgf + evzx + sigma.m2

  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sige2 / n
  variance.naive = 4 * sig.epsi2 / n
  estimate = unname(hat.t[1] - hat.t[2])
  prob = eb_reject(estimate, variance)
  prob.n = eb_reject(estimate, variance.naive)
  if (details) return(list(rejection = c(prob, prob.n), estimate = estimate,
                           variance = variance, variance_naive = variance.naive,
                           sigma_E2 = evzx, sigma_m2 = sigma.m2,
                           noise_component = sig.epsi2 - vgf,
                           fitted_m = efz, bootf = bootf))
  return(c(prob, prob.n))
}

#' Simulate the proposed linear working-model test in Scenario (i)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Fit original-data nuisance models once, resample observed donor rows, and restart CAR
#' to estimate the variance of the fitted conditional-mean imbalance.
#'
#' @return By default, a length-2 numeric vector c(adjusted, naive), with 1 for
#'   rejection, 0 otherwise, and possible NA for invalid variances. With details = TRUE,
#'   a list containing rejection, estimate (treatment contrast), variance and
#'   variance_naive (contrast variances), sigma_E2 (conditional residual variance),
#'   sigma_m2 (bootstrap imbalance variance), noise_component, fitted_m (length-n fixed
#'   fitted means), and bootf (length-b imbalance statistics). The linear-model list also
#'   contains working_residual_variance (SSE/residual df), signal_variance (Vgf),
#'   projection_correction (Q), working_covariates (covariate names),
#'   working_coefficients (fitted slopes), conditional_response (fitted f(Z) minus the
#'   working Z slope times Z), and residual_df (working-model residual degrees of
#'   freedom).
test2.dis <- function(i, n, u, b, fun, details = FALSE){
  
  eb_check_bootstrap(b)
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
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  # Vgf is the fitted prognostic variance to subtract from the working residual variance.
  vgf = var(hats)
  mnew = cbind(x, z)
  cm = cov(mnew)
  v = cov(mnew, hats)
  # Q adds back the prognostic variation already removed by the linear working model.
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

  bootf = rep(0, b)
  for (j in seq_len(b)) {
    resample = sample.int(n, size = n, replace = TRUE)
    # The same donor row selects X, its FIXED features, and fitted m(X).
    # With two covariates this is joint row resampling, not independent margins.
    treatnew = eb_assign_car(design, index = resample)
    efznew = efz[resample]
    bootf[j] = sum((2 * treatnew - 1) * efznew) / sqrt(n)
  }
  # Sample variance over b independent CAR bootstrap paths estimates sigma_m^2.
  sigma.m2 = var(bootf)
  # Exact plug-in decomposition specified in the manuscript table.
  noise.component = sig.epsi2 - vgf + e1
  sige2 = noise.component + evzx + sigma.m2

  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sige2 / n
  variance.naive = 4 * sig.epsi2 / n
  estimate = unname(hat.t[1] - hat.t[2])
  prob = eb_reject(estimate, variance)
  prob.n = eb_reject(estimate, variance.naive)
  if (details) return(list(rejection = c(prob, prob.n), estimate = estimate,
                           variance = variance, variance_naive = variance.naive,
                           sigma_E2 = evzx, sigma_m2 = sigma.m2,
                           noise_component = noise.component,
                           working_residual_variance = sig.epsi2,
                           signal_variance = vgf,
                           projection_correction = e1,
                           working_covariates = colnames(mnew),
                           working_coefficients = beta.work,
                           conditional_response = fzminus,
                           residual_df = df.residual(md.l),
                           fitted_m = efz, bootf = bootf))
  return(c(prob, prob.n))
}

#' Simulate the proposed partially linear test in Scenario (i)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Estimate the treatment contrast directly from the outcome GAM and use its residual
#' mean square SSE/n.
#'
#' @return Scalar numeric rejection indicator (0/1, or NA for invalid variance). With
#'   details = TRUE, a list with rejection, estimate (treatment contrast), variance
#'   (contrast variance), and noise_variance (SSE/n).
test3.dis <- function(i, n, u, b, fun, details = FALSE){
  
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
  
  md = gam(Y ~ 0 + treat + treat2 + x + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  hat.t = coef(md)[c("treat", "treat2")]
  sig.epsi2 = sum(residuals(md)^2) / n  # manuscript table
  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sig.epsi2 / n
  prob = eb_reject(hat.t[1] - hat.t[2], variance)
  if (details) return(list(rejection = prob, estimate = unname(hat.t[1] - hat.t[2]),
                           variance = variance, noise_variance = sig.epsi2))
  return(prob)
}

#' Simulate the stratified comparison test, variant 1, in Scenario (i)
#'
#' Allocate with carat::StrPBR(). Use a treatment-only contrast and the within-/between-
#' stratum variance decomposition.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest1.p.dis <- function(i, n, u, fun){
  eb_require_carat()
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
  
  bs = matrix(rep(0, 4 * n), ncol = 4)
  for (i in 1: n) {
    if(x[i] == -2){
      bs[i, 1] = 1
    }
    if(x[i] == -1){
      bs[i, 2] = 1
    }
    if(x[i] == 1){
      bs[i, 3] = 1
    }
    if(x[i] == 2){
      bs[i, 4] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0)
  X.d = data.frame(xd = x, zd = rep(1, n))
  
  trt <- ifelse(StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  # Do not discard a populated stratum with an unobserved treatment arm.
  used = colSums(bs) > 0
  arm1 = colSums(bs[trt == 1, , drop = FALSE])
  arm0 = colSums(bs[trt == 0, , drop = FALSE])
  if (any(used & (arm1 == 0 | arm0 == 0))) {
    warning("A populated stratum has an empty treatment arm; returning NA.")
    return(NA_real_)
  }
  
  yk1 = rep(0, 4)
  for (i in 1: 4) {
    yk1[i] = mean(Y.f[trt == 1 & bs[,i] == 1])
  }
  
  yk0 = rep(0, 4)
  for (i in 1: 4) {
    yk0[i] = mean(Y.f[trt == 0 & bs[,i] == 1])
  }
  
  yk1.square = rep(0, 4)
  for (i in 1: 4) {
    yk1.square[i] = sum((Y.f[trt == 1 & bs[,i] == 1] - yk1[i])^2) * (sum( bs[,i] == 1) / n) / sum(trt == 1 & bs[,i] == 1)
  }
  
  yk0.square = rep(0, 4)
  for (i in 1: 4) {
    yk0.square[i] = sum((Y.f[trt == 0 & bs[,i] == 1] - yk0[i])^2) * (sum( bs[,i] == 1) / n) / sum(trt == 0 & bs[,i] == 1)
  }
  
  # Aggregate within-stratum variances from the two treatment groups.
  sig1.square = 2 * sum(yk1.square, na.rm = T) + 2 * sum(yk0.square, na.rm = T)
  
  ay = rep(0, 4)
  for (i in 1: 4) {
    ay[i] = (yk1[i] - mean(Y.f[trt == 1]) - (yk0[i] - mean(Y.f[trt == 0])))^2 * (sum( bs[,i] == 1) / n)
  }
  
  # Add the between-stratum contribution to the treatment-contrast variance.
  sig2.square = sum(ay, na.rm = T)
  md.l.f <- lm(Y.f ~ trt)
  
  hat.t.f = md.l.f$coefficients[2]
  
  
  s.f <- hat.t.f / sqrt((sig1.square + sig2.square) / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratified comparison test, variant 2, in Scenario (i)
#'
#' Allocate with carat::StrBCD(). Regress on treatment and stratum indicators, plus Z
#' where included by this scenario, and use the treatment coefficient standard error.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest2.p.dis <- function(i, n, u, fun){
  eb_require_carat()
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
  
  bs = matrix(rep(0, 4 * n), ncol = 4)
  for (i in 1: n) {
    if(x[i] == -2){
      bs[i, 1] = 1
    }
    if(x[i] == -1){
      bs[i, 2] = 1
    }
    if(x[i] == 1){
      bs[i, 3] = 1
    }
    if(x[i] == 2){
      bs[i, 4] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0)
  X.d = data.frame(xd = x, e = rep(1, n))
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  
  md.l.f <- lm(Y.f ~ trt + bs + z)
  hat.t.f = md.l.f$coefficients[2]
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the variance-ratio-corrected comparator, variant 1, in Scenario (i)
#'
#' Allocate with carat::StrBCD() and compute the variant-1 treatment contrast using the
#' tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest.bnz1.dis <- function(i, n, u, fun){
  eb_require_carat()
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
  
  bs = matrix(rep(0, 5 * n), ncol = 5)
  for (i in 1: n) {
    if(x[i] == -2){
      bs[i, 1] = 1
    }
    if(x[i] == -1){
      bs[i, 2] = 1
    }
    if(x[i] == 0){
      bs[i, 3] = 1
    }
    if(x[i] == 1){
      bs[i, 4] = 1
    }
    if(x[i] == 2){
      bs[i, 5] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0)
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5])
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  dis2 = (z > 0); dis1 = (x > 0)
  
  e2.f = beta2[4]^2 * (var(z[dis2]) * sum(dis2) / n + 
                         var(z[!dis2] ) * (n - sum(dis2)) / n) 
  e1.f = beta2[3]^2 * var(x) + beta2[4]^2 * var(z)
  
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  M = matrix(c(trt, 1 - trt), ncol = 2)
  L = matrix(c(1, -1), ncol = 2)
  # L contrasts the two treatment means; this is its regression design variance factor.
  temp = L %*% solve(t(M) %*% M) %*% t(L)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * temp * tau2.f)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the variance-ratio-corrected comparator, variant 2, in Scenario (i)
#'
#' Allocate with carat::StrBCD() and compute the variant-2 treatment contrast using the
#' tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA
fftest.bnz2.dis <- function(i, n, u, fun){
  
  eb_require_carat()
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
  
  bs = matrix(rep(0, 5 * n), ncol = 5)
  for (i in 1: n) {
    if(x[i] == -2){
      bs[i, 1] = 1
    }
    if(x[i] == -1){
      bs[i, 2] = 1
    }
    if(x[i] == 0){
      bs[i, 3] = 1
    }
    if(x[i] == 1){
      bs[i, 4] = 1
    }
    if(x[i] == 2){
      bs[i, 5] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0)
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5])
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx + z)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  dis1 = (x > 0)
  
  e2.f = 0
  e1.f = beta2[3]^2 * var(x)
  
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * tau2.f * 4 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratum-residual variance comparator, variant 1, in Scenario (i)
#'
#' Allocate with carat::AdjBCD(), fit the variant-1 outcome regression, and aggregate
#' residual variances within the comparison strata.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
testye.bnz1.dis <- function(i, n, u, fun){
  eb_require_carat()
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
  
  bs = matrix(rep(0, 5 * n), ncol = 5)
  for (i in 1: n) {
    if(x[i] == -2){
      bs[i, 1] = 1
    }
    if(x[i] == -1){
      bs[i, 2] = 1
    }
    if(x[i] == 0){
      bs[i, 3] = 1
    }
    if(x[i] == 1){
      bs[i, 4] = 1
    }
    if(x[i] == 2){
      bs[i, 5] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0); xd1 = bs[,1]; xd2 = bs[,2]; xd3 = bs[,3]; xd4 = bs[, 4]; xd5 = bs[,5]
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5])
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2, xd3, xd4, xd5))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratum-residual variance comparator, variant 2, in Scenario (i)
#'
#' Allocate with carat::AdjBCD(), fit the variant-2 outcome regression, and aggregate
#' residual variances within the comparison strata.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
testye.bnz2.dis <- function(i, n, u, fun){
  eb_require_carat()
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
  
  bs = matrix(rep(0, 5 * n), ncol = 5)
  for (i in 1: n) {
    if(x[i] == -2){
      bs[i, 1] = 1
    }
    if(x[i] == -1){
      bs[i, 2] = 1
    }
    if(x[i] == 0){
      bs[i, 3] = 1
    }
    if(x[i] == 1){
      bs[i, 4] = 1
    }
    if(x[i] == 2){
      bs[i, 5] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0); xd1 = bs[,1]; xd2 = bs[,2]; xd3 = bs[,3]; xd4 = bs[, 4]; xd5 = bs[,5]
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5])
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2, xd3, xd4, xd5))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the proposed two-sample test in Scenario (ii)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Fit original-data nuisance models once, resample observed donor rows, and restart CAR
#' to estimate the variance of the fitted conditional-mean imbalance.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param b Bootstrap count, integer at least 2, for estimating sigma_m^2.
#' @param fun Vectorized prognostic mean function fun(x, z), returning one g(X) + f(Z)
#'   value per participant.
#' @param details Logical; FALSE returns rejection decisions only, TRUE returns the
#'   diagnostic list described below.
#' @return By default, a length-2 numeric vector c(adjusted, naive), with 1 for
#'   rejection, 0 otherwise, and possible NA for invalid variances. With details = TRUE,
#'   a list containing rejection, estimate (treatment contrast), variance and
#'   variance_naive (contrast variances), sigma_E2 (conditional residual variance),
#'   sigma_m2 (bootstrap imbalance variance), noise_component, fitted_m (length-n fixed
#'   fitted means), and bootf (length-b imbalance statistics).
test1.ind <- function(i, n, u, b, fun, details = FALSE){
  
  eb_check_bootstrap(b)
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  design = eb_car_design(x, "ind")
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
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  # Vgf is the fitted prognostic variance to subtract from the working residual variance.
  vgf = var(hats)
  fzmd = gam(fz ~ s(x))
  if (!isTRUE(fzmd$converged)) stop("Conditional-mean GAM did not converge.")
  efz = as.vector(fitted(fzmd))

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fz - efz)

  bootf = rep(0, b)
  for (j in seq_len(b)) {
    resample = sample.int(n, size = n, replace = TRUE)
    # The same donor row selects X, its FIXED features, and fitted m(X).
    # With two covariates this is joint row resampling, not independent margins.
    treatnew = eb_assign_car(design, index = resample)
    efznew = efz[resample]
    bootf[j] = sum((2 * treatnew - 1) * efznew) / sqrt(n)
  }
  # Sample variance over b independent CAR bootstrap paths estimates sigma_m^2.
  sigma.m2 = var(bootf)
  sige2 = sig.epsi2 - vgf + evzx + sigma.m2

  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sige2 / n
  variance.naive = 4 * sig.epsi2 / n
  estimate = unname(hat.t[1] - hat.t[2])
  prob = eb_reject(estimate, variance)
  prob.n = eb_reject(estimate, variance.naive)
  if (details) return(list(rejection = c(prob, prob.n), estimate = estimate,
                           variance = variance, variance_naive = variance.naive,
                           sigma_E2 = evzx, sigma_m2 = sigma.m2,
                           noise_component = sig.epsi2 - vgf,
                           fitted_m = efz, bootf = bootf))
  return(c(prob, prob.n))
}

#' Simulate the proposed linear working-model test in Scenario (ii)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Fit original-data nuisance models once, resample observed donor rows, and restart CAR
#' to estimate the variance of the fitted conditional-mean imbalance.
#'
#' @return By default, a length-2 numeric vector c(adjusted, naive), with 1 for
#'   rejection, 0 otherwise, and possible NA for invalid variances. With details = TRUE,
#'   a list containing rejection, estimate (treatment contrast), variance and
#'   variance_naive (contrast variances), sigma_E2 (conditional residual variance),
#'   sigma_m2 (bootstrap imbalance variance), noise_component, fitted_m (length-n fixed
#'   fitted means), and bootf (length-b imbalance statistics). The linear-model list also
#'   contains working_residual_variance (SSE/residual df), signal_variance (Vgf),
#'   projection_correction (Q), working_covariates (covariate names),
#'   working_coefficients (fitted slopes), conditional_response (fitted f(Z) minus the
#'   working Z slope times Z), and residual_df (working-model residual degrees of
#'   freedom).
test2.ind <- function(i, n, u, b, fun, details = FALSE){
  
  eb_check_bootstrap(b)
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  design = eb_car_design(x, "ind")
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
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  # Vgf is the fitted prognostic variance to subtract from the working residual variance.
  vgf = var(hats)
  mnew = cbind(x, z)
  cm = cov(mnew)
  v = cov(mnew, hats)
  # Q adds back the prognostic variation already removed by the linear working model.
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

  bootf = rep(0, b)
  for (j in seq_len(b)) {
    resample = sample.int(n, size = n, replace = TRUE)
    # The same donor row selects X, its FIXED features, and fitted m(X).
    # With two covariates this is joint row resampling, not independent margins.
    treatnew = eb_assign_car(design, index = resample)
    efznew = efz[resample]
    bootf[j] = sum((2 * treatnew - 1) * efznew) / sqrt(n)
  }
  # Sample variance over b independent CAR bootstrap paths estimates sigma_m^2.
  sigma.m2 = var(bootf)
  # Exact plug-in decomposition specified in the manuscript table.
  noise.component = sig.epsi2 - vgf + e1
  sige2 = noise.component + evzx + sigma.m2

  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sige2 / n
  variance.naive = 4 * sig.epsi2 / n
  estimate = unname(hat.t[1] - hat.t[2])
  prob = eb_reject(estimate, variance)
  prob.n = eb_reject(estimate, variance.naive)
  if (details) return(list(rejection = c(prob, prob.n), estimate = estimate,
                           variance = variance, variance_naive = variance.naive,
                           sigma_E2 = evzx, sigma_m2 = sigma.m2,
                           noise_component = noise.component,
                           working_residual_variance = sig.epsi2,
                           signal_variance = vgf,
                           projection_correction = e1,
                           working_covariates = colnames(mnew),
                           working_coefficients = beta.work,
                           conditional_response = fzminus,
                           residual_df = df.residual(md.l),
                           fitted_m = efz, bootf = bootf))
  return(c(prob, prob.n))
}

#' Simulate the proposed partially linear test in Scenario (ii)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Estimate the treatment contrast directly from the outcome GAM and use its residual
#' mean square SSE/n.
#'
#' @return Scalar numeric rejection indicator (0/1, or NA for invalid variance). With
#'   details = TRUE, a list with rejection, estimate (treatment contrast), variance
#'   (contrast variance), and noise_variance (SSE/n).
test3.ind <- function(i, n, u, b, fun, details = FALSE){
  
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  
  design = eb_car_design(x, "ind")
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  epsi <- rnorm(n)
  
  Y = treat.mx %*% beta + g + epsi
  
  md = gam(Y ~ 0 + treat + treat2 + s(x) + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  hat.t = coef(md)[c("treat", "treat2")]
  sig.epsi2 = sum(residuals(md)^2) / n  # manuscript table
  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sig.epsi2 / n
  prob = eb_reject(hat.t[1] - hat.t[2], variance)
  if (details) return(list(rejection = prob, estimate = unname(hat.t[1] - hat.t[2]),
                           variance = variance, noise_variance = sig.epsi2))
  return(prob)
}

#' Simulate the variance-ratio-corrected comparator, variant 1, in Scenario (ii)
#'
#' Allocate with carat::StrBCD() and compute the variant-1 treatment contrast using the
#' tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest.bnz1 <- function(i, n, u, b, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  xd = disc2(x, b); zd = disc2(z, b)
  X.d = data.frame(xd = xd, extra = rep(1, n))
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  strx = model.matrix( ~ factor(xd) - 1); strz = model.matrix( ~ factor(zd) - 1)
  varx = 0; varz = 0
  for (k in 1: ncol(strx)) {
    varx = varx + var(x[strx[,k] == 1]) * sum(strx[,k] == 1) / n
  }
  for (k in 1: ncol(strz)) {
    varz = varz + var(z[strz[,k] == 1]) * sum(strz[,k] == 1) / n
  }
  
  # Z is independent of X and is NOT balanced by this design.
  varz = var(z)

  e2.f = beta2[3]^2 * varx + beta2[4]^2 * varz
  e1.f = beta2[3]^2 * var(x) + beta2[4]^2 * var(z)
  
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  M = matrix(c(trt, 1 - trt), ncol = 2)
  L = matrix(c(1, -1), ncol = 2)
  # L contrasts the two treatment means; this is its regression design variance factor.
  temp = L %*% solve(t(M) %*% M) %*% t(L)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * temp * tau2.f)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the variance-ratio-corrected comparator, variant 2, in Scenario (ii)
#'
#' Allocate with carat::StrBCD() and compute the variant-2 treatment contrast using the
#' tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest.bnz2 <- function(i, n, u, b, fun){
  
  eb_require_carat()
  set.seed(i)
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  xd = disc2(x, b); zd = disc2(z, b)
  X.d = data.frame(xd = xd, extra = rep(1, n))
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx + z)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  strx = model.matrix( ~ factor(xd) - 1)
  varx = 0
  for (k in 1: ncol(strx)) {
    varx = varx + var(x[strx[,k] == 1]) * sum(strx[,k] == 1) / n
  }
  
  e2.f = beta2[3]^2 * varx 
  e1.f = beta2[3]^2 * var(x)
  
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * tau2.f * 4 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratified comparison test, variant 1, in Scenario (ii)
#'
#' Allocate with carat::StrPBR(). Use a treatment-only contrast and the within-/between-
#' stratum variance decomposition.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest1.p.c2 <- function(i, n, u, b, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b)
  bs =  model.matrix( ~ factor(dis1) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  
  X.d = data.frame(xd = dis1, zd = rep(1, n))        
  
  trt <- ifelse(StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  # Do not discard a populated stratum with an unobserved treatment arm.
  used = colSums(bs) > 0
  arm1 = colSums(bs[trt == 1, , drop = FALSE])
  arm0 = colSums(bs[trt == 0, , drop = FALSE])
  if (any(used & (arm1 == 0 | arm0 == 0))) {
    warning("A populated stratum has an empty treatment arm; returning NA.")
    return(NA_real_)
  }
  
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
  
  # Aggregate within-stratum variances from the two treatment groups.
  sig1.square = 2 * sum(yk1.square, na.rm = T) + 2 * sum(yk0.square, na.rm = T)
  
  ay = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    ay[i] = (yk1[i] - mean(Y.f[trt == 1]) - (yk0[i] - mean(Y.f[trt == 0])))^2 * (sum( bs[,i] == 1) / n)
  }
  
  # Add the between-stratum contribution to the treatment-contrast variance.
  sig2.square = sum(ay, na.rm = T)
  md.l.f <- lm(Y.f ~ trt)
  
  hat.t.f = md.l.f$coefficients[2]
  
  
  s.f <- hat.t.f / sqrt((sig1.square + sig2.square) / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratified comparison test, variant 2, in Scenario (ii)
#'
#' Allocate with carat::StrBCD(). Regress on treatment and stratum indicators, plus Z
#' where included by this scenario, and use the treatment coefficient standard error.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest2.p.c2 <- function(i, n, u, b, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b)
  bs =  model.matrix(~ factor(dis1) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  xd = dis1
  X.d = data.frame(xd = xd, e = rep(1, n))
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  
  md.l.f <- lm(Y.f ~ trt + bs + z)
  hat.t.f = md.l.f$coefficients[2]
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratum-residual variance comparator, variant 1, in Scenario (ii)
#'
#' Allocate with carat::AdjBCD(), fit the variant-1 outcome regression, and aggregate
#' residual variances within the comparison strata.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
testye.bnz1.c2 <- function(i, n, u, b, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b)
  bs =  model.matrix(~ factor(dis1) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  xd = dis1
  X.d = data.frame(xd = xd, zd = rep(1, n))
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(dis1))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratum-residual variance comparator, variant 2, in Scenario (ii)
#'
#' Allocate with carat::AdjBCD(), fit the variant-2 outcome regression, and aggregate
#' residual variances within the comparison strata.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
testye.bnz2.c2 <- function(i, n, u, b, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b)
  bs =  model.matrix(~ factor(dis1) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  xd = dis1
  X.d = data.frame(xd = xd, zd = rep(1, n))
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  md.l2 <- lm(Y.f ~ 0 + trt.mx + xd)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(dis1))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the proposed two-sample test in Scenario (iii)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Fit original-data nuisance models once, resample observed donor rows, and restart CAR
#' to estimate the variance of the fitted conditional-mean imbalance.
#'
#' @return By default, a length-2 numeric vector c(adjusted, naive), with 1 for
#'   rejection, 0 otherwise, and possible NA for invalid variances. With details = TRUE,
#'   a list containing rejection, estimate (treatment contrast), variance and
#'   variance_naive (contrast variances), sigma_E2 (conditional residual variance),
#'   sigma_m2 (bootstrap imbalance variance), noise_component, fitted_m (length-n fixed
#'   fitted means), and bootf (length-b imbalance statistics).
test1.ts <- function(i, n, u, b, fun, details = FALSE){
  
  eb_check_bootstrap(b)
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
  
  # Row callback: input i is a participant index; output is the full numeric
  # tensor row formed from the two marginal spline rows, before any later removal.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  # Remove the redundant constant direction in the OUTCOME regression.
  xs = xs[, -1, drop = FALSE]
  g = fun(x1, x2, z)
  
  design = eb_car_design(x, "ts", k = 2)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  Y = treat.mx %*% beta + g + epsi
  
  md.l <- lm(Y ~ 0 + treat.mx)
  hat.t = coef(md.l)
  sig.epsi2 = sum(residuals(md.l)^2) / df.residual(md.l)

  # Fit the outcome and conditional-mean models ONCE, on the original data.
  md = gam(Y ~ 0 + treat + treat2 + xs + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  fz = as.vector(predict(md, type = "terms", terms = "s(z)"))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  # Vgf is the fitted prognostic variance to subtract from the working residual variance.
  vgf = var(hats)
  fzmd = eb_fit_mixed_mean(x, fz)
  efz = fzmd$fitted.values

  # sigma_E^2: original-data conditional residual variance, not bootstrap noise.
  evzx = var(fz - efz)

  bootf = rep(0, b)
  for (j in seq_len(b)) {
    resample = sample.int(n, size = n, replace = TRUE)
    # The same donor row selects X, its FIXED features, and fitted m(X).
    # With two covariates this is joint row resampling, not independent margins.
    treatnew = eb_assign_car(design, index = resample)
    efznew = efz[resample]
    bootf[j] = sum((2 * treatnew - 1) * efznew) / sqrt(n)
  }
  # Sample variance over b independent CAR bootstrap paths estimates sigma_m^2.
  sigma.m2 = var(bootf)
  sige2 = sig.epsi2 - vgf + evzx + sigma.m2

  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sige2 / n
  variance.naive = 4 * sig.epsi2 / n
  estimate = unname(hat.t[1] - hat.t[2])
  prob = eb_reject(estimate, variance)
  prob.n = eb_reject(estimate, variance.naive)
  if (details) return(list(rejection = c(prob, prob.n), estimate = estimate,
                           variance = variance, variance_naive = variance.naive,
                           sigma_E2 = evzx, sigma_m2 = sigma.m2,
                           noise_component = sig.epsi2 - vgf,
                           fitted_m = efz, bootf = bootf))
  return(c(prob, prob.n))
}

#' Simulate the proposed linear working-model test in Scenario (iii)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Fit original-data nuisance models once, resample observed donor rows, and restart CAR
#' to estimate the variance of the fitted conditional-mean imbalance.
#'
#' @return By default, a length-2 numeric vector c(adjusted, naive), with 1 for
#'   rejection, 0 otherwise, and possible NA for invalid variances. With details = TRUE,
#'   a list containing rejection, estimate (treatment contrast), variance and
#'   variance_naive (contrast variances), sigma_E2 (conditional residual variance),
#'   sigma_m2 (bootstrap imbalance variance), noise_component, fitted_m (length-n fixed
#'   fitted means), and bootf (length-b imbalance statistics). The linear-model list also
#'   contains working_residual_variance (SSE/residual df), signal_variance (Vgf),
#'   projection_correction (Q), working_covariates (covariate names),
#'   working_coefficients (fitted slopes), conditional_response (fitted f(Z) minus the
#'   working Z slope times Z), and residual_df (working-model residual degrees of
#'   freedom).
test2.ts <- function(i, n, u, b, fun, details = FALSE){
  
  eb_check_bootstrap(b)
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
  
  # Row callback: input i is a participant index; output is the full numeric
  # tensor row formed from the two marginal spline rows, before any later removal.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  # Remove the redundant constant direction in the OUTCOME regression.
  xs = xs[, -1, drop = FALSE]
  
  g = fun(x1, x2, z)
  
  design = eb_car_design(x, "ts", k = 2)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
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
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = as.vector(fitted(md) - coef(md)["treat"] * treat - coef(md)["treat2"] * treat2)
  # Vgf is the fitted prognostic variance to subtract from the working residual variance.
  vgf = var(hats)
  mnew = cbind(x1, x2, z)
  cm = cov(mnew)
  v = cov(mnew, hats)
  # Q adds back the prognostic variation already removed by the linear working model.
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

  bootf = rep(0, b)
  for (j in seq_len(b)) {
    resample = sample.int(n, size = n, replace = TRUE)
    # The same donor row selects X, its FIXED features, and fitted m(X).
    # With two covariates this is joint row resampling, not independent margins.
    treatnew = eb_assign_car(design, index = resample)
    efznew = efz[resample]
    bootf[j] = sum((2 * treatnew - 1) * efznew) / sqrt(n)
  }
  # Sample variance over b independent CAR bootstrap paths estimates sigma_m^2.
  sigma.m2 = var(bootf)
  # Exact plug-in decomposition specified in the manuscript table.
  noise.component = sig.epsi2 - vgf + e1
  sige2 = noise.component + evzx + sigma.m2

  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sige2 / n
  variance.naive = 4 * sig.epsi2 / n
  estimate = unname(hat.t[1] - hat.t[2])
  prob = eb_reject(estimate, variance)
  prob.n = eb_reject(estimate, variance.naive)
  if (details) return(list(rejection = c(prob, prob.n), estimate = estimate,
                           variance = variance, variance_naive = variance.naive,
                           sigma_E2 = evzx, sigma_m2 = sigma.m2,
                           noise_component = noise.component,
                           working_residual_variance = sig.epsi2,
                           signal_variance = vgf,
                           projection_correction = e1,
                           working_covariates = colnames(mnew),
                           working_coefficients = beta.work,
                           conditional_response = fzminus,
                           residual_df = df.residual(md.l),
                           fitted_m = efz, bootf = bootf))
  return(c(prob, prob.n))
}

#' Simulate the proposed partially linear test in Scenario (iii)
#'
#' Balance X using the fixed discrete, univariate, or tensor design for this scenario.
#' Estimate the treatment contrast directly from the outcome GAM and use its residual
#' mean square SSE/n.
#'
#' @return Scalar numeric rejection indicator (0/1, or NA for invalid variance). With
#'   details = TRUE, a list with rejection, estimate (treatment contrast), variance
#'   (contrast variance), and noise_variance (SSE/n).
test3.ts <- function(i, n, u, b, fun, details = FALSE){
  
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
  
  # Row callback: input i is a participant index; output is the full numeric
  # tensor row formed from the two marginal spline rows, before any later removal.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  # Remove the redundant constant direction in the OUTCOME regression.
  xs = xs[, -1, drop = FALSE]
  
  g = fun(x1, x2, z)
  
  design = eb_car_design(x, "ts", k = 2)
  treat <- eb_assign_car(design)
  treat2 = 1 - treat
  treat.mx <- matrix(c(treat, 1 - treat), ncol = 2)
  beta <- matrix(c(1 + u , 1), ncol = 1)
  
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  Y = treat.mx %*% beta + g + epsi
  
  md = gam(Y ~ 0 + treat + treat2 + xs + s(z))
  if (!isTRUE(md$converged)) stop("Outcome GAM did not converge.")
  hat.t = coef(md)[c("treat", "treat2")]
  sig.epsi2 = sum(residuals(md)^2) / n  # manuscript table
  # Convert the component scale to the variance of the treatment contrast.
  variance = 4 * sig.epsi2 / n
  prob = eb_reject(hat.t[1] - hat.t[2], variance)
  if (details) return(list(rejection = prob, estimate = unname(hat.t[1] - hat.t[2]),
                           variance = variance, noise_variance = sig.epsi2))
  return(prob)
}

#' Simulate the variance-ratio-corrected comparator, variant 1, in Scenario (iii)
#'
#' Allocate with carat::StrBCD() and compute the variant-1 treatment contrast using the
#' tau2.f variance correction.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param bx Integer number of standard-normal quantile bins for each X coordinate and
#'   for the computed Z bin codes; only X bins enter the empirical comparator allocation.
#' @param fun Vectorized prognostic mean function fun(x1, x2, z), returning one g(X) +
#'   f(Z) value per participant.
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest.bnz1.c3 <- function(i, n, u, bx, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n, sd = sqrt(log(1 + 2*x1^2 + 1.5*z^2)))
  
  dis1 = factor(disc2(x1, bx)); dis2 = factor(disc2(x2, bx)); dis3 = factor(disc2(z, bx))
  xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2)
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x1 + x2 + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  strx1 = model.matrix( ~ dis1 - 1); strx2 = model.matrix( ~ dis2 - 1); strz = model.matrix( ~ dis3 - 1)
  varx1 = 0; varx2 = 0; varz = 0
  for (k in 1: ncol(strx1)) {
    varx1 = varx1 + var(x1[strx1[,k] == 1]) * sum(strx1[,k] == 1) / n
  }
  for (k in 1: ncol(strx2)) {
    varx2 = varx2 + var(x2[strx2[,k] == 1]) * sum(strx2[,k] == 1) / n
  }
  for (k in 1: ncol(strz)) {
    varz = varz + var(z[strz[,k] == 1]) * sum(strz[,k] == 1) / n
  }
  
  e2.f = beta2[3]^2 * varx1 + beta2[4]^2 * varx2 + beta2[5]^2 * varz
  e1.f = beta2[3]^2 * var(x1) + beta2[4]^2 * var(x2) + beta2[5]^2 * var(z)
  
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  M = matrix(c(trt, 1 - trt), ncol = 2)
  L = matrix(c(1, -1), ncol = 2)
  # L contrasts the two treatment means; this is its regression design variance factor.
  temp = L %*% solve(t(M) %*% M) %*% t(L)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * temp * tau2.f)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the variance-ratio-corrected comparator, variant 2, in Scenario (iii)
#'
#' Allocate with carat::StrBCD() and compute the variant-2 treatment contrast using the
#' tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest.bnz2.c3 <- function(i, n, u, bx, fun){
  
  eb_require_carat()
  set.seed(i)
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n, sd = sqrt(log(1 + 2*x1^2 + 1.5*z^2)))
  
  dis1 = factor(disc2(x1, bx)); dis2 = factor(disc2(x2, bx)); dis3 = factor(disc2(z, bx))
  xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2)
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x1 + x2 + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx + z)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  strx1 = model.matrix( ~ dis1 - 1); strx2 = model.matrix( ~ dis2 - 1)
  varx1 = 0; varx2 = 0
  for (k in 1: ncol(strx1)) {
    varx1 = varx1 + var(x1[strx1[,k] == 1]) * sum(strx1[,k] == 1) / n
  }
  for (k in 1: ncol(strx2)) {
    varx2 = varx2 + var(x2[strx2[,k] == 1]) * sum(strx2[,k] == 1) / n
  }
  
  e2.f = beta2[3]^2 * varx1 + beta2[4]^2 * varx2
  e1.f = beta2[3]^2 * var(x1) + beta2[4]^2 * var(x2) 
  
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * tau2.f * 4 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratified comparison test, variant 1, in Scenario (iii)
#'
#' Allocate with carat::StrPBR(). Use a treatment-only contrast and the within-/between-
#' stratum variance decomposition.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest1.p.c3 <- function(i, n, u, bx, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = disc2(x1, bx); dis2 = disc2(x2, bx)
  bs = model.matrix(~interaction(dis1, dis2) - 1)
  
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n, sd = sqrt(log(1 + 2*x1^2 + 1.5*z^2)))
  
  xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2)
  
  trt <- ifelse(StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  # Do not discard a populated stratum with an unobserved treatment arm.
  used = colSums(bs) > 0
  arm1 = colSums(bs[trt == 1, , drop = FALSE])
  arm0 = colSums(bs[trt == 0, , drop = FALSE])
  if (any(used & (arm1 == 0 | arm0 == 0))) {
    warning("A populated stratum has an empty treatment arm; returning NA.")
    return(NA_real_)
  }
  
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
  
  # Aggregate within-stratum variances from the two treatment groups.
  sig1.square = 2 * sum(yk1.square, na.rm = T) + 2 * sum(yk0.square, na.rm = T)
  
  ay = rep(0, ncol(bs))
  for (i in 1: ncol(bs)) {
    ay[i] = (yk1[i] - mean(Y.f[trt == 1]) - (yk0[i] - mean(Y.f[trt == 0])))^2 * (sum( bs[,i] == 1) / n)
  }
  
  # Add the between-stratum contribution to the treatment-contrast variance.
  sig2.square = sum(ay, na.rm = T)
  md.l.f <- lm(Y.f ~ trt)
  
  hat.t.f = md.l.f$coefficients[2]
  
  
  s.f <- hat.t.f / sqrt((sig1.square + sig2.square) / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratified comparison test, variant 2, in Scenario (iii)
#'
#' Allocate with carat::StrBCD(). Regress on treatment and stratum indicators, plus Z
#' where included by this scenario, and use the treatment coefficient standard error.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest2.p.c3 <- function(i, n, u, bx, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = disc2(x1, bx); dis2 = disc2(x2, bx)
  bs = model.matrix(~interaction(dis1, dis2) - 1)
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n, sd = sqrt(log(1 + 2*x1^2 + 1.5*z^2)))
  
  xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2)
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  
  md.l.f <- lm(Y.f ~ trt + bs + z)
  hat.t.f = md.l.f$coefficients[2]
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratum-residual variance comparator, variant 1, in Scenario (iii)
#'
#' Allocate with carat::AdjBCD(), fit the variant-1 outcome regression, and aggregate
#' residual variances within the comparison strata.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
testye.bnz1.c3 <- function(i, n, u, bx, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = disc2(x1, bx); dis2 = disc2(x2, bx)
  bs = model.matrix(~interaction(dis1, dis2) - 1)
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n, sd = sqrt(log(1 + 2*x1^2 + 1.5*z^2)))
  
  xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the stratum-residual variance comparator, variant 2, in Scenario (iii)
#'
#' Allocate with carat::AdjBCD(), fit the variant-2 outcome regression, and aggregate
#' residual variances within the comparison strata.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
testye.bnz2.c3 <- function(i, n, u, bx, fun){
  eb_require_carat()
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = disc2(x1, bx); dis2 = disc2(x2, bx)
  bs = model.matrix(~interaction(dis1, dis2) - 1)
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n, sd = sqrt(log(1 + 2*x1^2 + 1.5*z^2)))
  
  xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  md.l2 <- lm(Y.f ~ 0 + trt.mx + xd1 + xd2)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) <= 1){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Run one outer replicate across all scenarios and methods
#'
#' Simulate the three proposed methods and six comparator variants at every effect size,
#' preserving the return order expected by the companion plotting report.
#'
#' @param i Integer seed for one outer simulation replicate; individual tests reset to
#'   this seed.
#' @param n Integer trial sample size; default 1000.
#' @param b Integer empirical bootstrap count for test1/test2, at least 2; default 500.
#'   Comparator bin counts are separately fixed in this function.
#' @param u Numeric vector of treatment effects in saved order; default c(0, 0.1, 0.2,
#'   0.3).
#' @return Numeric vector of length 33 * length(u), concatenated in Scenario (i), (ii),
#'   (iii) order; length 132 at defaults. Values are individual 0/1 rejection decisions
#'   or NA, not averaged power estimates.
#' @details In each 44-entry default scenario block, positions 1:8 alternate
#'   adjusted/naive two-sample decisions, 9:16 alternate adjusted/naive linear-model
#'   decisions, and 17:20 contain partially linear decisions. Comparator blocks are 21:24
#'   (fftest variant 1), 25:28 (fftest variant 2), 29:32 (ftest variant 1), 33:36 (ftest
#'   variant 2), 37:40 (testye variant 1), and 41:44 (testye variant 2). For a different
#'   length(u), each 4-entry block has length(u) entries and each 8-entry block has 2 *
#'   length(u) entries. Calls require the loaded simulation dependencies, including
#'   carat. This function returns results to the caller; the job script is responsible
#'   for saving res1 to a seed file.
delta <- function(i, n = 1000, b = 500, u = c(0, 0.1, 0.2, 0.3)){
  eb_require_carat()
  
  # Scenario 1 callbacks: each sapply function takes scalar effect h and
  # returns the rejection decision(s) from the named trial test. The fun
  # callback takes raw x, z vectors and returns g(X)+f(Z) per row.
  # Column-major flattening keeps effects ordered; two-decision tests alternate
  # adjusted and naive decisions within each effect.
  # case 1
  power.bnz1 = sapply(u, function(h) test1.dis(i, n, h, b = b, fun = function(x, z) 2*x + exp(z)/2.5))
  
  power.bnz2 = sapply(u, function(h) test2.dis(i, n, h, b = b, fun = function(x, z) 2*x + exp(z)/2.5))
  
  power.bnz3 = sapply(u, function(h) test3.dis(i, n, h, b = b, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerf1.p = sapply(u, function(h) ftest1.p.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerf2.p = sapply(u, function(h) ftest2.p.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerye.bnz1 = sapply(u, function(h) testye.bnz1.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerye.bnz2 = sapply(u, function(h) testye.bnz2.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  
  powerff.bnz1 = sapply(u, function(h) fftest.bnz1.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerff.bnz2 = sapply(u, function(h) fftest.bnz2.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  # Scenario 2 callbacks: each sapply function takes scalar effect h and
  # returns the rejection decision(s) from the named trial test. The fun
  # callback takes raw x, z vectors and returns g(X)+f(Z) per row.
  # Column-major flattening keeps effects ordered; two-decision tests alternate
  # adjusted and naive decisions within each effect.
  # case 2
  power.bnz1.ind = sapply(u, function(h) test1.ind(i, n, h, b, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  power.bnz2.ind = sapply(u, function(h) test2.ind(i, n, h, b, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  power.bnz3.ind = sapply(u, function(h) test3.ind(i, n, h, b, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  f1.p.c2 = sapply(u, function(h) fftest.bnz1(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  f2.p.c2 = sapply(u, function(h) fftest.bnz2(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerf1.p.c2 = sapply(u, function(h) ftest1.p.c2(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerf2.p.c2 = sapply(u, function(h) ftest2.p.c2(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powery1.p.c2 = sapply(u, function(h) testye.bnz1.c2(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powery2.p.c2 = sapply(u, function(h) testye.bnz2.c2(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  # Scenario 3 callbacks: each sapply function takes scalar effect h and
  # returns the rejection decision(s) from the named trial test. The fun
  # callback takes raw x1, x2, z vectors and returns g(X)+f(Z) per row.
  # Column-major flattening keeps effects ordered; two-decision tests alternate
  # adjusted and naive decisions within each effect.
  # case 3
  power.bnz1.ts = sapply(u, function(h) test1.ts(i, n, h, b, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  power.bnz2.ts = sapply(u, function(h) test2.ts(i, n, h, b, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  power.bnz3.ts = sapply(u, function(h) test3.ts(i, n, h, b, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  f1.p.c3 = sapply(u, function(h) fftest.bnz1.c3(i, n, h, 10, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  f2.p.c3 = sapply(u, function(h) fftest.bnz2.c3(i, n, h, 10, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerf1.p.c3 = sapply(u, function(h) ftest1.p.c3(i, n, h, 10, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerf2.p.c3 = sapply(u, function(h) ftest2.p.c3(i, n, h, 10, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powery1.p.c3 = sapply(u, function(h) testye.bnz1.c3(i, n, h, 10, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powery2.p.c3 = sapply(u, function(h) testye.bnz2.c3(i, n, h, 10, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  
  return(c(power.bnz1, power.bnz2, power.bnz3, powerff.bnz1, powerff.bnz2, powerf1.p, powerf2.p, powerye.bnz1, powerye.bnz2, 
           power.bnz1.ind, power.bnz2.ind, power.bnz3.ind, f1.p.c2, f2.p.c2, powerf1.p.c2, powerf2.p.c2, powery1.p.c2, powery2.p.c2,
           power.bnz1.ts, power.bnz2.ts, power.bnz3.ts, f1.p.c3, f2.p.c3, powerf1.p.c3, powerf2.p.c3, powery1.p.c3, powery2.p.c3))
  
}
