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

# Corrected all-covariate balancing simulation; analytic variance, no bootstrap.
# Spline order = degree + 1; no bootstrap is used in the testbs* methods.
library(mgcv)
library(MASS)
library(splines2)
library(carat)
# carat is required for comparison methods and delta(), but not testbs*.


#' Convert a variance component to a treatment-effect standard error
#'
#' Use the study scaling Var(treatment contrast) = 4 * sigma / n.
#'
#' @param sigma Positive finite scalar variance component before the 4/n scaling.
#' @param n Positive trial sample size used in the variance scaling.
#' @return Scalar standard error sqrt(4 * sigma / n); NA_real_ with a warning when sigma
#'   is nonpositive or nonfinite.
ball_se <- function(sigma, n) {
  if (length(sigma) != 1L || !is.finite(sigma) || sigma <= 0) {
    warning("Nonpositive or nonfinite adjusted variance; returning NA, not a rejection decision.")
    return(NA_real_)
  }
  sqrt(4 * sigma / n)
}

##############
#randomization
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

#' Assign treatments using local additive B-spline balance in X and Z
#'
#' Identify the knot cell for each participant, concatenate the active marginal splines
#' for X and Z, and update that cell alone.
#'
#' @param x Numeric vector of observed X values in arrival order.
#' @param z Numeric vector of observed Z values, of the same length as x.
#' @param internal_knots_x Interior X knots, or NULL for no interior knot.
#' @param internal_knots_z Interior Z knots, or NULL for no interior knot.
#' @param order_x X spline order, equal to degree + 1; default 3.
#' @param order_z Z spline order, equal to degree + 1; default 3.
#' @param rho Biased-coin probability of reducing local imbalance; default 0.9.
#' @return Numeric vector of length length(x) with treatment labels 0 or 1.
covRcpp_additive_local <- function(x, z, 
                                   internal_knots_x = NULL, 
                                   internal_knots_z = NULL, 
                                   order_x = 3, order_z = 3, 
                                   rho = 0.9) {
  n <- length(x)
  
  # 1. Convert Order to Degree (Degree = Order - 1)
  degree_x <- order_x - 1
  degree_z <- order_z - 1
  
  # 2. Setup Full Knots for findInterval
  # We need the boundary knots (min/max) combined with internal knots
  #' Construct the knot boundaries used to locate a local cell
  #'
  #' Combine the observed endpoints with the sorted interior knots inside
  #' covRcpp_additive_local().
  #'
  #' @param var_data Numeric observations for one marginal covariate.
  #' @param int_knots Interior knots for that covariate, or NULL.
  #' @return Numeric vector consisting of the minimum, sorted interior knots, and
  #'   maximum; just range(var_data) when no interior knots are supplied.
  get_full_knots <- function(var_data, int_knots) {
    if (is.null(int_knots) || length(int_knots) == 0) {
      return(range(var_data))
    } else {
      return(c(min(var_data), sort(int_knots), max(var_data)))
    }
  }
  
  knots_x_full <- get_full_knots(x, internal_knots_x)
  knots_z_full <- get_full_knots(z, internal_knots_z)
  
  # 3. Generate Full B-spline Matrices
  bs_x <- bSpline(x, knots = internal_knots_x, degree = degree_x, intercept = TRUE)
  bs_z <- bSpline(z, knots = internal_knots_z, degree = degree_z, intercept = TRUE)
  
  # 4. Initialization
  T_ls <- numeric(n)
  Lambda_cubics <- new.env() 
  
  # The dimension of the local phi vector is simply the sum of the orders
  local_dim <- order_x + order_z 
  
  # 5. Main Loop
  for (i in 1:n) {
    # A. Identify the interval for X and Z
    jx <- findInterval(x[i], knots_x_full, rightmost.closed = TRUE, all.inside = TRUE)
    jz <- findInterval(z[i], knots_z_full, rightmost.closed = TRUE, all.inside = TRUE)
    
    cubic_id <- paste(jx, jz, sep = "-")
    
    # B. Retrieve local Lambda^i for this cell
    if (exists(cubic_id, envir = Lambda_cubics)) {
      Lambda_j <- get(cubic_id, envir = Lambda_cubics)
    } else {
      Lambda_j <- rep(0, local_dim) 
    }
    
    # C. Extract ACTIVE B-splines only
    # The active indices span from j to j + degree
    active_x <- bs_x[i, jx:(jx + degree_x)]
    active_z <- bs_z[i, jz:(jz + degree_z)]
    
    # Concatenate to form the additive local basis vector phi^i
    phi_i <- c(active_x, active_z)
    
    # D. Calculate Inner Product
    test_val <- sum(Lambda_j * phi_i)
    
    # E. Assignment Decision 
    if (i == 1) {
      T_ls[i] <- rbinom(1, 1, 0.5)
    } else {
      if (test_val < 0) {
        T_ls[i] <- rbinom(1, 1, rho)       
      } else if (test_val > 0) {
        T_ls[i] <- rbinom(1, 1, 1 - rho)   
      } else {
        T_ls[i] <- rbinom(1, 1, 0.5)       
      }
    }
    
    # F. Update Local Lambda
    Lambda_new <- Lambda_j + (2 * T_ls[i] - 1) * phi_i
    assign(cubic_id, Lambda_new, envir = Lambda_cubics)
  }
  
  return(T_ls)
}

# pocsim
#' Choose a sign-dependent allocation probability
#'
#' Map a scalar imbalance sign to the supplied positive- and negative-sign probabilities.
#'
#' @param x Scalar imbalance statistic.
#' @param p Probability returned when x < 0.
#' @param q Probability returned when x > 0.
#' @return Scalar q for positive x, p for negative x, or 0.5 for a tie.
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

############
#discrete case 1

#' Simulate the proposed two-sample test in Scenario (i)
#'
#' Use local additive spline balancing of X and Z. Estimate the difference in group means
#' and remove the fitted prognostic variance from the residual mean square.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param o2 Order of the Z balancing splines, equal to degree + 1.
#' @param knots2 Interior Z knots, or NULL.
#' @param fun Vectorized prognostic mean function fun(x, z), returning one g(X) + f(Z)
#'   value per participant.
#' @param details Logical; FALSE returns rejection decisions only, TRUE returns the
#'   diagnostic list described below.
#' @return By default, a length-2 numeric vector c(adjusted, naive). Decisions are 1 for
#'   rejection, 0 otherwise, and may be NA for invalid variances. With details = TRUE, a
#'   list with reject, estimate (treatment contrast), sigma (component before 4/n
#'   scaling), outcome and treatment (length-n vectors), nuisance_fit, and working_fit
#'   (fitted model objects).
testbs1.dis = function(i, n, u, o2, knots2, fun, details = FALSE){
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
  
  
  # treatment assignment
  treat = covRcpp_additive_local(x, z, internal_knots_x = 0, internal_knots_z = knots2, order_x = 3, order_z = o2)
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  y = treat.mx %*% beta + g + rnorm(n)
  
  # t-test
  md.l = lm(y ~ 0 + treat.mx)
  hat.t = md.l$coefficients
  sig.t = sum((md.l$residuals)^2) / df.residual(md.l)
  
  # glm estimate
  md.g = gam(y ~ 0 + treat.mx + x + s(z))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = md.g$fitted.values - treat.mx %*% matrix(md.g$coefficients[c("treat.mx1", "treat.mx2")], ncol = 1)
  # sig.g estimates Var{g(X)+f(Z)} after removing the treatment contribution.
  sig.g = var(hats)
  
  sig = sig.t - sig.g
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  s.n = (hat.t[1] - hat.t[2]) / ball_se(sig.t, n)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  
  if (details) return(list(reject = c(prob, prob.n), estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.l))
  return(c(prob, prob.n))
}

#' Simulate the proposed linear working-model test in Scenario (i)
#'
#' Use local additive spline balancing of X and Z. Fit a linear outcome regression on raw
#' X and Z, remove the fitted prognostic variance, and add back its linear covariance
#' projection.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param o2 Order of the Z balancing splines, equal to degree + 1.
#' @param knots2 Interior Z knots, or NULL.
#' @param fun Vectorized prognostic mean function fun(x, z), returning one g(X) + f(Z)
#'   value per participant.
#' @param details Logical; FALSE returns rejection decisions only, TRUE returns the
#'   diagnostic list described below.
#' @return By default, a length-2 numeric vector c(adjusted, naive). Decisions are 1 for
#'   rejection, 0 otherwise, and may be NA for invalid variances. With details = TRUE, a
#'   list with reject, estimate (treatment contrast), sigma (component before 4/n
#'   scaling), outcome and treatment (length-n vectors), nuisance_fit, and working_fit
#'   (fitted model objects).
testbs2.dis = function(i, n, u, o2, knots2, fun, details = FALSE){
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
  
  
  # treatment assignment
  treat = covRcpp_additive_local(x, z, internal_knots_x = 0, internal_knots_z = knots2, order_x = 3, order_z = o2)
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  y = treat.mx %*% beta + g + rnorm(n)
  
  # t-test
  md.l = lm(y ~ 0 + treat.mx + x + z)
  hat.t = md.l$coefficients
  sig.t = sum((md.l$residuals)^2) / df.residual(md.l)
  
  # glm estimate
  md.g = gam(y ~ 0 + treat.mx + x + s(z))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = md.g$fitted.values - treat.mx %*% matrix(md.g$coefficients[c("treat.mx1", "treat.mx2")], ncol = 1)
  # sig.g estimates Var{g(X)+f(Z)} after removing the treatment contribution.
  sig.g = var(hats)
  
  vv = cbind(x, z)
  v = cov(vv, as.vector(hats))
  # Q = v^T Cov(U)^(-1) v is the linear projection contribution for raw U.
  sig.c = t(v) %*% solve(cov(vv)) %*% v
  
  sig = sig.t - sig.g + sig.c
  
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  s.n = (hat.t[1] - hat.t[2]) / ball_se(sig.t, n)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  
  if (details) return(list(reject = c(prob, prob.n), estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.l))
  return(c(prob, prob.n))
}

#' Simulate the proposed partially linear test in Scenario (i)
#'
#' Use local additive spline balancing of X and Z. Estimate the contrast from the outcome
#' GAM and use the mean squared GAM residuals.
#' @return By default, one numeric rejection indicator. Decisions are 1 for rejection, 0
#'   otherwise, and may be NA for invalid variances. With details = TRUE, a list with
#'   reject, estimate (treatment contrast), sigma (component before 4/n scaling), outcome
#'   and treatment (length-n vectors), nuisance_fit, and working_fit (fitted model
#'   objects).
testbs3.dis = function(i, n, u, o2, knots2, fun, details = FALSE){
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
  
  # treatment assignment
  treat = covRcpp_additive_local(x, z, internal_knots_x = 0, internal_knots_z = knots2, order_x = 3, order_z = o2)
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  y = treat.mx %*% beta + g + rnorm(n)
  
  # glm test
  md.g = gam(y ~ 0 + treat.mx + x + s(z))
  hat.t = md.g$coefficients[c("treat.mx1", "treat.mx2")]
  sig = mean((md.g$residuals)^2)
  
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  if (details) return(list(reject = prob, estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.g))
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
ftest1.dis <- function(i, n, u, fun){
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
  
  bs = matrix(rep(0, 8 * n), ncol = 8)
  for (i in 1: n) {
    if(x[i] == -2 & z[i] > 0){
      bs[i, 1] = 1
    }
    if(x[i] == -1 & z[i] > 0){
      bs[i, 2] = 1
    }
    if(x[i] == 1 & z[i] > 0){
      bs[i, 3] = 1
    }
    if(x[i] == 2 & z[i] > 0){
      bs[i, 4] = 1
    }
    if(x[i] == -2 & z[i] < 0){
      bs[i, 5] = 1
    }
    if(x[i] == -1 & z[i] < 0){
      bs[i, 6] = 1
    }
    if(x[i] == 1 & z[i] < 0){
      bs[i, 7] = 1
    }
    if(x[i] == 2 & z[i] < 0){
      bs[i, 8] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0)
  X.d = data.frame(xd = x, zd = zd)
  
  trt <- ifelse(StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  yk1 = rep(0, 8)
  for (i in 1: 8) {
    yk1[i] = mean(Y.f[trt == 1 & bs[,i] == 1])
  }
  
  yk0 = rep(0, 8)
  for (i in 1: 8) {
    yk0[i] = mean(Y.f[trt == 0 & bs[,i] == 1])
  }
  
  yk1.square = rep(0, 8)
  for (i in 1: 8) {
    yk1.square[i] = sum((Y.f[trt == 1 & bs[,i] == 1] - yk1[i])^2) * (sum( bs[,i] == 1) / n) / sum(trt == 1 & bs[,i] == 1)
  }
  
  yk0.square = rep(0, 8)
  for (i in 1: 8) {
    yk0.square[i] = sum((Y.f[trt == 0 & bs[,i] == 1] - yk0[i])^2) * (sum( bs[,i] == 1) / n) / sum(trt == 0 & bs[,i] == 1)
  }
  
  # Aggregate within-stratum variances from the two treatment groups.
  sig1.square = 2 * sum(yk1.square, na.rm = T) + 2 * sum(yk0.square, na.rm = T)
  
  ay = rep(0, 8)
  for (i in 1: 8) {
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
ftest2.dis <- function(i, n, u, fun){
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
  
  bs = matrix(rep(0, 7 * n), ncol = 7)
  for (i in 1: n) {
    if(x[i] == -2 & z[i] > 0){
      bs[i, 1] = 1
    }
    if(x[i] == -1 & z[i] > 0){
      bs[i, 2] = 1
    }
    if(x[i] == 1 & z[i] > 0){
      bs[i, 3] = 1
    }
    if(x[i] == 2 & z[i] > 0){
      bs[i, 4] = 1
    }
    if(x[i] == -2 & z[i] < 0){
      bs[i, 5] = 1
    }
    if(x[i] == -1 & z[i] < 0){
      bs[i, 6] = 1
    }
    if(x[i] == 1 & z[i] < 0){
      bs[i, 7] = 1
    }
  }
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  zd = as.numeric(z > 0)
  X.d = data.frame(xd = x, zd = zd)
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  
  md.l.f <- lm(Y.f ~ trt + bs)
  hat.t.f = md.l.f$coefficients[2]
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
} 

#' Simulate the variance-ratio-corrected comparator, variant 1, in Scenario (i)
#'
#' Allocate with carat::PocSimMIN() and compute the variant-1 treatment contrast using
#' the tau2.f variance correction.
#'
fftest1.dis <- function(i, n, u, fun){
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
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5], zd = zd)
  
  trt <- ifelse(PocSimMIN(X.d)$assignments == 'A', 1, 0)
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
  
  # Legacy comparator formula: assumes unit error variance and omits cross-covariances.
  # Retained for method identity; not validated for heteroscedastic/nonlinear scenarios.
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
#' Allocate with carat::PocSimMIN() and compute the variant-2 treatment contrast using
#' the tau2.f variance correction.
#'
fftest2.dis <- function(i, n, u, fun){
  
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
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5], zd = zd)
  
  trt <- ifelse(PocSimMIN(X.d)$assignments == 'A', 1, 0)
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
  
  # Legacy comparator formula: assumes unit error variance and omits cross-covariances.
  # Retained for method identity; not validated for heteroscedastic/nonlinear scenarios.
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
testye1.dis <- function(i, n, u, fun){
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
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5], zd = zd)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2, xd3, xd4, xd5, zd))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) == 1){
      warning("Within-stratum variance is undefined for a singleton; returning NA.")
      return(NA_real_)
    }
    else if(length(stratas[[i]]) == 0){ sig2 = sig2 }
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

testye2.dis <- function(i, n, u, fun){
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
  X.d = data.frame(xd1 = bs[,1], xd2 = bs[,2], xd3 = bs[,3], xd4 = bs[, 4], xd5 = bs[,5], zd = zd)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + zd)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2, xd3, xd4, xd5, zd))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) == 1){
      warning("Within-stratum variance is undefined for a singleton; returning NA.")
      return(NA_real_)
    }
    else if(length(stratas[[i]]) == 0){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
} 

############
# ind case 2
#' Simulate the proposed two-sample test in Scenario (ii)
#'
#' Use local additive spline balancing of X and Z. Estimate the difference in group means
#' and remove the fitted prognostic variance from the residual mean square.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param o1 Order of the X balancing splines, equal to degree + 1.
#' @param o2 Order of the Z balancing splines, equal to degree + 1.
#' @param knots1 Interior X knots (used in both X coordinates for Scenario (iii)), or
#'   NULL.
#' @param knots2 Interior Z knots, or NULL.
#' @param fun Vectorized prognostic mean function fun(x, z), returning one g(X) + f(Z)
#'   value per participant.
#' @param details Logical; FALSE returns rejection decisions only, TRUE returns the
#'   diagnostic list described below.
#' @return By default, a length-2 numeric vector c(adjusted, naive). Decisions are 1 for
#'   rejection, 0 otherwise, and may be NA for invalid variances. With details = TRUE, a
#'   list with reject, estimate (treatment contrast), sigma (component before 4/n
#'   scaling), outcome and treatment (length-n vectors), nuisance_fit, and working_fit
#'   (fitted model objects).
testbs1.ind = function(i, n, u, o1, o2, knots1, knots2, fun, details = FALSE){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  
  # treatment assignment
  treat = covRcpp_additive_local(x, z, knots1, knots2, o1, o2)
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  y = treat.mx %*% beta + g + rnorm(n)
  
  # t-test
  md.l = lm(y ~ 0 + treat.mx)
  hat.t = md.l$coefficients
  sig.t = sum((md.l$residuals)^2) / df.residual(md.l)
  
  # glm estimate
  md.g = gam(y ~ 0 + treat.mx + s(x) + s(z))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = md.g$fitted.values - treat.mx %*% matrix(md.g$coefficients[c("treat.mx1", "treat.mx2")], ncol = 1)
  # sig.g estimates Var{g(X)+f(Z)} after removing the treatment contribution.
  sig.g = var(hats)
  
  sig = sig.t - sig.g
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  s.n = (hat.t[1] - hat.t[2]) / ball_se(sig.t, n)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  
  if (details) return(list(reject = c(prob, prob.n), estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.l))
  return(c(prob, prob.n))
}

#' Simulate the proposed linear working-model test in Scenario (ii)
#'
#' Use local additive spline balancing of X and Z. Fit a linear outcome regression on raw
#' X and Z, remove the fitted prognostic variance, and add back its linear covariance
#' projection.
#'

testbs2.ind = function(i, n, u, o1, o2, knots1, knots2, fun, details = FALSE){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  
  # treatment assignment
  treat = covRcpp_additive_local(x, z, knots1, knots2, o1, o2)
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  y = treat.mx %*% beta + g + rnorm(n)
  
  # t-test
  md.l = lm(y ~ 0 + treat.mx + x + z)
  hat.t = md.l$coefficients
  sig.t = sum((md.l$residuals)^2) / df.residual(md.l)
  
  # glm estimate
  md.g = gam(y ~ 0 + treat.mx + s(x) + s(z))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = md.g$fitted.values - treat.mx %*% matrix(md.g$coefficients[c("treat.mx1", "treat.mx2")], ncol = 1)
  # sig.g estimates Var{g(X)+f(Z)} after removing the treatment contribution.
  sig.g = var(hats)
  
  vv = cbind(x, z)
  v = cov(vv, as.vector(hats))
  # Q = v^T Cov(U)^(-1) v is the linear projection contribution for raw U.
  sig.c = t(v) %*% solve(cov(vv)) %*% v
  
  sig = sig.t - sig.g + sig.c
  
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  s.n = (hat.t[1] - hat.t[2]) / ball_se(sig.t, n)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  
  if (details) return(list(reject = c(prob, prob.n), estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.l))
  return(c(prob, prob.n))
}

#' Simulate the proposed partially linear test in Scenario (ii)
#'
#' Use local additive spline balancing of X and Z. Estimate the contrast from the outcome
#' GAM and use the mean squared GAM residuals.
#'

testbs3.ind = function(i, n, u, o1, o2, knots1, knots2, fun, details = FALSE){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  
  # treatment assignment
  treat = covRcpp_additive_local(x, z, knots1, knots2, o1, o2)
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  y = treat.mx %*% beta + g + rnorm(n)
  
  # glm test
  md.g = gam(y ~ 0 + treat.mx + s(x) + s(z))
  hat.t = md.g$coefficients[c("treat.mx1", "treat.mx2")]
  sig = mean((md.g$residuals)^2)
  
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  if (details) return(list(reject = prob, estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.g))
  return(prob)
}

#' Simulate the stratified comparison test, variant 1, in Scenario (ii)
#'
#' Allocate with carat::StrPBR(). Use a treatment-only contrast and the within-/between-
#' stratum variance decomposition.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param b Integer number of standard-normal quantile bins for continuous covariates in
#'   the comparison method; this is not a bootstrap count.
#' @param fun Vectorized prognostic mean function fun(x, z), returning one g(X) + f(Z)
#'   value per participant.
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest1.mf.ind <- function(i, n, u, b, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b); dis2 <- disc2(z, b)
  bs =  model.matrix( ~ interaction(dis1, dis2) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  X.d = data.frame(xd = dis1, zd = dis2)
  
  trt <- ifelse(StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  occupied = colSums(bs) > 0
  if (any(occupied & (colSums(bs * trt) == 0 | colSums(bs * (1-trt)) == 0))) {
    warning("A populated stratum has no observations in one arm; returning NA.")
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

fftest2.mf.ind <- function(i, n, u, b, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b); dis2 <- disc2(z, b)
  bs =  model.matrix( ~ interaction(dis1, dis2) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  X.d = data.frame(xd = dis1, zd = dis2)
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  
  md.l.f <- lm(Y.f ~ trt + bs)
  hat.t.f = md.l.f$coefficients[2]
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

#' Simulate the variance-ratio-corrected comparator, variant 1, in Scenario (ii)
#'
#' Allocate with carat::PocSimMIN() and compute the variant-1 treatment contrast using
#' the tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest1.ind <- function(i, n, u, b, fun){
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
  
  dis1 <- disc2(x, b); dis2 <- disc2(z, b)
  bs =  model.matrix( ~ interaction(dis1, dis2) - 1)
  
  xd = dis1; zd = dis2
  X.d = data.frame(xd = xd, zd = zd)
  
  trt <- ifelse(PocSimMIN(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx + x + z)
  beta2 <- md.l2$coefficients
  
  md.l.f <- lm(Y.f ~ 0 + trt.mx)
  hat.t.f = md.l.f$coefficients
  sig.epsi2 <- sum(md.l.f$residuals^2) / df.residual(md.l.f)
  
  
  strx = model.matrix( ~ factor(dis1) - 1); strz = model.matrix( ~ factor(dis2) - 1)
  
  varx = 0; varz = 0
  for (k in 1: ncol(strx)) {
    varx = varx + var(x[strx[,k] == 1]) * sum(strx[,k] == 1) / n
  }
  for (k in 1: ncol(strz)) {
    varz = varz + var(z[strz[,k] == 1]) * sum(strz[,k] == 1) / n
  }
  
  e2.f = beta2[3]^2 * varx + beta2[4]^2 * varz
  e1.f = beta2[3]^2 * var(x) + beta2[4]^2 * var(z)
  
  # Legacy comparator formula: assumes unit error variance and omits cross-covariances.
  # Retained for method identity; not validated for heteroscedastic/nonlinear scenarios.
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
#' Allocate with carat::PocSimMIN() and compute the variant-2 treatment contrast using
#' the tau2.f variance correction.
#'

fftest2.ind <- function(i, n, u, b, fun){
  
  set.seed(i)
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  xd <- disc2(x, b); zd <- disc2(z, b)
  X.d = data.frame(xd = xd, zd = zd)
  
  trt <- ifelse(PocSimMIN(X.d)$assignments == 'A', 1, 0)
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
  
  # Legacy comparator formula: assumes unit error variance and omits cross-covariances.
  # Retained for method identity; not validated for heteroscedastic/nonlinear scenarios.
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * tau2.f * 4 / n)
  
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
testye1.ind <- function(i, n, u, b, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b); dis2 <- disc2(z, b)
  bs =  model.matrix( ~ interaction(dis1, dis2) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  X.d = data.frame(xd1 = dis1, zd = dis2)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(dis1, dis2))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) == 1){
      warning("Within-stratum variance is undefined for a singleton; returning NA.")
      return(NA_real_)
    }
    else if(length(stratas[[i]]) == 0){ sig2 = sig2 }
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
testye2.ind <- function(i, n, u, b, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0), matrix(c(1, 0, 0, 1), ncol = 2))
  x = s2_conditional_x(data[,1])  # N(0,1) conditional on [-2,2]
  z = data[,2]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 <- disc2(x, b); dis2 <- disc2(z, b)
  bs =  model.matrix( ~ interaction(dis1, dis2) - 1)
  
  g = fun(x, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  epsi = rnorm(n)
  
  X.d = data.frame(xd = dis1, zd = dis2)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  md.l2 <- lm(Y.f ~ 0 + trt.mx + dis1 + dis2)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(dis1, dis2))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) == 1){
      warning("Within-stratum variance is undefined for a singleton; returning NA.")
      return(NA_real_)
    }
    else if(length(stratas[[i]]) == 0){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
}

########
# tensor case 3

#' Simulate the proposed two-sample test in Scenario (iii)
#'
#' Use global balancing of the full X tensor and Z splines, with the weighted intercept
#' in covRcpp11(). Estimate the difference in group means and remove the fitted
#' prognostic variance from the residual mean square.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param d1 Degree for both X-coordinate marginal splines; order is d1 + 1.
#' @param d2 Degree for the Z balancing splines; order is d2 + 1.
#' @param knots1 Interior X knots (used in both X coordinates for Scenario (iii)), or
#'   NULL.
#' @param knots2 Interior Z knots, or NULL.
#' @param fun Vectorized prognostic mean function fun(x1, x2, z), returning one g(X) +
#'   f(Z) value per participant.
#' @param details Logical; FALSE returns rejection decisions only, TRUE returns the
#'   diagnostic list described below.
#' @return By default, a length-2 numeric vector c(adjusted, naive). Decisions are 1 for
#'   rejection, 0 otherwise, and may be NA for invalid variances. With details = TRUE, a
#'   list with reject, estimate (treatment contrast), sigma (component before 4/n
#'   scaling), outcome and treatment (length-n vectors), nuisance_fit, and working_fit
#'   (fitted model objects).
testbs1.ts = function(i, n, u, d1, d2, knots1, knots2, fun, details = FALSE){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  # x below uses the clipped covariates
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
  xs1 = bSpline(x1, degree = d1, knots = knots1, Boundary.knots = c(-1, 1), intercept = TRUE)
  xs2 = bSpline(x2, degree = d1, knots = knots1, Boundary.knots = c(-1, 1), intercept = TRUE)
  # Row callback: input i is a participant index; output is the full numeric
  # tensor row formed from the two marginal spline rows, before any later removal.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  bs = bSpline(z, degree = d2, knots = knots2, Boundary.knots = c(-3, 3), intercept = TRUE)
  bs = cbind(bs, xs)
  
  # treatment assignment
  treat = covRcpp11(bs)
  # Remove one FULL tensor column only for the outcome regression.
  # Allocation above retains the original full-basis weighting.
  xs = xs[, -1, drop = FALSE]
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x1, x2, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  y = treat.mx %*% beta + g + epsi
  
  # t-test
  md.l = lm(y ~ 0 + treat.mx)
  hat.t = md.l$coefficients
  sig.t = sum((md.l$residuals)^2) / df.residual(md.l)
  
  # glm estimate
  md.g = gam(y ~ 0 + treat.mx + xs + s(z))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = md.g$fitted.values - treat.mx %*% matrix(md.g$coefficients[c("treat.mx1", "treat.mx2")], ncol = 1)
  # sig.g estimates Var{g(X)+f(Z)} after removing the treatment contribution.
  sig.g = var(hats)
  
  sig = sig.t - sig.g
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  s.n = (hat.t[1] - hat.t[2]) / ball_se(sig.t, n)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  
  if (details) return(list(reject = c(prob, prob.n), estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.l))
  return(c(prob, prob.n))
}

#' Simulate the proposed linear working-model test in Scenario (iii)
#'
#' Use global balancing of the full X tensor and Z splines, with the weighted intercept
#' in covRcpp11(). Fit a linear outcome regression on raw X and Z, remove the fitted
#' prognostic variance, and add back its linear covariance projection.
#'

testbs2.ts = function(i, n, u, d1, d2, knots1, knots2, fun, details = FALSE){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  # x below uses the clipped covariates
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
  xs1 = bSpline(x1, degree = d1, knots = knots1, Boundary.knots = c(-1, 1), intercept = TRUE)
  xs2 = bSpline(x2, degree = d1, knots = knots1, Boundary.knots = c(-1, 1), intercept = TRUE)
  # Row callback: input i is a participant index; output is the full numeric
  # tensor row formed from the two marginal spline rows, before any later removal.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  bs = bSpline(z, degree = d2, knots = knots2, Boundary.knots = c(-3, 3), intercept = TRUE)
  bs = cbind(bs, xs)
  
  # treatment assignment
  treat = covRcpp11(bs)
  # Remove one FULL tensor column only for the outcome regression.
  # Allocation above retains the original full-basis weighting.
  xs = xs[, -1, drop = FALSE]
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x1, x2, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  y = treat.mx %*% beta + g + epsi
  
  # t-test
  md.l = lm(y ~ 0 + treat.mx + x1 + x2 + z)
  hat.t = md.l$coefficients
  sig.t = sum((md.l$residuals)^2) / df.residual(md.l)
  
  # glm estimate
  md.g = gam(y ~ 0 + treat.mx + xs + s(z))
  # Remove both fitted treatment means to isolate the prognostic signal g(X)+f(Z).
  hats = md.g$fitted.values - treat.mx %*% matrix(md.g$coefficients[c("treat.mx1", "treat.mx2")], ncol = 1)
  # sig.g estimates Var{g(X)+f(Z)} after removing the treatment contribution.
  sig.g = var(hats)
  
  vv = cbind(x, z)
  v = cov(vv, as.vector(hats))
  # Q = v^T Cov(U)^(-1) v is the linear projection contribution for raw U.
  sig.c = t(v) %*% solve(cov(vv)) %*% v
  
  sig = sig.t - sig.g + sig.c
  
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  s.n = (hat.t[1] - hat.t[2]) / ball_se(sig.t, n)
  prob.n = ifelse(abs(s.n) < qnorm(0.975), 0, 1)
  
  
  if (details) return(list(reject = c(prob, prob.n), estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.l))
  return(c(prob, prob.n))
}

#' Simulate the proposed partially linear test in Scenario (iii)
#'
#' Use global balancing of the full X tensor and Z splines, with the weighted intercept
#' in covRcpp11(). Estimate the contrast from the outcome GAM and use the mean squared
#' GAM residuals.
#'

testbs3.ts = function(i, n, u, d1, d2, knots1, knots2, fun, details = FALSE){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  # x below uses the clipped covariates
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
  xs1 = bSpline(x1, degree = d1, knots = knots1, Boundary.knots = c(-1, 1), intercept = TRUE)
  xs2 = bSpline(x2, degree = d1, knots = knots1, Boundary.knots = c(-1, 1), intercept = TRUE)
  # Row callback: input i is a participant index; output is the full numeric
  # tensor row formed from the two marginal spline rows, before any later removal.
  xs = t(sapply(1: n, function(i) kronecker(xs1[i,], xs2[i,])))
  bs = bSpline(z, degree = d2, knots = knots2, Boundary.knots = c(-3, 3), intercept = TRUE)
  bs = cbind(bs, xs)
  
  # treatment assignment
  treat = covRcpp11(bs)
  # Remove one FULL tensor column only for the outcome regression.
  # Allocation above retains the original full-basis weighting.
  xs = xs[, -1, drop = FALSE]
  treat.mx = matrix(c(treat, 1 - treat), ncol = 2)
  # imbalance check
  g = fun(x1, x2, z)
  
  # Internal diagnostic: scaled signed imbalance in the supplied prognostic mean.
  imb = sum((2 * treat - 1) * g) / sqrt(n)
  
  beta = matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  y = treat.mx %*% beta + g + epsi
  
  # glm test
  md.g = gam(y ~ 0 + treat.mx + xs + s(z))
  hat.t = md.g$coefficients[c("treat.mx1", "treat.mx2")]
  sig = mean((md.g$residuals)^2)
  
  s = (hat.t[1] - hat.t[2]) / ball_se(sig, n)
  prob = ifelse(abs(s) < qnorm(0.975), 0, 1)
  
  if (details) return(list(reject = prob, estimate = unname(hat.t[1] - hat.t[2]), sigma = as.numeric(sig),
                          outcome = as.vector(y), treatment = treat, nuisance_fit = md.g,
                          working_fit = md.g))
  return(prob)
}


#' Simulate the stratified comparison test, variant 1, in Scenario (iii)
#'
#' Allocate with carat::StrPBR(). Use a treatment-only contrast and the within-/between-
#' stratum variance decomposition.
#'
#' @param i Integer simulation seed passed to set.seed(); every call resets the RNG.
#' @param n Integer number of participants in this simulated trial.
#' @param u Scalar treatment effect: treatment 1 has mean 1 + u and treatment 0 has mean
#'   1.
#' @param bx Integer number of standard-normal quantile bins for each X coordinate.
#' @param bz Integer number of standard-normal quantile bins for Z.
#' @param fun Vectorized prognostic mean function fun(x1, x2, z), returning one g(X) +
#'   f(Z) value per participant.
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
ftest1.c3 <- function(i, n, u, bx, bz, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x = data[,1:2]
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = factor(disc2(x1, bx)); dis2 = factor(disc2(x2, bx)); dis3 = factor(disc2(z, bz))
  bs = model.matrix(~interaction(dis1, dis2, dis3) - 1)
  
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  zd = dis3; xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2, zd = zd)
  
  trt <- ifelse(StrPBR(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  occupied = colSums(bs) > 0
  if (any(occupied & (colSums(bs * trt) == 0 | colSums(bs * (1-trt)) == 0))) {
    warning("A populated stratum has no observations in one arm; returning NA.")
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

ftest2.c3 <- function(i, n, u, bx, bz, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x = data[,1:2]
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = factor(disc2(x1, bx)); dis2 = factor(disc2(x2, bx)); dis3 = factor(disc2(z, bz))
  bs = model.matrix(~interaction(dis1, dis2, dis3) - 1)
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  zd = dis3; xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2, zd = zd)
  
  trt <- ifelse(StrBCD(X.d)$assignments == 'A', 1, 0)
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  
  md.l.f <- lm(Y.f ~ trt + bs)
  hat.t.f = md.l.f$coefficients[2]
  sig.epsi2 <- diag(vcov(md.l.f))[2]
  
  
  s.f <- hat.t.f / sqrt(sig.epsi2)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
} 

#' Simulate the variance-ratio-corrected comparator, variant 1, in Scenario (iii)
#'
#' Allocate with carat::PocSimMIN() and compute the variant-1 treatment contrast using
#' the tau2.f variance correction.
#'
#' @return Scalar numeric indicator: 1 if the two-sided 0.05 normal test rejects, 0
#'   otherwise. Existing checks or degenerate strata/variance calculations can produce
#'   NA.
fftest1.ts <- function(i, n, u, bx, bz, fun){
  set.seed(i)
  
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x = data[,1:2]
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
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  dis1 = factor(disc2(x1, bx)); dis2 = factor(disc2(x2, bx)); dis3 = factor(disc2(z, bz))
  bs = model.matrix(~interaction(dis1, dis2, dis3) - 1)
  xd1 = dis1; xd2 = dis2; zd = dis3
  X.d = data.frame(xd1 = xd1, xd2 = xd2, zd = zd)
  
  
  trt <- ifelse(PocSimMIN(X.d)$assignments == 'A', 1, 0)
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
  
  # Legacy comparator formula: assumes unit error variance and omits cross-covariances.
  # Retained for method identity; not validated for heteroscedastic/nonlinear scenarios.
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
#' Allocate with carat::PocSimMIN() and compute the variant-2 treatment contrast using
#' the tau2.f variance correction.
#'

fftest2.ts <- function(i, n, u, bx, bz, fun){
  
  set.seed(i)
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x = data[,1:2]
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
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  dis1 = factor(disc2(x1, bx)); dis2 = factor(disc2(x2, bx)); dis3 = factor(disc2(z, bz))
  xd1 = dis1; xd2 = dis2; zd = dis3
  X.d = data.frame(xd1 = xd1, xd2 = xd2, zd = zd)
  
  trt <- ifelse(PocSimMIN(X.d)$assignments == 'A', 1, 0)
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
  
  # Legacy comparator formula: assumes unit error variance and omits cross-covariances.
  # Retained for method identity; not validated for heteroscedastic/nonlinear scenarios.
  # Apply the comparator ratio of remaining to total linear prognostic variation.
  tau2.f <- (1 + e2.f) / (1 + e1.f)
  
  s.f <- (hat.t.f[1] - hat.t.f[2]) / sqrt(sig.epsi2 * tau2.f * 4 / n)
  
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
testye1.ts <- function(i, n, u, bx, bz, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x = data[,1:2]
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = disc2(x1, bx); dis2 = disc2(x2, bx); dis3 = disc2(z, bz)
  bs = model.matrix(~interaction(dis1, dis2, dis3) - 1)
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  zd = dis3; xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2, zd = zd)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  
  md.l2 <- lm(Y.f ~ 0 + trt.mx)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2, zd))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) == 1){
      warning("Within-stratum variance is undefined for a singleton; returning NA.")
      return(NA_real_)
    }
    else if(length(stratas[[i]]) == 0){ sig2 = sig2 }
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
testye2.ts <- function(i, n, u, bx, bz, fun){
  set.seed(i)
  # data generation
  data = mvrnorm(n, c(0, 0, 0), matrix(c(1, 0.5, 0.7, 0.5, 1, 0, 0.7, 0, 1), ncol = 3))
  x = data[,1:2]
  x1 = data[,1]
  x2 = data[,2]
  x1[x1 < -1] = -1
  x2[x2 < -1] = -1
  x1[x1 > 1] = 1
  x2[x2 > 1] = 1
  
  z = data[,3]
  z[z < -3] = -3
  z[z > 3] = 3
  
  dis1 = disc2(x1, bx); dis2 = disc2(x2, bx); dis3 = disc2(z, bz)
  bs = model.matrix(~interaction(dis1, dis2, dis3) - 1)
  
  g = fun(x1, x2, z)
  
  beta <- matrix(c(1 + u, 1), ncol = 1)
  # Noise callback: input i is a participant index; output is one zero-mean
  # normal draw with the Scenario (iii) conditional variance for that row.
  epsi = sapply(1: n, function(i) rnorm(1, 0, sqrt(log(1 + 2*x1[i]^2 + 1.5*z[i]^2))))
  
  zd = dis3; xd1 = dis1; xd2 = dis2
  X.d = data.frame(xd1 = xd1, xd2 = xd2, zd = zd)
  
  trt = ifelse(AdjBCD(X.d)$assignments == 'A', 1, 0)
  
  trt.mx <- matrix(c(trt, 1 - trt), ncol = 2)
  
  Y.f = trt.mx %*% beta + g + epsi
  md.l2 <- lm(Y.f ~ 0 + trt.mx + xd1 + xd2 + zd)
  beta.hat <- md.l2$coefficients
  
  stratas = split(Y.f, list(xd1, xd2, zd))
  
  sig2 = 0
  
  for (i in 1: length(stratas)) {
    if(length(stratas[[i]]) == 1){
      warning("Within-stratum variance is undefined for a singleton; returning NA.")
      return(NA_real_)
    }
    else if(length(stratas[[i]]) == 0){ sig2 = sig2 }
    else {sig2 = sig2 + var(stratas[[i]]) * length(stratas[[i]])}
  }
  
  # Normalize the size-weighted residual variance by the full trial sample size.
  sig2 = sig2 / n
  
  s.f <- (beta.hat[1] - beta.hat[2]) / sqrt(4 * sig2 / n)
  
  prob.f = ifelse(abs(s.f) < qnorm(0.975), 0, 1)
  
  return(prob.f)
} 

##########

#' Run one outer replicate across all scenarios and methods
#'
#' Simulate the three proposed methods and six comparator variants at every effect size,
#' preserving the return order expected by the companion plotting report.
#'
#' @param i Integer seed for one outer simulation replicate; individual tests reset to
#'   this seed.
#' @param n Integer trial sample size; default 1000.
#' @return Length-132 numeric vector of individual 0/1 rejection decisions or NA,
#'   concatenated in Scenario (i), (ii), (iii) order. Effects are fixed at c(0, 0.1, 0.2,
#'   0.3); the output is not an averaged power estimate.
#' @details In each 44-entry default scenario block, positions 1:8 alternate
#'   adjusted/naive two-sample decisions, 9:16 alternate adjusted/naive linear-model
#'   decisions, and 17:20 contain partially linear decisions. In Scenarios (i)/(iii),
#'   comparator blocks are 21:24 (ftest variant 1), 25:28 (ftest variant 2), 29:32
#'   (testye variant 1), 33:36 (testye variant 2), 37:40 (fftest variant 1), and 41:44
#'   (fftest variant 2). In Scenario (ii), fftest1.ind/fftest2.ind occupy 21:24/25:28,
#'   testye occupies 29:32/33:36, and fftest1.mf.ind/fftest2.mf.ind occupy 37:40/41:44.
#'   Calls require the loaded simulation dependencies, including carat. This function
#'   returns results to the caller; the job script is responsible for saving res1 to a
#'   seed file.
delta <- function(i, n = 1000){
  if (!requireNamespace("carat", quietly = TRUE)) stop("Install carat and its dependencies before running delta().")
  
  # Scenario 1 callbacks: each sapply function takes scalar effect h and
  # returns the rejection decision(s) from the named trial test. The fun
  # callback takes raw x, z vectors and returns g(X)+f(Z) per row.
  # Column-major flattening keeps effects ordered; two-decision tests alternate
  # adjusted and naive decisions within each effect.
  # case 1
  power1 = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs1.dis(i, n, h, o2 = 3, 
                                                               knots2 = c(0), fun = function(x, z) 2*x + exp(z)/2.5))
  power2 = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs2.dis(i, n, h, o2 = 3, 
                                                               knots2 = c(0), fun = function(x, z) 2*x + exp(z)/2.5))
  power3 = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs3.dis(i, n, h, o2 = 3, 
                                                               knots2 = c(0), fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerye1 = sapply(c(0, 0.1, 0.2, 0.3), function(h) testye1.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerye2 = sapply(c(0, 0.1, 0.2, 0.3), function(h) testye2.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerf1 = sapply(c(0, 0.1, 0.2, 0.3), function(h) ftest1.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerf2 = sapply(c(0, 0.1, 0.2, 0.3), function(h) ftest2.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerff1 = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest1.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  powerff2 = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest2.dis(i, n, h, fun = function(x, z) 2*x + exp(z)/2.5))
  
  # Scenario 2 callbacks: each sapply function takes scalar effect h and
  # returns the rejection decision(s) from the named trial test. The fun
  # callback takes raw x, z vectors and returns g(X)+f(Z) per row.
  # Column-major flattening keeps effects ordered; two-decision tests alternate
  # adjusted and naive decisions within each effect.
  # case 2
  power1.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs1.ind(i, n, h, o1 = 3, o2 = 3, knots1 = c(0), 
                                                                   knots2 = c(0), fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  power2.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs2.ind(i, n, h, o1 = 3, o2 = 3, knots1 = c(0), 
                                                                   knots2 = c(0), fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  power3.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs3.ind(i, n, h, o1 = 3, o2 = 3, knots1 = c(0), 
                                                                   knots2 = c(0), fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerye1.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) testye1.ind(i, n, h, 10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerye2.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) testye2.ind(i, n, h, 10,fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerff1.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest1.ind(i, n, h,10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerff2.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest2.ind(i, n, h,10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerff1.mf.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest1.mf.ind(i, n, h,10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
  
  powerff2.mf.ind = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest2.mf.ind(i, n, h,10, fun = function(x, z) exp(x) + 4 * log(z + 4)))
 
  # Scenario 3 callbacks: each sapply function takes scalar effect h and
  # returns the rejection decision(s) from the named trial test. The fun
  # callback takes raw x1, x2, z vectors and returns g(X)+f(Z) per row.
  # Column-major flattening keeps effects ordered; two-decision tests alternate
  # adjusted and naive decisions within each effect.
  # case 3
  power1.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs1.ts(i, n, h, d1 = 3, d2 = 3, knots1 = NULL, 
                                                                 knots2 = NULL, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  power2.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs2.ts(i, n, h, d1 = 3, d2 = 3, knots1 = NULL, 
                                                                 knots2 = NULL, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  power3.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) testbs3.ts(i, n, h, d1 = 3, d2 = 3, knots1 = NULL, 
                                                                 knots2 = NULL, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerf1.c3 = sapply(c(0, 0.1, 0.2, 0.3), function(h) ftest1.c3(i, n, h, 3, 4, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerf2.c3 = sapply(c(0, 0.1, 0.2, 0.3), function(h) ftest2.c3(i, n, h, 3, 4, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerye1.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) testye1.ts(i, n, h,3,4, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerye2.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) testye2.ts(i, n, h,3,4, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerff1.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest1.ts(i, n, h,3,4, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  powerff2.ts = sapply(c(0, 0.1, 0.2, 0.3), function(h) fftest2.ts(i, n, h, 3,4, fun = function(x1, x2, z) 2 * x1 * x2^2 + exp(z) / 2.5))
  
  return(c(power1, power2, power3, powerf1, powerf2, powerye1, powerye2, powerff1, powerff2,
           power1.ind, power2.ind, power3.ind, powerff1.ind, powerff2.ind, powerye1.ind, powerye2.ind, powerff1.mf.ind, powerff2.mf.ind,
           power1.ts, power2.ts, power3.ts, powerf1.c3, powerf2.c3, powerye1.ts, powerye2.ts, powerff1.ts, powerff2.ts))
  
}
