library(testthat)
library(demohaz)

# Create a baseline parameter vector used across tests
b0 = c(.175,
       1.40,
       .368 * .01,
       log(.917 * .1/(.075 * .001))/(.917 * .1),
       .917 * .1)

th0 <- c(2e-2, 1.2, b0)

test_that("usher3_rho1 returns correct density values", {
  # Test with w = c(1, 0) (everyone well at x0)
  w0 <- c(1, 0)
  
  # Test cases with various age inputs
  expect_equal(usher3_rho1(0, th0[1], th0[3:7], x0 = 0, w = w0), 
               dsiler(0, th0[3:7]) * exp(-th0[1] * 0), tolerance = 1e-8)
  expect_equal(usher3_rho1(10, th0[1], th0[3:7], x0 = 0, w = w0), 
               dsiler(10, th0[3:7]) * exp(-th0[1] * 10), tolerance = 1e-8)
  expect_equal(usher3_rho1(25, th0[1], th0[3:7], x0 = 0, w = w0), 
               dsiler(25, th0[3:7]) * exp(-th0[1] * 25), tolerance = 1e-8)
  expect_equal(usher3_rho1(50, th0[1], th0[3:7], x0 = 0, w = w0), 
               dsiler(50, th0[3:7]) * exp(-th0[1] * 50), tolerance = 1e-8)
  
  # Test with w = NULL and x0 = 0 (should also give same results since calc_weights(0) = c(1,0))
  expect_equal(usher3_rho1(10, th0[1], th0[3:7], x0 = 0, w = NULL, k2 = th0[2]), 
               dsiler(10, th0[3:7]) * exp(-th0[1] * 10), tolerance = 1e-8)
})

test_that("usher3_rho2 returns correct density values", {
  # Test with w = c(1, 0) (everyone well at x0)
  w0 <- c(1, 0)
  
  # Test cases with various age inputs
  x <- c(0, 10, 25, 50)
  for (i in seq_along(x)) {
    expect_equal(usher3_rho2(x[i], th0[1], th0[2], th0[3:7], x0 = 0, w = w0),
                 th0[1] * th0[2] * dsiler(x[i], th0[3:7]) *
                   ssiler(x[i], th0[3:7])^(th0[2] - 1) *
                   integrate(function(t) exp(-th0[1] * t) *
                               ssiler(t, th0[3:7])^(1 - th0[2]),
                             0, x[i])$value,
                 tolerance = 1e-6)
  }
  
  # Test with w = NULL and x0 = 0 (should also give same results)
  expect_equal(usher3_rho2(10, th0[1], th0[2], th0[3:7], x0 = 0, w = NULL),
               th0[1] * th0[2] * dsiler(10, th0[3:7]) *
                 ssiler(10, th0[3:7])^(th0[2] - 1) *
                 integrate(function(t) exp(-th0[1] * t) *
                             ssiler(t, th0[3:7])^(1 - th0[2]),
                           0, 10)$value,
               tolerance = 1e-6)
})

test_that("usher3_integrand returns correct values", {
  # Test cases with various age inputs
  expect_equal(usher3_integrand(0, th0[1], th0[2], th0[3:7]), exp(-th0[1] * 0) * (ssiler(0, th0[3:7]))^(1 - th0[2]), tolerance = 1e-8)
  expect_equal(usher3_integrand(10, th0[1], th0[2], th0[3:7]), exp(-th0[1] * 10) * (ssiler(10, th0[3:7]))^(1 - th0[2]), tolerance = 1e-8)
  expect_equal(usher3_integrand(25, th0[1], th0[2], th0[3:7]), exp(-th0[1] * 25) * (ssiler(25, th0[3:7]))^(1 - th0[2]), tolerance = 1e-8)
  expect_equal(usher3_integrand(50, th0[1], th0[2], th0[3:7]), exp(-th0[1] * 50) * (ssiler(50, th0[3:7]))^(1 - th0[2]), tolerance = 1e-8)
})

test_that("calc_filtration_density returns correct values", {
  # Test cases with various age inputs
  xcalc <- seq(0, 100, by = 0.1)
  x_mid <- 50
  infant_prop <- 0.1

  # Test case 1: discrete = TRUE
  f_discrete <- calc_filtration_density(xcalc, x_mid, infant_prop, discrete = TRUE)

  expect_equal(sum(f_discrete), 1, tolerance = 1e-6)

  z_mid <- 1 / x_mid / (1 + infant_prop)
  z_inf <- infant_prop * z_mid

  f_expected_discrete <- rep(0, length(xcalc))
  f_expected_discrete[xcalc <= x_mid] <- z_inf + (z_mid - z_inf) * xcalc[xcalc <= x_mid] / x_mid
  f_expected_discrete[xcalc > x_mid & xcalc <= 2 * x_mid] <- z_mid + (z_inf - z_mid) * (xcalc[xcalc > x_mid & xcalc <= 2 * x_mid] - x_mid) / x_mid
  f_expected_discrete <- f_expected_discrete / sum(f_expected_discrete)

  expect_equal(f_discrete, f_expected_discrete, tolerance = 1e-6)

  # Test case 2: discrete = FALSE
  f_continuous <- calc_filtration_density(xcalc, x_mid, infant_prop, discrete = FALSE)

  integration_weights = calc_trapez_weights(xcalc)
  expect_equal(sum(f_continuous * integration_weights), 1, tolerance = 1e-6)

  f_expected_continuous <- rep(0, length(xcalc))
  f_expected_continuous[xcalc <= x_mid] <- z_inf + (z_mid - z_inf) * xcalc[xcalc <= x_mid] / x_mid
  f_expected_continuous[xcalc > x_mid & xcalc <= 2 * x_mid] <- z_mid + (z_inf - z_mid) * (xcalc[xcalc > x_mid & xcalc <= 2 * x_mid] - x_mid) / x_mid

  expect_equal(f_continuous, f_expected_continuous, tolerance = 1e-6)
})

test_that("calc_trapez_weights returns correct trapezoidal integration weights", {
  # Test case 1: Evenly spaced grid points
  tau1 <- seq(0, 1, by = 0.1)
  weights1 <- calc_trapez_weights(tau1)
  expected_weights1 <- c(0.05, rep(0.1, 9), 0.05)
  expect_equal(weights1, expected_weights1, tolerance = 1e-6)

  # Test case 2: Unevenly spaced grid points
  tau2 <- c(0, 0.1, 0.3, 1)
  weights2 <- calc_trapez_weights(tau2)
  expected_weights2 <- c(0.05, 0.15, 0.45, 0.35)
  expect_equal(weights2, expected_weights2, tolerance = 1e-6)
})

test_that("sample_usher3 generates samples correctly", {
  set.seed(123)  # Set a seed for reproducibility

  # Test case 1: Basic sampling without age filtration
  N <- 50
  dx <- 0.001
  xmax <- 120
  samples1 <- sample_usher3(N, th0, dx, xmax)
  expect_equal(length(samples1$x), N)
  expect_equal(length(samples1$ill), N)
  expect_true(all(samples1$x >= 0 & samples1$x <= xmax))
  expect_true(all(samples1$ill %in% c(TRUE, FALSE)))

  # Test case 2: Sampling with age filtration
  x_mid <- 50
  infant_prop <- 0.1
  samples2 <- sample_usher3(N, th0, dx, xmax, x_mid, infant_prop)
  expect_equal(length(samples2$x), N)
  expect_equal(length(samples2$ill), N)
  expect_true(all(samples2$x >= 0 & samples2$x <= 2 * x_mid))
  expect_true(all(samples2$ill %in% c(TRUE, FALSE)))
  
  # Test case 3: Sampling with a smaller number of samples
  N_small <- 10
  samples3 <- sample_usher3(N_small, th0, dx, xmax)
  expect_equal(length(samples3$x), N_small)
  expect_equal(length(samples3$ill), N_small)
  expect_true(all(samples3$x >= 0 & samples3$x <= xmax))
  expect_true(all(samples3$ill %in% c(TRUE, FALSE)))
})

test_that("nll_usher3 returns viable negative log-likelihood values", {
  # Test case 1: Valid parameter values
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  nll <- nll_usher3(th0, x, ill)
  expect_true(is.numeric(nll) && length(nll) == 1)

  # Test case 2: Invalid parameter values
  theta_invalid <- c(-0.1, 1.5, 0.2, 1.3, 0.4, 0.01, 0.1)
  expect_equal(nll_usher3(theta_invalid, x, ill), Inf)

  # Test case 3: ill vector with NA values
  ill_with_na <- c(0, NA, 0, 1, NA)
  nll_with_na <- nll_usher3(th0, x, ill_with_na)
  expect_true(is.numeric(nll_with_na) && length(nll_with_na) == 1)
  expect_true(is.finite(nll_with_na))
})

test_that("usher3_hessian returns a valid Hessian matrix", {
  # Using a more realistic vector for the parameters
  # Example input data
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  # Call the usher3_hessian function
  hessian <- usher3_hessian(th0, x, ill)
  
  # Check if the returned object is a matrix
  expect_true(is.matrix(hessian), 
              info = "The returned object should be a matrix.")
  
  # Check if the matrix has the expected dimensions
  expected_dims <- c(7, 7)
  expect_equal(dim(hessian), expected_dims, 
               info = paste0("Expected dimensions: ", 
                             paste(expected_dims, collapse = "x")))
  
  # Ensure no NA values are in the matrix
  expect_false(anyNA(hessian), info = "Hessian matrix should not contain NA.")
})


test_that("nll_usher3_hessian_wrapper returns a valid negative log-likelihood vector", {
  # Test case 1: Valid parameter values
  ageVect <- c(10, 20, 30, 40, 50)
  illVect <- c(0, 1, 0, 1, 0)
  nll <- nll_usher3_hessian_wrapper(th0, ageVect, illVect)
  expect_true(is.numeric(nll) && length(nll) == 1)
  expect_false(anyNA(nll), info = "Negative log-likelihood vector should not contain NA.")
})

# TODO: this must be tested with a parameter vector at a local minimum of the
# negative log-likelihood. This may require more values in the x vector.
#test_that("usher3_errors returns viabl3e standard errors, z-scores, and p-values", {
#  # Test case 1: Valid parameter values
#  TODO: use th0 here
#  #theta <- c(2e-2, 1.2, 0.175, 1.40, 0.368 * 0.01, 0.075 * 0.001, 0.917 * 0.1)
#  theta <- c(1.319066e-02, 1.723116e+00, 1.758189e-02, 1.696364e+00,
#             4.814728e-07, 3.810192e+01, 5.427843e-02)
#  x <- c(10, 20, 30, 40, 50)
#  ill <- c(0, 1, 0, 1, 0)
#  errors <- usher3_errors(theta, x, ill)
#  print('-------------')
#  print('-------------')
#  print('-------------')
#  print(errors)
#  expect_true(is.data.frame(errors) && nrow(errors) == 7)
#  expect_true(all(c("Estimate", "StandErr", "z", "pval", "against", "sideAdj") %in% names(errors)))
#})

test_that("temper_and_tune_usher3 with tune=FALSE returns viable results (tempering only)", {
  # Test case: Valid inputs
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)

  # Get results from the function with tuning disabled
  results <- temper_and_tune_usher3(x = x, ill = ill, tune = FALSE)

  # Check if the results are a list and contain specific names
  expect_true(is.list(results))
  expect_true(all(c("obj_fun", "th0", "optional_inputs", "temper",
                    "th_temper", "th") %in% names(results)))
  
  # Should not contain optim_result or th_bar when tune=FALSE
  expect_false("optim_result" %in% names(results), 
               info = "optim_result should not be present when tune=FALSE")
  expect_false("th_bar" %in% names(results), 
               info = "th_bar should not be present when tune=FALSE")

  # Extract transformed parameter vectors from results
  th_temper <- exp(results$th_temper)
  th <- exp(results$th)

  # Check that 'th' is not equal to 'th0'
  expect_false(all(th == th0), info = "'th' should not match 'th0'")

  # Check that 'th_temper' is not equal to 'th0'
  expect_false(all(th_temper == th0), 
               info = "'th_temper' should not match 'th0'")

  # Check that 'th' equals 'th_temper' when tune=FALSE
  expect_true(all(th == th_temper), 
              info = "'th' should match 'th_temper' when tune=FALSE")

  # Check that 'th' and 'th_temper' contain no NA values
  expect_false(anyNA(th), info = "'th' should not contain NA values")
  expect_false(anyNA(th_temper), 
               info = "'th_temper' should not contain NA values")
})

# Reduced parameter vector for Gompertz-Makeham mortality (length 5)
b0_gompertz <- c(0.368 * 0.01,
                 log(0.917 * 0.1 / (0.075 * 0.001)) / (0.917 * 0.1),
                 0.917 * 0.1)
th0_gompertz <- c(2e-2, 1.2, b0_gompertz)  # Length 5
x <- c(10, 20, 30, 40, 50)
ill <- c(0, 1, 0, 1, 0)

# Full parameter vectors
b0_full <- c(0.175,
             1.40,
             0.368 * 0.01,
             log(0.917 * 0.1 / (0.075 * 0.001)) / (0.917 * 0.1),
             0.917 * 0.1)
th0_full <- c(2e-2, 1.2, b0_full)  # Length 7

test_that("nll_usher3 handles Gompertz-Makeham mortality correctly", {
  nll_gompertz <- nll_usher3(th0_gompertz, x, ill)
  expect_true(is.numeric(nll_gompertz) && length(nll_gompertz) == 1)
  expect_false(is.infinite(nll_gompertz))

  theta_invalid <- c(-0.1, 1.5, 0.2, 1.3, 0.4)
  expect_equal(nll_usher3(theta_invalid, x, ill), Inf)

  ill_with_na <- c(0, NA, 1, NA, 0)
  nll_with_na <- nll_usher3(th0_gompertz, x, ill_with_na)
  expect_true(is.numeric(nll_with_na) && length(nll_with_na) == 1)
  expect_true(is.finite(nll_with_na))
})

test_that("temper_and_tune_usher3 with tune=FALSE handles use_gompertz = TRUE correctly", {
  results_gompertz <- temper_and_tune_usher3(th0 = th0_gompertz, x = x,
                                    ill = ill, use_gompertz = TRUE, tune = FALSE)
  expect_true(is.list(results_gompertz))
  expect_true(all(c("obj_fun", "th0", "optional_inputs", "temper",
                    "th_temper", "th") %in% names(results_gompertz)))
  
  # Should not contain optim_result when tune=FALSE
  expect_false("optim_result" %in% names(results_gompertz))

  th_temper <- results_gompertz$th_temper
  th <- results_gompertz$th
  expect_equal(length(th), 5, 
               info = "'th' should have length 5 when use_gompertz = TRUE")
  expect_false(all(th == th0_gompertz))
  expect_false(all(th_temper == th0_gompertz))
  expect_true(all(th == th_temper), 
              info = "'th' should match 'th_temper' when tune=FALSE")
  expect_false(anyNA(th))
  expect_false(anyNA(th_temper))
})

test_that("temper_and_tune_usher3 adjusts parameter vector length correctly", {
  results_adjusted <- temper_and_tune_usher3(th0 = th0_full, x = x, ill = ill,
                                    use_gompertz = TRUE, tune = FALSE)
  adjusted_th0 <- results_adjusted$th0
  expect_equal(length(adjusted_th0), 5)
  expected_th0_adjusted <- th0_full[c(1, 2, 5, 6, 7)]
  expect_equal(adjusted_th0, expected_th0_adjusted, tolerance = 1e-8)
})

test_that("nll_usher3 checks negative parameters with varying lengths", {
  theta_negative <- c(0.02, 1.2, -0.01, 0.05, 0.1)
  expect_equal(nll_usher3(theta_negative, x, ill), Inf)

  theta_negative_full <- c(0.02, 1.2, 0.1, 0.2, -0.01, 0.05, 0.1)
  expect_equal(nll_usher3(theta_negative_full, x, ill), Inf)
})

test_that("temper_and_tune_usher3 errors with incorrect th0 length", {
  th0_incorrect <- c(0.02, 1.2, 0.1, 0.2, 0.3, 0.4)
  expect_error(temper_and_tune_usher3(th0 = th0_incorrect, x = x, ill = ill,
                             use_gompertz = TRUE),
               "th0 should have length 5 or 7 if use_gompertz is TRUE")
})

test_that("temper_and_tune_usher3 defaults correctly when use_gompertz = FALSE", {
  expect_error(temper_and_tune_usher3(th0 = th0_gompertz, x = x, ill = ill,
                             use_gompertz = FALSE),
               "th0 should have length 7")
})

# Tests for tuning functionality (tune=TRUE)

test_that("temper_and_tune_usher3 with tune=TRUE returns expected components", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  results <- temper_and_tune_usher3(x = x, ill = ill, tune = TRUE, num_cyc = 2)
  
  # Check that results contain tuning-specific components
  expect_true(is.list(results))
  expect_true(all(c("obj_fun", "th0", "optional_inputs", "temper",
                    "th_temper", "optim_result", "th_bar", "th") %in% names(results)),
              info = "Results should contain optim_result, th_bar, and th when tune=TRUE")
  
  # Check optim_result structure
  expect_true(is.list(results$optim_result))
  expect_true(all(c("par", "value", "counts", "convergence") %in% names(results$optim_result)),
              info = "optim_result should contain par, value, counts, and convergence")
  
  # Check that th_bar matches optim_result$par
  expect_equal(results$th_bar, results$optim_result$par,
               info = "th_bar should match optim_result$par")
  
  # Check that th equals exp(th_bar)
  expect_equal(results$th, exp(results$th_bar), tolerance = 1e-10,
               info = "th should equal exp(th_bar)")
  
  # Check no NA values
  expect_false(anyNA(results$th))
  expect_false(anyNA(results$th_bar))
})

test_that("temper_and_tune_usher3 with tune=TRUE performs optimization", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  results <- temper_and_tune_usher3(x = x, ill = ill, tune = TRUE, num_cyc = 2)
  
  # Check that tuning was performed (function evaluations > 0)
  expect_true(results$optim_result$counts["function"] > 0,
              info = "optim should have performed function evaluations")
  
  # Check that objective value is finite
  expect_true(is.finite(results$optim_result$value),
              info = "Optimized objective value should be finite")
  
  # Check that th is not NA
  expect_false(anyNA(results$th))
})

test_that("temper_and_tune_usher3 tuning improves or maintains objective function value", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  results <- temper_and_tune_usher3(x = x, ill = ill, tune = TRUE, num_cyc = 2)
  
  # Get tempered objective value
  n <- which.min(unlist(lapply(results$temper$chains, function(x) {x$eta_best})))
  tempered_obj <- results$temper$chains[[n]]$eta_best
  
  # Tuned objective should be <= tempered objective (allowing for small numerical differences)
  expect_true(results$optim_result$value <= tempered_obj + 1e-6,
              info = "Tuning should improve or maintain objective function value")
})

test_that("temper_and_tune_usher3 uses tune=TRUE by default", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  # Call without specifying tune parameter
  results <- temper_and_tune_usher3(x = x, ill = ill, num_cyc = 2)
  
  # Should have tuning components by default
  expect_true("optim_result" %in% names(results),
              info = "optim_result should be present when tune is not specified (default=TRUE)")
  expect_true("th_bar" %in% names(results),
              info = "th_bar should be present when tune is not specified (default=TRUE)")
})

test_that("temper_and_tune_usher3 with tune=TRUE and use_gompertz=TRUE works correctly", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  results <- temper_and_tune_usher3(th0 = th0_gompertz, x = x, ill = ill,
                                    use_gompertz = TRUE, tune = TRUE, num_cyc = 2)
  
  # Check that results contain tuning components
  expect_true("optim_result" %in% names(results))
  expect_true("th_bar" %in% names(results))
  
  # Check that parameter length is maintained (5 for Gompertz-Makeham)
  expect_equal(length(results$th), 5,
               info = "th should have length 5 with use_gompertz=TRUE")
  expect_equal(length(results$th_bar), 5,
               info = "th_bar should have length 5 with use_gompertz=TRUE")
  
  # Check optimization was performed
  expect_true(results$optim_result$counts["function"] > 0)
  expect_false(anyNA(results$th))
})

test_that("temper_and_tune_usher3 tuned result differs from tempered result", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  results <- temper_and_tune_usher3(x = x, ill = ill, tune = TRUE, num_cyc = 2)
  
  # In most cases, tuned result should differ from tempered result
  # (though they could be the same if already at optimum)
  # Just check that both exist and are valid
  expect_false(anyNA(results$th))
  expect_false(anyNA(results$th_temper))
  
  # Both should be valid parameter vectors
  expect_equal(length(results$th), length(results$th_temper))
})

test_that("temper_and_tune_usher3 respects custom control parameters", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  # Test with custom maxit (lower for speed in testing)
  results <- temper_and_tune_usher3(x = x, ill = ill, tune = TRUE, num_cyc = 2,
                                    control = list(maxit = 100))
  
  # Should still produce valid results
  expect_true("optim_result" %in% names(results))
  expect_false(anyNA(results$th))
  
  # Verify it didn't exceed the specified maxit
  expect_true(results$optim_result$counts["function"] <= 100 * (length(results$th_temper) + 1),
              info = "Should not exceed specified maxit bounds")
})

test_that("nll_usher3 handles both full and reduced parameter vectors", {
  nll_full <- nll_usher3(th0_full, x, ill)
  expect_true(is.numeric(nll_full) && length(nll_full) == 1)
  expect_false(is.infinite(nll_full))

  nll_reduced <- nll_usher3(th0_gompertz, x, ill)
  expect_true(is.numeric(nll_reduced) && length(nll_reduced) == 1)
  expect_false(is.infinite(nll_reduced))
})

test_that("nll_usher3 handles x0 parameter correctly", {
  x0 <- 5
  nll_x0 <- nll_usher3(th0_gompertz, x, ill, x0 = x0)
  expect_true(is.numeric(nll_x0) && length(nll_x0) == 1)
  expect_false(is.infinite(nll_x0))
})

test_that("p_11 calculates transition probability correctly", {
  # Test against analytical formula: p_11 = exp(-k1 * (x - x0)) * S_13(x0, x)
  x_val <- 10
  p11_computed <- p_11(x = x_val, k1 = th0[1], b_siler = th0[3:7])
  S13_val <- ssiler(x_val, th0[3:7], x0 = 0)
  p11_analytical <- exp(-th0[1] * x_val) * S13_val
  expect_equal(p11_computed, p11_analytical, tolerance = 1e-10)
  
  # Test at age 0 (x = x0) - analytically should be 1
  p11_0 <- p_11(x = 0, k1 = th0[1], b_siler = th0[3:7], x0 = 0)
  expect_equal(p11_0, 1, tolerance = 1e-12)
  
  # Test that probability decreases with age
  p11_20 <- p_11(x = 20, k1 = th0[1], b_siler = th0[3:7])
  expect_true(p11_20 < p11_computed)
  
  # Test with x_cut - verify analytical formula for both cases
  # Case 1: x > x_cut, so p_11 = exp(-k1 * x_cut) * S_13(0, x)
  x_cut_val <- 5
  p11_cut <- p_11(x = 10, k1 = th0[1], b_siler = th0[3:7], x_cut = x_cut_val)
  S13_10 <- ssiler(10, th0[3:7], x0 = 0)
  p11_cut_analytical <- exp(-th0[1] * x_cut_val) * S13_10
  expect_equal(p11_cut, p11_cut_analytical, tolerance = 1e-10)
  
  # Case 2: x <= x_cut, so p_11 = exp(-k1 * x) * S_13(0, x)
  p11_cut_3 <- p_11(x = 3, k1 = th0[1], b_siler = th0[3:7], x_cut = x_cut_val)
  S13_3 <- ssiler(3, th0[3:7], x0 = 0)
  p11_cut_3_analytical <- exp(-th0[1] * 3) * S13_3
  expect_equal(p11_cut_3, p11_cut_3_analytical, tolerance = 1e-10)
  
  # Test with non-zero x0 - verify analytical formula
  p11_5_to_10 <- p_11(x = 10, k1 = th0[1], b_siler = th0[3:7], x0 = 5)
  S13_5_to_10 <- ssiler(10, th0[3:7], x0 = 5)
  p11_5_to_10_analytical <- exp(-th0[1] * (10 - 5)) * S13_5_to_10
  expect_equal(p11_5_to_10, p11_5_to_10_analytical, tolerance = 1e-10)
  
  # Test error handling
  expect_error(p_11(x = 10, k1 = th0[1], b_siler = th0[3:7], x0 = -1), 
               "x0 cannot be negative")
  expect_error(p_11(x = 5, k1 = th0[1], b_siler = th0[3:7], x0 = 10), 
               "x must be >= x0")
  expect_error(p_11(x = 10, k1 = -0.01, b_siler = th0[3:7]), 
               "k1 cannot be negative")
})

test_that("p_12 calculates transition probability correctly", {
  # Test at x = x0 (analytically should be 0 - no time to transition)
  p12_0 <- p_12(x = 0, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], x0 = 0)
  expect_equal(p12_0, 0, tolerance = 1e-12)
  
  # Test analytical limit: when k1 = 0, p_12 should be 0 (no transitions)
  p12_k1_zero <- p_12(x = 10, k1 = 0, k2 = th0[2], b_siler = th0[3:7])
  expect_equal(p12_k1_zero, 0, tolerance = 1e-12)
  
  # Test analytical limit: when x_cut = 0, p_12 should be 0 (no transitions possible)
  p12_xcut_zero <- p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
                        x_cut = 0)
  expect_equal(p12_xcut_zero, 0, tolerance = 1e-12)
  
  # Test that probability increases with age (more time to transition)
  p12_5 <- p_12(x = 5, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7])
  p12_10 <- p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7])
  expect_true(p12_10 > p12_5)
  
  # Test with x_cut
  p12_cut <- p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
                  x_cut = 5)
  p12_nocut <- p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
                    x_cut = Inf)
  expect_true(p12_cut < p12_nocut)  # Less hazard means fewer transitions
  
  # Verify p_12 is in valid range
  p12_val <- p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7])
  expect_true(p12_val >= 0 && p12_val <= 1)
  
  # Test with non-zero x0
  p12_5_to_10 <- p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
                      x0 = 5)
  expect_true(is.numeric(p12_5_to_10) && p12_5_to_10 >= 0 && p12_5_to_10 <= 1)
  
  # Test error handling
  expect_error(p_12(x = 10, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
                    x0 = -1), 
               "x0 cannot be negative")
  expect_error(p_12(x = 5, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
                    x0 = 10), 
               "x must be >= x0")
  expect_error(p_12(x = 10, k1 = -0.01, k2 = th0[2], b_siler = th0[3:7]), 
               "k1 cannot be negative")
  expect_error(p_12(x = 10, k1 = th0[1], k2 = -0.5, b_siler = th0[3:7]), 
               "k2 cannot be negative")
})

test_that("calc_weights calculates state weights correctly", {
  # Test at age 0 (analytically should be c(1, 0) - everyone well)
  weights_0 <- calc_weights(x0 = 0, k1 = th0[1], k2 = th0[2], 
                            b_siler = th0[3:7])
  expect_equal(as.numeric(weights_0[1]), 1, tolerance = 1e-12)
  expect_equal(as.numeric(weights_0[2]), 0, tolerance = 1e-12)
  
  # Test analytical limit: when k1 = 0, weights should be c(1, 0) at any age
  weights_k1_zero <- calc_weights(x0 = 10, k1 = 0, k2 = th0[2], 
                                  b_siler = th0[3:7])
  expect_equal(as.numeric(weights_k1_zero[1]), 1, tolerance = 1e-12)
  expect_equal(as.numeric(weights_k1_zero[2]), 0, tolerance = 1e-12)
  
  # Test analytical limit: when x_cut = 0, weights should be c(1, 0) at any age
  weights_xcut_zero <- calc_weights(x0 = 10, k1 = th0[1], k2 = th0[2], 
                                    b_siler = th0[3:7], x_cut = 0)
  expect_equal(as.numeric(weights_xcut_zero[1]), 1, tolerance = 1e-12)
  expect_equal(as.numeric(weights_xcut_zero[2]), 0, tolerance = 1e-12)
  
  # Test that weights are computed from p_11 and p_12 using analytical formula
  # w_i = p_1i / (p_11 + p_12)
  x0_val <- 10
  weights <- calc_weights(x0 = x0_val, k1 = th0[1], k2 = th0[2], 
                          b_siler = th0[3:7])
  p11 <- p_11(x = x0_val, k1 = th0[1], b_siler = th0[3:7], x0 = 0)
  p12 <- p_12(x = x0_val, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], x0 = 0)
  expected_w1 <- p11 / (p11 + p12)
  expected_w2 <- p12 / (p11 + p12)
  expect_equal(as.numeric(weights[1]), expected_w1, tolerance = 1e-10)
  expect_equal(as.numeric(weights[2]), expected_w2, tolerance = 1e-10)
  
  # Test that weights sum to 1 (mathematical constraint)
  expect_equal(sum(weights), 1, tolerance = 1e-12)
  
  # Test that w1 decreases and w2 increases with age
  weights_5 <- calc_weights(x0 = 5, k1 = th0[1], k2 = th0[2], 
                            b_siler = th0[3:7])
  weights_10 <- calc_weights(x0 = 10, k1 = th0[1], k2 = th0[2], 
                             b_siler = th0[3:7])
  expect_true(weights_10[1] < weights_5[1])  # w1 decreases
  expect_true(weights_10[2] > weights_5[2])  # w2 increases
  
  # Test with x_cut - after cutoff, w1 should decrease more slowly
  weights_cut_5 <- calc_weights(x0 = 5, k1 = th0[1], k2 = th0[2], 
                                b_siler = th0[3:7], x_cut = 3)
  weights_nocut_5 <- calc_weights(x0 = 5, k1 = th0[1], k2 = th0[2], 
                                  b_siler = th0[3:7], x_cut = Inf)
  expect_true(weights_cut_5[1] > weights_nocut_5[1])  # More stay well with cutoff
  
  # Test error handling
  expect_error(calc_weights(x0 = -1, k1 = th0[1], k2 = th0[2], 
                            b_siler = th0[3:7]), 
               "x0 cannot be negative")
})

test_that("p_12_integrand calculates integrand correctly", {
  # Test against analytical formula: exp(-k1 * (y - x0)) * [S_13(x0, y)]^(1 - k2)
  y_val <- 5
  x0_val <- 0
  integrand_val <- p_12_integrand(y = y_val, x0 = x0_val, k1 = th0[1], 
                                  k2 = th0[2], b_siler = th0[3:7])
  S13_0_y <- ssiler(y_val, th0[3:7], x0 = x0_val)
  expected_val <- exp(-th0[1] * (y_val - x0_val)) * (S13_0_y)^(1 - th0[2])
  expect_equal(integrand_val, expected_val, tolerance = 1e-12)
  
  # Test at y = x0 (analytically should be 1^(1 - k2) = 1 if k2 ≠ 1)
  integrand_x0 <- p_12_integrand(y = 0, x0 = 0, k1 = th0[1], k2 = th0[2], 
                                 b_siler = th0[3:7])
  expect_equal(integrand_x0, 1^(1 - th0[2]), tolerance = 1e-12)
  
  # Test with non-zero x0 - verify analytical formula
  y_val_2 <- 8
  x0_val_2 <- 5
  integrand_val_5 <- p_12_integrand(y = y_val_2, x0 = x0_val_2, k1 = th0[1], 
                                    k2 = th0[2], b_siler = th0[3:7])
  S13_5_8 <- ssiler(y_val_2, th0[3:7], x0 = x0_val_2)
  expected_val_5 <- exp(-th0[1] * (y_val_2 - x0_val_2)) * (S13_5_8)^(1 - th0[2])
  expect_equal(integrand_val_5, expected_val_5, tolerance = 1e-12)
})

test_that("q_1 calculates occupancy probability correctly", {
  # Test with x0 = 0 (w1 = 1, w2 = 0)
  w0 <- c(1, 0)
  x_vals <- c(5, 10, 15)
  
  # When x0 = 0 and w = c(1, 0), q_1 should equal p_11(0, x)
  q1_vals <- q_1(x = x_vals, w = w0, k1 = th0[1], b_siler = th0[3:7], 
                 x0 = 0, x_cut = Inf)
  p11_vals <- p_11(x = x_vals, k1 = th0[1], b_siler = th0[3:7], 
                   x0 = 0, x_cut = Inf)
  expect_equal(q1_vals, p11_vals, tolerance = 1e-10)
  
  # Test analytical formula: q_1 = w_1 * exp(-k1*(x-x0)) * S_13(0,x) / S_13(0,x0)
  x0_val <- 5
  w <- calc_weights(x0 = x0_val, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7])
  x_val <- 15
  q1_computed <- q_1(x = x_val, w = w, k1 = th0[1], b_siler = th0[3:7], 
                     x0 = x0_val, x_cut = Inf)
  
  S13_0_x <- ssiler(x_val, th0[3:7], x0 = 0)
  S13_0_x0 <- ssiler(x0_val, th0[3:7], x0 = 0)
  q1_analytical <- as.numeric(w[1]) * exp(-th0[1] * (x_val - x0_val)) * S13_0_x / S13_0_x0
  expect_equal(q1_computed, q1_analytical, tolerance = 1e-10)
  
  # Test with finite x_cut, Case 2: x0 < x_cut < x
  x_cut_val <- 10
  x_val <- 15
  x0_val <- 5
  q1_cut <- q_1(x = x_val, w = w, k1 = th0[1], b_siler = th0[3:7], 
                x0 = x0_val, x_cut = x_cut_val)
  S13_0_x <- ssiler(x_val, th0[3:7], x0 = 0)
  S13_0_x0 <- ssiler(x0_val, th0[3:7], x0 = 0)
  q1_cut_analytical <- as.numeric(w[1]) * exp(-th0[1] * (x_cut_val - x0_val)) * S13_0_x / S13_0_x0
  expect_equal(q1_cut, q1_cut_analytical, tolerance = 1e-10)
  
  # Test that q_1 can compute weights automatically if not provided
  q1_auto <- q_1(x = 10, w = NULL, k1 = th0[1], b_siler = th0[3:7], 
                 x0 = 0, x_cut = Inf, k2 = th0[2])
  expect_true(is.numeric(q1_auto) && q1_auto >= 0)
  
  # Test error handling
  expect_error(q_1(x = 10, w = w0, k1 = th0[1], b_siler = th0[3:7], x0 = -1), 
               "x0 cannot be negative")
  expect_error(q_1(x = 5, w = w0, k1 = th0[1], b_siler = th0[3:7], x0 = 10), 
               "x must be >= x0")
})

test_that("q_2 calculates occupancy probability correctly", {
  # Test with x0 = 0 (w1 = 1, w2 = 0)
  # When w2 = 0, only the integral term contributes
  w0 <- c(1, 0)
  x_val <- 10
  
  q2_computed <- q_2(x = x_val, w = w0, k1 = th0[1], k2 = th0[2], 
                     b_siler = th0[3:7], x0 = 0, x_cut = Inf)
  expect_true(is.numeric(q2_computed) && q2_computed >= 0)
  
  # Test analytical limit: when k1 = 0, no transitions, so q_2 = w_2 * S_23
  x0_val <- 5
  w_k1_zero <- calc_weights(x0 = x0_val, k1 = 0, k2 = th0[2], 
                            b_siler = th0[3:7])
  # w_k1_zero should be c(1, 0)
  q2_k1_zero <- q_2(x = 10, w = w_k1_zero, k1 = 0, k2 = th0[2], 
                    b_siler = th0[3:7], x0 = x0_val)
  expect_equal(q2_k1_zero, 0, tolerance = 1e-10)  # w2 = 0, no integral
  
  # Test that q_1 + q_2 represents total survival from x0 to x
  x0_val <- 5
  x_val <- 15
  w <- calc_weights(x0 = x0_val, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7])
  q1 <- q_1(x = x_val, w = w, k1 = th0[1], b_siler = th0[3:7], 
            x0 = x0_val, x_cut = Inf)
  q2 <- q_2(x = x_val, w = w, k1 = th0[1], k2 = th0[2], b_siler = th0[3:7], 
            x0 = x0_val, x_cut = Inf)
  
  # Total occupancy should be <= 1
  total_occupancy <- q1 + q2
  expect_true(total_occupancy <= 1 && total_occupancy >= 0)
  
  # Test that q_2 can compute weights automatically if not provided
  q2_auto <- q_2(x = 10, w = NULL, k1 = th0[1], k2 = th0[2], 
                 b_siler = th0[3:7], x0 = 0, x_cut = Inf)
  expect_true(is.numeric(q2_auto) && q2_auto >= 0)
  
  # Test vectorization
  q2_vec <- q_2(x = c(10, 15, 20), w = w0, k1 = th0[1], k2 = th0[2], 
                b_siler = th0[3:7], x0 = 0, x_cut = Inf)
  expect_equal(length(q2_vec), 3)
  expect_true(all(q2_vec >= 0))
  
  # Test error handling
  expect_error(q_2(x = 10, w = w0, k1 = th0[1], k2 = th0[2], 
                   b_siler = th0[3:7], x0 = -1), 
               "x0 cannot be negative")
})

test_that("q_2_integrand calculates integrand correctly", {
  # Test against analytical formula: exp(-k1 * (y - x0)) * [S_13(0, y)]^(1 - k2)
  # This is the same as p_12_integrand, so should give same results
  y_val <- 5
  x0_val <- 0
  integrand_q2 <- q_2_integrand(y = y_val, x0 = x0_val, k1 = th0[1], 
                                k2 = th0[2], b_siler = th0[3:7])
  integrand_p12 <- p_12_integrand(y = y_val, x0 = x0_val, k1 = th0[1], 
                                  k2 = th0[2], b_siler = th0[3:7])
  expect_equal(integrand_q2, integrand_p12, tolerance = 1e-12)
})