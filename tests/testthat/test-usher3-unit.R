library(testthat)
library(demohaz)

# Create a baseline parameter vector used across tests
th0 <- c(2e-2, 1.2, 0.175, 1.40, 0.368 * 0.01, 0.075 * 0.001, 0.917 * 0.1)

test_that("usher3_rho1 returns correct density values", {
  # Test cases with various age inputs
  expect_equal(usher3_rho1(0, th0[1], th0[3:7]), dsiler(0, th0[3:7]) * exp(-th0[1] * 0), tolerance = 1e-8)
  expect_equal(usher3_rho1(10, th0[1], th0[3:7]), dsiler(10, th0[3:7]) * exp(-th0[1] * 10), tolerance = 1e-8)
  expect_equal(usher3_rho1(25, th0[1], th0[3:7]), dsiler(25, th0[3:7]) * exp(-th0[1] * 25), tolerance = 1e-8)
  expect_equal(usher3_rho1(50, th0[1], th0[3:7]), dsiler(50, th0[3:7]) * exp(-th0[1] * 50), tolerance = 1e-8)
})

test_that("usher3_rho2 returns correct density values", {
  # Test cases with various age inputs
  x <- c(0, 10, 25, 50)
  for (i in seq_along(x)) {
    expect_equal(usher3_rho2(x[i], th0[1], th0[2], th0[3:7]),
                 th0[1] * th0[2] * dsiler(x[i], th0[3:7]) *
                   ssiler(x[i], th0[3:7])^(th0[2] - 1) *
                   integrate(function(t) exp(-th0[1] * t) *
                               ssiler(t, th0[3:7])^(1 - th0[2]),
                             0, x[i])$value,
                 tolerance = 1e-6)
  }
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
  th <- c(2e-2, 1.2, 0.175, 1.40, 0.368, 0.01, 0.075, 0.001, 0.917, 0.1)
  dx <- 0.001
  xmax <- 120
  samples1 <- sample_usher3(N, th, dx, xmax)
  expect_equal(length(samples1$x), N)
  expect_equal(length(samples1$ill), N)
  expect_true(all(samples1$x >= 0 & samples1$x <= xmax))
  expect_true(all(samples1$ill %in% c(TRUE, FALSE)))

  # Test case 2: Sampling with age filtration
  x_mid <- 50
  infant_prop <- 0.1
  samples2 <- sample_usher3(N, th, dx, xmax, x_mid, infant_prop)
  expect_equal(length(samples2$x), N)
  expect_equal(length(samples2$ill), N)
  expect_true(all(samples2$x >= 0 & samples2$x <= 2 * x_mid))
  expect_true(all(samples2$ill %in% c(TRUE, FALSE)))
  
  # Test case 3: Sampling with a smaller number of samples
  N_small <- 10
  samples3 <- sample_usher3(N_small, th, dx, xmax)
  expect_equal(length(samples3$x), N_small)
  expect_equal(length(samples3$ill), N_small)
  expect_true(all(samples3$x >= 0 & samples3$x <= xmax))
  expect_true(all(samples3$ill %in% c(TRUE, FALSE)))
})

test_that("nll_usher3 returns the correct negative log-likelihood value", {
  # Test case 1: Valid parameter values
  theta <- c(0.1, 1.5, 0.2, 1.3, 0.4, 0.01, 0.1)
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  nll <- nll_usher3(theta, x, ill)
  expect_true(is.numeric(nll) && length(nll) == 1)

  # Test case 2: Invalid parameter values
  theta_invalid <- c(-0.1, 1.5, 0.2, 1.3, 0.4, 0.01, 0.1)
  expect_equal(nll_usher3(theta_invalid, x, ill), Inf)
})

test_that("usher3_hessian returns a valid Hessian matrix", {
  # Using a more realistic vector for the parameters
  th0 <- c(2e-2, 1.2, 0.175, 1.40, 0.368 * 0.01, 0.075 * 0.001, 0.917 * 0.1)
  
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
  paramVect <- c(2e-2, 1.2, 0.175, 1.40, 0.368 * 0.01, 0.075 * 0.001, 0.917 * 0.1)
  ageVect <- c(10, 20, 30, 40, 50)
  illVect <- c(0, 1, 0, 1, 0)
  nll <- nll_usher3_hessian_wrapper(paramVect, ageVect, illVect)
  expect_true(is.numeric(nll) && length(nll) == 1)
  expect_false(anyNA(nll), info = "Negative log-likelihood vector should not contain NA.")
})

# TODO: this must be tested with a parameter vector at a local minimum of the
# negative log-likelihood. This may require more values in the x vector.
#test_that("usher3_errors returns viabl3e standard errors, z-scores, and p-values", {
#  # Test case 1: Valid parameter values
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

test_that("temper_and_tune_usher3 returns viable results", {
  # Test case: Valid inputs
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)

  # Get results from the function
  results <- temper_and_tune_usher3(x = x, ill = ill)

  # Check if the results are a list and contain specific names
  expect_true(is.list(results))
  expect_true(all(c("obj_fun", "th0", "optional_inputs", "temper",
                    "th_temper", "th") %in% names(results)))

  # Extract transformed parameter vectors from results
  th_temper <- exp(results$th_temper)
  th <- exp(results$th)

  # Check that 'th' is not equal to 'th0'
  expect_false(all(th == th0), info = "'th' should not match 'th0'")

  # Check that 'th_temper' is not equal to 'th0'
  expect_false(all(th_temper == th0), 
               info = "'th_temper' should not match 'th0'")

  # Check that 'th' equals 'th_temper'
  expect_true(all(th == th_temper), 
              info = "'th' should match 'th_temper'")

  # Check that 'th' and 'th_temper' contain no NA values
  # TODO: actually implement the tuning step
  expect_false(anyNA(th), info = "'th' should not contain NA values")
  expect_false(anyNA(th_temper), 
               info = "'th_temper' should not contain NA values")
})
