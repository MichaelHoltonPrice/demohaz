library(testthat)
library(demohaz)

# Create a baseline parameter vector used across tests
# a0 is the baseline parameter vector used to test Siler functions. It is from
# Gage and Dyke 1986, Table 2, Level 15.
a0 <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)
# Convert traditional parameterization to demohaz parameterization
b0 <- trad_to_demohaz_siler_param(a0)

test_that("hsiler function returns correct hazard values", {
  # Test cases with various age inputs
  expect_equal(hsiler(0, b0), a0[1] * exp(-a0[2] * 0) + a0[3] + a0[4] * exp(a0[5] * 0), tolerance = 1e-8)
  expect_equal(hsiler(10, b0), a0[1] * exp(-a0[2] * 10) + a0[3] + a0[4] * exp(a0[5] * 10), tolerance = 1e-8)
  expect_equal(hsiler(25, b0), a0[1] * exp(-a0[2] * 25) + a0[3] + a0[4] * exp(a0[5] * 25), tolerance = 1e-8)
  expect_equal(hsiler(50, b0), a0[1] * exp(-a0[2] * 50) + a0[3] + a0[4] * exp(a0[5] * 50), tolerance = 1e-8)
})

test_that("chsiler function returns correct cumulative hazard values", {
  # Test cases with various age inputs and optional starting ages
  expect_equal(chsiler(0, b0), 0, tolerance = 1e-8)
  expect_equal(chsiler(10, b0), -b0[1] / b0[2] * (exp(-b0[2] * 10) - 1) +
                 b0[3] * 10 + exp(b0[5] * (10 - b0[4])) - exp(-b0[5]*b0[4]), tolerance = 1e-8)
  expect_equal(chsiler(25, b0, x0 = 10), chsiler(25, b0) - chsiler(10, b0), tolerance = 1e-8)
  expect_equal(chsiler(50, b0, x0 = 25), chsiler(50, b0) - chsiler(25, b0), tolerance = 1e-8)
})

test_that("ssiler function returns correct survival values", {
  # Test cases with various age inputs and optional starting ages
  expect_equal(ssiler(0, b0), 1, tolerance = 1e-8)
  expect_equal(ssiler(10, b0), exp(-chsiler(10, b0)), tolerance = 1e-8)
  expect_equal(ssiler(25, b0, x0 = 10), exp(-chsiler(25, b0, x0 = 10)), tolerance = 1e-8)
  expect_equal(ssiler(50, b0, x0 = 25), exp(-chsiler(50, b0, x0 = 25)), tolerance = 1e-8)
})

test_that("dsiler function returns correct probability density values", {
  # Test cases with various age inputs and optional starting ages
  expect_equal(dsiler(0, b0), hsiler(0, b0), tolerance = 1e-8)
  expect_equal(dsiler(10, b0), hsiler(10, b0) * ssiler(10, b0), tolerance = 1e-8)
  expect_equal(dsiler(25, b0, x0 = 10), hsiler(25, b0) * ssiler(25, b0, x0 = 10), tolerance = 1e-8)
  expect_equal(dsiler(50, b0, x0 = 25), hsiler(50, b0) * ssiler(50, b0, x0 = 25), tolerance = 1e-8)
})

test_that("psiler function returns correct cumulative distribution values", {
  # Test cases with various age inputs and optional starting ages
  expect_equal(psiler(0, b0), 1 - ssiler(0, b0), tolerance = 1e-8)
  expect_equal(psiler(10, b0), 1 - ssiler(10, b0), tolerance = 1e-8)
  expect_equal(psiler(25, b0, x0 = 10), 1 - ssiler(25, b0, x0 = 10), tolerance = 1e-8)
  expect_equal(psiler(50, b0, x0 = 25), 1 - ssiler(50, b0, x0 = 25), tolerance = 1e-8)
})

# functional tests for qsiler are in test-siler-function.R
test_that("qsiler does not throw errors for various inputs", {
  # Test case 1: quantiles in the middle of the range
  expect_error(qsiler(c(0.25, 0.5, 0.75), b0), NA)

  # Test case 2: extreme quantiles
  expect_error(qsiler(c(0.001, 0.999), b0), NA)

  # Test case 3: handling x0 (starting age) correctly
  expect_error(qsiler(0.5, b0, 20), NA)
})

test_that("qsiler throws errors for invalid inputs", {
  # Test case 1: negative quantiles
  expect_error(qsiler(c(-0.1, 0.5, 0.75), b0),
                      "Quantiles should be between 0 and 1",
                      fixed=TRUE)

  # Test case 2: quantiles greater than 1
  expect_error(qsiler(c(0.5, 1.1), b0),
                      "Quantiles should be between 0 and 1",
                      fixed=TRUE)
})

test_that("rsiler returns expected types", {
  set.seed(123)
  N <- 100
  b <- c(0.1, 0.1, 0.1, 0.1, 0.1)
  x0 <- 0

  result <- rsiler(N, b, x0)

  expect_type(result, "double")
  expect_length(result, N)
})

test_that("rsiler works with different x0", {
  set.seed(123)
  N <- 100
  b <- c(0.1, 0.1, 0.1, 0.1, 0.1)
  x0 <- 10

  result <- rsiler(N, b, x0)

  expect_type(result, "double")
  expect_length(result, N)
})

test_that("rsiler returns different results with different seeds", {
  N <- 100
  b <- c(0.1, 0.1, 0.1, 0.1, 0.1)
  x0 <- 0

  set.seed(123)
  result1 <- rsiler(N, b, x0)

  set.seed(456)
  result2 <- rsiler(N, b, x0)

  expect_false(identical(result1, result2))
})

test_that("hesiler returns a matrix with the correct dimensions", {
  b <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x <- seq(1, 5)
  x0 <- 0

  result <- hesiler(b, x, x0)

  expect_type(result, "double")
  expect_equal(dim(result), c(5, 5))
})

test_that("hesiler returns a symmetric matrix", {
  b <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x <- seq(1, 5)
  x0 <- 0

  result <- hesiler(b, x, x0)

  expect_true(all(result == t(result)))
})

test_that("hesiler works with different x0", {
  b <- c(0.1, 0.2, 0.3, 0.4, 0.5)
  x <- seq(1, 5)
  x0 <- 10

  result <- hesiler(b, x, x0)

  expect_type(result, "double")
  expect_equal(dim(result), c(5, 5))
})

test_that("grad_hsiler_vect returns the correct gradient matrix", {
  set.seed(123)
  b <- c(1, 0.5, 0.3, 2, 0.1)
  x <- runif(100, 0, 10)
  
  grad_matrix <- grad_hsiler_vect(x, b)
  
  # Check the output dimensions
  expect_is(grad_matrix, "matrix")
  expect_equal(dim(grad_matrix), c(length(x), 5))

  # Calculate the gradients using analytical expressions
  exp1 <- exp(-b[2]*x)
  exp3 <- exp(-b[5]*(x-b[4]))
  grad_analytical <- matrix(NA, length(x), 5)
  grad_analytical[,1] <- exp1
  grad_analytical[,2] <- -b[1]*x*exp1
  grad_analytical[,3] <- 1
  grad_analytical[,4] <- -b[5]^2 * exp3
  grad_analytical[,5] <- exp3 + b[5]*(x-b[4])*exp3

  # Compare the result with the analytical expressions
  expect_equal(grad_matrix, grad_analytical, tolerance = 1e-6)
})

test_that("hsiler_vect returns the correct hazard matrix", {
  set.seed(123)
  b <- c(1, 0.5, 0.3, 2, 0.1)
  x <- runif(100, 0, 10)
  
  hazard_matrix <- hsiler_vect(x, b)
  
  # Check the output dimensions
  expect_is(hazard_matrix, "matrix")
  expect_equal(dim(hazard_matrix), c(length(x), 5))

  # Calculate the hazard using hsiler function
  hazard <- hsiler(x, b)

  # Compare the result with the replicated hazards
  for (i in 1:5) {
    expect_equal(hazard_matrix[, i], hazard, tolerance = 1e-6)
  }
})

test_that("grad_chsiler_vect returns the correct gradient matrix", {
  set.seed(123)
  b <- c(1, 0.5, 0.3, 2, 0.1)
  x <- runif(100, 0, 10)
  x0 <- .25

  # Calculate the gradient matrix
  grad_matrix <- grad_chsiler_vect(x, b, x0)

  # Check the output dimensions
  expect_is(grad_matrix, "matrix")
  expect_equal(dim(grad_matrix), c(length(x), 5))

  # Calculate the gradient matrix using a numerical method
  num_grad_matrix <- matrix(NA, length(x), 5)
  delta <- 1e-6
  for (i in 1:5) {
    b_up <- b
    b_down <- b
    b_up[i] <- b[i] + delta
    b_down[i] <- b[i] - delta

    ch_up <- chsiler(x, b_up) - chsiler(x0, b_up)
    ch_down <- chsiler(x, b_down) - chsiler(x0, b_down)

    num_grad_matrix[, i] <- (ch_up - ch_down) / (2 * delta)
  }

  # Compare the result with the numerical gradient matrix
  for (i in 1:5) {
    expect_equal(grad_matrix[, i], num_grad_matrix[, i], tolerance = 1e-6)
  }
})

test_grads <- function(b, x, x0 = 0, tol = 1e-6) {
  xtable <- table(x)
  xvalues <- as.numeric(names(xtable))
  xcounts <- as.numeric(xtable)

  analytical_grad <- gradnllsiler(b, x, x0)
  fast_grad <- fast_gradnllsiler(b, x, xcounts, x0)

  max_diff <- max(abs(analytical_grad - fast_grad))
  return(max_diff <= tol)
}

# TODO: this is more a functional test. Move? It really doesn't matter.
test_that("fast and analytical gradient calculations match", {
  # Generate random data for testing
  set.seed(10000)
  x <- runif(100, 0, 100)
  b_for_grad <- c(0.01, 0.02, 0.03, 10, 0.05)
  
  # Test the gradient implementations
  expect_true(test_grads(b_for_grad, x))
  
  # Test with different starting age
  x0 <- 5
  expect_true(test_grads(b_for_grad, x, x0))

})

test_that("gradnllsiler does not throw errors for valid inputs", {
  # Test case 1: valid b vector and x vector with x0 = 0
  expect_error(gradnllsiler(c(1, 2, 3, 4, 5), c(1, 2, 3, 4, 5), 0), NA)

  # Test case 2: valid b vector and x vector with x0 > 0
  expect_error(gradnllsiler(c(1, 2, 3, 4, 5), c(1, 2, 3, 4, 5), 1), NA)
})

# functional tests are in test-siler-function.R
test_that("fit_siler does not throw errors for valid inputs", {
  set.seed(100000)
  x <- rsiler(20, b0)
  x <- c(0,x)
  # Test case 1: Default parameters
  expect_error(fit_siler(x), NA)

  # Test case 2: Non-default x0
  expect_error(fit_siler(x, x0 = 1), NA)

  # Test case 3: Non-default b0
  expect_error(fit_siler(x, b0 = c(0.1, 0.2, 0.3, 0.4, 0.5)), NA)

  # Test case 4: Setting calc_hessian to TRUE
  expect_error(fit_siler(x, calc_hessian = TRUE), NA)

  # Test case 5: Setting verbose to TRUE
  # Send output to a temporary file so that test progress is not obscured
  sink(tempfile())
  expect_error(fit_siler(x, verbose = TRUE), NA)

  # Test case 6: Setting lr to a non-default value
  expect_error(fit_siler(x, lr = 1e-4), NA)

  # Test case 7: Setting lr to a non-default value
  expect_error(fit_siler(x, lr = 1e-4), NA)

  # Test case 8: Setting show_plot to TRUE
  expect_error(fit_siler(x, verbose=TRUE, show_plot = TRUE), NA)
  sink()
})

# Test trad_to_demohaz_siler_param function
test_that("trad_to_demohaz_siler_param is correct", {
  # Test case 1
  a1 <- c(1, 2, 3, 4, 5)
  b1 <- c(1, 2, 3, log(5/4) / 5, 5)
  expect_equal(trad_to_demohaz_siler_param(a1), b1)
  
  # Test case 2
  a2 <- c(2, 1, 5, 3, 4)
  b2 <- c(2, 1, 5, log(4/3) / 4, 4)
  expect_equal(trad_to_demohaz_siler_param(a2), b2)
})

# Test demohaz_to_trad_siler_param function
test_that("demohaz_to_trad_siler_param is correct", {
  # Test case 1
  b1 <- c(1, 2, 3, 1, 5)
  a1 <- c(1, 2, 3, 5 * exp(-1 * 5), 5)
  expect_equal(demohaz_to_trad_siler_param(b1), a1)
  
  # Test case 2
  b2 <- c(2, 1, 5, 0.5, 4)
  a2 <- c(2, 1, 5, 4 * exp(-0.5 * 4), 4)
  expect_equal(demohaz_to_trad_siler_param(b2), a2)
})

# Test conversion functions are inverses of each other
test_that("trad_to_demohaz_siler_param and demohaz_to_trad_siler_param are inverses", {
  # Test case 1
  a1 <- c(1, 2, 3, 4, 5)
  expect_equal(demohaz_to_trad_siler_param(trad_to_demohaz_siler_param(a1)), a1)
  
  # Test case 2
  b1 <- c(2, 1, 5, 0.5, 4)
  expect_equal(trad_to_demohaz_siler_param(demohaz_to_trad_siler_param(b1)), b1)
})

test_that("check_siler_inputs does not throw errors for valid inputs", {
  # Test case 1: valid b vector and x0
  expect_error(check_siler_inputs(c(1, 2, 3, 4, 5), 0), NA)

  # Test case 2: b vector with 0
  expect_error(check_siler_inputs(c(0, 2, 3, 4, 5), 0), NA)
})

test_that("check_siler_inputs throws errors for invalid inputs", {
  # Test case 1: b vector with negative numbers
  expect_error(check_siler_inputs(c(-1, 2, 3, 4, 5), 0),
               "No element of b can be negative",
               fixed=TRUE)

  # Test case 2: x0 is negative
  expect_error(check_siler_inputs(c(1, 2, 3, 4, 5), -1),
               "x0 cannot be negative",
               fixed=TRUE)
})
