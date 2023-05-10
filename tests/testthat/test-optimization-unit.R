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
})# Load necessary library
library(testthat)

# Generate a smaller normal distribution sample
set.seed(123)
x <- rnorm(50, mean = 5, sd = 2)

# Define negative log-likelihood function and its gradient for Gaussian
fn0 <- function(par, x) {
  mu <- par[1]
  sigma <- par[2]
  n <- length(x)
  return(0.5 * n * log(2 * pi * sigma^2) + sum((x - mu)^2) / (2 * sigma^2))
}

gr0 <- function(par, x) {
  mu <- par[1]
  sigma <- par[2]
  n <- length(x)
  dmu <- -sum(x - mu) / sigma^2
  dsigma <- n / sigma - sum((x - mu)^2) / sigma^3
  return(c(dmu, dsigma))
}

# Unit tests for gradient descent. Functionality is tested in
# test-optimization-functional.R. The units tests are representative but not
# comprehensive.
test_that("gradient descent returns a list", {
  result <- gradient_descent(th0 = c(mu = 0, sigma = 1), fn0 = fn0, gr0 = gr0, x = x)
  expect_is(result, "list")
})

test_that("gradient descent returns the correct keys", {
  result <- gradient_descent(th0 = c(mu = 0, sigma = 1), fn0 = fn0, gr0 = gr0, x = x)
  expect_setequal(names(result), c("par", "value", "feval", "inputs"))
})

test_that("gradient descent works with rescale=TRUE", {
  result <- gradient_descent(th0 = c(mu = 0, sigma = 1), fn0 = fn0, gr0 = gr0, rescale = TRUE, x = x)
  expect_is(result, "list")
})

test_that("gradient descent works with verbose=TRUE", {
  result <- gradient_descent(th0 = c(mu = 0, sigma = 1), fn0 = fn0, gr0 = gr0, rescale = TRUE, x = x)
  expect_is(result, "list")
})

test_that("gradient descent halts at maxiter", {
  result <- gradient_descent(th0 = c(mu = 0, sigma = 1), fn0 = fn0, gr0 = gr0, maxiter = 10, x = x)
  expect_equal(result$feval, 10)
})