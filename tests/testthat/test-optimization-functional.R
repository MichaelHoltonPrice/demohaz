library(testthat)
library(demohaz)

# Generate a larger normal distribution sample
set.seed(123)
x <- rnorm(10000, mean = 5, sd = 2)

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

# Define the functional tests
test_that("gradient descent recovers correct parameters with rescale = FALSE", {
  result <- gradient_descent(th0 = c(mu = 0.5, sigma = 1), fn0 = fn0, gr0 = gr0, rescale = FALSE, x = x)
  expect_equal(as.numeric(result$par[1]), 5, tolerance = .01)
  expect_equal(as.numeric(result$par[2]), 2, tolerance = .01)
})

# Define the functional tests
test_that("gradient descent recovers correct parameters with rescale = TRUE", {
  result <- gradient_descent(th0 = c(mu = 0.5, sigma = 1), fn0 = fn0, gr0 = gr0, rescale = TRUE, x = x)
  expect_equal(as.numeric(result$par[1])*0.5, 5, tolerance = .01)
  expect_equal(as.numeric(result$par[2]), 2, tolerance = .01)
})