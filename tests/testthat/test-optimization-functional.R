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

# A check on the statistical identifiability behavior could, conceptually, be
# placed in either the unit tests or functional tests. I have chosen to place
# it here. The following data yield b[5] equal to about 1, which fails the
# check.
test_that("Warning message when b[5] > 0.5", {
  xvalues <- c(0.5, 1.5, 5.5, 10.5, 15.5, 20.5, 25.5, 30.5, 35.5, 40.5, 45.5,
               50.5, 55.5, 60.5, 65.5, 70.5)
  xcounts <- c(948, 759, 213, 394, 204, 284, 326, 202, 255, 120, 100, 169, 198,
               279, 88, 883)
  # Create x from these values and counts
  x <- c()
  for (n in 1:length(xvalues)) {
    x <- c(x, rep(xvalues[n], xcounts[n]))
  }
  set.seed(1234)
  expect_warning(output <- fit_siler(x), 
                 "b[5] is greater than 0.5, which indicates a potential issue with statistical identifiability. Not doing gradient descent. The user should investigate this further.",
                 fixed=TRUE)
})

