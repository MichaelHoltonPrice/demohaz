library(testthat)
library(demohaz)

# a0 is the baseline parameter vector used to test Siler functions. It is from
# Gage and Dyke 1986, Table 2, Level 15.
a0 <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)
# Convert traditional parameterization to demohaz parameterization
b0 <- trad_to_demohaz_siler_param(a0)

test_that("qsiler followed by psiler yields the same quantiles", {

  
  # Check the functionality of the Siler code. These tests should be run prior to
  # pushing back any code to github, but only the unit tests need to be run when
  # actively doing development (they are much quicker).
  
  # Check that one recovers the quantiles when going from quantiles to x then back
  # to quantiles. Do this for both x0=0 (the default) and x0=10.
  qin <- seq(0, .999, by = .001) # quantiles in
  x <- qsiler(qin, a0) # x from quantiles
  qout <- psiler(x, a0) # quantiles from x
  expect_equal(
    qin,
    qout
  )
  
  x0 <- 10
  qin <- seq(0, .999, by = .001) # quantiles in
  x <- qsiler(qin, a0, x0) # x from quantiles
  qout <- psiler(x, a0, x0) # quantiles from x
  expect_equal(
    qin,
    qout
  )
})

test_that("fit_siler works on simulated data", {
  # Check that when we sample then fit a Siler model the fit gives the original
  # parameter vector (to a tolerance of 5e-3). Do this for both x0=0 (the default)
  # and x0=1 (use x0=0.1 rather than x0=10 so that the first two parameters can be
  # accurately fit).
  set.seed(291542) # from random.org
  N <- 5000
  xsamp <- rsiler(N, b0)
  # jitter starting vector for fit so that we don't start already at the optimum.
  b1 <- b0 * runif(5, min = .9, max = 1.1)
  silerFit <- fit_siler(xsamp, b0=b1)
  expect_equal(
    b0,
    silerFit$b, tol = 5e-3
  )
  
  set.seed(795075) # from random.org
  x0 <- 0.1
  N <- 5000
  xsamp <- rsiler(N, b0, x0)
  # jitter starting vector for fit so that we don't start already at the optimum.
  b1 <- b0 * runif(5, min = .9, max = 1.1)
  silerFit <- fit_siler(xsamp,b0=b1, x0=x0)
  expect_equal(
    b0,
    silerFit$b, tol = 5e-3
  )
})

# This test takes quite long
test_that("rsamp yields the correct empirical CDF", {
  # Check that random samples made by rsiler match the cumulative distribution
  # function calcualted by psiler; to achieve a tolerance of 1e-3 on the check,
  # one million samples are needed. Do this for both x0=0 (the default) and
  # x0=10.
  set.seed(825201) # from random.org
  N <- 1000000
  x <- seq(0, 120, by = .01)
  cdf <- psiler(x, a0)
  xsamp <- rsiler(N, a0)
  expect_equal(
    cdf,
    ecdf(xsamp)(x), tolerance = 1e-3
  )
  
  set.seed(400532) # from random.org
  x0 <- 10
  N <- 1000000
  x <- seq(x0, 120, by = .01)
  cdf <- psiler(x, a0, x0)
  xsamp <- rsiler(N, a0, x0)
  expect_equal(
    cdf,
    ecdf(xsamp)(x), tolerance = 1e-3
  )
})