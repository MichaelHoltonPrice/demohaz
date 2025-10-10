library(testthat)
library(demohaz)

test_that("H12 matches constant hazard formula when x_cutoff = Inf", {
  k1 <- 0.02
  x0 <- 0
  x <- c(0, 5, 10, 20, 50)
  
  # When x_cutoff = Inf, should match k1 * (x - x0)
  expected <- k1 * (x - x0)
  result <- H12(x, k1, x0, x_cutoff = Inf)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("H12 stops accumulating after x_cutoff", {
  k1 <- 0.02
  x0 <- 0
  x_cutoff <- 6
  
  # Before cutoff: normal accumulation
  x_before <- c(0, 3, 5, 6)
  expected_before <- k1 * (x_before - x0)
  result_before <- H12(x_before, k1, x0, x_cutoff)
  expect_equal(result_before, expected_before, tolerance = 1e-10)
  
  # After cutoff: hazard stays constant at cutoff value
  x_after <- c(7, 10, 20, 50)
  expected_after <- rep(k1 * (x_cutoff - x0), length(x_after))
  result_after <- H12(x_after, k1, x0, x_cutoff)
  expect_equal(result_after, expected_after, tolerance = 1e-10)
})

test_that("H12 handles x0 > 0 correctly", {
  k1 <- 0.02
  x0 <- 5
  x_cutoff <- 10
  
  # Before x0: hazard should be zero
  x_before_x0 <- c(0, 2, 4)
  expected_before_x0 <- rep(0, length(x_before_x0))
  result_before_x0 <- H12(x_before_x0, k1, x0, x_cutoff)
  expect_equal(result_before_x0, expected_before_x0, tolerance = 1e-10)
  
  # Between x0 and cutoff
  x_between <- c(5, 7, 10)
  expected_between <- k1 * (x_between - x0)
  result_between <- H12(x_between, k1, x0, x_cutoff)
  expect_equal(result_between, expected_between, tolerance = 1e-10)
  
  # After cutoff
  x_after <- c(15, 20, 50)
  expected_after <- rep(k1 * (x_cutoff - x0), length(x_after))
  result_after <- H12(x_after, k1, x0, x_cutoff)
  expect_equal(result_after, expected_after, tolerance = 1e-10)
})

test_that("H12 validates input parameters", {
  expect_error(H12(10, k1 = -0.01, x0 = 0), "k1 cannot be negative")
  expect_error(H12(10, k1 = 0.02, x0 = -1), "x0 cannot be negative")
  expect_error(H12(10, k1 = 0.02, x0 = 5, x_cutoff = 3), 
               "x_cutoff must be greater than or equal to x0")
})

test_that("H12 is vectorized over x", {
  k1 <- 0.02
  x0 <- 0
  x_cutoff <- 6
  x <- c(0, 2, 5, 6, 8, 10, 20)
  
  result <- H12(x, k1, x0, x_cutoff)
  
  expect_equal(length(result), length(x))
  expect_true(is.numeric(result))
  
  # Manually calculate expected values
  expected <- k1 * pmin(x - x0, x_cutoff - x0)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("S12 matches exponential of constant hazard when x_cutoff = Inf", {
  k1 <- 0.02
  x0 <- 0
  x <- c(0, 5, 10, 20, 50)
  
  # When x_cutoff = Inf, should match exp(-k1 * (x - x0))
  expected <- exp(-k1 * (x - x0))
  result <- S12(x, k1, x0, x_cutoff = Inf)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("S12 stays constant after x_cutoff", {
  k1 <- 0.02
  x0 <- 0
  x_cutoff <- 6
  
  # Before cutoff: normal exponential decay
  x_before <- c(0, 3, 6)
  expected_before <- exp(-k1 * (x_before - x0))
  result_before <- S12(x_before, k1, x0, x_cutoff)
  expect_equal(result_before, expected_before, tolerance = 1e-10)
  
  # After cutoff: survival stays constant (no more transitions possible)
  x_after <- c(7, 10, 20, 50)
  expected_after <- rep(exp(-k1 * (x_cutoff - x0)), length(x_after))
  result_after <- S12(x_after, k1, x0, x_cutoff)
  expect_equal(result_after, expected_after, tolerance = 1e-10)
})

test_that("S12 is between 0 and 1", {
  k1 <- 0.02
  x0 <- 0
  x_cutoff <- 6
  x <- seq(0, 50, by = 1)
  
  result <- S12(x, k1, x0, x_cutoff)
  
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("S12 is non-increasing", {
  k1 <- 0.02
  x0 <- 0
  x_cutoff <- 6
  x <- seq(0, 50, by = 0.1)
  
  result <- S12(x, k1, x0, x_cutoff)
  
  # Check that differences are non-positive (non-increasing)
  diffs <- diff(result)
  expect_true(all(diffs <= 1e-10))  # Allow for numerical precision
})

test_that("LEH example: 6-year cutoff", {
  # Linear Enamel Hypoplasia can only occur before age 6
  k1 <- 0.02  # Hazard rate for developing LEH before age 6
  x0 <- 0
  x_cutoff <- 6
  
  # At age 3: can still develop LEH
  H_age3 <- H12(3, k1, x0, x_cutoff)
  expect_equal(H_age3, 0.02 * 3, tolerance = 1e-10)
  
  S_age3 <- S12(3, k1, x0, x_cutoff)
  expect_equal(S_age3, exp(-0.02 * 3), tolerance = 1e-10)
  
  # At age 6: last moment LEH can occur
  H_age6 <- H12(6, k1, x0, x_cutoff)
  expect_equal(H_age6, 0.02 * 6, tolerance = 1e-10)
  
  S_age6 <- S12(6, k1, x0, x_cutoff)
  expect_equal(S_age6, exp(-0.02 * 6), tolerance = 1e-10)
  
  # At age 10: LEH can no longer occur (hazard stopped at age 6)
  H_age10 <- H12(10, k1, x0, x_cutoff)
  expect_equal(H_age10, 0.02 * 6, tolerance = 1e-10)  # Same as age 6
  
  S_age10 <- S12(10, k1, x0, x_cutoff)
  expect_equal(S_age10, exp(-0.02 * 6), tolerance = 1e-10)  # Same as age 6
  
  # At age 50: still same as age 6
  H_age50 <- H12(50, k1, x0, x_cutoff)
  S_age50 <- S12(50, k1, x0, x_cutoff)
  expect_equal(H_age50, H_age6, tolerance = 1e-10)
  expect_equal(S_age50, S_age6, tolerance = 1e-10)
})

