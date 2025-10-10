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

# Tests for integration with usher3 functions

# Create a baseline parameter vector
b0 <- c(.175, 1.40, .368 * .01, 
        log(.917 * .1/(.075 * .001))/(.917 * .1), 
        .917 * .1)
th0 <- c(2e-2, 1.2, b0)

test_that("usher3_rho1 with x_cutoff = Inf matches original behavior", {
  k1 <- th0[1]
  b_siler <- th0[3:7]
  x <- c(0, 10, 25, 50)
  x0 <- 0
  
  # Calculate with x_cutoff = Inf
  rho1_inf <- usher3_rho1(x, k1, b_siler, x0, x_cutoff = Inf)
  
  # Calculate expected (original formula)
  f13 <- dsiler(x, b_siler, x0)
  S12_expected <- exp(-k1 * (x - x0))
  rho1_expected <- f13 * S12_expected
  
  expect_equal(rho1_inf, rho1_expected, tolerance = 1e-10)
})

test_that("usher3_rho1 with x_cutoff behaves correctly", {
  k1 <- 0.02
  b_siler <- th0[3:7]
  x0 <- 0
  x_cutoff <- 6
  
  # Before cutoff
  x_before <- c(3, 5, 6)
  rho1_before <- usher3_rho1(x_before, k1, b_siler, x0, x_cutoff)
  
  # After cutoff: rho1 should still be calculated, but with constant S12
  x_after <- c(10, 20, 50)
  rho1_after <- usher3_rho1(x_after, k1, b_siler, x0, x_cutoff)
  
  # Verify all are positive (valid density values)
  expect_true(all(rho1_before > 0))
  expect_true(all(rho1_after > 0))
  
  # Manual calculation for one point after cutoff
  x_test <- 20
  f13_test <- dsiler(x_test, b_siler, x0)
  S12_test <- exp(-k1 * (x_cutoff - x0))  # Should be constant after cutoff
  rho1_expected <- f13_test * S12_test
  rho1_actual <- usher3_rho1(x_test, k1, b_siler, x0, x_cutoff)
  
  expect_equal(rho1_actual, rho1_expected, tolerance = 1e-10)
})

test_that("usher3_integrand with x_cutoff = Inf matches original behavior", {
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  y <- c(0, 5, 10, 20)
  
  # Calculate with x_cutoff = Inf
  integrand_inf <- usher3_integrand(y, k1, k2, b_siler, x_cutoff = Inf)
  
  # Calculate expected (original formula)
  S12_expected <- exp(-k1 * y)
  S13_expected <- ssiler(y, b_siler)
  integrand_expected <- S12_expected * (S13_expected)^(1 - k2)
  
  expect_equal(integrand_inf, integrand_expected, tolerance = 1e-10)
})

test_that("usher3_integrand with x_cutoff behaves correctly", {
  k1 <- 0.02
  k2 <- 1.2
  b_siler <- th0[3:7]
  x_cutoff <- 6
  
  # Before cutoff: normal behavior
  y_before <- c(0, 3, 6)
  integrand_before <- usher3_integrand(y_before, k1, k2, b_siler, x_cutoff)
  
  # After cutoff: S12 should be constant
  y_after <- c(10, 20, 50)
  integrand_after <- usher3_integrand(y_after, k1, k2, b_siler, x_cutoff)
  
  # Verify all are positive
  expect_true(all(integrand_before > 0))
  expect_true(all(integrand_after > 0))
  
  # Manual check for point after cutoff
  y_test <- 20
  S12_test <- exp(-k1 * x_cutoff)  # Constant after cutoff
  S13_test <- ssiler(y_test, b_siler)
  integrand_expected <- S12_test * (S13_test)^(1 - k2)
  integrand_actual <- usher3_integrand(y_test, k1, k2, b_siler, x_cutoff)
  
  expect_equal(integrand_actual, integrand_expected, tolerance = 1e-10)
})

test_that("usher3_rho2 with x_cutoff = Inf matches original behavior", {
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  x <- c(10, 25, 50)
  x0 <- 0
  
  # Calculate with x_cutoff = Inf
  rho2_inf <- usher3_rho2(x, k1, k2, b_siler, x0, x_cutoff = Inf)
  
  # Calculate expected using manual integration
  # (This is a complex integral, so we just verify it's finite and positive)
  expect_true(all(is.finite(rho2_inf)))
  expect_true(all(rho2_inf > 0))
  
  # Also verify it matches rho2 calculated with the old approach by
  # checking the integral term calculation
  f13_0_x <- dsiler(x, b_siler)
  S13_0_x <- ssiler(x, b_siler)
  S12_0_x0 <- exp(-k1 * x0)
  S13_0_x0 <- ssiler(x0, b_siler)
  
  integralTerm <- rep(NA, length(x))
  for (ii in 1:length(x)) {
    integralTerm[ii] <- integrate(
      usher3_integrand, x0, x[ii],
      k1 = k1, k2 = k2, b_siler = b_siler, x_cutoff = Inf
    )$value
  }
  
  rho2_expected <- k1 * k2 * f13_0_x * S13_0_x^(k2 - 1) * 
                   integralTerm / S12_0_x0 / S13_0_x0
  
  expect_equal(rho2_inf, rho2_expected, tolerance = 1e-8)
})

test_that("usher3_rho2 with x_cutoff produces valid densities", {
  k1 <- 0.02
  k2 <- 1.2
  b_siler <- th0[3:7]
  x0 <- 0
  x_cutoff <- 6
  
  # Test at various ages
  x <- c(10, 20, 30, 50)
  rho2 <- usher3_rho2(x, k1, k2, b_siler, x0, x_cutoff)
  
  # Should all be finite and positive
  expect_true(all(is.finite(rho2)))
  expect_true(all(rho2 > 0))
})

# Tests for nll_usher3 with x_cutoff

test_that("nll_usher3 with x_cutoff = Inf matches original behavior", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  
  # Calculate with x_cutoff = Inf
  nll_inf <- nll_usher3(th0, x, ill, x0, x_cutoff = Inf)
  
  # Should be finite and positive
  expect_true(is.finite(nll_inf))
  expect_true(nll_inf > 0)
  
  # Compare with manually calculated likelihood components
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]
  
  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0, x_cutoff = Inf)
  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0, x_cutoff = Inf)
  
  ll_expected <- sum(log(rho1_wll)) + sum(log(rho2_ill))
  nll_expected <- -ll_expected
  
  expect_equal(nll_inf, nll_expected, tolerance = 1e-8)
})

test_that("nll_usher3 with x_cutoff produces valid likelihood", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cutoff <- 6
  
  # Calculate with finite x_cutoff
  nll_cutoff <- nll_usher3(th0, x, ill, x0, x_cutoff)
  
  # Should be finite and positive
  expect_true(is.finite(nll_cutoff))
  expect_true(nll_cutoff > 0)
  
  # Manual verification
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]
  
  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0, x_cutoff)
  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0, x_cutoff)
  
  ll_expected <- sum(log(rho1_wll)) + sum(log(rho2_ill))
  nll_expected <- -ll_expected
  
  expect_equal(nll_cutoff, nll_expected, tolerance = 1e-8)
})

test_that("nll_usher3 with x_cutoff handles NA illness indicators", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, NA, 0, 1, NA)
  x0 <- 0
  x_cutoff <- 6
  
  # Should handle NA values correctly
  nll_cutoff <- nll_usher3(th0, x, ill, x0, x_cutoff)
  
  expect_true(is.finite(nll_cutoff))
  expect_true(nll_cutoff > 0)
})

test_that("nll_usher3 with x_cutoff rejects negative parameters", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cutoff <- 6
  
  # Negative k1
  theta_neg_k1 <- th0
  theta_neg_k1[1] <- -0.01
  expect_equal(nll_usher3(theta_neg_k1, x, ill, x0, x_cutoff), Inf)
  
  # Negative k2
  theta_neg_k2 <- th0
  theta_neg_k2[2] <- -0.5
  expect_equal(nll_usher3(theta_neg_k2, x, ill, x0, x_cutoff), Inf)
  
  # Negative Siler parameter
  theta_neg_siler <- th0
  theta_neg_siler[3] <- -0.1
  expect_equal(nll_usher3(theta_neg_siler, x, ill, x0, x_cutoff), Inf)
})

test_that("nll_usher3 with x_cutoff supports Gompertz-Makeham", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cutoff <- 6
  
  # Use 5-parameter version (Gompertz-Makeham)
  th0_gompertz <- c(2e-2, 1.2, 0.368 * 0.01,
                    log(0.917 * 0.1 / (0.075 * 0.001)) / (0.917 * 0.1),
                    0.917 * 0.1)
  
  nll_gm <- nll_usher3(th0_gompertz, x, ill, x0, x_cutoff)
  
  expect_true(is.finite(nll_gm))
  expect_true(nll_gm > 0)
})

test_that("nll_usher3 likelihood changes with x_cutoff", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  
  # Calculate likelihoods with different cutoffs
  nll_inf <- nll_usher3(th0, x, ill, x0, x_cutoff = Inf)
  nll_6 <- nll_usher3(th0, x, ill, x0, x_cutoff = 6)
  nll_12 <- nll_usher3(th0, x, ill, x0, x_cutoff = 12)
  
  # All should be finite and positive
  expect_true(all(is.finite(c(nll_inf, nll_6, nll_12))))
  expect_true(all(c(nll_inf, nll_6, nll_12) > 0))
  
  # They should be different (different models)
  expect_false(isTRUE(all.equal(nll_inf, nll_6)))
  expect_false(isTRUE(all.equal(nll_inf, nll_12)))
  expect_false(isTRUE(all.equal(nll_6, nll_12)))
})

test_that("usher3_hessian with x_cutoff produces valid Hessian", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cutoff <- 6
  
  H <- usher3_hessian(th0, x, ill, x0, x_cutoff)
  
  # Should be a matrix
  expect_true(is.matrix(H))
  
  # Should have correct dimensions
  expect_equal(dim(H), c(7, 7))
  
  # Should not contain NA values
  expect_false(anyNA(H))
})

test_that("usher3_errors with x_cutoff runs without error", {
  x <- c(10, 20, 30, 40, 50, 15, 25, 35, 45)
  ill <- c(0, 1, 0, 1, 0, 1, 0, 1, 0)
  x0 <- 0
  x_cutoff <- 6
  
  # The main test is that this runs without throwing an error
  # Note: With small sample sizes, the Hessian may be near-singular
  # leading to numerical issues in error estimation
  expect_error(
    errors <- usher3_errors(th0, x, ill, x0, x_cutoff),
    NA  # NA means we expect NO error
  )
  
  # If it runs, it should return a data frame with the right structure
  errors <- usher3_errors(th0, x, ill, x0, x_cutoff)
  expect_true(is.data.frame(errors))
  expect_equal(nrow(errors), 7)
  expect_true(all(c("Estimate", "StandErr", "z", "pval", "against", "sideAdj") %in% names(errors)))
})
