library(testthat)
library(demohaz)

test_that("H12 matches constant hazard formula when x_cut = Inf", {
  k1 <- 0.02
  x0 <- 0
  x <- c(0, 5, 10, 20, 50)
  
  # When x_cut = Inf, should match k1 * (x - x0)
  expected <- k1 * (x - x0)
  result <- H12(x, k1, x0, x_cut = Inf)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("H12 stops accumulating after x_cut", {
  k1 <- 0.02
  x0 <- 0
  x_cut <- 6
  
  # Before cutoff: normal accumulation
  x_before <- c(0, 3, 5, 6)
  expected_before <- k1 * (x_before - x0)
  result_before <- H12(x_before, k1, x0, x_cut)
  expect_equal(result_before, expected_before, tolerance = 1e-10)
  
  # After cutoff: hazard stays constant at cutoff value
  x_after <- c(7, 10, 20, 50)
  expected_after <- rep(k1 * (x_cut - x0), length(x_after))
  result_after <- H12(x_after, k1, x0, x_cut)
  expect_equal(result_after, expected_after, tolerance = 1e-10)
})

test_that("H12 handles x0 > 0 correctly", {
  k1 <- 0.02
  x0 <- 5
  x_cut <- 10
  
  # Before x0: hazard should be zero
  x_before_x0 <- c(0, 2, 4)
  expected_before_x0 <- rep(0, length(x_before_x0))
  result_before_x0 <- H12(x_before_x0, k1, x0, x_cut)
  expect_equal(result_before_x0, expected_before_x0, tolerance = 1e-10)
  
  # Between x0 and cutoff
  x_between <- c(5, 7, 10)
  expected_between <- k1 * (x_between - x0)
  result_between <- H12(x_between, k1, x0, x_cut)
  expect_equal(result_between, expected_between, tolerance = 1e-10)
  
  # After cutoff
  x_after <- c(15, 20, 50)
  expected_after <- rep(k1 * (x_cut - x0), length(x_after))
  result_after <- H12(x_after, k1, x0, x_cut)
  expect_equal(result_after, expected_after, tolerance = 1e-10)
})

test_that("H12 validates input parameters", {
  expect_error(H12(10, k1 = -0.01, x0 = 0), "k1 cannot be negative")
  expect_error(H12(10, k1 = 0.02, x0 = -1), "x0 cannot be negative")
  expect_error(H12(10, k1 = 0.02, x0 = 0, x_cut = 0), 
               "x_cut must be positive")
  expect_error(H12(10, k1 = 0.02, x0 = 0, x_cut = -1), 
               "x_cut must be positive")
})

test_that("H12 returns 0 when x0 >= x_cut", {
  k1 <- 0.02
  x_cut <- 6
  x0 <- 7  # Past the cutoff
  x <- c(8, 10, 20, 50)
  
  # All hazard should be 0 when starting past the cutoff
  result <- H12(x, k1, x0, x_cut)
  expect_equal(result, rep(0, length(x)), tolerance = 1e-10)
  
  # Also test when x0 exactly equals x_cut
  x0 <- 6
  result2 <- H12(x, k1, x0, x_cut)
  expect_equal(result2, rep(0, length(x)), tolerance = 1e-10)
})

test_that("H12 is vectorized over x", {
  k1 <- 0.02
  x0 <- 0
  x_cut <- 6
  x <- c(0, 2, 5, 6, 8, 10, 20)
  
  result <- H12(x, k1, x0, x_cut)
  
  expect_equal(length(result), length(x))
  expect_true(is.numeric(result))
  
  # Manually calculate expected values
  expected <- k1 * pmin(x - x0, x_cut - x0)
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("S12 matches exponential of constant hazard when x_cut = Inf", {
  k1 <- 0.02
  x0 <- 0
  x <- c(0, 5, 10, 20, 50)
  
  # When x_cut = Inf, should match exp(-k1 * (x - x0))
  expected <- exp(-k1 * (x - x0))
  result <- S12(x, k1, x0, x_cut = Inf)
  
  expect_equal(result, expected, tolerance = 1e-10)
})

test_that("S12 stays constant after x_cut", {
  k1 <- 0.02
  x0 <- 0
  x_cut <- 6
  
  # Before cutoff: normal exponential decay
  x_before <- c(0, 3, 6)
  expected_before <- exp(-k1 * (x_before - x0))
  result_before <- S12(x_before, k1, x0, x_cut)
  expect_equal(result_before, expected_before, tolerance = 1e-10)
  
  # After cutoff: survival stays constant (no more transitions possible)
  x_after <- c(7, 10, 20, 50)
  expected_after <- rep(exp(-k1 * (x_cut - x0)), length(x_after))
  result_after <- S12(x_after, k1, x0, x_cut)
  expect_equal(result_after, expected_after, tolerance = 1e-10)
})

test_that("S12 is between 0 and 1", {
  k1 <- 0.02
  x0 <- 0
  x_cut <- 6
  x <- seq(0, 50, by = 1)
  
  result <- S12(x, k1, x0, x_cut)
  
  expect_true(all(result >= 0))
  expect_true(all(result <= 1))
})

test_that("S12 is non-increasing", {
  k1 <- 0.02
  x0 <- 0
  x_cut <- 6
  x <- seq(0, 50, by = 0.1)
  
  result <- S12(x, k1, x0, x_cut)
  
  # Check that differences are non-positive (non-increasing)
  diffs <- diff(result)
  expect_true(all(diffs <= 1e-10))  # Allow for numerical precision
})

test_that("LEH example: 6-year cutoff", {
  # Linear Enamel Hypoplasia can only occur before age 6
  k1 <- 0.02  # Hazard rate for developing LEH before age 6
  x0 <- 0
  x_cut <- 6
  
  # At age 3: can still develop LEH
  H_age3 <- H12(3, k1, x0, x_cut)
  expect_equal(H_age3, 0.02 * 3, tolerance = 1e-10)
  
  S_age3 <- S12(3, k1, x0, x_cut)
  expect_equal(S_age3, exp(-0.02 * 3), tolerance = 1e-10)
  
  # At age 6: last moment LEH can occur
  H_age6 <- H12(6, k1, x0, x_cut)
  expect_equal(H_age6, 0.02 * 6, tolerance = 1e-10)
  
  S_age6 <- S12(6, k1, x0, x_cut)
  expect_equal(S_age6, exp(-0.02 * 6), tolerance = 1e-10)
  
  # At age 10: LEH can no longer occur (hazard stopped at age 6)
  H_age10 <- H12(10, k1, x0, x_cut)
  expect_equal(H_age10, 0.02 * 6, tolerance = 1e-10)  # Same as age 6
  
  S_age10 <- S12(10, k1, x0, x_cut)
  expect_equal(S_age10, exp(-0.02 * 6), tolerance = 1e-10)  # Same as age 6
  
  # At age 50: still same as age 6
  H_age50 <- H12(50, k1, x0, x_cut)
  S_age50 <- S12(50, k1, x0, x_cut)
  expect_equal(H_age50, H_age6, tolerance = 1e-10)
  expect_equal(S_age50, S_age6, tolerance = 1e-10)
})

# Tests for integration with usher3 functions

# Create a baseline parameter vector
b0 <- c(.175, 1.40, .368 * .01, 
        log(.917 * .1/(.075 * .001))/(.917 * .1), 
        .917 * .1)
th0 <- c(2e-2, 1.2, b0)

test_that("usher3_rho1 with x_cut = Inf is calculated correctly", {
  k1 <- th0[1]
  b_siler <- th0[3:7]
  x <- c(0, 10, 25, 50)
  x0 <- 0
  w0 <- c(1, 0)  # Everyone well at x0
  
  # Calculate with x_cut = Inf
  rho1_inf <- usher3_rho1(x, k1, b_siler, x0, x_cut = Inf, w = w0)
  
  # Calculate expected
  f13 <- dsiler(x, b_siler, x0)
  S12_expected <- exp(-k1 * (x - x0))
  rho1_expected <- f13 * S12_expected
  
  expect_equal(rho1_inf, rho1_expected, tolerance = 1e-10)
})

test_that("usher3_rho1 with x_cut behaves correctly", {
  k1 <- 0.02
  b_siler <- th0[3:7]
  x0 <- 0
  x_cut <- 6
  w0 <- c(1, 0)  # Everyone well at x0
  
  # Before cutoff
  x_before <- c(3, 5, 6)
  rho1_before <- usher3_rho1(x_before, k1, b_siler, x0, x_cut, w = w0)
  
  # After cutoff: rho1 should still be calculated, but with constant S12
  x_after <- c(10, 20, 50)
  rho1_after <- usher3_rho1(x_after, k1, b_siler, x0, x_cut, w = w0)
  
  # Verify all are positive (valid density values)
  expect_true(all(rho1_before > 0))
  expect_true(all(rho1_after > 0))
  
  # Manual calculation for one point after cutoff
  x_test <- 20
  f13_test <- dsiler(x_test, b_siler, x0)
  S12_test <- exp(-k1 * (x_cut - x0))  # Should be constant after cutoff
  rho1_expected <- f13_test * S12_test
  rho1_actual <- usher3_rho1(x_test, k1, b_siler, x0, x_cut, w = w0)
  
  expect_equal(rho1_actual, rho1_expected, tolerance = 1e-10)
})

test_that("usher3_integrand with x_cut = Inf is calculated correctly", {
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  y <- c(0, 5, 10, 20)
  
  # Calculate with x_cut = Inf
  integrand_inf <- usher3_integrand(y, k1, k2, b_siler, x_cut = Inf)
  
  # Calculate expected
  S12_expected <- exp(-k1 * y)
  S13_expected <- ssiler(y, b_siler)
  integrand_expected <- S12_expected * (S13_expected)^(1 - k2)
  
  expect_equal(integrand_inf, integrand_expected, tolerance = 1e-10)
})

test_that("usher3_integrand with x_cut behaves correctly", {
  k1 <- 0.02
  k2 <- 1.2
  b_siler <- th0[3:7]
  x_cut <- 6
  
  # Before cutoff: normal behavior
  y_before <- c(0, 3, 6)
  integrand_before <- usher3_integrand(y_before, k1, k2, b_siler, x_cut)
  
  # After cutoff: S12 should be constant
  y_after <- c(10, 20, 50)
  integrand_after <- usher3_integrand(y_after, k1, k2, b_siler, x_cut)
  
  # Verify all are positive
  expect_true(all(integrand_before > 0))
  expect_true(all(integrand_after > 0))
  
  # Manual check for point after cutoff
  y_test <- 20
  S12_test <- exp(-k1 * x_cut)  # Constant after cutoff
  S13_test <- ssiler(y_test, b_siler)
  integrand_expected <- S12_test * (S13_test)^(1 - k2)
  integrand_actual <- usher3_integrand(y_test, k1, k2, b_siler, x_cut)
  
  expect_equal(integrand_actual, integrand_expected, tolerance = 1e-10)
})

test_that("usher3_rho2 with x_cut = Inf is calculated correctly", {
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  x <- c(10, 25, 50)
  x0 <- 0
  w0 <- c(1, 0)  # Everyone well at x0
  
  # Calculate with x_cut = Inf
  rho2_inf <- usher3_rho2(x, k1, k2, b_siler, x0, x_cut = Inf, w = w0)
  
  # Calculate expected using manual integration
  # (This is a complex integral, so we just verify it's finite and positive)
  expect_true(all(is.finite(rho2_inf)))
  expect_true(all(rho2_inf > 0))
  
  # Also verify using direct calculation of the integral term
  f13_0_x <- dsiler(x, b_siler)
  S13_0_x <- ssiler(x, b_siler)
  S12_0_x0 <- exp(-k1 * x0)
  S13_0_x0 <- ssiler(x0, b_siler)
  
  integralTerm <- rep(NA, length(x))
  for (ii in 1:length(x)) {
    integralTerm[ii] <- integrate(
      usher3_integrand, x0, x[ii],
      k1 = k1, k2 = k2, b_siler = b_siler, x_cut = Inf
    )$value
  }
  
  rho2_expected <- k1 * k2 * f13_0_x * S13_0_x^(k2 - 1) * 
                   integralTerm / S12_0_x0 / S13_0_x0
  
  expect_equal(rho2_inf, rho2_expected, tolerance = 1e-8)
})

test_that("usher3_rho2 with x_cut produces valid densities", {
  k1 <- 0.02
  k2 <- 1.2
  b_siler <- th0[3:7]
  x0 <- 0
  x_cut <- 6
  w0 <- c(1, 0)  # Everyone well at x0
  
  # Test at various ages
  x <- c(10, 20, 30, 50)
  rho2 <- usher3_rho2(x, k1, k2, b_siler, x0, x_cut, w = w0)
  
  # Should all be finite and positive
  expect_true(all(is.finite(rho2)))
  expect_true(all(rho2 > 0))
})

# Tests for nll_usher3 with x_cut

test_that("nll_usher3 with x_cut = Inf is calculated correctly", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  
  # Calculate with x_cut = Inf
  nll_inf <- nll_usher3(th0, x, ill, x0, x_cut = Inf)
  
  # Should be finite and positive
  expect_true(is.finite(nll_inf))
  expect_true(nll_inf > 0)
  
  # Compare with manually calculated likelihood components
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]
  
  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0, x_cut = Inf, k2 = k2)
  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0, x_cut = Inf)
  
  ll_expected <- sum(log(rho1_wll)) + sum(log(rho2_ill))
  nll_expected <- -ll_expected
  
  expect_equal(nll_inf, nll_expected, tolerance = 1e-8)
})

test_that("nll_usher3 with x_cut produces valid likelihood", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cut <- 6
  
  # Calculate with finite x_cut
  nll_cutoff <- nll_usher3(th0, x, ill, x0, x_cut)
  
  # Should be finite and positive
  expect_true(is.finite(nll_cutoff))
  expect_true(nll_cutoff > 0)
  
  # Manual verification
  k1 <- th0[1]
  k2 <- th0[2]
  b_siler <- th0[3:7]
  
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]
  
  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0, x_cut, k2 = k2)
  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0, x_cut)
  
  ll_expected <- sum(log(rho1_wll)) + sum(log(rho2_ill))
  nll_expected <- -ll_expected
  
  expect_equal(nll_cutoff, nll_expected, tolerance = 1e-8)
})

test_that("nll_usher3 with x_cut handles NA illness indicators", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, NA, 0, 1, NA)
  x0 <- 0
  x_cut <- 6
  
  # Should handle NA values correctly
  nll_cutoff <- nll_usher3(th0, x, ill, x0, x_cut)
  
  expect_true(is.finite(nll_cutoff))
  expect_true(nll_cutoff > 0)
})

test_that("nll_usher3 with x_cut rejects negative parameters", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cut <- 6
  
  # Negative k1
  theta_neg_k1 <- th0
  theta_neg_k1[1] <- -0.01
  expect_equal(nll_usher3(theta_neg_k1, x, ill, x0, x_cut), Inf)
  
  # Negative k2
  theta_neg_k2 <- th0
  theta_neg_k2[2] <- -0.5
  expect_equal(nll_usher3(theta_neg_k2, x, ill, x0, x_cut), Inf)
  
  # Negative Siler parameter
  theta_neg_siler <- th0
  theta_neg_siler[3] <- -0.1
  expect_equal(nll_usher3(theta_neg_siler, x, ill, x0, x_cut), Inf)
})

test_that("nll_usher3 with x_cut supports Gompertz-Makeham", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cut <- 6
  
  # Use 5-parameter version (Gompertz-Makeham)
  th0_gompertz <- c(2e-2, 1.2, 0.368 * 0.01,
                    log(0.917 * 0.1 / (0.075 * 0.001)) / (0.917 * 0.1),
                    0.917 * 0.1)
  
  nll_gm <- nll_usher3(th0_gompertz, x, ill, x0, x_cut)
  
  expect_true(is.finite(nll_gm))
  expect_true(nll_gm > 0)
})

test_that("nll_usher3 likelihood changes with x_cut", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  
  # Calculate likelihoods with different cutoffs
  nll_inf <- nll_usher3(th0, x, ill, x0, x_cut = Inf)
  nll_6 <- nll_usher3(th0, x, ill, x0, x_cut = 6)
  nll_12 <- nll_usher3(th0, x, ill, x0, x_cut = 12)
  
  # All should be finite and positive
  expect_true(all(is.finite(c(nll_inf, nll_6, nll_12))))
  expect_true(all(c(nll_inf, nll_6, nll_12) > 0))
  
  # They should be different (different models)
  expect_false(isTRUE(all.equal(nll_inf, nll_6)))
  expect_false(isTRUE(all.equal(nll_inf, nll_12)))
  expect_false(isTRUE(all.equal(nll_6, nll_12)))
})

test_that("usher3_hessian raises error with finite x_cut", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cut <- 6
  
  # Should raise an error because finite x_cut is not yet supported
  expect_error(
    usher3_hessian(th0, x, ill, x0, x_cut),
    "Finite x_cut is not yet supported"
  )
})

test_that("usher3_hessian works with x_cut = Inf", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  
  # With x_cut = Inf (default), should work
  H <- usher3_hessian(th0, x, ill, x0, x_cut = Inf)
  
  # Should be a matrix
  expect_true(is.matrix(H))
  
  # Should have correct dimensions
  expect_equal(dim(H), c(7, 7))
  
  # Should not contain NA values
  expect_false(anyNA(H))
})

test_that("usher3_errors raises error with finite x_cut", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cut <- 6
  
  # Should raise an error because finite x_cut is not yet supported
  expect_error(
    usher3_errors(th0, x, ill, x0, x_cut),
    "Finite x_cut is not yet supported"
  )
})

test_that("usher3_errors works with x_cut = Inf", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  
  # With x_cut = Inf (default), should work
  # Note: th0 is not at a local minimum for this data, so the Hessian
  # may have negative diagonal values (expected behavior, not a bug)
  # We suppress warnings and only check structural properties
  suppressWarnings({
    errors <- usher3_errors(th0, x, ill, x0, x_cut = Inf)
  })
  
  expect_true(is.data.frame(errors))
  expect_equal(nrow(errors), 7)
  expect_true(all(c("Estimate", "StandErr", "z", "pval", "against", "sideAdj") %in% names(errors)))
})

# Tests for sample_usher3 with x_cut

test_that("sample_usher3 with x_cut = Inf works correctly", {
  set.seed(123)
  N <- 50
  dx <- 0.001  # Smaller dx for better integration
  xmax <- 120
  
  samples <- sample_usher3(N, th0, dx, xmax, x_cut = Inf)
  
  expect_equal(length(samples$x), N)
  expect_equal(length(samples$ill), N)
  expect_true(all(samples$x >= 0 & samples$x <= xmax))
  expect_true(all(samples$ill %in% c(TRUE, FALSE)))
})

test_that("sample_usher3 raises error with finite x_cut", {
  set.seed(456)
  N <- 50
  dx <- 0.001
  xmax <- 120
  x_cut <- 6
  
  expect_error(
    sample_usher3(N, th0, dx, xmax, x0 = 0, x_cut = x_cut),
    "Finite x_cut is not yet supported"
  )
  
  # Also test with age filtration
  x_mid <- 50
  infant_prop <- 0.1
  expect_error(
    sample_usher3(N, th0, dx, xmax, x_mid, infant_prop, x_cut = x_cut),
    "Finite x_cut is not yet supported"
  )
})

# Tests for temper_usher3 with x_cut

test_that("temper_usher3 with x_cut runs without error", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x0 <- 0
  x_cut <- 6
  
  # Test with very few cycles for speed
  results <- temper_usher3(th0, x = x, ill = ill, num_cyc = 2, 
                           x0 = x0, x_cut = x_cut)
  
  # Should return a list with expected components
  expect_true(is.list(results))
  expect_true(all(c("obj_fun", "th0", "optional_inputs", "temper",
                    "th_temper", "th") %in% names(results)))
})

test_that("temper_usher3 with x_cut = Inf works", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  
  # Test with very few cycles and x_cut = Inf (default behavior)
  results <- temper_usher3(th0, x = x, ill = ill, num_cyc = 2, 
                           x_cut = Inf)
  
  expect_true(is.list(results))
})

test_that("temper_usher3 with x_cut and Gompertz-Makeham works", {
  x <- c(10, 20, 30, 40, 50)
  ill <- c(0, 1, 0, 1, 0)
  x_cut <- 6
  
  th0_gompertz <- c(2e-2, 1.2, 0.368 * 0.01,
                    log(0.917 * 0.1 / (0.075 * 0.001)) / (0.917 * 0.1),
                    0.917 * 0.1)
  
  # Test with very few cycles
  results <- temper_usher3(th0_gompertz, x = x, ill = ill, num_cyc = 2,
                           use_gompertz = TRUE, x_cut = x_cut)
  
  expect_true(is.list(results))
  expect_equal(length(results$th), 5)
})
