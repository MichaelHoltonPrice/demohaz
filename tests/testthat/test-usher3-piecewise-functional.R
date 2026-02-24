library(testthat)
library(demohaz)

# Baseline parameter vector (LEH example)
b0 <- c(.175, 1.40, .368 * .01,
        log(.917 * .1/(.075 * .001))/(.917 * .1),
        .917 * .1)
th0 <- c(2e-2, 1.2, b0)

test_that("sample_usher3 with finite x_cut passes round-trip NLL validation", {
  set.seed(42)
  N <- 2000
  dx <- 0.001
  xmax <- 120
  x_cut <- 6

  samples <- sample_usher3(N, th0, dx, xmax, x0 = 0, x_cut = x_cut)

  # NLL at the true parameters
  nll_true <- nll_usher3(th0, samples$x, samples$ill, x0 = 0, x_cut = x_cut)

  # NLL at perturbed parameters (multiply each by 1.5)
  th_perturbed <- th0 * 1.5
  nll_perturbed <- nll_usher3(th_perturbed, samples$x, samples$ill,
                               x0 = 0, x_cut = x_cut)

  # True parameters should yield a better (lower) NLL than wrong parameters
  expect_true(nll_true < nll_perturbed)
})
