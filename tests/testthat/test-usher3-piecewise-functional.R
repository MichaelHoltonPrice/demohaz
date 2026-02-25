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

test_that("sample then fit recovers parameters with x_cut = Inf", {
  set.seed(314159)
  N <- 5000
  dx <- 0.001
  xmax <- 120

  samples <- sample_usher3(N, th0, dx, xmax, x_cut = Inf)

  # Jitter starting vector so we don't start at the optimum
  th_start <- th0 * runif(7, min = 0.8, max = 1.2)

  # Light tempering since we start close to the solution; Nelder-Mead does the rest
  result <- temper_and_tune_usher3(th_start, x = samples$x, ill = samples$ill,
                                    x0 = 0, x_cut = Inf,
                                    num_cyc = 5, samps_per_cyc = 5, tune = TRUE,
                                    control = list(maxit = 1000, reltol = 1e-12, abstol = 1e-12))

  expect_equal(th0, result$th, tolerance = 5e-3)
})

test_that("sample then fit recovers parameters with finite x_cut", {
  set.seed(271828)
  N <- 5000
  dx <- 0.001
  xmax <- 120
  x_cut <- 6

  samples <- sample_usher3(N, th0, dx, xmax, x0 = 0, x_cut = x_cut)

  # Jitter starting vector so we don't start at the optimum
  th_start <- th0 * runif(7, min = 0.8, max = 1.2)

  # Light tempering since we start close to the solution; Nelder-Mead does the rest
  result <- temper_and_tune_usher3(th_start, x = samples$x, ill = samples$ill,
                                    x0 = 0, x_cut = x_cut,
                                    num_cyc = 5, samps_per_cyc = 5, tune = TRUE,
                                    control = list(maxit = 1000, reltol = 1e-12, abstol = 1e-12))

  expect_equal(th0, result$th, tolerance = 5e-3)
})
