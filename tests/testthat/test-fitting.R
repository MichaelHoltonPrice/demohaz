# Since fitting only has one function, do all tests in a single script (rather
# than dividing unit and functional tests).
test_that("Doing fitting unit and functional tests", {
  # Test temper_and_tune on a univariate Gaussian fit
  mu0  <- -0.5
  sig0 <- 1.5
  
  th0 <- c(mu0, sig0)
  th_bar0 <- c(mu0, log(sig0))
  
  uni_gauss_obj_fun <- function(th_bar, xsamp) {
    mu  <- th_bar[1]
    sig <- exp(th_bar[2]) # Force sig to be positive
    N <- length(xsamp)
    eta <- N*log(sig) + 0.5 * sum((xsamp - mu)^2)/sig^2
    return(eta)
  }
  
  # Do unit tests for a variety of input patterns to temper_and_tune
  N <- 100
  set.seed(1000)
  xsamp <- mu0 + rnorm(N) * sig0
  
  # Case 1
  expect_error(
    fit_obj <- temper_and_tune(uni_gauss_obj_fun, c(0,0), xsamp=xsamp),
    NA
  )
  
  expect_equal(
    names(fit_obj),
    c("obj_fun", "th0", "temper", "th_temper", "hjk_output", "th")
  )
  
  expect_equal(
    length(fit_obj$th0),
    2
  )
  
  expect_false(
    any(is.na(fit_obj$th))
  )
  
  expect_equal(
    length(fit_obj$temper$chains),
    9
  )
  
  expect_equal(
    length(fit_obj$temper$chains[[1]]$eta_vect),
    100*20
  )
  
  expect_equal(
    fit_obj$temper$inputs$prop_scale,
    t(replicate(2, rev(seq(0.001,.1,len=9))))
  )
  
  # Case 2
  prop_scale_mat <- t(replicate(2, rev(seq(0.01,.1,len=5))))
  expect_error(
    fit_obj <- temper_and_tune(uni_gauss_obj_fun, c(0,0), xsamp=xsamp,
                               num_cyc=40,
                               samps_per_cyc=15,
                               temp_vect = 10^(rev(seq(-1,1,by=.5))),
                               prop_scale_mat=prop_scale_mat),
    NA
  )
  
  expect_equal(
    names(fit_obj),
    c("obj_fun", "th0", "temper", "th_temper", "hjk_output", "th")
  )
  
  expect_equal(
    length(fit_obj$th0),
    2
  )
  
  expect_false(
    any(is.na(fit_obj$th))
  )
  
  expect_equal(
    length(fit_obj$temper$chains),
    5
  )
  
  expect_equal(
    length(fit_obj$temper$chains[[1]]$eta_vect),
    40*15
  )
  
  expect_equal(
    fit_obj$temper$inputs$prop_scale,
    prop_scale_mat
  )
  
  # Check the functioning of temper_and_tune
  N <- 1000000
  set.seed(1001)
  xsamp <- mu0 + rnorm(N) * sig0
  
  expect_error(
    fit_obj <- temper_and_tune(uni_gauss_obj_fun, c(0,0), xsamp=xsamp),
    NA
  )
  
  expect_equal(
    fit_obj$th,
    th_bar0,
    tol=1e-3
  )
})