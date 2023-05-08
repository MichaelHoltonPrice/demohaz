test_that("Doing Siler distribution unit tests", {
  # a0 is the baseline parameter vector used to test Siler functions. It is from
  # Gage and D 1986, Table 2, Level 15.
  a0 <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)
  
  xtest <- seq(0,120,by=.1)
  
  x0 <- 10
  
  # Unit test hsiler, the hazard
  expect_equal(
    hsiler(xtest,a0),
    a0[1] * exp(-a0[2] * xtest) + a0[3] + a0[4] * exp(a0[5] * xtest)
  )
  
  # Unit test chsiler, the cumulative hazard
  expect_equal(
    chsiler(xtest,a0),
    -a0[1] / a0[2] * (exp(-a0[2] * xtest) - 1) +
               a0[3] * xtest + a0[4] / a0[5] * (exp(a0[5] * xtest) - 1)
  )
  
  expect_equal(
    chsiler(xtest,a0,x0),
    chsiler(xtest,a0) - chsiler(x0,a0)
  )
  
  # Unit test ssiler, the survival
  expect_equal(
    ssiler(xtest,a0),
    exp(-chsiler(xtest,a0))
  )
  
  expect_equal(
    ssiler(xtest,a0,x0),
    exp(-chsiler(xtest,a0,x0))
  )
  
  # Unit test dsiler, the density
  expect_equal(
    dsiler(xtest,a0),
    hsiler(xtest,a0) * ssiler(xtest,a0)
  )
  
  expect_equal(
    dsiler(xtest,a0,x0),
    hsiler(xtest,a0) * ssiler(xtest,a0,x0)
  )
  
  # Unit test psiler, the cumulative distribution function (CDF)
  expect_equal(
    psiler(xtest,a0),
    1 - ssiler(xtest,a0)
  )
  
  expect_equal(
    psiler(xtest,a0,x0),
    1 - ssiler(xtest,a0,x0)
  )
  
  # Unit test qsiler, the quantile.
  qtest <- seq(0,1-.01,by=.01)
  expect_error(
    x_of_q <- qsiler(qtest,a0),
    NA
  )
  
  expect_equal(
    length(x_of_q),
    length(qtest)
  )
  
  expect_equal(
    any(is.na(x_of_q)),
    FALSE
  )
  
  expect_error(
    qsiler(1,a0)
  )
  
  expect_error(
    x_of_q <- qsiler(qtest,a0,x0),
    NA
  )
  
  expect_equal(
    length(x_of_q),
    length(qtest)
  )
  
  expect_equal(
    any(is.na(x_of_q)),
    FALSE
  )
  
  expect_error(
    qsiler(1,a0,x0)
  )
  
  # Unit test rsiler, the random sample generator
  set.seed(147144) # from random.org
  N <- 200
  
  expect_error(
    xsamp <- rsiler(N,a0),
    NA
  )
  
  expect_equal(
    length(xsamp),
    N
  )
  
  expect_equal(
    any(is.na(xsamp)),
    FALSE
  )
  
  expect_error(
    xsamp <- rsiler(N,a0,x0),
    NA
  )
  
  expect_equal(
    length(xsamp),
    N
  )
  
  expect_equal(
    any(is.na(xsamp)),
    FALSE
  )
  
  
  # Unit test nllsiler, the negative log-likelihood
  expect_equal(
    nllsiler(a0,xtest),
    -sum(log(dsiler(xtest, a0)))
  )
  
  expect_equal(
    nllsiler(a0,xtest,x0),
    -sum(log(dsiler(xtest, a0,x0)))
  )
  
  # Unit test hesiler, the Hessian
  expect_error(
    hess <- hesiler(a0,xtest),
    NA
  )
  
  expect_equal(
    dim(hess),
    c(5,5)
  )
  
  expect_equal(
    any(is.na(hess)),
    FALSE
  )
  
  expect_error(
    hess <- hesiler(a0,xtest,x0),
    NA
  )
  
  expect_equal(
    dim(hess),
    c(5,5)
  )
  
  expect_equal(
    any(is.na(hess)),
    FALSE
  )
  
  # Unit test fit_siler, the maximum likelihood fitting
  expect_error(
    fit_obj <- fit_siler(xsamp),
    NA
  )
  
  expect_equal(
    names(fit_obj),
    c("a", "fit")
  )
  
  expect_equal(
    length(fit_obj$a),
    5
  )
  
  expect_equal(
    any(is.na(fit_obj$a)),
    FALSE
  )
  
  expect_equal(
    names(fit_obj$fit),
    c("obj_fun","th0","temper","th_temper","hjk_output","th")
  )
  
  expect_equal(
    length(fit_obj$fit$th),
    5
  )
  
  expect_equal(
    any(is.na(fit_obj$fit$th)),
    FALSE
  )
  
  # 0 for convergence means successful completion
  expect_equal(
    fit_obj$fit$hjk_output$convergence,
    0
  )
  
  expect_error(
    fit_obj <- fit_siler(xsamp,x0=x0,calc_hessian=TRUE),
    NA
  )
  
  expect_equal(
    names(fit_obj),
    c("a","fit","hessian")
  )
  
  expect_equal(
    length(fit_obj$a),
    5
  )
  
  expect_equal(
    any(is.na(fit_obj$a)),
    FALSE
  )
  
  expect_equal(
    names(fit_obj$fit),
    c("obj_fun","th0","temper","th_temper","hjk_output","th")
  )
  
  expect_equal(
    length(fit_obj$fit$th),
    5
  )
  
  expect_equal(
    any(is.na(fit_obj$fit$th)),
    FALSE
  )
  
  # 0 for convergence means successful completion
  expect_equal(
    fit_obj$fit$hjk_output$convergence,
    0
  )
  
  expect_equal(
    dim(fit_obj$hessian),
    c(5,5)
  )
  
  expect_equal(
    any(is.na(fit_obj$hessian)),
    FALSE
  )
  
  # Unit test fast_transformed_nllsiler, the function for quickly and robustly
  # calculating the negative log-likelihood
  xmax <- 120
  alpha <- a0
  beta <- alpha
  beta[4]
  beta[4] <- alpha[4] * exp(alpha[5] * xmax)
  x <- c(2,2,20,30,110)
  xtable <- table(x)
  xvalues <- as.numeric(names(xtable))
  xcounts <- as.numeric(xtable)
  
  bbar <- rep(0.05,5)
  b <- beta * exp(bbar)
  a <- b
  a[4] <- b[4] / exp(b[5] * xmax)
  
  expect_equal(
    fast_transformed_nllsiler(bbar,xvalues,xcounts,beta,xmax=xmax),
    nllsiler(a, x)
  )
  
  expect_equal(
    fast_transformed_nllsiler(bbar,xvalues,xcounts,beta,xmax=xmax,x0=1),
    nllsiler(a, x, x0=1)
  )
})