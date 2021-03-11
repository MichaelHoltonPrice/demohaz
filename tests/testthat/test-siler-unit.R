# a0 is the baseline parameter vector used to test Siler functions. It is from
# Gage and Dyke 1986, Table 2, Level 15.
a0 <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)

xtest <- seq(0,120,by=.1)

x0 <- 10

# TODO: add unit tests for chsiler
# TODO: add unit tests for siler
# TODO: add unit tests for dsiler

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
  c("fit","a")
)

expect_equal(
  length(fit_obj$a),
  5
)

expect_equal(
  any(is.na(fit_obj$a)),
  FALSE
)

# 0 for convergence means successful completion
expect_equal(
  fit_obj$fit$convergence,
  0
)

expect_error(
  fit_obj <- fit_siler(xsamp,x0=x0,calc_hessian=TRUE),
  NA
)

expect_equal(
  names(fit_obj),
  c("fit","a","hessian")
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
  dim(fit_obj$hessian),
  c(5,5)
)

expect_equal(
  any(is.na(fit_obj$hessian)),
  FALSE
)
# 0 for convergence means successful completion
expect_equal(
  fit_obj$fit$convergence,
  0
)
