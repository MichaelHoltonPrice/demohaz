#' @title Siler hazard model
#'
#' @description
#' \code{dsiler} gives the probability density function, \code{psiler} the
#' cumulative density function, \code{qsiler} the quantile function,
#' \code{rsiler} random draws, \code{hsiler} the hazard, \code{chsiler} the
#' cumulative hazard, \code{ssiler} the survival, \code{nllsiler} the negative
#' log-likelihood, and \code{hesiler} the Hessian.
#'
#' @details The Siler hazard is
#'
#' \deqn{\lambda(x) = a_1 \exp(-a_2 \, x) + a_3 + a_4 \exp(a_5 \, x)}
#'
#' where \eqn{a_i} (index notation) is the parameter vector. The cumulative
#' hazard is found by integrating the hazard from some initial point \eqn{x_0}
#' to \eqn{x},
#'
#' \deqn{\Lambda(x_0,x) = - \frac{a_1}{a_2} (e^{-a_2 x}
#'                        + e^{-a_2 x_0}) + a_3 (x-x_0)
#'                        + \frac{a_4}{a_5} (e^{a_5 x} - e^{a_5 x_0})}.
#'
#' If \eqn{x_0} (optional) is not input it is set equal to 0. The survival and
#' cumulative density function are, respectively,
#'
#' \deqn{S(x_0,x) = \exp(-\Lambda(x_0,x))
#'                = \frac{S(0,x)}{S(0,x_0)} = \frac{S(x)}{S(x_0)}}
#'
#' and
#'
#' \deqn{F(x_0,x) = 1 - S(x_0,x)}.
#'
#' The probability density function is
#'
#' \deqn{f(x_0,x) = \frac{d}{dx} F(x_0,x) = \lambda(x) S(x_0,x)}.
#'
#' @param x the ending point
#' @param x0 the initial point
#' @param a parameterization vector
#' @param qvect vector of quantiles
#' @param N number of observations
#'
#' @keywords Siler, hazard, survival, demography
#'
#' @return A numeric vector of length equal to the input
#
#' @author Michael Holton Price <MichaelHoltonPrice@gmail.com>

#' Calculate the Siler hazard
#'
#' @param x The parameter vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of hazards (same length as x)
#' @export
hsiler <- function(x, a) {
  return(a[1] * exp(-a[2] * x) + a[3] + a[4] * exp(a[5] * x))
}

#' Calculate the Siler cumulative hazard
#'
#' @param x The parameter vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of cumulative hazards (same length as x)
#' @export
chsiler <- function(x, a, x0 = 0) {
  if (x0 == 0) {
    return(-a[1] / a[2] * (exp(-a[2] * x) - 1) +
             a[3] * x + a[4] / a[5] * (exp(a[5] * x) - 1))
  } else {
    return(chsiler(x, a) - chsiler(x0, a))
  }
}

#' Calculate the Siler survival
#'
#' @param x The parameter vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of cumulative hazards (same length as x)
#' @export
ssiler <- function(x, a, x0 = 0) {
  return(exp(-chsiler(x, a, x0)))
}

#' Calculate the Siler probability density
#'
#' @param x The parameter vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of probability densities (same length as x)
#' @export
dsiler <- function(x, a, x0 = 0) {
  return(hsiler(x, a) * ssiler(x, a, x0))
}

#' Calculate the Siler cumulative distribution function (CDF) value
#'
#' @param x The parameter vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of cumulative distribution function values (same length as
#'   x)
#' @export
psiler <- function(x, a, x0 = 0) {
  return(1 - ssiler(x, a, x0))
}

#' Calculate x given input quantile(s)
#'
#' Solve the inverse problem for the cumulative density with a change of
#' variable
#'
#' xbar = a5*(x - x0) / (a5*x + 1)
#'
#' or, equivalently,
#'
#' x = ((xbar/a5) + x0)/(1-xbar).
#'
#' xbar lies between 0 to 1, which allows R's root-finder to be invoked on a
#' bounded interval. Normalizing by a5 regularizes the root finding for large
#' ages where the positive exponential term dominates the hazard.
#'
#' @param x The parameter vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of quantiles (same length as  x)
#'
#' @export
qsiler <- function(qvect, a, x0 = 0) {
  # TODO: consider returning Inf if the quantile is 1.
  N <- length(qvect)
  xvect <- rep(NA, N)
  for (n in 1:N) {
    Sn <- 1 - qvect[n]
    froot <- function(xbar) {
      log(Sn) + chsiler(((xbar / a[5]) + x0) / (1 - xbar), a, x0)
    }
    xbar_n <- uniroot(froot, c(0, 1), tol = .Machine$double.eps / 2)$root
    xvect[n] <- ((xbar_n / a[5]) + x0) / (1 - xbar_n)
  }
  return(xvect)
}

#' Make random samples from the Siler distribution
#'
#' @param N The number of samples to make
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @return The vector of samples (length N)
#'
#' @export
rsiler <- function(N, a, x0 = 0) {
  cdf <- runif(N)
  return(qsiler(cdf, a, x0))
}

#' Calculate the negative log-likelihood for the Siler distribution
#'
#' @param a	The parameter vector
#' @param x The parameter vector of ages
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The vector of negative log-likelihoods (same length as x)
#'
#' @export
nllsiler <- function(a, x, x0 = 0) {
  # a is the first input as expected by most R optimization routines
  return(-sum(log(dsiler(x, a, x0))))
}

#' Calculate the Hessian for the Siler distribution
#'
#' @param a	The parameter vector
#' @param x The parameter vector of ages
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The Hessian matrix (5 by 5)
#'
#' @export
hesiler <- function(a, x, x0 = 0) {
  # a is the first input as expected by most R optimization routines
  aStr <- c('a1','a2','a3','a4','a5')

  H <- matrix(NA,5,5) # Initialize the matrix

  # eta is the log-likelihood
  eta <- function(x,a1,a2,a3,a4,a5,x0) {
    log(a1*exp(-a2*x)+a3+a4*exp(a5*x))
    + a1*(exp(-a2*x)-exp(-a2*x0))/a2
    - a3*(x-x0)
    - a4*(exp(a5*x)-exp(a5*x0))/a5
  }

  # Analytically calculate Hessian diagonals
  for(n in 1:5) {
    H[n,n] <- sum(
      Deriv::Deriv(
        Deriv::Deriv(eta,aStr[n]),
        aStr[n])(x,a[1],a[2],a[3],a[4],a[5],x0)
    )
  }

  # Analytically calculate Hessian off-diagonals
  for(n1 in 1:4) {
    for(n2 in (n1+1):5) {
      H[n1,n2] <- sum(
        Deriv::Deriv(Deriv::Deriv(eta,aStr[n1]),aStr[n2])
         (x,a[1],a[2],a[3],a[4],a[5],x0)
      )
      H[n2,n1] <- H[n1,n2]
    }
  }
  return(H)
}

#' Maximum likelihood fit for Siler distribution
#'
#' Do a maximum likelihood Siler fit for the input data vector x. To allow
#' unconstrained optimization to be used, optimization is done on the
#' transformed variable abar = log(a/a0).
#'
#' @param x Locations at which to evaluate probability density function
#' @param a0 Initial value for optimization. Default from Gage and Dyke 1986,
#'   Table 2, Level 15
#' @param x0 x-offset. The default is 0 (no offset).
#' @param calc_hessian Whether to calculate the Hessian (default FALSE)
#'
#' @return A list consisting of the fit (on the transformed variable abar) and
#'   maximum likelihood estimate of a. Optionally, the Hessian of the
#'   log-likelihood is returned, which allows estimation of the standard errors
#'   of the maximum likelihood estimate via the observed Fisher information
#'   matrix.
#'
#' @export
fit_siler <- function(x,
                      a0 = c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1),
                      x0 = 0,
                      calc_hessian = FALSE) {

  # Create function to calculate the negative log-likelihood for the barred
  # parameter vector
  nllbar <- function(abar, x, a0) {
    a <- a0 * exp(abar)
    -sum(log(dsiler(x, a, x0)))
  }

  # Call optim to do the fit
  # TODO: consider adding support for parallel tempering
  fit <- optim(rep(0, 5), nllbar, x = x, a0 = a0)

  afit <- a0 * exp(fit$par)
  if (calc_hessian) {
    H <- hesiler(afit,x,x0)
    return(list(fit = fit, a = afit, hessian = H))
  } else {
    return(list(fit = fit, a = afit))
  }
}