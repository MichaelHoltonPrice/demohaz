#' @title
#' Calculate the Siler hazard
#'
#' @description
#' Calculate the Siler hazard at the locations in the input vector x given
#' the parameter vector a. See details for a definition of the Siler hazard
#' and related quantities.
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
#'                        - e^{-a_2 x_0}) + a_3 (x-x_0)
#'                        + \frac{a_4}{a_5} (e^{a_5 x} - e^{a_5 x_0})}.
#'
#' For all functions, if \eqn{x_0} (optional) is not input it is set equal to
#' 0. The survival and cumulative density function are, respectively,
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
#' \code{dsiler} gives the probability density function, \code{psiler} the
#' cumulative distribution function, \code{qsiler} the quantile function,
#' \code{rsiler} random draws, \code{hsiler} the hazard, \code{chsiler} the
#' cumulative hazard, \code{ssiler} the survival, \code{nllsiler} the negative
#' log-likelihood, and \code{hesiler} the Hessian. fit_siler performs a robust,
#' maximum likelihood fit to yield the best-fit parameter vector a.
#'
#' @param x The vector of ages
#' @param a The parameter vector
#'
#' @return The hazard evaluated at the locations in the input vector x
#'
#' @export
hsiler <- function(x, a) {
  return(a[1] * exp(-a[2] * x) + a[3] + a[4] * exp(a[5] * x))
}

#' @param qvect The vector of quantiles
#' @param N The number of observations
#'
#' @title
#' Calculate the Siler cumulative hazard
#'
#' @description
#' Calculate the Siler cumulative hazard at the locations in the input vector x
#' given the parameter vector a. Optionally, a conditional starting age can be
#' input. See details in the documentation for hsiler for a definition of the
#' Siler hazard and related quantities, including the cumulative hazard.
#'
#' @param x The vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The cumulative hazard evaluated at the locations in the input vector
#'   x
#'
#' @export
chsiler <- function(x, a, x0 = 0) {
  if (x0 == 0) {
    return(-a[1] / a[2] * (exp(-a[2] * x) - 1) +
             a[3] * x + a[4] / a[5] * (exp(a[5] * x) - 1))
  } else {
    return(chsiler(x, a) - chsiler(x0, a))
  }
}

#' @title
#' Calculate the Siler survival
#'
#' @description
#' Calculate the Siler survival function at the locations in the input vector x
#' given the parameter vector a. Optionally, a conditional starting age can be
#' input. See details in the documentation for hsiler for a definition of the
#' Siler hazard and related quantities, including the survival function.
#
#'
#' @param x The vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The survival evaluated at the locations in the input vector x
#'
#' @export
ssiler <- function(x, a, x0 = 0) {
  return(exp(-chsiler(x, a, x0)))
}

#' @title
#' Calculate the Siler probability density
#'
#' @description
#' Calculate the Siler probablity density at the locations in the input vector
#' x given the parameter vector a. Optionally, a conditional starting age can
#' be input. See details in the documentation for hsiler for a definition of
#' the Siler hazard and related quantities, including the probability density.
#'
#' @param x The vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The probability density evaluated at the locations in the input
#'   vector x
#'
#' @export
dsiler <- function(x, a, x0 = 0) {
  return(hsiler(x, a) * ssiler(x, a, x0))
}

#' @title
#' Calculate the Siler cumulative distribution function (CDF) value
#'
#' @description
#' Calculate the Siler cumulative distribuition function (CDF) at the locations
#' in the input vector x given the parameter vector a. Optionally, a
#' conditional starting age can be input. See details in the documentation for
#' hsiler for a definition of the Siler hazard and related quantities,
#' including the CDF.
#'
#' @param x The vector of ages
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The cumulative distribution function evaluated at the locations in
#' the input vector x
#'
#' @export
psiler <- function(x, a, x0 = 0) {
  return(1 - ssiler(x, a, x0))
}

#' @title
#' Calculate the the inverse of the cumulative distribution function (CDF)
#'
#' @description
#' Calculate x given the input quantiles qvect. That is, solve the inverse
#' problem for the cumulative distribution. This is accomplished using the
#' following change of variables:
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
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param qvect The vector of quantiles
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The age, x, evaluated at the locations in the input vector qvect
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

#' @title
#' Make random samples from the Siler distribution
#'
#' @description
#' Make N random samples from the Siler distribution given the input parameter
#' vector a and a conditional starting age x0.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param N The number of samples to make
#' @param a	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return A vector of N random samples
#'
#' @export
rsiler <- function(N, a, x0 = 0) {
  cdf <- runif(N)
  return(qsiler(cdf, a, x0))
}

#' @title
#' Calculate the negative log-likelihood for the Siler distribution
#'
#' @description
#' Given the parameter vector a and the vector of ages x, calculate the Siler
#' negative log-likelihood.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities, including the negative log-likelihood.
#'
#' @param a	The parameter vector
#' @param x The vector of ages
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The negative log-likelihood (a scalar, not a vector)
#'
#' @export
nllsiler <- function(a, x, x0 = 0) {
  # a is the first input as expected by most R optimization routines
  return(-sum(log(dsiler(x, a, x0))))
}

#' @title
#' Calculate the Hessian for the Siler distribution
#'
#' @description
#' Given the parameter vector a and the vector of ages x, calculate the Hessian
#' for a negative log-likelihood fit.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param a	The parameter vector
#' @param x The vector of ages
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

#' @title
#' Maximum likelihood fit for the Siler distribution
#'
#' @description
#' Do a maximum likelihood Siler fit given the input vector of ages, x. To
#' allow unconstrained optimization to be used, optimization is done on the
#' transformed variable bbar (see the description for the function
#' fast_transformed_nllsiler).
#'
#' @param x Locations at which to evaluate probability density function
#' @param x0 The conditional starting age. The default is 0 (no offset).
#' @param alpha Initial value for optimization. Default from Gage and Dyke 1986,
#'   Table 2, Level 15. This is also used for the baseline age in the
#'   tranformed negative log-likelihood calculation.
#' @param calc_hessian Whether to calculate the Hessian (default FALSE)
#' @param temper Whether to use parallel tempering for the optimization (as
#'   opposed to base R's optim; default FALSE)
#'
#' @return A list consisting of the fit (on the transformed variable abar) and
#'   maximum likelihood estimate of a. Optionally, the Hessian of the
#'   log-likelihood is returned, which allows estimation of the standard errors
#'   of the maximum likelihood estimate via the observed Fisher information
#'   matrix.
#'
#' @export
fit_siler <- function(x,
                      x0 = 0,
                      alpha = c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1),
                      xmax = 120,
                      calc_hessian = FALSE,
                      temper=FALSE) {

  # Create function to calculate the negative log-likelihood for the barred
  # parameter vector
  # alpha is the baseline parameter vector in the baseline representation that
  # follows (a = alpha):
  # hazard = a1 * exp(-a2*x) + a3 + a4*exp(a5*x)
  #
  # beta is the baseline parameter vector in the following representation:
  # hazard = b1 * exp(-b2*x) + b3 + (b4 / exp(b5*xmax))*exp(b5*x)
  beta <- alpha
  beta[4] <- alpha[4] * exp(alpha[5] * xmax)

  xtable <- table(x)
  xvalues <- as.numeric(names(xtable))
  xcounts <- as.numeric(xtable)
  obj_fun <- function(bbar) {
    return(fast_transformed_nllsiler(bbar,xvalues,xcounts,beta,x0,xmax))
  }

  # Call optim to do the fit
  if(temper) {
    prop_scale_mat <- t(replicate(5,rev(seq(0.001,.1,len=21))))
    fit <- enneal::par_temper(
      rep(0,5),
      obj_fun,
      prop_scale=prop_scale_mat,
      num_cyc=1000,
      samps_per_cyc=20)
    n <- which.min(unlist(lapply(fit$chains,function(chain){chain$eta_best})))
    bfit <- beta * exp(fit$chains[[n]]$theta_best)
  } else {
    fit <- optim(rep(0, 5), obj_fun,control=list(maxit=10000))
    bfit <- beta * exp(fit$par)
  }
  afit <- bfit
  afit[4] <- bfit[4] / exp(bfit[5]*xmax)
#  return(list(fit=fit, afit=afit))

  if (calc_hessian) {
    H <- hesiler(afit,x,x0)
    return(list(fit = fit, a = afit, hessian = H))
  } else {
    return(list(fit = fit, a = afit))
  }
}

#' @title
#' Calculate the negative log-likelihood for the Siler distribution using counts
#' (hence fast) and assuming a transformed parameter vector.
#'
#' @description
#' The vector of ages may contain repeats. If there are many repeats, the
#' negative log-likelihood calculation will be slow if it is repeated for each
#' repeated value in x. To avoid this, the inputs to this function are a
#' vector of xvalues and corresponding counts. These two vectors can be created
#' from the vector of ages, x, using the following code:
#'
#' xtable <- table(x)
#' xvalues <- as.numeric(names(xtable))
#' xcounts <- as.numeric(xtable)
#'
#' To improve the robustness of the negative log-likelihood fit, a transformed
#' vector is used for the actual optimization and the fourth parameter is
#' adjusted so that it is measured with respect to a notional maximum age,
#' xmax. To make this precise, define the vector b to be:
#'
#' b1 = a1
#' b2 = a2
#' b3 = a3
#' b4 = a4 * exp(a5 * xmax)
#' b5 = a5
#'
#' The Siler hazard given these definitions is:
#'
#' lambda(x) = b1 * exp(-b_2*x) + b_3 + b4 * exp(b5*(x-xmax))
#'
#' In addition, let beta be a baseline parameter vector and define the
#' transformed parameter vector bbar to be
#'
#' b1 = beta1 * exp(bbar1)
#' b2 = beta2 * exp(bbar2)
#' b3 = beta3 * exp(bbar3)
#' b4 = beta4 * exp(bbar4)
#' b5 = beta5 * exp(bbar5)
#'
#' This second transformation allows the parameter vector bbar to be
#' unconstrained (positive or negative) while the parameter vector b is
#' positive (assuming that the baseline vector beta is positive, which should
#' be ensured in its construction). In principle, the choices of beta and xmax
#' are arbitrary as far as the maximum likelihood optimization is concerned.
#' However, appropriate choices can make the fitting far more robust. If
#' fitting human demographic data, we suggest that xmax = 120 (years) and
#' derive beta from beta from Gage and Dyke 1986, Table 2, Level 15.
#'
#' @param bbar The transformed parameter vector, used to increase robustness of
#'   fitting
#' @param xvalues Locations at which to evaluate probability density function
#' @param xcounts The number of observations for each entry in xvalues
#' @param beta The baseline parameter vector
#' @param x0 The conditional starting age [default: 0].
#' @param xmax The baseline maximum age (default: 120)
#'
#' @return The negative log-likelihood
#'
#' @export
fast_transformed_nllsiler <- function(bbar,
                                      xvalues,
                                      xcounts,
                                      beta,
                                      x0=0,
                                      xmax=120) {
  b <- beta * exp(bbar)
  x <- xvalues
  haz <- b[1] * exp(-b[2]*x) + b[3] + (b[4]/exp(b[5]*xmax))*exp(b[5]*x)
  cumhaz <- (b[1]/b[2])*(-exp(-b[2]*x)+exp(-b[2]*x0)) + b[3]*(x-x0) +
     (b[4]/exp(b[5]*xmax)/b[5])*(exp(b[5]*x)-exp(b[5]*x0))
  return(sum(xcounts*(-log(haz) + cumhaz)))
}
