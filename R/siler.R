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
#' \deqn{\lambda(x) = b_1 \exp(-b_2 \, x) + b_3 + b_5 \exp(b_5 \, (x - b_4))}
#'
#' where \eqn{b_i} (index notation) is the parameter vector. The cumulative
#' hazard is found by integrating the hazard from some initial point \eqn{x_0}
#' to \eqn{x},
#'
#' \deqn{\Lambda(x_0,x) = - \frac{b_1}{b_2} (e^{-b_2 x}
#'                        - e^{-b_2 x_0}) + b_3 (x-x_0)
#'                        + (e^{b_5 (x - b_4)} - e^{b_5 (x_0 - b_4)})}.
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
#' @param b The parameter vector
#'
#' @return The hazard evaluated at the locations in the input vector x
#'
#' @export
hsiler <- function(x, b) {
  if (any(b < 0)) {
    stop('No element of b can be negative')
  }
  return(b[1] * exp(-b[2] * x) + b[3] + b[5] * exp(b[5] * (x - b[4])))
}

#' @param qvect The vector of quantiles
#' @param N The number of observations
#'
#' @title
#' Calculate the Siler cumulative hazard
#'
#' @description
#' Calculate the Siler cumulative hazard at the locations in the input vector x
#' given the parameter vector b. Optionally, a conditional starting age can be
#' input. See details in the documentation for hsiler for a definition of the
#' Siler hazard and related quantities, including the cumulative hazard.
#'
#' @param x The vector of ages
#' @param b	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The cumulative hazard evaluated at the locations in the input vector
#'   x
#'
#' @export
chsiler <- function(x, b, x0 = 0) {
  check_siler_inputs(b, x0)

  if (x0 == 0) {
    return(-b[1] / b[2] * (exp(-b[2] * x) - 1) +
             b[3] * x + exp(b[5] * (x - b[4])) - exp(-b[5]*b[4]))
  } else {
    return(chsiler(x, b) - chsiler(x0, b))
  }
}

#' @title
#' Calculate the Siler survival
#'
#' @description
#' Calculate the Siler survival function at the locations in the input vector x
#' given the parameter vector b. Optionally, a conditional starting age can be
#' input. See details in the documentation for hsiler for a definition of the
#' Siler hazard and related quantities, including the survival function.
#
#'
#' @param x The vector of ages
#' @param b	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The survival evaluated at the locations in the input vector x
#'
#' @export
ssiler <- function(x, b, x0 = 0) {
  check_siler_inputs(b, x0)
  return(exp(-chsiler(x, b, x0)))
}

#' @title
#' Calculate the Siler probability density
#'
#' @description
#' Calculate the Siler probablity density at the locations in the input vector
#' x given the parameter vector b. Optionally, a conditional starting age can
#' be input. See details in the documentation for hsiler for a definition of
#' the Siler hazard and related quantities, including the probability density.
#'
#' @param x The vector of ages
#' @param b	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The probability density evaluated at the locations in the input
#'   vector x
#'
#' @export
dsiler <- function(x, b, x0 = 0) {
  check_siler_inputs(b, x0)
  return(hsiler(x, b) * ssiler(x, b, x0))
}

#' @title
#' Calculate the Siler cumulative distribution function (CDF) value
#'
#' @description
#' Calculate the Siler cumulative distribuition function (CDF) at the locations
#' in the input vector x given the parameter vector b. Optionally, a
#' conditional starting age can be input. See details in the documentation for
#' hsiler for a definition of the Siler hazard and related quantities,
#' including the CDF.
#'
#' @param x The vector of ages
#' @param b	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The cumulative distribution function evaluated at the locations in
#' the input vector x
#'
#' @export
psiler <- function(x, b, x0 = 0) {
  check_siler_inputs(b, x0)
  return(1 - ssiler(x, b, x0))
}

#' @title
#' Calculate the the inverse of the cumulative distribution function (CDF)
#'
#' @description
#' Calculate x given the input quantiles qvect. That is, solve the inverse
#' problem for the cumulative distribution. This is accomplished using the
#' following change of variables:
#'
#' xbar = b5*(x - x0) / (b5*x + 1)
#'
#' or, equivalently,
#'
#' x = ((xbar/b5) + x0)/(1-xbar).
#'
#' xbar lies between 0 to 1, which allows R's root-finder to be invoked on a
#' bounded interval. Normalizing by b5 regularizes the root finding for large
#' ages where the positive exponential term dominates the hazard.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param qvect The vector of quantiles
#' @param b	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The age, x, evaluated at the locations in the input vector qvect
#'
#' @export
qsiler <- function(qvect, b, x0 = 0) {
  # TODO: consider returning Inf if the quantile is 1.
  check_siler_inputs(b, x0)

  if (any(qvect < 0) | any(qvect > 1)) {
      stop("Quantiles should be between 0 and 1")
  }

  N <- length(qvect)
  xvect <- rep(NA, N)
  for (n in 1:N) {
    Sn <- 1 - qvect[n]
    froot <- function(xbar) {
      log(Sn) + chsiler(((xbar / b[5]) + x0) / (1 - xbar), b, x0)
    }
    xbar_n <- uniroot(froot, c(0, 1), tol = .Machine$double.eps / 2)$root
    xvect[n] <- ((xbar_n / b[5]) + x0) / (1 - xbar_n)
  }
  return(xvect)
}

#' @title
#' Make random samples from the Siler distribution
#'
#' @description
#' Make N random samples from the Siler distribution given the input parameter
#' vector b and a conditional starting age x0.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param N The number of samples to make
#' @param b	The parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return A vector of N random samples
#'
#' @export
rsiler <- function(N, b, x0 = 0) {
  check_siler_inputs(b, x0)
  cdf <- runif(N)
  return(qsiler(cdf, b, x0))
}

#' @title
#' Calculate the Hessian for the Siler distribution
#'
#' @description
#' Given the parameter vector b and the vector of ages x, calculate the Hessian
#' for a negative log-likelihood fit.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param b	The parameter vector
#' @param x The vector of ages
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The Hessian matrix (5 by 5)
#'
#' @export
hesiler <- function(b, x, x0 = 0) {
  # b is the first input as expected by most R optimization routines
  check_siler_inputs(b, x0)
  bStr <- c('b1','b2','b3','b4','b5')

  H <- matrix(NA,5,5) # Initialize the matrix

  # eta is the log-likelihood
  eta <- function(x,b1,b2,b3,b4,b5,x0) {
    -log(b1*exp(-b2*x)+b3+b5*exp(b5*(x-b4))) - b1*(exp(-b2*x)-exp(-b2*x0))/b2 + b3*(x-x0) + exp(b5*(x-b4))-exp(b5*(x0-b4))
  }

  # Analytically calculate Hessian diagonals
  for(n in 1:5) {
    H[n,n] <- sum(
      Deriv::Deriv(
        Deriv::Deriv(eta,bStr[n]),
        bStr[n])(x,b[1],b[2],b[3],b[4],b[5],x0)
    )
  }

  # Analytically calculate Hessian off-diagonals
  for(n1 in 1:4) {
    for(n2 in (n1+1):5) {
      H[n1,n2] <- sum(
        Deriv::Deriv(Deriv::Deriv(eta,bStr[n1]),bStr[n2])
         (x,b[1],b[2],b[3],b[4],b[5],x0)
      )
      H[n2,n1] <- H[n1,n2]
    }
  }
  return(H)
}

#' @title
#' Calculate the gradient for the Siler distribution
#'
#' @description
#' Given the parameter vector b and the vector of ages x, calculate the gradient
#' for a negative log-likelihood fit. Do this symbolically using the Deriv
#' package.
#'
#' See details in the documentation for hsiler for a definition of the Siler
#' hazard and related quantities.
#'
#' @param b The parameter vector
#' @param x The vector of ages
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The gradient vector (length 5)
#'
#' @export
gradnllsiler <- function(b, x, x0 = 0) {
  check_siler_inputs(b, x0)
  bStr <- c('b1','b2','b3','b4','b5')
  bStr <- c('b1', 'b2', 'b3', 'b4', 'b5')

  # eta is the log-likelihood
  eta <- function(x,b1,b2,b3,b4,b5,x0) {
    -log(b1*exp(-b2*x)+b3+b5*exp(b5*(x-b4))) - (b1/b2)*(exp(-b2*x)-exp(-b2*x0)) + b3*(x-x0) + exp(b5*(x-b4))-exp(b5*(x0-b4))
  }

  # Initialize the gradient vector
  grad <- numeric(length = 5)

  # Analytically calculate the gradient components
  for (n in 1:5) {
    grad[n] <- sum(Deriv::Deriv(eta, bStr[n])(x, b[1], b[2], b[3], b[4], b[5], x0))
  }

  return(grad)
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
#' @param b0 Initial value for optimization. Default from Gage and Dyke 1986,
#'   Table 2, Level 15. This is also used for the baseline age in the
#'   tranformed negative log-likelihood calculation.
#' @param calc_hessian Whether to calculate the Hessian (default FALSE)
#' @param lr The learning rate to use for fine tuning using gradient descent
#'   (default: 1e-5)
#' @param verbose Whether to print out optimization information
#'   (default: FALSE)
#'
#' @return A list consisting of the fit (on the transformed variable bbar) and
#'   maximum likelihood estimate of b. Optionally, the Hessian of the
#'   log-likelihood is returned, which allows estimation of the standard errors
#'   of the maximum likelihood estimate via the observed Fisher information
#'   matrix.
#'
#' @export
fit_siler <- function(x,
                      x0 = 0,
                      b0 = c(.175,
                             1.40,
                             .368 * .01,
                             log(.917 * .1/(.075 * .001))/(.917 * .1),
                             .917 * .1),
                      calc_hessian = FALSE,
                      lr=1e-5,
                      verbose=FALSE) {
  # The traditional parameterization of the Siler hazard is
  # a1 * exp(-a2*x) + a3 + a4*exp(-a5*x)
  # The demohaz parameterization is related to this one per
  # b1 = a1
  # b2 = a2
  # b3 = a3
  # b4 = log(a5/a4)/a5
  # b5 = a5
  # The Gage and Dyke vector in this traditional parameterization is
  # a0 = c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)
  # Hence, the default intial value for the demohaz parameterization is:
  # b0 = c(.175, 1.40, .368 * .01, log(.917 * .1/(.075 * .001))/(.917 * .1), .917 * .1)
  check_siler_inputs(b0, x0)
  bStr <- c('b1','b2','b3','b4','b5')

  xtable <- table(x)
  xvalues <- as.numeric(names(xtable))
  xcounts <- as.numeric(xtable)

  # fast_transformed_nllsiler gives the negative log-likelihood for the
  # parameterization exp(bbar) = 
  fit <- temper_and_tune(fast_nllsiler,
                         b0,
                         verbose=verbose,
                         xvalues = xvalues,
                         xcounts = xcounts,
                         lr=lr,
                         num_cyc=1000,
                         samps_per_cyc=20,
                         x0 = x0)
  
  b <- fit$th
  if (calc_hessian) {
    H <- hesiler(b,x,x0)
    return(list(b=b, fit=fit, hessian=H))
  } else {
    return(list(b=b, fit=fit))
  }
}

#' @title
#' Use counts to quickly calculate the Siler negative log-likelihood
#'
#' @description
#' Calculate the negative log-likelihood of the Siler probability density.
#' This can often be done much faster by tabulating repeated x-values If there
#' are no (or few) repeats, there is little speed loss from using tabulated
#' values. The inputs are a vector of x-values, xvalues, and corresponding
#' counts, xcounts, both of which are pre-computed. These two vectors can be
#' created from the vector of ages, x, using the following code:
#'
#' xtable <- table(x)
#' xvalues <- as.numeric(names(xtable))
#' xcounts <- as.numeric(xtable)
#'
#' The Siler hazard is:
#'
#' lambda(x) = b1 * exp(-b_2*x) + b_3 + b5 * exp(b5*(x-b4))
#'
#' Let b0 be a baseline parameter vector and define the transformed parameter
#' vector bbar to be (the notation is modified slightly from, e.g., b1 to b_1
#' to accommodate the initial parameter vector being b0)
#'
#' b_1 = b0_1 * exp(bbar1)
#' b_2 = b0_2 * exp(bbar2)
#' b_3 = b0_3 * exp(bbar3)
#' b_4 = b0_4 * exp(bbar4)
#' b_5 = b0_5 * exp(bbar5)
#'
#' This second transformation allows the parameter vector bbar to be
#' unconstrained (positive or negative) while the parameter vector b is
#' positive (assuming that the baseline vector beta is positive, which should
#' be ensured in its construction).
#'
#' @param bbar The transformed parameter vector, used to increase robustness of
#'   fitting
#' @param xvalues Locations at which to evaluate probability density function
#' @param xcounts The number of observations for each entry in xvalues
#' @param b0 The baseline parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return The negative log-likelihood
#'
#' @export
fast_nllsiler <- function(b,
                                      xvalues,
                                      xcounts,
                                      x0=0) {
  if (any(b < 0)) {
    return(Inf)
  }

  x <- xvalues
  haz <- b[1] * exp(-b[2]*x) + b[3] + b[5]*exp(b[5]*(x-b[4]))
  cumhaz <- (-b[1] / b[2]) * (exp(-b[2] * x) - exp(-b[2]*x0)) +
            b[3] * (x-x0) +
            exp(b[5] * (x - b[4])) - exp(b[5] * (x0 - b[4]))

  eta <- sum(xcounts*(-log(haz) + cumhaz))
  if (is.na(eta)) {
    eta <- Inf # Use infinity since this allows the hjk algorithm to handle it
  }
  return(eta)
}

#' @title
#' Vectorized calculation of the gradient of the Siler hazard, lambda
#'
#' @description
#' Calculate the gradient of the Siler hazard, lambda, for each age in the
#' input vector x. This is a vectorized calculation since no sum is done over
#' x. Hence, the return object is a matrix with dimensions [N x 5]
#' rather than a vector of length 5 (N = length(x)). This allows fast
#' calculation of the gradient when there are repeated entries in x, as in
#' gradsiler_fast.
#'
#' @param x Vector of ages
#' @param b Siler parameter vector
#'
#' @return A matrix of gradients with dimensions [length(x) x 5]
#'
#' @export
grad_hsiler_vect <- function(x, b) {
  if (any(b < 0)) {
    stop('No element of b can be negative')
  }
  # Initialize a matrix with dimensions [N x 5] to hold the gradients
  grad_matrix <- matrix(NA, length(x), 5)

  # Precalculate the exponential for the first and third Siler term
  exp1 <- exp(-b[2]*x)
  exp3 <- exp(-b[5]*(x-b[4]))

  # Calculate each gradient
  grad_matrix[,1] <- exp1
  grad_matrix[,2] <- -b[1]*x*exp1
  grad_matrix[,3] <- 1
  grad_matrix[,4] <- -b[5]^2 * exp3
  grad_matrix[,5] <- exp3 + b[5]*(x-b[4])*exp3
  return(grad_matrix)
}

#' @title
#' Vectorized calculation of the Siler hazard, lambda
#'
#' @description
#' Calculate the the Siler hazard, lambda, for each age in the input vector x.
#' This is a vectorized calculation that is used by, e.g., fast_gradnllsiler.
#' The return object is a matrix with dimensions [N x 5] rather than a vector
#' of length N (N = length(x)). The hazard is replicated 5 times for each
#' column to achieve these dimensions.
#'
#' @param x Vector of ages
#' @param b Siler parameter vector
#'
#' @return A matrix of hazards with dimensions [length(x) x 5]
#'
#' @export
hsiler_vect <- function(x, b) {
  if (any(b < 0)) {
    stop('No element of b can be negative')
  }
  hazard <- hsiler(x, b)
  return(replicate(5, hazard))
}

#' @title
#' Vectorized calculation of the gradient of the Siler cumulative hazard
#'
#' @description
#' Calculate the gradient of the Siler cumulative hazard for each age in the
#' input vector x. This is a vectorized calculation since no sum is done over
#' x. Hence, the return object is a matrix with dimensions [N x 5] rather than
#' a vector of length 5 (N = length(x)). This allows fast calculation of the
#' gradient when there are repeated entries in x, as in gradsiler_fast.
#'
#' @param x Vector of ages
#' @param b Siler parameter vector
#' @param x0 The conditional starting age [default: 0].
#'
#' @return A matrix of gradients with dimensions [length(x) x 5]
#'
#' @export
grad_chsiler_vect <- function(x, b, x0=0) {
  if (any(b < 0)) {
    stop('No element of b can be negative')
  }

  # Initialize a matrix with dimensions [N x 5] to hold the gradients
  grad_matrix <- matrix(NA, length(x), 5)

  # Precalculate the exponential for the first and third Siler term
  exp1 <- exp(-b[2]*x)
  exp3 <- exp(b[5]*(x-b[4]))

  # Calculate each gradient
  grad_matrix[,1] <- (-1/b[2]) * (exp1 - exp(-b[2]*x0))
  grad_matrix[,2] <- (b[1]/b[2]^2)*(exp1 - exp(-b[2]*x0)) + (b[1]/b[2])*(x*exp1 - x0*exp(-b[2]*x0))
  grad_matrix[,3] <- x-x0
  grad_matrix[,4] <- -b[5]*(exp3 - exp(b[5]*(x0-b[4])))
  grad_matrix[,5] <- (x-b[4])*exp3 - (x0-b[4])*exp(b[5]*(x0-b[4]))
  return(grad_matrix)
}

#' @export
fast_gradnllsiler <- function(b,
                              xvalues,
                              xcounts,
                              x0=0) {

  x <- xvalues
  grad_matrix <- matrix(NA, length(x), 5)

  # Precalculate the exponential for the first and third Siler term
  exp1 <- exp(-b[2]*x)
  exp3 <- exp(b[5]*(x-b[4]))

  # Directly calculate the hazard (it's faster do do it directly here since
  # exp1 and exp3 are alread calculated)
  lambda <- b[1]*exp1 + b[3] + b[5]*exp3

  grad_matrix[,1] <- -exp1 / lambda + (-1/b[2]) * (exp1 - exp(-b[2]*x0))
  grad_matrix[,2] <- b[1]*x*exp1/lambda + (b[1]/b[2]^2)*(exp1 - exp(-b[2]*x0)) + (b[1]/b[2])*(x*exp1 - x0*exp(-b[2]*x0))
  grad_matrix[,3] <- -1/lambda + (x-x0)
  grad_matrix[,4] <- b[5]^2*exp3/lambda -b[5]*(exp3 - exp(b[5]*(x0-b[4])))
  grad_matrix[,5] <- -exp3*(1+b[5]*(x-b[4]))/lambda + (x-b[4])*exp3 - (x0-b[4])*exp(b[5]*(x0-b[4]))

  # Create the counts matrix by replicated across columns five times, then
  # multiply by the counts
  count_matrix <- replicate(5, xcounts)

  return(colSums(grad_matrix * count_matrix))
}

#' @title
#' Convert from the traditional to demohaz Siler parameterization
#'
#' @description
#' For the input vector a in the traditional parameterization, convert to the
#' demohaz parameterization. The traditional parameterization of the hazard is
#' a[1]*exp(a[2]*x) + a[3] + a[4]*exp(a[5]*x) whereas the demohaz
#' parameterization is b[1]*exp(b[2]*x) + b[3] + b[5]*exp(b[5]*(x-b[4])). They
#' are related via
#'
#' b[1] = a[1]
#' b[2] = a[2]
#' b[3] = a[3]
#' b[4] = log(a[5]/a[4])/a[5]
#' b[5] = a[5]
#'
#' @param a The traditional pramaterization of the Siler hazard
#'
#' @return The demohaz pramaterization of the Siler hazard
#'
#' @export
trad_to_demohaz_siler_param <- function(a) {
  if (any(a < 0)) {
    stop('No element of a can be negative')
  }

  b <- a
  b[4] <- log(a[5]/a[4]) / a[5]
  return(b)
}

#' @title
#' Convert from the demohaz to traditional Siler parameterization
#'
#' @description
#' For the input vector a in the traditional parameterization, convert to the
#' demohaz parameterization. The traditional parameterization of the hazard is
#' a[1]*exp(a[2]*x) + a[3] + a[4]*exp(a[5]*x) whereas the demohaz
#' parameterization is b[1]*exp(b[2]*x) + b[3] + b[5]*exp(b[5]*(x-b[4])). They
#' are related via
#'
#' a[1] = b[1]
#' a[2] = b[2]
#' a[3] = b[3]
#' a[4] = b[5]*exp(-b[4]*b[5])
#' a[5] = b[5]
#'
#' @param b The demohaz pramaterization of the Siler hazard
#'
#' @return The traditional parameterization of the Siler hazard
#'
#' @export
demohaz_to_trad_siler_param <- function(b) {
  if (any(b < 0)) {
    stop('No element of b can be negative')
  }

  a <- b
  a[4] <- b[5]*exp(-b[4]*b[5])
  return(a)
}

#' @title
#' A utility function to check Siler inputs (b and x0)
#'
#' @description
#' Ensure that b is positive. Ensure that x0>=0.
#'
#' @param b	The parameter vector
#' @param x0 The conditional starting age
#'
#' @return None. An error is thrown if the inputs fail the checks.
#'
#' @export
check_siler_inputs<- function(b, x0) {
  # b cannot be negative (0 is okay, though possibly should not be allowed)
  if (any(b < 0)) {
    stop('No element of b can be negative')
  }

  # x0 cannot be negative
  if(x0 < 0) {
    stop("x0 cannot be negative")
  }
}