#' @title Usher 3 Illness-Death Model Functions
#'
#' @description This file contains functions related to the Usher 3
#' illness-death model,  including the calculation of the pseudo-densities rho1
#' and rho2, and the integrand used in the calculation of rho2.
#'
#' @details The Usher 3 illness-death model is a multi-state model that
#' describes the  transitions of individuals among three states: healthy, ill,
#' and dead. The model is parameterized by transition rates between the states,
#' and the Siler hazard model is used to describe the baseline mortality.
#'
#'   The density rho1 represents the probability density of observing an
#'   individual who died in the healthy state at a given age, while rho2
#'   represents the probability density of observing an individual who died in
#'   the ill state at a given age. They are not proper densities since it is
#'   the some of the two, rho1 + rho2, that integrates to 1; hence, we
#'   sometimes call them pseudo-densities.
#'
#'   The calculation of rho2 involves an integral term, which is computed using
#'   the \code{integrate} function in R. The integrand for this term is given
#'   by the \code{usher3_integrand} function.
#'
#' @name usher3
NULL

#' @rdname usher3
#' @param x The vector of ages
#' @param k1 The transition rate from the healthy state to the ill state
#' @param b_siler The parameter vector for the Siler hazard model
#' @param x0 The conditional starting age [default: 0]
#'
#' @return The density rho1 evaluated at the locations in the input vector x
#'
#' @export
usher3_rho1 <- function(x, k1, b_siler, x0 = 0) {
  f13 <- dsiler(x, b_siler, x0)
  S12 <- exp(-k1 * (x - x0))
  rho1 <- f13 * S12
  return(rho1)
}

#' @rdname usher3
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#'
#' @return The density rho2 evaluated at the locations in the input vector x
#'
#' @export
usher3_rho2 <- function(x, k1, k2, b_siler, x0 = 0) {
  f13_0_x <- dsiler(x, b_siler)
  S13_0_x <- ssiler(x, b_siler)
  S12_0_x0 <- exp(-k1 * x0)
  S13_0_x0 <- ssiler(x0, b_siler)

  integralTerm <- rep(NA, length(x))
  for (ii in 1:length(x)) {
    integralTerm[ii] <- tryCatch(
      integrate(
        usher3_integrand, x0, x[ii],
        k1 = k1, k2 = k2, b_siler = b_siler
      )$value,
      error = function(e) { NA }
    )
  }

  rho2 <- k1 * k2 * f13_0_x * S13_0_x^(k2 - 1) * integralTerm / S12_0_x0 / S13_0_x0
  return(rho2)
}

#' @rdname usher3
#' @param y The variable of integration
#'
#' @return The integrand evaluated at y
#'
#' @export
usher3_integrand <- function(y, k1, k2, b_siler) {
  # _0_y indicates that the relevant quantity is relative to the interval 0 to y
  S12_0_y <- exp(-k1 * y)
  S13_0_y <- ssiler(y, b_siler)

  return(S12_0_y * (S13_0_y)^(1 - k2))
}

#' Calculate the negative log-likelihood for the Usher 3 illness-death model
#'
#' @param theta The parameter vector for the Usher illness-death model with
#'   ordering [k1, k2, mort]. The mortality parameter vector, mort, is
#'   [b1, b2, b3, b4, b5].
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators, which can have NA entries
#' @param x0 Conditional starting age [default: 0]
#'
#' @return The negative log-likelihood value
#'
#' @export
nll_usher3 <- function(theta, x, ill, x0 = 0) {
  k1 <- theta[1]
  if (k1 < 0) {
    return(Inf)
  }

  k2 <- theta[2]
  if (k2 < 0) {
    return(Inf)
  }

  if (any(theta[3:7] < 0)) {
    return(Inf)
  }
  b_siler <- theta[3:7]

  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]

  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0)
  if (any(is.na(rho2_ill))) {
    return(Inf)
  }

  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0)

  ll <- sum(log(rho1_wll)) + sum(log(rho2_ill))
  return(-ll)
}

#' Calculate the Hessian matrix for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#'
#' @return The Hessian matrix
#'
#' @export
usher3_hessian <- function(theta, x, ill, x0 = 0) {
  H <- numDeriv::hessian(usher3_nll, theta, method.args = list(eps = 1e-12),
                         ageVect = x, illVect = ill, x0 = x0)
  return(H)
}

#' Wrapper function for the negative log-likelihood of the Usher 3 model
#'
#' @param paramVect The parameter vector for the Usher illness-death model
#' @param ageVect The vector of ages-at-death
#' @param illVect The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#'
#' @return The negative log-likelihood value
#'
#' @export
usher3_nll <- function(paramVect, ageVect, illVect, x0 = 0) {
  return(nll_usher3(paramVect, ageVect, illVect, x0))
}

#' Calculate standard errors, z-scores, and p-values for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#'
#' @return A data frame with standard errors, z-scores, and p-values
#'
#' @export
usher3_errors <- function(theta, x, ill, x0 = 0) {
  H <- usher3_hessian(theta, x, ill, x0)
  against <- c(0, 1, 0, 0, 0, 0, 0)
  sideAdjustment <- c(1, 2, 1, 1, 1, 1, 1)
  varName <- c('k1', 'k2', 'a1', 'b1', 'a2', 'a3', 'b3')

  seVect <- sqrt(diag(solve(H)))
  zVect <- (theta - against) / seVect
  pvalVect <- sideAdjustment * pnorm(-abs(zVect))
  outputDf <- data.frame(Estimate = theta, StandErr = seVect, z = zVect,
                         pval = pvalVect, against = against,
                         sideAdj = sideAdjustment)
  rownames(outputDf) <- varName
  return(outputDf)
}

#' Calculate the age filtration density
#'
#' @param xcalc The vector of ages for calculation
#' @param x_mid The midpoint age for the filtration density
#' @param infant_prop The proportion of infants in the filtration density
#' @param discrete Whether to return a discrete density [default: TRUE]
#'
#' @return The age filtration density
#'
#' @export
calc_filtration_density <- function(xcalc, x_mid, infant_prop, discrete = T) {
  z_mid <- 1 / x_mid / (1 + infant_prop)
  z_inf <- infant_prop * z_mid

  f <- rep(0, length(xcalc))
  ind_lo <- which(xcalc <= x_mid)
  f[ind_lo] <- z_inf + (z_mid - z_inf) * xcalc[ind_lo] / x_mid

  ind_hi <- which(xcalc > x_mid & (xcalc <= 2 * x_mid))
  f[ind_hi] <- z_mid + (z_inf - z_mid) * (xcalc[ind_hi] - x_mid) / x_mid

  if (discrete) {
    f <- f / sum(f)
  }
  return(f)
}

#' Sample from the Usher 3 distribution using importance sampling
#'
#' This function generates samples from the Usher 3 distribution using
#' importance sampling. It allows for optional age filtration based on a
#' specified midpoint age and infant proportion.
#'
#' @param N The number of samples to generate.
#' @param th A vector of parameters for the Usher model.
#' @param dx The spacing for the discretized sampling.
#' @param xmax The maximum age for the discretized sampling.
#' @param x_mid The midpoint age for age filtration (optional).
#' @param infant_prop The proportion of infants for age filtration (optional).
#' @param area_tol The tolerance for checking the area under the density
#'   curves (default: 1e-6).
#'
#' @return A list containing the sampled ages (x), illness indicators (ill),
#'   and the calculated densities (rho1 and rho2).
#'
#' @examples
#' th <- c(2e-2, 1.2, 0.175, 1.40, 0.368, 0.01, 0.075, 0.001, 0.917, 0.1)
#' sample_usher3(100, th, 0.01, 120)
#'
#' @export
sample_usher3 <- function(N, th, dx, xmax, x_mid = NA, infant_prop = NA,
                          area_tol = 1e-6) {
  # Set up the sampling grid and calculate the densities
  filter_by_age <- !is.na(x_mid)
  k1 <- th[1]
  k2 <- th[2]
  b_siler <- th[3:7]

  # Create a sampling grid, xcalc, from xmax and dx
  # We use these to check that dx is small enough and/or xmax is large enough
  # We do so by ensuring that the unfiltered density integrates to 1 and,
  # if necessary, that the age filtration density integrates to 1
  # We also use xmax and the maximum value of the evaluated density to set
  # the limits of the uniform rectangular sampling used in importance sampling
  xcalc <- seq(0, xmax, by = dx)
  rho1 <- usher3_rho1(xcalc, k1, b_siler)
  rho2 <- usher3_rho2(xcalc, k1, k2, b_siler)

  # Treat NA as 0 in rho1 and rho2
  rho1[is.na(rho1)] <- 0
  rho2[is.na(rho2)] <- 0

  if (filter_by_age) {
    f_filt <- calc_filtration_density(xcalc,
                                      x_mid,
                                      infant_prop,
                                      discrete=FALSE)
  }

  integration_weights <- calc_trapez_weights(xcalc)
  area1 <- rho1 %*% integration_weights
  area2 <- rho2 %*% integration_weights
  total_area <- area1 + area2

  if (abs(total_area - 1) >= area_tol) {
    stop("total_area of rho1 + rho2 is not 1 to acceptable tolerance (",
         "total_area = ", total_area, "). Decrease dx or increase xmax.")
  }
  
  if (filter_by_age) {
    area_filt <- f_filt %*% integration_weights
    if (abs(area_filt - 1) >= area_tol) {
      stop("Age filtration density does not integrate to 1 to acceptable ",
           "tolerance (area_filt = ", area_filt, "). Decrease dx or increase ",
           "xmax.")
    }
  }

  # Importance sampling
  #
  # For the importance sampling, we first sample uniformly on the rectangle
  # [0, x_lim] x [0, y_lim]. To ensure the accuracy of the sampling, we set
  # x_lim to (1+alpha)*xmax and y_lim to (1+alpha)*max(y), where alpha=0.1 and
  # y is the total value of the density, possibly normalized by the age
  # filtration density (we do not rescale this normalized density since the
  # importance sampling is impervious to the overall y-scale of the constituent
  # curves).
  alpha <- 0.1
  x_lim <- (1 + alpha) * xmax
  y <- rho1 + rho2
  if (filter_by_age) {
    y <- y * f_filt
  }
  y_lim <- (1 + alpha) * max(y)

  # Initialize objects to store the samples
  x <- numeric(N)
  ill <- logical(N)
  n_sampled <- 0

  # Enter a while loop to do the actual importance sampling
  while (n_sampled < N) {
    x_samp <- runif(1, min = 0, max = x_lim)
    y_samp <- runif(1, min = 0, max = y_lim)
    y1 <- usher3_rho1(x_samp, k1, b_siler)
    y2 <- usher3_rho2(x_samp, k1, k2, b_siler)
    # Treat NA as 0 in y1 and y2
    if (is.na(y1)) {
      y1 <- 0
    }
    if (is.na(y2)) {
      y2 <- 0
    }

    if (filter_by_age) {
      f_filt <- calc_filtration_density(x_samp,
                                        x_mid,
                                        infant_prop,
                                        discrete=FALSE)
      y1 <- y1 * f_filt
      y2 <- y2 * f_filt
    }

    if (y_samp <= y1) {
      n_sampled <- n_sampled + 1
      x[n_sampled] <- x_samp
      ill[n_sampled] <- FALSE
    } else if (y_samp <= y1 + y2) {
      n_sampled <- n_sampled + 1
      x[n_sampled] <- x_samp
      ill[n_sampled] <- TRUE
    }
  }

  return(list(x = x, ill = ill, rho1 = rho1, rho2 = rho2))
}

usher3_nll_transform <- function(th_bar, x, ill) {
  # A wrapper function with a transformation to ensure that parameters are
# positive
  th <- exp(th_bar)
  return(nll_usher3(th, x, ill))
}

#' Temper and tune the Usher 3 model
#'
#' @param obj_fun The objective function
#' @param th0 The initial parameter vector
#' @param verbose Whether to print verbose output [default: FALSE]
#' @param fn_plot The function for plotting [default: NULL]
#' @param num_cyc The number of cycles for tempering [default: 100]
#' @param samps_per_cyc The number of samples per cycle [default: 20]
#' @param temp_vect The temperature vector for tempering
#'   [default: 10^(rev(seq(-1, 1, by = 0.25)))]
#' @param prop_scale_mat The proposal scaling matrix [default: NULL]
#' @param lr The learning rate for gradient descent [default: 1e-5]
#' @param func_tol The tolerance for the objective function [default: 1e-6]
#' @param miniter The minimum number of iterations [default: 1]
#' @param maxiter The maximum number of iterations [default: 1000]
#' @param report_period The reporting period for verbose output [default: 50]
#' @param ... Additional arguments passed to obj_fun
#'
#' @return A list containing the results of tempering and tuning
#'
#' @export
temper_and_tune_usher3 <- function(obj_fun, th0, verbose = FALSE,
                                   fn_plot = NULL, num_cyc = 100,
                                   samps_per_cyc = 20,
                                   temp_vect = 10^(rev(seq(-1, 1, by = 0.25))),
                                   prop_scale_mat = NULL, lr = 1e-5,
                                   func_tol = 1e-6, miniter = 1,
                                   maxiter = 1000, report_period = 50, ...) {
  num_param <- length(th0)

  optional_inputs <- list(verbose = verbose, fn_plot = fn_plot,
                          num_cyc = num_cyc, samps_per_cyc = samps_per_cyc,
                          temp_vect = temp_vect, prop_scale_mat = prop_scale_mat,
                          lr = lr, func_tol = func_tol, miniter = miniter,
                          maxiter = maxiter, report_period = report_period)

  if (is.null(prop_scale_mat)) {
    just_scale <- t(replicate(num_param,
                              rev(seq(0.001, 0.1, len = length(temp_vect)))))
    just_param <- replicate(length(temp_vect), th0)
    prop_scale_mat <- just_scale * just_param
  }
  if (verbose) {
    print('Starting parameter vector th0 = ')
    print(th0)
    print('Starting objective function value eta = ')
    print(obj_fun(th0, ...))
  }
  temper <- par_temper(th0, obj_fun, temp_vect = temp_vect,
                       prop_scale = prop_scale_mat, num_cyc = num_cyc,
                       samps_per_cyc = samps_per_cyc, fn_plot = fn_plot,
                       verbose = verbose, report_period = report_period, ...)

  n <- which.min(unlist(lapply(temper$chains, function(x) {x$eta_best})))
  th_temper <- temper$chains[[n]]$theta_best

  if (verbose) {
    print('Best parameter vector from initial tempering:')
    print(th_temper)
    print('The corresponding best value of the objective function:')
    print(temper$chains[[n]]$eta_best)
  }

  # TODO: for now return just the result of the tempering;
  #       need to identify a good approach for tuning
  #       (see TODO above regarding the gradient).
  return(list(obj_fun = obj_fun, th0 = th0, optional_inputs = optional_inputs,
              temper = temper, th_temper = th_temper, th = th_temper))
}

#' Calculate trapezoidal integration weights
#'
#' This function calculates the weights for trapezoidal integration given a
#' vector of grid points.
#'
#' @param tau A vector of grid points.
#'
#' @return A vector of trapezoidal integration weights.
#'
#' @details The function computes the weights for trapezoidal integration
#'   based on the provided grid points. The weights for the central points
#'   are calculated as the average of the differences between the adjacent
#'   grid points. The weights for the first and last points are calculated
#'   as half the difference between the first two and last two grid points,
#'   respectively.
#'
#' @examples
#' tau <- seq(0, 1, by = 0.1)
#' weights <- calc_trapez_weights(tau)
#'
#' @export
calc_trapez_weights <- function(tau) {
  G <- length(tau)
  weight_vect <- rep(NA, length(tau))
  
  # Calculate weights for central points
  ind_cent <- 2:(G - 1)
  weight_vect[ind_cent] <- (tau[ind_cent + 1] - tau[ind_cent - 1]) / 2
  
  # Calculate weight for the first point
  weight_vect[1] <- (tau[2] - tau[1]) / 2
  
  # Calculate weight for the last point
  weight_vect[G] <- (tau[G] - tau[G - 1]) / 2
  
  return(weight_vect)
}