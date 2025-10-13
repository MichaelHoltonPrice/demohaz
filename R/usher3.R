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

#' Calculate the cumulative hazard for the well-to-ill transition
#'
#' This function computes the cumulative hazard for the transition from the
#' healthy (well) state to the ill state over the age interval [x0, x]. It
#' supports a single cutoff age x_cut, after which the hazard becomes zero.
#' This is useful for modeling conditions like Linear Enamel Hypoplasia (LEH)
#' that can only occur during specific age windows (e.g., 0-6 years for LEH).
#'
#' @param x The age(s) at which to evaluate the cumulative hazard
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the hazard becomes zero [default: Inf].
#'   When x_cut = Inf, this reduces to the standard constant hazard model.
#'
#' @return The cumulative hazard evaluated at the location(s) in x
#'
#' @details
#' For each age x, the cumulative hazard is calculated as:
#' - If x0 >= x_cut: H12 = 0 (cutoff already passed)
#' - If x0 < x_cut and x <= x_cut: H12 = k1 * (x - x0)
#' - If x0 < x_cut and x > x_cut: H12 = k1 * (x_cut - x0)
#'
#' The function is vectorized over x.
#'
#' Future extensions could support piecewise constant hazards with multiple
#' breakpoints and different hazard rates in each interval.
#'
#' @export
H12 <- function(x, k1, x0 = 0, x_cut = Inf) {
  if (k1 < 0) {
    stop('k1 cannot be negative')
  }
  if (x0 < 0) {
    stop('x0 cannot be negative')
  }
  if (x_cut <= 0) {
    stop('x_cut must be positive')
  }
  
  # If starting age is at or past the cutoff, no transitions can occur
  if (is.finite(x_cut) && x0 >= x_cut) {
    return(rep(0, length(x)))
  }
  
  # For ages before cutoff: hazard accumulates normally
  # For ages after cutoff: hazard stops accumulating at cutoff
  cumhaz <- k1 * pmin(pmax(x - x0, 0), x_cut - x0)
  
  return(cumhaz)
}

#' Calculate the survival function for the well-to-ill transition
#'
#' This function computes the survival function (probability of remaining in
#' the healthy state) for the transition from well to ill over the age interval
#' [x0, x]. It supports a single cutoff age x_cut, after which the hazard
#' becomes zero.
#'
#' @param x The age(s) at which to evaluate the survival function
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the hazard becomes zero [default: Inf].
#'   When x_cut = Inf, this reduces to the standard constant hazard model.
#'
#' @return The survival probability evaluated at the location(s) in x
#'
#' @details
#' The survival function is S12(x) = exp(-H12(x)), where H12 is the cumulative
#' hazard calculated by the H12 function.
#'
#' When x_cut = Inf, this gives S12(x) = exp(-k1 * (x - x0)), which matches
#' the standard constant hazard model currently used in usher3_rho1 and
#' usher3_integrand.
#'
#' Future extensions could support piecewise constant hazards with multiple
#' breakpoints and different hazard rates in each interval.
#'
#' @export
S12 <- function(x, k1, x0 = 0, x_cut = Inf) {
  return(exp(-H12(x, k1, x0, x_cut)))
}

#' Calculate the transition probability p_11(x0, x)
#'
#' This function computes the probability of remaining in the well state from
#' age x0 to age x, accounting for both the well-to-ill transition (with
#' cutoff) and the well-to-dead transition (baseline mortality).
#'
#' @param x The age at which to evaluate the transition probability
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param b_siler The parameter vector for the Siler hazard model (baseline
#'   mortality)
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @return The transition probability p_11(x0, x)
#'
#' @details
#' Computes the probability of staying in the well state from age x0 to age x.
#' When x0 = 0, this is: p_11(0, x) = exp(-k1 * min(x_cut, x)) * S_13(0, x)
#'
#' This is the probability of staying in the well state, which requires 
#' surviving both the well-to-ill transition and the baseline mortality hazard.
#'
#' @export
p_11 <- function(x, k1, b_siler, x0 = 0, x_cut = Inf) {
  if (x0 < 0) {
    stop('x0 cannot be negative')
  }
  if (any(x < x0)) {
    stop('x must be >= x0')
  }
  if (k1 < 0) {
    stop('k1 cannot be negative')
  }
  
  # S_13(x0, x) is the survival from baseline mortality
  S13_x0_x <- ssiler(x, b_siler, x0 = x0)
  
  # Survival from well-to-ill transition from x0 to x
  # Accumulate hazard only up to x_cut
  if (x_cut <= x0) {
    # No hazard after cutoff
    S12_x0_x <- rep(1, length(x))
  } else {
    # Hazard accumulates from x0 to min(x, x_cut)
    x_limit <- pmin(x_cut, x)
    S12_x0_x <- exp(-k1 * (x_limit - x0))
  }
  
  # Combined probability of staying in well state
  return(S12_x0_x * S13_x0_x)
}

#' Integrand for p_12 calculation
#'
#' Helper function for computing the integral in p_12(x0, x).
#'
#' @param y The integration variable (age of illness onset)
#' @param x0 The starting age
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#' @param b_siler The parameter vector for the Siler hazard model
#'
#' @return The integrand value at y
#'
#' @details
#' The integrand is: exp(-k1 * (y - x0)) * [S_13(x0, y)]^(1 - k2)
#'
#' @keywords internal
p_12_integrand <- function(y, x0, k1, k2, b_siler) {
  S13_x0_y <- ssiler(y, b_siler, x0 = x0)
  return(exp(-k1 * (y - x0)) * (S13_x0_y)^(1 - k2))
}

#' Calculate the transition probability p_12(x0, x)
#'
#' This function computes the probability of being in the ill state at age x,
#' given that the individual started in the well state at age x0. This accounts
#' for transitioning from well to ill at some intermediate age y and then
#' surviving in the ill state to age x.
#'
#' @param x The age at which to evaluate the transition probability
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#' @param b_siler The parameter vector for the Siler hazard model (baseline
#'   mortality)
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @return The transition probability p_12(x0, x)
#'
#' @details
#' Computes the probability of transitioning from well to ill between ages x0 
#' and x. When x0 = 0: p_12(0, x) = k1 * [S_13(0, x)]^k2 * integral from 0 to 
#' min(x_cut, x) of exp(-k1 * y) * [S_13(0, y)]^(1 - k2) dy
#'
#' The integral is computed numerically using R's integrate() function.
#'
#' @export
p_12 <- function(x, k1, k2, b_siler, x0 = 0, x_cut = Inf) {
  if (x0 < 0) {
    stop('x0 cannot be negative')
  }
  if (any(x < x0)) {
    stop('x must be >= x0')
  }
  if (k1 < 0) {
    stop('k1 cannot be negative')
  }
  if (k2 < 0) {
    stop('k2 cannot be negative')
  }
  
  # Initialize result vector
  result <- numeric(length(x))
  
  # Vectorized check: where x == x0, result is 0 (no time to transition)
  is_at_x0 <- (x == x0)
  result[is_at_x0] <- 0
  
  # Vectorized check: if starting after cutoff, no transitions can occur
  if (x0 >= x_cut) {
    return(result)  # All zeros
  }
  
  # Vectorized check: where upper_limit <= x0, result is 0
  upper_limits <- pmin(x_cut, x)
  no_transition <- (upper_limits <= x0) | is_at_x0
  
  # Only process x values where transitions can occur
  x_to_process <- x[!no_transition]
  if (length(x_to_process) == 0) {
    return(result)
  }
  
  # Vectorized calculation of S_13(x0, x)^k2 for all x values to process
  S13_x0_x <- ssiler(x_to_process, b_siler, x0 = x0)
  S13_x0_x_k2 <- S13_x0_x^k2
  
  # Compute integrals (this part must be done iteratively since 
  # each x value has a different upper limit)
  indices_to_process <- which(!no_transition)
  for (i in seq_along(indices_to_process)) {
    idx <- indices_to_process[i]
    upper_limit <- upper_limits[idx]
    
    # Compute the integral from x0 to upper_limit
    integral_result <- tryCatch(
      integrate(
        p_12_integrand,
        lower = x0,
        upper = upper_limit,
        x0 = x0,
        k1 = k1,
        k2 = k2,
        b_siler = b_siler
      )$value,
      error = function(e) {
        warning(paste("Integration failed in p_12 at x =", x[idx], ":", e$message))
        return(NA)
      }
    )
    
    if (is.na(integral_result)) {
      result[idx] <- NA
    } else {
      result[idx] <- k1 * S13_x0_x_k2[i] * integral_result
    }
  }
  
  return(result)
}

#' Calculate the state weights at age x0
#'
#' This function computes the proportion of living individuals who are in the
#' well state (w_1) and ill state (w_2) at age x0, given that all individuals
#' start in the well state at age 0. These weights reflect disease etiology
#' from birth to age x0.
#'
#' @param x0 The age at which to calculate the weights
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#' @param b_siler The parameter vector for the Siler hazard model (baseline
#'   mortality)
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @return A vector of length 2 containing [w_1(x0), w_2(x0)], where:
#'   w_1(x0) is the proportion of living individuals who are well at age x0
#'   w_2(x0) is the proportion of living individuals who are ill at age x0
#'
#' @details
#' Computes the state weights from birth (age 0) to age x0:
#' w_1(x0) = p_11(0, x0) / [p_11(0, x0) + p_12(0, x0)]
#' w_2(x0) = p_12(0, x0) / [p_11(0, x0) + p_12(0, x0)]
#'
#' These weights sum to 1 and represent the state distribution at age x0 for
#' individuals who survive to that age, assuming all individuals start in the 
#' well state at birth. They are used in calculating the weighted age-at-death 
#' densities for analyses restricted to individuals who survive to age x0.
#'
#' At age 0, returns c(1, 0) since everyone starts in the well state.
#'
#' @examples
#' # Example parameter values
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' k1 <- 0.02
#' k2 <- 1.2
#' 
#' # Calculate weights at age 10 (from birth)
#' calc_weights(x0 = 10, k1 = k1, k2 = k2, b_siler = b_siler)
#'
#' @export
calc_weights <- function(x0, k1, k2, b_siler, x_cut = Inf) {
  if (x0 < 0) {
    stop('x0 cannot be negative')
  }
  
  # At age 0, everyone is in the well state
  if (x0 == 0) {
    return(c(w1 = 1, w2 = 0))
  }
  
  # Calculate transition probabilities from age 0 to x0
  p11_val <- p_11(x0, k1, b_siler, x0 = 0, x_cut)
  p12_val <- p_12(x0, k1, k2, b_siler, x0 = 0, x_cut)
  
  # Check for integration failure
  if (is.na(p12_val)) {
    warning("Failed to compute p_12, returning NA weights")
    return(c(w1 = NA, w2 = NA))
  }
  
  # Total survival probability from age 0 to x0
  total_survival <- p11_val + p12_val
  
  # Check for numerical issues
  if (total_survival <= 0) {
    warning("Total survival probability is non-positive, returning NA weights")
    return(c(w1 = NA, w2 = NA))
  }
  
  # Calculate weights (proportions)
  w1 <- p11_val / total_survival
  w2 <- p12_val / total_survival
  
  return(c(w1 = w1, w2 = w2))
}

#' Calculate the occupancy probability q_1(x0, x; w)
#'
#' This function computes the probability of being in state 1 (well) at age x,
#' given the state weights at age x0. This accounts for survival from x0 to x
#' while remaining in the well state.
#'
#' @param x The age at which to evaluate the occupancy probability (can be vector)
#' @param w A vector of length 2 containing [w_1(x0), w_2(x0)], the state 
#'   weights at age x0. If NULL, computed using calc_weights.
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param b_siler The parameter vector for the Siler hazard model (baseline
#'   mortality)
#' @param x0 The starting age (default: 0)
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#' @param k2 Optional k2 parameter (only needed if w is NULL)
#'
#' @return The occupancy probability q_1(x0, x; w) evaluated at x
#'
#' @details
#' Computes the occupancy probability based on three cases:
#' - Case x_cut <= x0 <= x: q_1 = w_1(x0) * S_13(0,x) / S_13(0,x0)
#' - Case x0 <= x_cut <= x: q_1 = w_1(x0) * exp(-k1*(x_cut - x0)) * S_13(0,x) / S_13(0,x0)
#' - Case x0 <= x <= x_cut: q_1 = w_1(x0) * exp(-k1*(x - x0)) * S_13(0,x) / S_13(0,x0)
#'
#' @export
q_1 <- function(x, w, k1, b_siler, x0 = 0, x_cut = Inf, k2 = NULL) {
  if (x0 < 0) {
    stop('x0 cannot be negative')
  }
  if (any(x < x0)) {
    stop('x must be >= x0')
  }
  if (k1 < 0) {
    stop('k1 cannot be negative')
  }
  
  # If w is not provided, calculate it
  if (is.null(w)) {
    if (is.null(k2)) {
      stop('k2 must be provided if w is NULL')
    }
    w <- calc_weights(x0, k1, k2, b_siler, x_cut)
  }
  
  w1 <- w[1]
  
  # Vectorized calculation
  result <- numeric(length(x))
  
  # S_13(0, x) / S_13(0, x0) term (vectorized)
  S13_0_x <- ssiler(x, b_siler, x0 = 0)
  S13_0_x0 <- ssiler(x0, b_siler, x0 = 0)
  S13_ratio <- S13_0_x / S13_0_x0
  
  # Determine which case applies for each x value
  # Case 1: x_cut <= x0 <= x (no hazard from x0 to x)
  case1 <- (x_cut <= x0)
  
  # Case 2: x0 <= x_cut <= x (hazard from x0 to x_cut, then none)
  case2 <- (x0 < x_cut) & (x >= x_cut)
  
  # Case 3: x0 <= x <= x_cut (hazard throughout)
  case3 <- (x0 < x_cut) & (x < x_cut)
  
  # Case 1: no transitions possible
  result[case1] <- w1 * S13_ratio[case1]
  
  # Case 2: transitions from x0 to x_cut
  result[case2] <- w1 * exp(-k1 * (x_cut - x0)) * S13_ratio[case2]
  
  # Case 3: transitions from x0 to x
  result[case3] <- w1 * exp(-k1 * (x[case3] - x0)) * S13_ratio[case3]
  
  return(result)
}

#' Integrand for q_2 calculation
#'
#' Helper function for computing the integral in q_2(x0, x; w).
#'
#' @param y The integration variable (age of illness onset)
#' @param x0 The starting age
#' @param k1 The constant transition rate from the healthy state to the ill state
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#' @param b_siler The parameter vector for the Siler hazard model
#'
#' @return The integrand value at y
#'
#' @keywords internal
q_2_integrand <- function(y, x0, k1, k2, b_siler) {
  S13_0_y <- ssiler(y, b_siler, x0 = 0)
  return(exp(-k1 * (y - x0)) * (S13_0_y)^(1 - k2))
}

#' Calculate the occupancy probability q_2(x0, x; w)
#'
#' This function computes the probability of being in state 2 (ill) at age x,
#' given the state weights at age x0. This accounts for individuals who were
#' already ill at x0 and survived, plus individuals who transitioned from well
#' to ill between x0 and x.
#'
#' @param x The age at which to evaluate the occupancy probability (can be vector)
#' @param w A vector of length 2 containing [w_1(x0), w_2(x0)], the state 
#'   weights at age x0. If NULL, computed using calc_weights.
#' @param k1 The constant transition rate from the healthy state to the ill
#'   state (before x_cut)
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#' @param b_siler The parameter vector for the Siler hazard model (baseline
#'   mortality)
#' @param x0 The starting age (default: 0)
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#'
#' @return The occupancy probability q_2(x0, x; w) evaluated at x
#'
#' @details
#' Computes the occupancy probability based on three cases with integral terms
#' for individuals transitioning from well to ill between x0 and min(x, x_cut).
#'
#' @export
q_2 <- function(x, w, k1, k2, b_siler, x0 = 0, x_cut = Inf) {
  if (x0 < 0) {
    stop('x0 cannot be negative')
  }
  if (any(x < x0)) {
    stop('x must be >= x0')
  }
  if (k1 < 0) {
    stop('k1 cannot be negative')
  }
  if (k2 < 0) {
    stop('k2 cannot be negative')
  }
  
  # If w is not provided, calculate it
  if (is.null(w)) {
    w <- calc_weights(x0, k1, k2, b_siler, x_cut)
  }
  
  w1 <- w[1]
  w2 <- w[2]
  
  # Initialize result vector
  result <- numeric(length(x))
  
  # Common terms
  S13_0_x0 <- ssiler(x0, b_siler, x0 = 0)
  
  # Process each x value
  for (i in seq_along(x)) {
    x_val <- x[i]
    
    # S_13 terms
    S13_0_x <- ssiler(x_val, b_siler, x0 = 0)
    S13_ratio_k2 <- (S13_0_x / S13_0_x0)^k2
    
    # First term: individuals already ill at x0
    term1 <- w2 * S13_ratio_k2
    
    # Second term: individuals transitioning from well to ill
    # Determine integration limits based on case
    if (x_cut <= x0) {
      # Case 1: no transitions possible
      upper_limit <- x0  # No integration
    } else if (x_val >= x_cut) {
      # Case 2: integrate from x0 to x_cut
      upper_limit <- x_cut
    } else {
      # Case 3: integrate from x0 to x
      upper_limit <- x_val
    }
    
    # Compute integral if there's a range to integrate over
    if (upper_limit > x0) {
      integral_result <- tryCatch(
        integrate(
          q_2_integrand,
          lower = x0,
          upper = upper_limit,
          x0 = x0,
          k1 = k1,
          k2 = k2,
          b_siler = b_siler
        )$value,
        error = function(e) {
          warning(paste("Integration failed in q_2 at x =", x_val, ":", e$message))
          return(NA)
        }
      )
      
      if (is.na(integral_result)) {
        result[i] <- NA
      } else {
        term2 <- w1 * k1 * (S13_0_x^k2) / S13_0_x0 * integral_result
        result[i] <- term1 + term2
      }
    } else {
      # No transitions, only term1
      result[i] <- term1
    }
  }
  
  return(result)
}

#' @rdname usher3
#' @param x The vector of ages
#' @param k1 The transition rate from the healthy state to the ill state
#' @param b_siler The parameter vector for the Siler hazard model
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @return The density rho1 evaluated at the locations in the input vector x
#'
#' @export
usher3_rho1 <- function(x, k1, b_siler, x0 = 0, x_cut = Inf) {
  f13 <- dsiler(x, b_siler, x0)
  S12_val <- S12(x, k1, x0, x_cut)
  rho1 <- f13 * S12_val
  return(rho1)
}

#' @rdname usher3
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#'
#' @return The density rho2 evaluated at the locations in the input vector x
#'
#' @export
usher3_rho2 <- function(x, k1, k2, b_siler, x0 = 0, x_cut = Inf) {
  f13_0_x <- dsiler(x, b_siler)
  S13_0_x <- ssiler(x, b_siler)
  S12_0_x0 <- S12(x0, k1, x0 = 0, x_cut)
  S13_0_x0 <- ssiler(x0, b_siler)

  integralTerm <- rep(NA, length(x))
  for (ii in 1:length(x)) {
    integralTerm[ii] <- tryCatch(
      integrate(
        usher3_integrand, x0, x[ii],
        k1 = k1, k2 = k2, b_siler = b_siler, x_cut = x_cut
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
usher3_integrand <- function(y, k1, k2, b_siler, x_cut = Inf) {
  # _0_y indicates that the relevant quantity is relative to the interval 0 to y
  S12_0_y <- S12(y, k1, x0 = 0, x_cut)
  S13_0_y <- ssiler(y, b_siler)

  return(S12_0_y * (S13_0_y)^(1 - k2))
}

#' Calculate the negative log-likelihood for the Usher 3 illness-death model.
#' This method also supports Gompertz-Makeham mortality, for which the infant
#' mortality term is zero; this is accomplished by detecting if the input
#' parameter vector, theta, has length 5, which implies no infant mortality
#' term.
#'
#' @param theta The parameter vector for the Usher illness-death model with
#'   ordering [k1, k2, mort]. The mortality parameter vector is
#'   [b1, b2, b3, b4, b5], or [b3, b4, b5] for Gompertz-Makeham mortality.
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators, which can have NA entries
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @return The negative log-likelihood value
#'
#' @export
nll_usher3 <- function(theta, x, ill, x0 = 0, x_cut = Inf) {
  k1 <- theta[1]
  if (k1 < 0) {
    return(Inf)
  }

  k2 <- theta[2]
  if (k2 < 0) {
    return(Inf)
  }

  if (length(theta) == 5) {
    # Then assume Gompertz-Makeham mortality and add an infant mortality term
    # that equals zero
    b_siler <- c(0, 1, theta[3:5])
  } else {
    # Then we do have an infant mortality term
    b_siler <- theta[3:7]
  }

  if (any(b_siler < 0)) {
    return(Inf)
  }

  # If the ill state is not known (NA) then we can still use the datapoint
  # by using rho1 + rho2 for the likelihood.
  ind_na <- is.na(ill)
  x_na <- x[ind_na]
  rho2_na_ill <- usher3_rho2(x_na, k1, k2, b_siler, x0, x_cut)
  if (any(is.na(rho2_na_ill))) {
    return(Inf)
  }
  rho1_na_wll <- usher3_rho1(x_na, k1, b_siler, x0, x_cut)
  # rho_na contributes to the  likelihood below
  rho_na <- rho2_na_ill + rho1_na_wll

  # Subet x with the non-missing values calculate rho1 and rho2
  x <- x[!ind_na]
  ill <- ill[!ind_na]
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]

  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0, x_cut)
  if (any(is.na(rho2_ill))) {
    return(Inf)
  }

  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0, x_cut)

  # Calculate and return the negative log-likelihood
  ll <- sum(log(rho_na)) + sum(log(rho1_wll)) + sum(log(rho2_ill))
  return(-ll)
}

#' Calculate the Hessian matrix for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#'
#' @return The Hessian matrix
#'
#' @export
usher3_hessian <- function(theta, x, ill, x0 = 0, x_cut = Inf) {
  # Finite x_cut not yet supported for Hessian calculation
  if (is.finite(x_cut)) {
    stop("Finite x_cut is not yet supported for usher3_hessian. ",
         "Use x_cut = Inf (the default).")
  }
  
  H <- numDeriv::hessian(nll_usher3_hessian_wrapper,
                         theta, method.args = list(eps = 1e-12),
                         ageVect = x,
                         illVect = ill, x0 = x0, x_cut = x_cut)
  return(H)
}

#' Wrapper function for the negative log-likelihood of the Usher 3 model
#'
#' @param paramVect The parameter vector for the Usher illness-death model
#' @param ageVect The vector of ages-at-death
#' @param illVect The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#'
#' @return The negative log-likelihood value
#'
nll_usher3_hessian_wrapper <- function(paramVect, ageVect, illVect, x0 = 0, 
                                       x_cut = Inf) {
  return(nll_usher3(paramVect, ageVect, illVect, x0, x_cut))
}

#' Calculate standard errors, z-scores, and p-values for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#'
#' @return A data frame with standard errors, z-scores, and p-values
#'
#' @export
usher3_errors <- function(theta, x, ill, x0 = 0, x_cut = Inf) {
  # Finite x_cut not yet supported for error estimation
  if (is.finite(x_cut)) {
    stop("Finite x_cut is not yet supported for usher3_errors. ",
         "Use x_cut = Inf (the default).")
  }
  
  H <- usher3_hessian(theta, x, ill, x0, x_cut)
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
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @examples
#' th <- c(2e-2, 1.2, 0.175, 1.40, 0.368, 0.01, 0.075, 0.001, 0.917, 0.1)
#' sample_usher3(100, th, 0.01, 120)
#'
#' @export
sample_usher3 <- function(N, th, dx, xmax, x_mid = NA, infant_prop = NA,
                          area_tol = 1e-6, x0 = 0, x_cut = Inf) {
  # Finite x_cut not yet supported for importance sampling
  if (is.finite(x_cut)) {
    stop("Finite x_cut is not yet supported for sample_usher3. ",
         "Use x_cut = Inf (the default).")
  }
  
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
  rho1 <- usher3_rho1(xcalc, k1, b_siler, x0, x_cut)
  rho2 <- usher3_rho2(xcalc, k1, k2, b_siler, x0, x_cut)

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
    y1 <- usher3_rho1(x_samp, k1, b_siler, x0, x_cut)
    y2 <- usher3_rho2(x_samp, k1, k2, b_siler, x0, x_cut)
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

nll_usher3_optim_wrapper <- function(th_bar, ...) {
  # A wrapper function use as the objective function for
  # temper_usher3
  th <- exp(th_bar)
  return(nll_usher3(th, ...))
}

#' Use parallel tempering to fit the usher3 model.
#'
#' If the flag use_gompertz is True, then the infant mortality hazard is
#' assumed to be zero and the parameter vector, th0, should only have five
# terms (k1, k2, and three mortality parameters). If necessary, th0 is subset
#' from 7 to 5 terms if use_gompertz is TRUE. If the scaling matrix,
#' prop_scale_mat, is specified when use_gommpertz is True its dimensions must
#' match the reduced length of the parameter vector.
#'
#' @param th0 The initial parameter vector with the ordering [k1, k2, b_siler]
#'   where b_siler uses the demohaz parameterization of the Siler hazard (see
#'   Siler documentation)
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
#' @param use_gompertz Whether to use Gompert-Makeham mortality [default: FAlSE]
#' @param ... Additional arguments passed to obj_fun, including x (ages), 
#'   ill (illness indicators), x0 (conditional starting age, default 0), and
#'   x_cut (age cutoff for transitions, default Inf)
#'
#' @return A list containing the results of tempering and tuning
#'
#' @export
temper_usher3 <- function(th0 = c(1e-2,
                                  1.0,
                                  0.175,
                                  1.40,
                                  .368 * .01,
                                  log(.917 * .1/(.075 * .001))/(.917 * .1),
                                  .917 * .1),
                          verbose = FALSE,
                          fn_plot = NULL, num_cyc = 100,
                          samps_per_cyc = 20,
                          temp_vect = 10^(rev(seq(-1, 1, by = 0.25))),
                          prop_scale_mat = NULL, lr = 1e-5,
                          func_tol = 1e-6, miniter = 1,
                          maxiter = 1000, report_period = 50,
                          use_gompertz=FALSE, ...) {
  obj_fun <- nll_usher3_optim_wrapper
  if (!use_gompertz) {
    if (length(th0) != 7) {
      stop('th0 should have length 7')
    }
  } else {
    if (!(length(th0) %in% c(5,7))) {
      stop('th0 should have length 5 or 7 if use_gompertz is TRUE')
    }
    if (length(th0) == 7) {
      # If necessary, remove the infant mortality term from th0
      th0 <- th0[c(1,2,5,6,7)]
    }
  }

  num_param <- length(th0)

  th0_bar <- log(th0)

  optional_inputs <- list(verbose = verbose, fn_plot = fn_plot,
                          num_cyc = num_cyc, samps_per_cyc = samps_per_cyc,
                          temp_vect = temp_vect, prop_scale_mat = prop_scale_mat,
                          lr = lr, func_tol = func_tol, miniter = miniter,
                          maxiter = maxiter, report_period = report_period)

  if (is.null(prop_scale_mat)) {
    just_scale <- t(replicate(num_param,
                              rev(seq(0.001, 0.1, len = length(temp_vect)))))
    just_param <- replicate(length(temp_vect), th0_bar)
    prop_scale_mat <- just_scale * just_param
  }
  if (verbose) {
    print('Starting parameter vector th0 = ')
    print(th0)
    print('Starting parameter vector th0_bar = ')
    print(th0_bar)
    print('Starting objective function value eta = ')
    print(obj_fun(th0_bar, ...))
  }
  temper <- par_temper(th0_bar, obj_fun, temp_vect = temp_vect,
                       prop_scale = prop_scale_mat, num_cyc = num_cyc,
                       samps_per_cyc = samps_per_cyc, fn_plot = fn_plot,
                       verbose = verbose, report_period = report_period, ...)

  n <- which.min(unlist(lapply(temper$chains, function(x) {x$eta_best})))
  th_bar_temper <- temper$chains[[n]]$theta_best
  th_temper <- exp(th_bar_temper)

  if (verbose) {
    print('Best parameter vector from initial tempering (th):')
    print(th_temper)
    print('Best parameter vector from initial tempering (th_bar):')
    print(th_bar_temper)
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