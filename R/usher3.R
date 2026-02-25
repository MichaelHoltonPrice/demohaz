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
#'   the sum of the two, rho1 + rho2, that integrates to 1; hence, we
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
#' @examples
#' # Standard constant hazard (x_cut = Inf)
#' H12(c(5, 10, 20), k1 = 0.02)
#'
#' # With cutoff at age 6 (e.g., LEH)
#' H12(c(3, 6, 10), k1 = 0.02, x_cut = 6)
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
#' @examples
#' # Probability of remaining healthy (not transitioning to ill)
#' S12(c(5, 10, 20), k1 = 0.02)
#'
#' # With cutoff: survival plateaus after x_cut
#' S12(c(3, 6, 10), k1 = 0.02, x_cut = 6)
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
#' surviving both the well-to-ill transition and the baseline mortality hazard
#' (S_13, the Siler survival function).
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' p_11(c(10, 30, 60), k1 = 0.02, b_siler = b_siler)
#'
#' # With cutoff: well-to-ill transitions stop after age 6
#' p_11(c(3, 6, 10), k1 = 0.02, b_siler = b_siler, x_cut = 6)
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
#' The integrand is: exp(-k1 * (y - x0)) * [S_13(x0, y)]^(1 - k2).
#' This arises from the product of the well-to-ill survival to age y and the
#' ratio of Siler survival functions raised to the (1 - k2) power, which
#' accounts for the differential mortality between the well and ill states.
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
#' and x. The formula is:
#'
#' p_12(x0, x) = k1 * [S_13(x0, x)]^k2 * integral from x0 to min(x_cut, x)
#' of exp(-k1 * (y - x0)) * [S_13(x0, y)]^(1 - k2) dy
#'
#' The integral is computed numerically using R's integrate() function. Each
#' element of x may have a different upper integration limit (min(x_cut, x[i])),
#' so the integral must be evaluated iteratively. If integrate() fails for a
#' given x value, NA is returned for that element.
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' p_12(c(10, 30, 60), k1 = 0.02, k2 = 1.2, b_siler = b_siler)
#'
#' # With cutoff: no new transitions after age 6
#' p_12(c(3, 6, 10), k1 = 0.02, k2 = 1.2, b_siler = b_siler, x_cut = 6)
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
        # Silently return NA - caller handles this
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
#' Returns c(NA, NA) if the numerical integration in p_12 fails or if the
#' total survival probability is zero or negative.
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
  
  # Check for integration failure - silently return NA
  if (is.na(p12_val)) {
    return(c(w1 = NA, w2 = NA))
  }
  
  # Total survival probability from age 0 to x0
  total_survival <- p11_val + p12_val
  
  # Check for numerical issues - silently return NA
  if (total_survival <= 0) {
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
#' @param w Optional weight vector [w_1(x0), w_2(x0)], the proportion of
#'   living individuals who are well vs ill at age x0. If NULL (default),
#'   weights are computed using calc_weights, which assumes disease etiology
#'   unfolds from birth. To assume everyone starts in the well state at x0,
#'   pass w = c(1, 0).
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
#' Computes the occupancy probability based on the relationship between x, x0,
#' and x_cut. All cases share the baseline mortality survival ratio
#' S_13(0, x) / S_13(0, x0), which captures the probability of surviving
#' background mortality from x0 to x. The cases differ in the well-to-ill
#' survival term:
#'
#' - Case x_cut <= x0 (cutoff already passed): q_1 = w_1 * S_13(0,x) / S_13(0,x0).
#'   No well-to-ill transitions possible, so only baseline mortality matters.
#' - Case x0 < x_cut <= x (cutoff within interval): q_1 = w_1 * exp(-k1*(x_cut - x0)) *
#'   S_13(0,x) / S_13(0,x0). Hazard accumulates from x0 to x_cut, then stops.
#' - Case x0 <= x < x_cut (entirely before cutoff): q_1 = w_1 * exp(-k1*(x - x0)) *
#'   S_13(0,x) / S_13(0,x0). Hazard accumulates over the full interval.
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' w <- c(1, 0)  # everyone starts well
#' q_1(c(10, 30, 60), w = w, k1 = 0.02, b_siler = b_siler)
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
#' @details
#' The integrand is: exp(-k1 * (y - x0)) * [S_13(0, y)]^(1 - k2).
#' This has the same mathematical form as p_12_integrand but uses S_13(0, y)
#' (survival from birth) rather than S_13(x0, y) (survival from x0), because
#' the q_2 formula factors the S_13 terms differently.
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
#' @param w Optional weight vector [w_1(x0), w_2(x0)], the proportion of
#'   living individuals who are well vs ill at age x0. If NULL (default),
#'   weights are computed using calc_weights, which assumes disease etiology
#'   unfolds from birth. To assume everyone starts in the well state at x0,
#'   pass w = c(1, 0).
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
#' The occupancy probability q_2 has two additive components:
#' - Term 1: w_2 * [S_13(0, x) / S_13(0, x0)]^k2 — individuals already ill at
#'   x0 who survive to age x under the elevated ill-state mortality (k2 times
#'   baseline).
#' - Term 2: w_1 * k1 * [S_13(0, x)]^k2 / S_13(0, x0) * integral — individuals
#'   who were well at x0 and transitioned to ill at some intermediate age y,
#'   then survived in the ill state to age x.
#'
#' The integration limits for Term 2 depend on three cases:
#' - x_cut <= x0: no transitions possible, Term 2 = 0.
#' - x0 < x_cut <= x: integrate from x0 to x_cut.
#' - x0 <= x < x_cut: integrate from x0 to x.
#'
#' Each element of x is processed iteratively because the upper integration
#' limit varies. If integrate() fails for a given x value, NA is returned.
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' w <- c(1, 0)  # everyone starts well
#' q_2(c(10, 30, 60), w = w, k1 = 0.02, k2 = 1.2, b_siler = b_siler)
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
          # Silently return NA - caller (e.g., nll_usher3) handles this
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
#'   (applies only before x_cut; zero after x_cut)
#' @param b_siler The parameter vector for the Siler hazard model
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#' @param w Optional weight vector [w_1(x0), w_2(x0)]. If NULL (default),
#'   weights are computed using calc_weights, which assumes disease etiology
#'   unfolds from birth. To assume everyone starts in the well state at x0,
#'   pass w = c(1, 0).
#' @param k2 Optional k2 parameter (only needed if w is NULL for weight calculation)
#'
#' @details
#' \code{usher3_rho1} computes the pseudo-density for dying in the well
#' (healthy) state at age x: rho1(x) = q_1(x) * lambda_13(x), where q_1 is
#' the occupancy probability of the well state and lambda_13 = hsiler(x, b)
#' is the baseline Siler mortality hazard.
#'
#' @return The density rho1 evaluated at the locations in the input vector x
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' x <- seq(0, 100, by = 1)
#' rho1 <- usher3_rho1(x, k1 = 0.02, b_siler = b_siler, k2 = 1.2)
#'
#' @export
usher3_rho1 <- function(x, k1, b_siler, x0 = 0, x_cut = Inf, w = NULL, k2 = NULL) {
  # If w is not provided, calculate it (requires k2)
  if (is.null(w)) {
    if (is.null(k2)) {
      stop('k2 must be provided if w is NULL')
    }
    w <- calc_weights(x0, k1, k2, b_siler, x_cut)
  }
  
  # Calculate occupancy probability
  q1_val <- q_1(x, w, k1, b_siler, x0, x_cut, k2)
  
  # Calculate hazard lambda_13(x)
  lambda_13 <- hsiler(x, b_siler)
  
  # rho1 = q_1 * lambda_13
  rho1 <- q1_val * lambda_13
  return(rho1)
}

#' @rdname usher3
#' @param k2 The factor by which the mortality hazard out of the ill state is
#'   larger than that out of the healthy state
#' @param w Optional weight vector [w_1(x0), w_2(x0)]. If NULL (default),
#'   weights are computed using calc_weights, which assumes disease etiology
#'   unfolds from birth. To assume everyone starts in the well state at x0,
#'   pass w = c(1, 0).
#'
#' @details
#' \code{usher3_rho2} computes the pseudo-density for dying in the ill state
#' at age x: rho2(x) = q_2(x) * lambda_23(x), where q_2 is the occupancy
#' probability of the ill state and lambda_23 = k2 * hsiler(x, b) is the
#' elevated mortality hazard for ill individuals.
#'
#' @return The density rho2 evaluated at the locations in the input vector x
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' x <- seq(0, 100, by = 1)
#' rho2 <- usher3_rho2(x, k1 = 0.02, k2 = 1.2, b_siler = b_siler)
#'
#' @export
usher3_rho2 <- function(x, k1, k2, b_siler, x0 = 0, x_cut = Inf, w = NULL) {
  # If w is not provided, calculate it
  if (is.null(w)) {
    w <- calc_weights(x0, k1, k2, b_siler, x_cut)
  }
  
  # Calculate occupancy probability
  q2_val <- q_2(x, w, k1, k2, b_siler, x0, x_cut)
  
  # Calculate hazard lambda_23(x) = k2 * lambda_13(x)
  lambda_13 <- hsiler(x, b_siler)
  lambda_23 <- k2 * lambda_13
  
  # rho2 = q_2 * lambda_23
  rho2 <- q2_val * lambda_23
  return(rho2)
}

#' @rdname usher3
#' @param y The variable of integration (age of illness onset)
#'
#' @details
#' \code{usher3_integrand} computes S12(0, y) * [S_13(0, y)]^(1 - k2), where
#' S12 is the well-to-ill survival function and S_13 is the Siler survival
#' function. The subscript notation "0_y" in the source code indicates that
#' the relevant survival quantities are evaluated over the interval [0, y].
#'
#' @return The integrand evaluated at y
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' usher3_integrand(5, k1 = 0.02, k2 = 1.2, b_siler = b_siler)
#'
#' @export
usher3_integrand <- function(y, k1, k2, b_siler, x_cut = Inf) {
  # _0_y indicates that the relevant quantity is relative to the interval 0 to y
  S12_0_y <- S12(y, k1, x0 = 0, x_cut)
  S13_0_y <- ssiler(y, b_siler)

  return(S12_0_y * (S13_0_y)^(1 - k2))
}

#' Calculate the negative log-likelihood for the Usher 3 illness-death model
#'
#' This method also supports Gompertz-Makeham mortality, for which the infant
#' mortality term is zero; this is accomplished by detecting if the input
#' parameter vector, theta, has length 5, which implies no infant mortality
#' term.
#'
#' @param theta The parameter vector for the Usher illness-death model with
#'   ordering [k1, k2, b_siler]. When length(theta) == 7, b_siler is the
#'   full five-parameter Siler mortality vector [b1, b2, b3, b4, b5] (see
#'   \code{\link{hsiler}}). When length(theta) == 5, the last three elements
#'   are treated as Gompertz-Makeham parameters [b3, b4, b5] and infant
#'   mortality is set to zero (b1 = 0, b2 = 1).
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators (0 = well, 1 = ill). Entries
#'   may be NA, in which case the corresponding observations contribute
#'   rho1 + rho2 to the likelihood (marginalizing over illness status).
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @details
#' The log-likelihood is the sum of log-densities for each observation:
#' - For observations with known illness status: log(rho1(x)) if well,
#'   log(rho2(x)) if ill.
#' - For observations with unknown illness status (ill = NA):
#'   log(rho1(x) + rho2(x)).
#'
#' Returns Inf (indicating an impossible parameter vector) if any parameter
#' is negative or if numerical integration fails.
#'
#' @return The negative log-likelihood value
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' th <- c(0.02, 1.2, b_siler)
#' x <- c(5, 30, 50, 70, 80)
#' ill <- c(FALSE, FALSE, TRUE, TRUE, FALSE)
#' nll_usher3(th, x, ill)
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
  rho1_na_wll <- usher3_rho1(x_na, k1, b_siler, x0, x_cut, k2 = k2)
  # rho_na contributes to the  likelihood below
  rho_na <- rho2_na_ill + rho1_na_wll

  # Subset x to the non-missing values and calculate rho1 and rho2
  x <- x[!ind_na]
  ill <- ill[!ind_na]
  x_wll <- x[ill == 0]
  x_ill <- x[ill == 1]

  rho2_ill <- usher3_rho2(x_ill, k1, k2, b_siler, x0, x_cut)
  if (any(is.na(rho2_ill))) {
    return(Inf)
  }

  rho1_wll <- usher3_rho1(x_wll, k1, b_siler, x0, x_cut, k2 = k2)

  # Calculate and return the negative log-likelihood
  ll <- sum(log(rho_na)) + sum(log(rho1_wll)) + sum(log(rho2_ill))
  return(-ll)
}

#' Calculate the Hessian matrix for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model.
#'   If log_transform=TRUE, this should be log(theta).
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#' @param log_transform If TRUE, compute Hessian with respect to log(theta)
#'   [default: FALSE]
#'
#' @details
#' Uses \code{numDeriv::hessian} to compute the matrix of second partial
#' derivatives of the negative log-likelihood with respect to theta (or
#' log(theta) when log_transform=TRUE). A small step size (eps = 1e-12) is
#' used for the finite-difference approximation. The log-transform option is
#' useful for computing standard errors on the log scale, where the
#' asymptotic normality assumption is more appropriate for positive parameters.
#'
#' @return The Hessian matrix (square, dimension = length(theta))
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' th <- c(0.02, 1.2, b_siler)
#' samp <- sample_usher3(200, th, 0.01, 120)
#' H <- usher3_hessian(th, samp$x, samp$ill)
#'
#' @export
usher3_hessian <- function(theta, x, ill, x0 = 0, x_cut = Inf, 
                           log_transform = FALSE) {
  if (log_transform) {
    H <- numDeriv::hessian(nll_usher3_log_hessian_wrapper,
                           theta, method.args = list(eps = 1e-12),
                           ageVect = x, illVect = ill, x0 = x0, x_cut = x_cut)
  } else {
    # Compute Hessian with respect to untransformed parameters
    H <- numDeriv::hessian(nll_usher3_hessian_wrapper,
                           theta, method.args = list(eps = 1e-12),
                           ageVect = x,
                           illVect = ill, x0 = x0, x_cut = x_cut)
  }
  return(H)
}

#' Calculate the gradient of the negative log-likelihood for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model.
#'   If log_transform=TRUE, this should be log(theta).
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#' @param log_transform If TRUE, compute gradient with respect to log(theta)
#'   [default: FALSE]
#'
#' @details
#' Uses \code{numDeriv::grad} to compute the vector of first partial
#' derivatives of the negative log-likelihood with respect to theta (or
#' log(theta) when log_transform=TRUE). A small step size (eps = 1e-12)
#' is used for the finite-difference approximation. Useful for checking
#' convergence of optimization (gradient should be near zero at the MLE).
#'
#' @return The gradient vector (length = length(theta))
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' th <- c(0.02, 1.2, b_siler)
#' samp <- sample_usher3(200, th, 0.01, 120)
#' g <- usher3_gradient(th, samp$x, samp$ill)
#'
#' @export
usher3_gradient <- function(theta, x, ill, x0 = 0, x_cut = Inf,
                            log_transform = FALSE) {
  if (log_transform) {
    grad_vec <- numDeriv::grad(nll_usher3_log_hessian_wrapper,
                               theta, method.args = list(eps = 1e-12),
                               ageVect = x, illVect = ill, x0 = x0, x_cut = x_cut)
  } else {
    # Compute gradient with respect to untransformed parameters
    grad_vec <- numDeriv::grad(nll_usher3_hessian_wrapper,
                               theta, method.args = list(eps = 1e-12),
                               ageVect = x,
                               illVect = ill, x0 = x0, x_cut = x_cut)
  }
  return(grad_vec)
}

#' Wrapper function for the negative log-likelihood of the Usher 3 model
#'
#' Thin wrapper around \code{\link{nll_usher3}} that uses parameter names
#' (paramVect, ageVect, illVect) that do not collide with \code{numDeriv}'s
#' internal use of \code{x} as the differentiation variable.
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
#' @keywords internal
nll_usher3_hessian_wrapper <- function(paramVect, ageVect, illVect, x0 = 0, 
                                       x_cut = Inf) {
  return(nll_usher3(paramVect, ageVect, illVect, x0, x_cut))
}

#' Wrapper for Hessian/gradient of the log-transformed Usher 3 NLL
#'
#' Like nll_usher3_hessian_wrapper but exponentiates the parameter vector
#' first. Uses ageVect/illVect to avoid name collision with numDeriv's x
#' parameter.
#'
#' @param paramVect The log-transformed parameter vector
#' @param ageVect The vector of ages-at-death
#' @param illVect The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#'
#' @return The negative log-likelihood value
#'
#' @keywords internal
nll_usher3_log_hessian_wrapper <- function(paramVect, ageVect, illVect, x0 = 0,
                                           x_cut = Inf) {
  return(nll_usher3(exp(paramVect), ageVect, illVect, x0, x_cut))
}

#' Calculate standard errors, z-scores, and p-values for the Usher 3 model
#'
#' @param theta The parameter vector for the Usher illness-death model.
#'   If log_transform=TRUE, this should be log(theta).
#' @param x The vector of ages-at-death
#' @param ill The vector of illness indicators
#' @param x0 Conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]
#' @param log_transform If TRUE, compute errors with respect to log(theta)
#'   [default: FALSE]
#'
#' @details
#' Standard errors are computed from the observed Fisher information matrix:
#' SE = sqrt(diag(H^{-1})), where H is the Hessian of the negative
#' log-likelihood. Z-scores test each parameter against a null hypothesis
#' value, and p-values are computed from the standard normal distribution.
#'
#' The null hypothesis values (the "against" column) depend on the
#' parameterization:
#' - Untransformed: k1 tested against 0 (one-sided, is there illness?),
#'   k2 tested against 1 (two-sided, does illness affect mortality?),
#'   Siler parameters tested against 0 (one-sided).
#' - Log-transformed: log(k2) tested against 0 (equivalent to k2 = 1);
#'   log(k1) and log(Siler params) are not tested (against = NA) because
#'   they are necessarily positive and a zero null is not meaningful on the
#'   log scale.
#'
#' The sideAdjustment vector controls one-sided (1) vs two-sided (2) tests.
#' k2 uses a two-sided test because deviations in either direction from 1
#' are meaningful; all other testable parameters use one-sided tests.
#'
#' @return A data frame with rows named by parameter and columns:
#'   \describe{
#'     \item{Estimate}{The parameter estimates (theta values)}
#'     \item{StandErr}{Standard errors from the inverse Hessian}
#'     \item{z}{Z-scores: (estimate - against) / SE. NA for untestable
#'       parameters.}
#'     \item{pval}{P-values from the standard normal. NA for untestable
#'       parameters.}
#'     \item{against}{Null hypothesis value for each parameter}
#'     \item{sideAdj}{Side adjustment (1 = one-sided, 2 = two-sided)}
#'   }
#'
#' @examples
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' th <- c(0.02, 1.2, b_siler)
#' samp <- sample_usher3(500, th, 0.01, 120)
#' usher3_errors(th, samp$x, samp$ill)
#'
#' @export
usher3_errors <- function(theta, x, ill, x0 = 0, x_cut = Inf,
                          log_transform = FALSE) {
  H <- usher3_hessian(theta, x, ill, x0, x_cut, log_transform = log_transform)
  
  if (log_transform) {
    against <- c(NA, 0, NA, NA, NA, NA, NA)
    varName <- c('log_k1', 'log_k2', 'log_a1', 'log_b1', 'log_a2', 'log_a3', 'log_b3')
  } else {
    against <- c(0, 1, 0, 0, 0, 0, 0)
    varName <- c('k1', 'k2', 'a1', 'b1', 'a2', 'a3', 'b3')
  }
  
  sideAdjustment <- c(1, 2, 1, 1, 1, 1, 1)

  seVect <- sqrt(diag(solve(H)))
  zVect <- (theta - against) / seVect
  pvalVect <- sideAdjustment * pnorm(-abs(zVect))
  
  # Where no null hypothesis is meaningful (against = NA), clear z and p
  untestable <- is.na(against)
  zVect[untestable] <- NA
  pvalVect[untestable] <- NA

  outputDf <- data.frame(Estimate = theta, StandErr = seVect, z = zVect,
                         pval = pvalVect, against = against,
                         sideAdj = sideAdjustment)
  rownames(outputDf) <- varName
  return(outputDf)
}

#' Calculate the age filtration density
#'
#' @param xcalc The vector of ages for calculation
#' @param x_mid The mode (peak) age of the triangular density
#' @param infant_prop The ratio of the density at age 0 to the density at
#'   x_mid. Controls the asymmetry of the triangle: infant_prop = 1 gives a
#'   symmetric triangle; infant_prop < 1 gives a density that rises to x_mid.
#' @param discrete Whether to return a discrete density that sums to 1
#'   [default: TRUE]. If FALSE, returns the continuous piecewise-linear density.
#'
#' @details
#' Constructs a piecewise-linear (triangular) density over the interval
#' [0, 2 * x_mid]. The density is zero outside this interval. The peak is at
#' x_mid with height z_mid = 1 / (x_mid * (1 + infant_prop)). At age 0 the
#' density is z_inf = infant_prop * z_mid. The density linearly interpolates
#' between z_inf at 0 and z_mid at x_mid, then between z_mid at x_mid and
#' z_inf at 2 * x_mid. This is used to model the age-at-death distribution
#' of an archaeological assemblage when filtering rejection-sampled data.
#'
#' @return The age filtration density evaluated at xcalc. If discrete=TRUE,
#'   the values sum to 1; if discrete=FALSE, they integrate to approximately 1
#'   (exactly 1 in the continuous limit).
#'
#' @examples
#' xcalc <- seq(0, 120, by = 0.1)
#' f <- calc_filtration_density(xcalc, x_mid = 30, infant_prop = 0.5)
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

#' Sample from the Usher 3 distribution using rejection sampling
#'
#' This function generates samples from the Usher 3 distribution using
#' rejection (acceptance-rejection) sampling. It allows for optional age
#' filtration based on a specified midpoint age and infant proportion.
#'
#' @param N The number of samples to generate.
#' @param th The parameter vector with ordering [k1, k2, b_siler], where
#'   b_siler is the five-parameter Siler mortality vector (see
#'   \code{\link{hsiler}}).
#' @param dx The grid spacing used to evaluate and validate the densities.
#'   Must be small enough that rho1 + rho2 integrates to 1 within area_tol.
#' @param xmax The maximum age for the sampling grid. Must be large enough
#'   that the density is negligible beyond xmax.
#' @param x_mid The midpoint age for age filtration (optional, default NA).
#' @param infant_prop The proportion of infants for age filtration (optional,
#'   default NA).
#' @param area_tol The tolerance for checking that rho1 + rho2 integrates
#'   to 1 [default: 1e-6]. An error is thrown if the check fails.
#' @param x0 The conditional starting age [default: 0]
#' @param x_cut The age at which the well-to-ill transition hazard becomes
#'   zero [default: Inf]. When x_cut = Inf, this reduces to the standard
#'   constant hazard model.
#'
#' @details
#' The algorithm proceeds in two phases:
#'
#' 1. **Validation**: Evaluates rho1 and rho2 on a grid from 0 to xmax with
#'    spacing dx, then checks via trapezoidal integration that rho1 + rho2
#'    integrates to 1 (within area_tol). This ensures that dx is fine enough
#'    and xmax is large enough to capture the full density.
#'
#' 2. **Rejection sampling**: Draws candidate points (x_samp, y_samp) uniformly
#'    from the rectangle [0, 1.1*xmax] x [0, 1.1*max(rho1+rho2)]. If
#'    y_samp <= rho1(x_samp), the point is accepted as a well-at-death sample;
#'    else if y_samp <= rho1(x_samp) + rho2(x_samp), it is accepted as an
#'    ill-at-death sample; otherwise it is rejected. When age filtration is
#'    active, the densities are multiplied by the filtration density before
#'    comparison. NA values in rho1/rho2 (from numerical integration failure)
#'    are treated as 0.
#'
#' @return A list with components:
#'   \describe{
#'     \item{x}{Numeric vector of length N: sampled ages at death}
#'     \item{ill}{Logical vector of length N: TRUE if ill at death}
#'     \item{rho1}{The well-at-death density evaluated on the grid}
#'     \item{rho2}{The ill-at-death density evaluated on the grid}
#'   }
#'
#' @examples
#' th <- c(2e-2, 1.2, 0.175, 1.40, 0.00368, 38.1, 0.0917)
#' samp <- sample_usher3(100, th, 0.01, 120)
#'
#' @export
sample_usher3 <- function(N, th, dx, xmax, x_mid = NA, infant_prop = NA,
                          area_tol = 1e-6, x0 = 0, x_cut = Inf) {
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
  # the limits of the uniform rectangular sampling used in rejection sampling
  xcalc <- seq(0, xmax, by = dx)
  rho1 <- usher3_rho1(xcalc, k1, b_siler, x0, x_cut, k2 = k2)
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

  # Rejection sampling
  #
  # For the rejection sampling, we first sample uniformly on the rectangle
  # [0, x_lim] x [0, y_lim]. To ensure the accuracy of the sampling, we set
  # x_lim to (1+alpha)*xmax and y_lim to (1+alpha)*max(y), where alpha=0.1 and
  # y is the total value of the density, possibly normalized by the age
  # filtration density (we do not rescale this normalized density since the
  # rejection sampling is impervious to the overall y-scale of the constituent
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

  # Enter a while loop to do the actual rejection sampling
  while (n_sampled < N) {
    x_samp <- runif(1, min = 0, max = x_lim)
    y_samp <- runif(1, min = 0, max = y_lim)
    y1 <- usher3_rho1(x_samp, k1, b_siler, x0, x_cut, k2 = k2)
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

#' Wrapper for log-transformed NLL used by temper_and_tune_usher3
#'
#' Exponentiates the parameter vector before passing it to
#' \code{\link{nll_usher3}}, enabling unconstrained optimization on the
#' log scale.
#'
#' @param th_bar The log-transformed parameter vector
#' @param ... Additional arguments passed to \code{\link{nll_usher3}}
#'
#' @return The negative log-likelihood value
#'
#' @keywords internal
nll_usher3_optim_wrapper <- function(th_bar, ...) {
  th <- exp(th_bar)
  return(nll_usher3(th, ...))
}

#' Use parallel tempering to fit the usher3 model.
#'
#' If the flag use_gompertz is True, then the infant mortality hazard is
#' assumed to be zero and the parameter vector, th0, should only have five
#' terms (k1, k2, and three mortality parameters). If necessary, th0 is subset
#' from 7 to 5 terms if use_gompertz is TRUE. If the scaling matrix,
#' prop_scale_mat, is specified when use_gompertz is True its dimensions must
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
#' @param func_tol The tolerance for the objective function [default: 1e-6]
#' @param miniter The minimum number of iterations [default: 1]
#' @param maxiter The maximum number of iterations [default: 1000]
#' @param report_period The reporting period for verbose output [default: 50]
#' @param use_gompertz Whether to use Gompertz-Makeham mortality [default: FALSE]
#' @param tune Whether to perform Nelder-Mead tuning after tempering [default: TRUE]
#' @param control Control parameters for optim() when tune=TRUE 
#'   [default: list(maxit = 10000, reltol = 1e-12, abstol = 1e-12)]. 
#'   See ?optim for available options.
#' @param ... Additional arguments passed to the objective function
#'   (\code{\link{nll_usher3}}), including x (ages-at-death vector),
#'   ill (illness indicators), x0 (conditional starting age, default 0),
#'   and x_cut (age cutoff for transitions, default Inf)
#'
#' @details
#' Fitting proceeds in two stages, both operating on the log-transformed
#' parameter space (th_bar = log(th)) so that all parameters are
#' unconstrained:
#'
#' 1. **Parallel tempering**: Runs \code{\link{par_temper}} with multiple
#'    temperature chains to explore the parameter space broadly and avoid
#'    local minima. The best parameter vector across all chains is selected.
#'
#' 2. **Nelder-Mead tuning** (when tune=TRUE): Refines the tempering result
#'    using \code{\link[stats]{optim}} with method "Nelder-Mead" and tight
#'    convergence tolerances. This polishes the solution to high precision.
#'
#' If prop_scale_mat is NULL, a default proposal scaling matrix is
#' constructed from the initial parameter vector, with scales ranging from
#' 0.001 to 0.1 across the temperature ladder.
#'
#' @return A list with components:
#'   \describe{
#'     \item{obj_fun}{The objective function used (nll_usher3_optim_wrapper)}
#'     \item{th0}{The initial parameter vector}
#'     \item{optional_inputs}{List of optional arguments used}
#'     \item{temper}{Full output from par_temper}
#'     \item{th_temper}{Best parameter vector from tempering (natural scale)}
#'     \item{optim_result}{Output from optim() (only when tune=TRUE)}
#'     \item{th_bar}{Best log-transformed parameter vector (only when tune=TRUE)}
#'     \item{th}{Best parameter vector on the natural scale (final result)}
#'   }
#'
#' @examples
#' \dontrun{
#' b_siler <- c(0.175, 1.40, 0.00368, 38.1, 0.0917)
#' th <- c(0.02, 1.2, b_siler)
#' samp <- sample_usher3(1000, th, 0.01, 120)
#' result <- temper_and_tune_usher3(th, x = samp$x, ill = samp$ill)
#' result$th  # recovered parameters
#' }
#'
#' @export
temper_and_tune_usher3 <- function(th0 = c(1e-2,
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
                          prop_scale_mat = NULL,
                          func_tol = 1e-6, miniter = 1,
                          maxiter = 1000, report_period = 50,
                          use_gompertz=FALSE, tune=TRUE, 
                          control = list(maxit = 10000, reltol = 1e-12, abstol = 1e-12), ...) {
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
                          func_tol = func_tol, miniter = miniter,
                          maxiter = maxiter, report_period = report_period,
                          tune = tune)

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

  # Tune using Nelder-Mead optimization on the unconstrained parameterization
  if (tune) {
    if (verbose) {
      print('Starting Nelder-Mead tuning...')
    }
    
    optim_result <- optim(par = th_bar_temper, fn = obj_fun, 
                          method = "Nelder-Mead", control = control, ...)
    
    th_bar <- optim_result$par
    th <- exp(th_bar)
    
    if (verbose) {
      print('Best parameter vector after tuning (th):')
      print(th)
      print('Best parameter vector after tuning (th_bar):')
      print(th_bar)
      print('The corresponding best value of the objective function:')
      print(optim_result$value)
      print('Number of function evaluations:')
      print(optim_result$counts)
      print('Convergence code:')
      print(optim_result$convergence)
    }
    
    return(list(obj_fun = obj_fun, th0 = th0, optional_inputs = optional_inputs,
                temper = temper, th_temper = th_temper, 
                optim_result = optim_result, th_bar = th_bar, th = th))
  } else {
    # Return tempered results without tuning
    return(list(obj_fun = obj_fun, th0 = th0, optional_inputs = optional_inputs,
                temper = temper, th_temper = th_temper, th = th_temper))
  }
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