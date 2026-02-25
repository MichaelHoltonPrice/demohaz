library(testthat)
library(demohaz)

# Baseline parameter vector (LEH example):
#   th0 = [k1, k2, b_siler]  where b_siler has 5 Siler hazard parameters.
#   k1  = 0.02  : well-to-ill transition rate
#   k2  = 1.2   : mortality multiplier for ill individuals
b0 <- c(.175, 1.40, .368 * .01,
        log(.917 * .1 / (.075 * .001)) / (.917 * .1),
        .917 * .1)
th0 <- c(2e-2, 1.2, b0)

# ---------------------------------------------------------------------------
# Forward simulation of the Usher3 model using only the raw hazard rates.
# This is an independent implementation that shares no code with the
# rejection sampler, making it a strong check.
#
# The simulation discretizes time into small steps of size dt. Within each
# step, every living individual faces independent competing risks whose
# transition probabilities are hazard * dt (valid when hazard * dt << 1).
#
# Args:
#   N     - number of individuals to simulate
#   th    - parameter vector [k1, k2, b_siler] (same format as sample_usher3)
#   dt    - time step size (smaller = more accurate; 0.01 keeps all
#           transition probabilities below ~0.05)
#   xmax  - maximum age; anyone still alive at xmax is recorded there
#   x0    - starting age (all individuals begin in the well state)
#   x_cut - age after which well-to-ill transitions stop (Inf = never stop)
#
# Returns:
#   list with $x (numeric: age at death) and $ill (logical: ill at death)
# ---------------------------------------------------------------------------
simulate_usher3_forward <- function(N, th, dt, xmax, x0 = 0, x_cut = Inf) {
  k1 <- th[1]
  k2 <- th[2]
  b_siler <- th[3:7]

  # Integer-coded state for each individual: 1 = well, 2 = ill, 3 = dead.
  # Using integer literals (1L, 2L, 3L) to match the integer storage type.
  state <- rep(1L, N)
  age <- rep(x0, N)

  # Default: if someone survives to xmax they are recorded as dying there.
  age_at_death <- rep(xmax, N)
  ill_at_death <- rep(FALSE, N)

  nsteps <- ceiling((xmax - x0) / dt)

  for (step in seq_len(nsteps)) {
    # Identify living individuals (state 1 or 2, i.e. < 3)
    alive <- state < 3L
    if (!any(alive)) break

    alive_idx <- which(alive)
    x <- age[alive_idx]        # current ages of living individuals
    s <- state[alive_idx]      # current states of living individuals
    n_alive <- length(alive_idx)

    # Baseline Siler mortality hazard at each individual's current age
    h <- hsiler(x, b_siler)

    # One uniform draw per individual to determine this step's outcome
    u <- runif(n_alive)

    # Boolean masks within the alive subset
    well_mask <- s == 1L
    ill_mask  <- s == 2L

    # --- Compute transition probabilities (discrete-time approximation) ---
    # For well individuals, two competing risks:
    #   die:  p = hsiler(x) * dt          (well -> dead)
    #   sick: p = k1 * dt  if x < x_cut   (well -> ill)
    # For ill individuals, one risk:
    #   die:  p = k2 * hsiler(x) * dt      (ill -> dead)
    p_die <- numeric(n_alive)
    p_ill_trans <- numeric(n_alive)

    p_die[well_mask] <- h[well_mask] * dt
    p_ill_trans[well_mask] <- ifelse(x[well_mask] < x_cut, k1 * dt, 0)
    p_die[ill_mask] <- k2 * h[ill_mask] * dt

    # --- Resolve outcomes using the single uniform draw ---
    # The [0, 1) interval is partitioned as:
    #   [0, p_die)                  -> individual dies
    #   [p_die, p_die+p_ill_trans)  -> individual becomes ill (well only)
    #   [p_die+p_ill_trans, 1)      -> nothing happens
    dies <- u < p_die
    becomes_ill <- !dies & well_mask & (u < p_die + p_ill_trans)

    # --- Record deaths ---
    dead_idx <- alive_idx[dies]
    state[dead_idx] <- 3L
    n_dead <- length(dead_idx)
    if (n_dead > 0) {
      # Place the death time within the step [age, age + dt] by sampling
      # from the truncated exponential distribution.  Given constant hazard
      # h_death over [0, dt] and conditional on death occurring in that
      # interval, the CDF is:
      #   F(t) = (1 - exp(-h_death * t)) / (1 - exp(-h_death * dt))
      # Inverting via the probability integral transform gives:
      #   t = -log(1 - U * (1 - exp(-h_death * dt))) / h_death
      # This avoids ties (which would invalidate the KS test) and is more
      # principled than a uniform jitter.
      h_death <- numeric(n_alive)
      h_death[well_mask] <- h[well_mask]
      h_death[ill_mask] <- k2 * h[ill_mask]
      h_dead <- h_death[dies]
      u_sub <- runif(n_dead)
      offset <- -log(1 - u_sub * (1 - exp(-h_dead * dt))) / h_dead
      age_at_death[dead_idx] <- age[dead_idx] + offset
    }
    ill_at_death[dead_idx] <- (s[dies] == 2L)

    # --- Record well -> ill transitions ---
    trans_idx <- alive_idx[becomes_ill]
    state[trans_idx] <- 2L

    # --- Advance age for all survivors ---
    survivors <- state < 3L
    age[survivors] <- age[survivors] + dt
  }

  # Anyone still alive at the end is recorded as dying at xmax
  still_alive <- state < 3L
  ill_at_death[still_alive] <- (state[still_alive] == 2L)

  list(x = age_at_death, ill = ill_at_death)
}

# ---------------------------------------------------------------------------
# Test 1: x_cut = Inf (standard model, well-to-ill transitions at all ages)
# ---------------------------------------------------------------------------
test_that("forward simulation matches rejection sampler (x_cut = Inf)", {
  set.seed(12345)
  N <- 10000       # sample size for each method
  dt <- 0.01       # forward simulation time step
  dx <- 0.001      # rejection sampler grid spacing
  xmax <- 120      # maximum age

  sim <- simulate_usher3_forward(N, th0, dt, xmax, x_cut = Inf)
  imp <- sample_usher3(N, th0, dx, xmax, x_cut = Inf)

  # Split each sample by illness status at death
  sim_well <- sim$x[!sim$ill]
  sim_ill  <- sim$x[sim$ill]
  imp_well <- imp$x[!imp$ill]
  imp_ill  <- imp$x[imp$ill]

  # Two-sample KS tests: if both methods sample from the same distribution,
  # p-values should be non-significant (> 0.05) with high probability.
  expect_gt(ks.test(sim_well, imp_well)$p.value, 0.05)
  expect_gt(ks.test(sim_ill, imp_ill)$p.value, 0.05)

  # The fraction of individuals who are ill at death should agree closely
  # between the two methods.
  prop_sim <- mean(sim$ill)
  prop_imp <- mean(imp$ill)
  expect_lt(abs(prop_sim - prop_imp), 0.03)
})

# ---------------------------------------------------------------------------
# Test 2: x_cut = 6 (well-to-ill transitions only before age 6, e.g. LEH)
# ---------------------------------------------------------------------------
test_that("forward simulation matches rejection sampler (x_cut = 6)", {
  set.seed(67890)
  N <- 10000
  dt <- 0.01
  dx <- 0.001
  xmax <- 120
  x_cut <- 6

  sim <- simulate_usher3_forward(N, th0, dt, xmax, x_cut = x_cut)
  imp <- sample_usher3(N, th0, dx, xmax, x_cut = x_cut)

  # Split each sample by illness status at death
  sim_well <- sim$x[!sim$ill]
  sim_ill  <- sim$x[sim$ill]
  imp_well <- imp$x[!imp$ill]
  imp_ill  <- imp$x[imp$ill]

  # Two-sample KS tests
  expect_gt(ks.test(sim_well, imp_well)$p.value, 0.05)
  expect_gt(ks.test(sim_ill, imp_ill)$p.value, 0.05)

  # Proportion ill at death
  prop_sim <- mean(sim$ill)
  prop_imp <- mean(imp$ill)
  expect_lt(abs(prop_sim - prop_imp), 0.03)
})
