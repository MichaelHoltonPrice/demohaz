% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization.R
\name{par_temper}
\alias{par_temper}
\title{Parallel tempering for function minimization}
\usage{
par_temper(
  theta0,
  neg_log_cost_func,
  samps_per_cyc = 200,
  temp_vect = 10^(rev(seq(-1, 1, by = 0.1))),
  prop_scale = 1,
  num_cyc = 100,
  num_cores = NA,
  verbose = FALSE,
  fn_plot = NULL,
  report_period = 50,
  ...
)
}
\arguments{
\item{theta0}{The starting point for sampling (used to initialize all chains)}

\item{neg_log_cost_func}{The negative log of the cost function}

\item{samps_per_cyc}{Number of within chain Metropolis samples per cycle}

\item{temp_vect}{Vector of temperatures for each chain (ordered highest to
lowest)}

\item{prop_scale}{A speficiation of the standard deviations to use for
proposing new parameter values. prop_scale should be a scalar, a vector of
length(theta), or a matrix with dimensions length(theta) by
length(temp_vect).}

\item{num_cyc}{Number of cycles (a swap is attempted after each cycle).}

\item{num_cores}{Number of cores to use in parallel for loop (default: NA, no
parallelization)}

\item{verbose}{Whether to print out optimization information
(default: FALSE)}

\item{fn_plot}{A function that can be used to plot over the histogram of
the data (xvalues/xcounts). If NULL, no plot is made
(default: FALSE)}

\item{report_period}{How often to update information, in steps. This is used
both for printing out information (if verbose is TRUE) and making a plot
(if fn_plot is not NULL) (default: 50)}

\item{...}{Variables required by neg_log_cost_func}
}
\value{
An object of class \code{par_temper} that consists of (a) chains (the
  sampled chains), (b) swap_mat, a matrix summarizing the results of the
  swap attempts, and (c) inputs, the original inputs to the function.
}
\description{
Minimize an input cost function using parallel tempering. The cost function
is often a negative log-likelihood, and is exponentiated at the accept/reject
steps so that each chain is sampling the likelihood function raised to the
temperature of the chain. The tempering consists of two steps:
}
\details{
(1) Do Metropolis-Hastings sampling of each chain samps_per_cyc (samples per
    cycle) times.
(2) Randomly choose two adjacent temperatures to attempt a swap.

These two steps are performed num_cyc (number of cycles) times. A precise
description of these steps follows.

Let k = 1...K index temperatures in the vector temp_vect. Each cycle of the
algorithm consists of a set of within-in chain Metropolis samples for each
temperature, followed by a single attempt to randomly swap samples across
chains for a single, randomly chosen pair of adjacent temperatures (k and
k+1). M within chain samples are made by calling the function
do_mh_sampling_at_temp. Following these within-chain samples, an
attempt is made to swap samples between one pair of adjacent temperatures,
also using a Metropolis criterion. The pair for which this is attempted is
randomly chosen. Let k be the lower temperature. The weighting accorded the
non-swapped configuration is exp(-eta_k/T_k)*exp(-eta_kp1/T_kp1), where kp1
stands for k+1 (k plus one) and eta_k (etc.) is the cost function evaluated
for the current value of the chain, theta_k. The weighting accorded the
swapped configuration is exp(-eta_kp1/T_k)*exp(-eta_k/T_kp1). The swap is
always accepted if the swapped configuration has a higher weighting;
otherwise, it is accepted with a probability equal to the proportion of these
weightings. To be precise, the Metropolis acceptance ratio is

a = min(1,exp(-(eta_kp1-eta_k)*(1/T_k-1/T_kp1)))

The swap is accepted with probability a, and otherwise rejected. The cycle of
within-chain sampling followed by a single swap attempt is repeated num_cyc
times.
}
