#' @title
#' Solve a minimization problem by doing an initial tempering then fine-tuning
#' the optimization usking the Hooke-Jeeves algorithm.
#'
#' @description
#' This method does robust function minimization by doing an initial
#' optimization using the parallel tempering algorithm in the this package,
#' then fine tuning that optimization using the Hooke-Jeeves algorithm from the
#' dfoptim package. The first input is the objective function to be minimized
#' and the second input is the starting parameter vector. Optionally, control
#' parameters for the tempering can be input (the number of cycles, number of
#' samples per cycle, temperature vector, and proposal scaling matrix).
#' Additional inputs needed by the objective function can also be input.
#'
#' @param obj_fun The objective function to minimize
#' @param grad_fun The gradient function for the objective
#' @param th0 The starting point for minimization
#' @param verbose Whether to print out optimization information
#'   (default: FALSE)
#' @param fn_plot A function that can be used to plot progress. If NULL, no
#'   plot is made
#' @param num_cyc Number of cycles for the tempering (default: 100)
#' @param samps_per_cyc Samples per cycle for the tempering (default: 20)
#' @param temp_vect Temperature vector for the tempering (default:
#'   10^(rev(seq(-1,1,by=.25)))
#' @param prop_scale_mat Proposal scale matrix for the tempering (default:
#'   ranges evenly from 0.1 for the hottest temperature to 0.001 for the
#'   coldest temperature, for all parameter values at each temperature).
#' @param rescale Whether to rescale the parameter th by dividing by th0
#'   (default: FALSE) [for gradient_descent]
#' @param lr The learning rate to use
#'   (default: 1e-5) [for gradient_descent]
#' @param func_tol A tolerance for changes in the objective function value. If
#'   the change in value is less than the tolerance, optimization can halt,
#'   though only if the grad_tol is condition is met (default: 1e-6).
#'   [for gradient_descent]
#' @param grad_tol A tolerance for the absolute value of the gradient. All
#'   gradients must be lower than this tolerance value. If rescale is TRUE,
#'   the condition is applied to the rescaled gradient (default: 1e-2). Both
#'   the func_tol and grad_tol conditions must be satisfied to halt.
#'   [for gradient_descent]
#' @param miniter The minimum number of iterations to use (default: 1)
#'   [for gradient_descent]
#' @param maxiter The maximum number of iterations to use (default: 1000)
#'   [for gradient_descent]
#' @param report_period How often to update information, in steps. This is used
#'   by both the tempering and gradient descent, but will not apply unless
#'   verbose is TRUE or fn_plot is not NULL. (default: 50)
#
#' @param ... Additional inputs to the objective function
#'
#' @return The best fit parameter vector
#'
#' @export
temper_and_tune <- function(obj_fun,
                            grad_fun,
                            th0,
                            verbose=FALSE,
                            fn_plot=NULL,
                            num_cyc=100,
                            samps_per_cyc=20,
                            temp_vect = 10^(rev(seq(-1,1,by=.25))),
                            prop_scale_mat = NULL,
                            lr=1e-5,
                            func_tol=1e-6,
                            grad_tol=1e-2,
                            miniter=1,
                            maxiter=1000,
                            report_period=50,
                            ...) {
  num_param <- length(th0)

  # For the scale of the proposal distribution, use 0.1 for the highest
  # temperature and 0.001 for the coldest temperature. This is multiplied by
  # each parameter value in b0 since they have different scales.
  if (is.null(prop_scale_mat)) {
    just_scale <- t(replicate(num_param,
                                  rev(seq(0.001,.1,len=length(temp_vect)))))
    just_param <- replicate(length(temp_vect),th0)
    prop_scale_mat <- just_scale * just_param
  }
  if(verbose) {
    print('Starting parameter vector th0 = ')
    print(th0)
    print('Starting objective function value eta = ')
    print(obj_fun(th0, ...))
  }
  # TODO: add verbosity to par_temper
  # TODO: support plotting in par_temper
  temper <- par_temper(th0,
                       obj_fun,
                       temp_vect=temp_vect,
                       prop_scale=prop_scale_mat,
                       num_cyc=num_cyc,
                       samps_per_cyc=samps_per_cyc,
                       fn_plot=fn_plot,
                       verbose=verbose,
                       ...)

  n <- which.min(unlist(lapply(temper$chains,function(x){x$eta_best})))
  th_temper <- temper$chains[[n]]$theta_best

  if(verbose) {
    print('Best parameter vector from initial tempering:')
    print(th_temper)
    print('The corresponding best value of the objective function:')
    print(temper$chains[[n]]$eta_best)
  }

  descent <- gradient_descent(th_temper,
                              obj_fun,
                              grad_fun,
                              rescale=TRUE,
                              lr=lr,
                              func_tol=func_tol,
                              grad_tol=grad_tol,
                              miniter=miniter,
                              maxiter=maxiter,
                              verbose=verbose,
                              fn_plot=fn_plot,
                              report_period=report_period,
                              ...)

  # TODO: support plotting in gradient_descent
  # Undo the scaling and, if necessary, print the result
  th <- descent$par*th_temper
  if(verbose) {
    print('Gradient descent best')
    print(th)
    print(descent$value)
  }

  return(list(obj_fun=obj_fun,
              th0=th0,
              temper=temper,
              th_temper=th_temper,
              descent=descent,
              th=th))
}

#' @title
#' Vanilla gradient descent for demohaz
#'
#' @description
#' Do vanilla gradient descent. Optionally, if rescale is TRUE, the gradient
#' is calculated after dividing by the starting parameter vector par. There are
#' two reasons why demohaz implements its own optimziation algorithms. First,
#' it reduces dependenceies on third-party packages. Second, it allows
#' precise control over how to provide diagnostic information (e.g., printing
#' optimization information if verbose is TRUE and plotting the fit if
#' show_plot is TRUE).
#'
#' @param th0 The starting point for minimization
#' @param fn0 The objective function to minimize
#' @param gr0 The gradient of the objective function to minimize
#' @param rescale Whether to rescale the parameter th by dividing by th0
#'   (default: FALSE)
#' @param lr The learning rate to use for fine tuning using gradient descent
#'   (default: 1e-5)
#' @param func_tol A tolerance for changes in the objective function value. If
#'   the change in value is less than the tolerance, optimization can halt,
#'   though only if the grad_tol is condition is met (default: 1e-6).
#' @param grad_tol A tolerance for the absolute value of the gradient. All
#'   gradients must be lower than this tolerance value. If rescale is TRUE,
#'   the condition is applied to the rescaled gradient (default: 1e-2). Both
#'   the func_tol and grad_tol conditions must be satisfied to halt.
#' @param miniter The minimum number of iterations to use (default: 1)
#' @param maxiter The maximum number of iterations to use (default: 10000)
#' @param verbose Whether to print out optimization information
#'   (default: FALSE)
#' @param fn_plot A function that can be used to plot over the histogram of
#'   the data (xvalues/xcounts). If NULL, no plot is made
#'   (default: FALSE)
#' @param report_period How often to update information, in steps. This is used
#'   both for printing out information (if verbose is TRUE) and making a plot
#'   (if fn_plot is not NULL) (default: 50)
#' @param ... Additional inputs to fn0, gr0, and fn_plot
#'
#' @return A list containing par (the best-fit parameter value), value (the
#'   the corresponding objective function value), feval (the number of
#'   iterations/steps), and inputs (a list with the input variable values).
#'
#' @export
gradient_descent <- function(th0,
                             fn0,
                             gr0,
                             rescale=FALSE,
                             lr=1e-5,
                             func_tol=1e-6,
                             grad_tol=1e-2,
                             miniter=1,
                             maxiter=10000,
                             verbose=FALSE,
                             fn_plot=NULL,
                             report_period=50,
                             ...) {
  # TODO: th0 could have zeros when rescale is TRUE
  if (miniter < 1) {
    stop('The minimum number of iterations, miniter, must be at least 1')
  }
  if (rescale) {
    th <- rep(1, length(th0))
    fn <- function(th,...) {
      return(fn0(th*th0, ...))
    }
    gr <- function(th,...) {
      grad <- gr0(th*th0,...)
      return(grad*th0)
    }
  } else {
    th <- th0
    fn <- function(th,...) {
      return(fn0(th, ...))
    }
    gr <- function(th,...) {
      return(gr0(th, ...))
    }
  }

  f <- fn(th, ...)
  grad <- gr(th, ...)
  if(verbose) {
    cat(sprintf("Iteration: %d, Objective function: %.3f, Gradient: %s\r",
                0, f, toString(round(grad, 3))))
    flush.console()
  }

  for (iter in 1:maxiter) {
    grad <- gr(th, ...)

    th_new <- th - lr * grad
    f_new <- fn(th_new, ...)


    th <- th_new

    if (iter %% report_period == 0) {
      if(verbose) {
        # TODO: gradient is calculated twice if report_period=1
        #       could store the calculation for when it is made here
        grad <- gr(th, ...)  # Update gradient
        print(paste0('Iteration: ', str(iter)))
        print(paste0('Obj. Func: ', str(f_new)))
        print(paste0('Gradient : ', toString(round(grad, 3))))
      }
      if (!is.null(fn_plot)) {
        if (rescale) {
          th_plot <- th * th0
        } else {
          th_plot <- th
        }
        fn_plot(th_plot, ...)
      }
    }
    
    if (iter > miniter) {
      func_condition <- abs(f_new - f) < func_tol
      grad_condition <- all(abs(grad) < grad_tol)
      if (func_condition & grad_condition) break
    }
    f <- f_new
  }

  inputs <- list(th0=th0,
                 fn0=fn0,
                 gr0=gr0,
                 rescale=rescale,
                 lr=lr,
                 func_tol=func_tol,
                 grad_tol=grad_tol,
                 miniter=miniter,
                 maxiter=maxiter,
                 verbose=verbose,
                 fn_plot=fn_plot,
                 report_period=report_period)

  return(list(par = th, value = f, feval = iter, inputs=inputs))
}

#' Do Metropolis-Hastings sampling at a given temperature
#'
#' Do Metropolis-Hastings sampling given an input negative log of some target
#' distribution (probably a negative log-likelihood) at a given temperature. The
#' chain can be initialized either with a starting parameter vector or with a
#' chain created by a previous call to do_mh_sampling_at_temp. For continuing
#' chains, some (but not all) of the inputs from the previously used values can
#' be input.
#'
#' @param init The initializing object (either a parameter vector or a chain
#'   created by a previous call to do_mh_sampling_at_temp).
#' @param neg_log_cost_func The function that calculates the negative log of the
#'   function being sampled.
#' @param num_samp The number of samples to make.
#' @param temp The temperature
#' @param prop_scale The standard deviation(s) to use for the proposal
#'   distribution (either a vector that is the same length as the parameter
#'   vector, theta, or a scalar that is used for all parameters).
#' @param save_theta Whether or not to save and return the accepted parameter
#'   values (assumed False if not input). Cannot be over-ridden for continuing
#'   chains.
#' @param ... Variables required by neg_log_cost_func.
#' @return A list-like object of class "mh_chain" with sampling information
#' @export
do_mh_sampling_at_temp <- function(init,
                           neg_log_cost_func=NA,
                           num_samp=NA,
                           temp=NA,
                           prop_scale=NA,
                           save_theta=NA,
                           ...) {

  # TODO: explore saving the ... parameters for use in continuing chains.

  # Is this a new chain, or are we continuining an existing chain?
  if ("mh_chain" %in% class(init)) {
    new_chain <- FALSE
  } else if (is.numeric(init)) {
    new_chain <- TRUE
  } else {
    stop(paste0("init should be a starting parameter vector or continued ",
                "chain of class mh_mcmc"))
  }

  if(new_chain) {
    # Do error checking on the inputs
    if (is.na(num_samp)) {
      stop("num_samp must be input for new chains")
    }

    if(class(neg_log_cost_func) != "function") {
      stop("neg_log_cost_func must be input for new chains")
    }

    if(is.na(temp)) {
      stop("temp must be input for new chains")
    }

    # prop_scale must always be input
    if (any(is.na(prop_scale))) {
      stop("prop_scale must be input for new chains")
    }

    if (length(prop_scale) != 1 && length(prop_scale) != length(init)) {
      stop("prop_scale must be a scalar or have the same length as theta")
    }

    if (length(prop_scale) == 1) {
      if (length(init) != 1) {
        prop_scale <- rep(prop_scale,length(init))
      }
    }

    # If not input, save_theta is set to False. The default in the function call
    # is NA so we can ensure that the user has not over-ridden the value for
    # continuing chains.
    if(is.na(save_theta)) {
      save_theta <- False
    }

    # For a new chain, calculate the starting negative log-likelihood and set
    # theta0 equal to the input parameter vector.
    theta0 <- init
    eta0 <- neg_log_cost_func(theta0,...)
    if (!is.finite(eta0)) {
      stop(paste0("The negative log cost function is not finite for the ",
                  "input initialization vector theta0"))
    }

    # Set eta_best and theta_best
    theta_best <- theta0
    eta_best   <- eta0

    # Set theta and eta, the currently accepted variables, to theta0 and eta0
    eta <- eta0
    theta <- theta0
  } else {
    # Currently, resetting the sampling function is not supported inside this
    # function. It can be rest by directly editing the function inside the
    # mh_chain object outside the function. Allowing it to be reset seems likely
    # to cause mis-uses more often than it is needed.
    if (class(neg_log_cost_func) != "function") {
       neg_log_cost_func <- init$neg_log_cost_func
    } else {
      stop("neg_log_cost_func should not be input for continuing chains.")
    }

    # For an existing chain, get eta0, theta0, eta_best, and theta_best from
    # init
    eta0       <- init$eta0
    theta0     <- init$theta0
    eta_best   <- init$eta_best
    theta_best <- init$theta_best

    # Set theta and eta, the currently accepted variables, using the last values
    # from the input chain.
    #num_prev_samples <- init$total_samples
    eta   <- init$eta
    theta <- init$theta

    # If necessary, set variables using the values stored in init for continuing
    # chains or, for new chains, throw errors if necessary variables are not
    # input.
    if (is.na(num_samp)) {
      num_samp <- init$num_samp
    }

    # The temperature can always be over-ridden
    if(is.na(temp)) {
      temp <- init$temp
    }

    # The proposal scale can always be over-ridden
    if(any(is.na(prop_scale))) {
      prop_scale <- init$prop_scale
    }

    # The number of samples can always be over-ridden
    if(is.na(num_samp)) {
      num_samp <- init$num_samp_vect[length(init$num_samp_vect)]
    }

    # Currently, at least, the value of save_theta cannot be overwritten.
    if(is.na(save_theta)) {
      save_theta <- init$save_theta
    } else {
      stop("save_theta should not be input for continuing chains.")
    }
  }

  # accept_vect is a vector of length num_samp that records whether each sample
  # was accepted or rejected
  accept_vect <- rep(NA,num_samp)

  # eta_vect is a vector of length num_samp that records the value of the
  # negative log-likelihood (or, if direct is TRUE, the sampled function).
  eta_vect <- rep(NA,num_samp)

  # If necessary, initialize the matrix of save parameter values, theta_mat,
  # which has dimensions Number of Parameters x Number of Samples.
  if(save_theta) {
    theta_mat <- matrix(NA,length(theta0),num_samp)
  }

  for(n in 1:num_samp) {
    # Create a new proposal parameter
    theta_prop <- theta + rnorm(length(theta))*prop_scale
    eta_prop <- neg_log_cost_func(theta_prop,...)

    if(!is.finite(eta_prop)) {
      accept <- F
    } else {
      alpha <- min(1,exp(-(eta_prop-eta)/temp))
      accept <- runif(1) < alpha
    }

    accept_vect[n] <- accept
    if(!accept) {
      theta_prop <- theta
      eta_prop   <- eta
    } else {
      if(eta_prop < eta_best) {
        theta_best <- theta_prop
        eta_best   <- eta_prop
      }
    }
    eta_vect[n] <- eta_prop

    # Get ready for next sample
    theta <- theta_prop
    eta   <- eta_prop

    if(save_theta) {
      theta_mat[,n] <- theta
    }
  }

  # If this is a new chain, create and return a new list-like object. If this is
  # a continuing chain, update the input chain and return it.
  if(new_chain) {
    output <- list(theta0=theta0,
                   eta0=eta0,
                   temp=temp,
                   prop_scale=prop_scale,
                   eta_best=eta_best,
                   theta_best=theta_best,
                   accept_vect=accept_vect,
                   eta_vect=eta_vect,
                   num_samp=num_samp,
                   eta=eta,
                   theta=theta,
                   num_samp_vect=num_samp,
                   neg_log_cost_func=neg_log_cost_func,
                   temp=temp,
                   save_theta=save_theta)
    if(save_theta) {
      output$theta_mat <- theta_mat
    }
    class(output) <- "mh_chain"
    return(output)
  } else {
    # Not a new chain

    # Update modified values
    init$eta <- eta
    init$theta <- theta
    init$eta_best <- eta_best
    init$theta_best <- theta_best
    init$accept_vect <- c(init$accept_vect,accept_vect)
    init$eta_vect <- c(init$eta_vect,eta_vect)
    init$num_samp_vect <- c(init$num_samp_vect,num_samp)
    return(init)
  }

  stop("This point should never be reached")
}

#' Parallel tempering for function minimization
#'
#' Minimize an input cost function using parallel tempering. The cost function
#' is often a negative log-likelihood, and is exponentiated at the accept/reject
#' steps so that each chain is sampling the likelihood function raised to the
#' temperature of the chain. The tempering consists of two steps:
#'
#' (1) Do Metropolis-Hastings sampling of each chain samps_per_cyc (samples per
#'     cycle) times.
#' (2) Randomly choose two adjacent temperatures to attempt a swap.
#'
#' These two steps are performed num_cyc (number of cycles) times. A precise
#' description of these steps follows.
#'
#' Let k = 1...K index temperatures in the vector temp_vect. Each cycle of the
#' algorithm consists of a set of within-in chain Metropolis samples for each
#' temperature, followed by a single attempt to randomly swap samples across
#' chains for a single, randomly chosen pair of adjacent temperatures (k and
#' k+1). M within chain samples are made by calling the function
#' do_mh_sampling_at_temp. Following these within-chain samples, an
#' attempt is made to swap samples between one pair of adjacent temperatures,
#' also using a Metropolis criterion. The pair for which this is attempted is
#' randomly chosen. Let k be the lower temperature. The weighting accorded the
#' non-swapped configuration is exp(-eta_k/T_k)*exp(-eta_kp1/T_kp1), where kp1
#' stands for k+1 (k plus one) and eta_k (etc.) is the cost function evaluated
#' for the current value of the chain, theta_k. The weighting accorded the
#' swapped configuration is exp(-eta_kp1/T_k)*exp(-eta_k/T_kp1). The swap is
#' always accepted if the swapped configuration has a higher weighting;
#' otherwise, it is accepted with a probability equal to the proportion of these
#' weightings. To be precise, the Metropolis acceptance ratio is
#'
#' a = min(1,exp(-(eta_kp1-eta_k)*(1/T_k-1/T_kp1)))
#'
#' The swap is accepted with probability a, and otherwise rejected. The cycle of
#' within-chain sampling followed by a single swap attempt is repeated num_cyc
#' times. The minimum value of this procedure can be extracted from the output
#' by calling the function par_temper_best.
#' TODO: implement par_temper_best
#'
#' @param theta0 The starting point for sampling (used to initialize all chains)
#' @param neg_log_cost_func The negative log of the cost function
#' @param samps_per_cyc Number of within chain Metropolis samples per cycle
#' @param temp_vect Vector of temperatures for each chain (ordered highest to
#'   lowest)
#' @param prop_scale A speficiation of the standard deviations to use for
#'   proposing new parameter values. prop_scale should be a scalar, a vector of
#'   length(theta), or a matrix with dimensions length(theta) by
#'   length(temp_vect).
#' @param num_cyc Number of cycles (a swap is attempted after each cycle).
#' @param num_cores Number of cores to use in parallel for loop (default: NA, no
#'   parallelization)
#' @param verbose Whether to print out optimization information
#'   (default: FALSE)
#' @param fn_plot A function that can be used to plot over the histogram of
#'   the data (xvalues/xcounts). If NULL, no plot is made
#'   (default: FALSE)
#' @param report_period How often to update information, in steps. This is used
#'   both for printing out information (if verbose is TRUE) and making a plot
#'   (if fn_plot is not NULL) (default: 50)

#' @param ... Variables required by neg_log_cost_func
#'
#' @return An object of class \code{par_temper} that consists of (a) chains (the
#'   sampled chains), (b) swap_mat, a matrix summarizing the results of the
#'   swap attempts, and (c) inputs, the original inputs to the function.

#' @import doParallel
#' @import foreach
#' @export
par_temper <- function(theta0,
                       neg_log_cost_func,
                       samps_per_cyc=200,
                       temp_vect=10^(rev(seq(-1,1,by=.1))),
                       prop_scale=1,
                       num_cyc=100,
                       num_cores=NA,
                       verbose=FALSE,
                       fn_plot=NULL,
                       report_period=50,
                       ...) {

  inputs <- list(theta0=theta0,
                 neg_log_cost_func=neg_log_cost_func,
                 samps_per_cyc=samps_per_cyc,
                 temp_vect=temp_vect,
                 prop_scale=prop_scale,
                 num_cyc=num_cyc)

  chains = list()
  # prop_scale should be a scalar, a vector of length theta, or a matrix with
  # dimensions length theta x length temp_vect.
  if(is.matrix(prop_scale)) {
    if(!all(dim(prop_scale) == c(length(theta0),length(temp_vect)))) {
      stop(paste0("If prop_scale is a matrix, it should have dimensions ",
                  "length(theta0) by length(temp_vect)"))
    }
  } else {
     if(length(prop_scale) == 1) {
      # A scalar
      prop_scale <- matrix(prop_scale,length(theta0),length(temp_vect))
    } else if (length(prop_scale) == length(theta0)) {
      prop_scale <- replicate(length(temp_vect),prop_scale)
    } else {
      stop(paste0("If prop_scale is a vector, it should be length 1 or the ",
                  "same length as theta0"))
    }
  }

  swap_mat <- matrix(NA,3,num_cyc)

  # If multiple cores should be used, register them.
  if (!is.na(num_cores)) {
    doParallel::registerDoParallel(num_cores)
  }

  if(verbose) {
    eta0 <- neg_log_cost_func(theta0,...)
    print(paste0('Cycle: ', 0))
    print(paste0('Obj. Func: ', eta0))

    if (!is.null(fn_plot)) {
      fn_plot(theta0, ...)
    }
  }


  # Iterate over number of cycles
  for (cc in 1:num_cyc) {
    if (cc == 1) {
      # Start new chains
      if (is.na(num_cores)) {
        # Start new chains using a conventional for loop
        for (k in 1:length(temp_vect)) {
          chains[[k]] <-
            do_mh_sampling_at_temp(theta0,
                                   num_samp=samps_per_cyc,
                                   neg_log_cost_func=neg_log_cost_func,
                                   temp=temp_vect[k],
                                   prop_scale=prop_scale[,k],
                                   save_theta=T,
                                   ...)
        }
      } else {
        # Start new chains using a parallel for loop
        chains <-
          foreach(k=1:length(temp_vect)) %dopar% {
            # Unclear whether, and when, demohaz:: is needed here
            demohaz::do_mh_sampling_at_temp(theta0,
                                            num_samp=samps_per_cyc,
                                            neg_log_cost_func=neg_log_cost_func,
                                            temp=temp_vect[k],
                                            prop_scale=prop_scale[,k],
                                            save_theta=T,
                                            ...)
          }
      }
    } else {
      # Extend chains
      if (is.na(num_cores)) {
        # Extend chains using a conventional for loop
        for (k in 1:length(temp_vect)) {
          chains[[k]] <- do_mh_sampling_at_temp(chains[[k]],...)
        }
      } else {
        # Extend chains using a parallel for loop
        chains <-
          foreach(k=1:length(temp_vect)) %dopar% {
            do_mh_sampling_at_temp(chains[[k]],...)
          }
      }
    }

    # Randomly choose an adjacent pair of temperatures to attempt a swap
    k <- sample(1:(length(temp_vect)-1),1)

    # Extract the current negative log cost from the two adjacent chains
    eta_k     <- chains[[k  ]]$eta
    eta_kp1   <- chains[[k+1]]$eta

    # Extract the temperatures of the two adjacent chains
    T_k     <- temp_vect[k]
    T_kp1   <- temp_vect[k+1]

    # Calculate the swap probability
    a_swap <- exp(-(eta_kp1-eta_k)*(1/T_k-1/T_kp1))
    if(a_swap >= 1) {
      accept <- T
    } else {
      accept <- runif(1) <= a_swap
    }

    if(accept) {
      theta_k     <- chains[[k  ]]$theta
      theta_kp1   <- chains[[k+1]]$theta

      # Swap eta and theta (there is no need to update eta_best and theta_best
      # here).
      chains[[k  ]]$eta <- eta_kp1
      chains[[k+1]]$eta <- eta_k
      chains[[k  ]]$theta <- theta_kp1
      chains[[k+1]]$theta <- theta_k
    }

    swap_mat[1,cc] <- samps_per_cyc*cc
    swap_mat[2,cc] <- k
    if(accept) {
      swap_mat[3,cc] <- 1
    } else {
      swap_mat[3,cc] <- 0
    }


    if (cc %% report_period == 0) {
      nbest <- which.min(unlist(lapply(chains,function(x){x$eta_best})))
      th_best <- chains[[nbest]]$theta_best
      if(verbose) {
        eta_best <- neg_log_cost_func(th_best,...)
        print(paste0('Cycle: ', cc))
        print(paste0('Obj. Func: ', eta_best))
      }

      if (!is.null(fn_plot)) {
        fn_plot(th_best, ...)
      }
    }
  }

  output <- list(chains=chains,swap_mat=swap_mat,inputs=inputs)
  class(output) <- "par_temper"
  return(output)
}