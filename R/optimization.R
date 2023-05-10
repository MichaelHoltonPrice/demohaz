#' @title
#' Solve a minimization problem by doing an initial tempering then fine-tuning
#' the optimization usking the Hooke-Jeeves algorithm.
#'
#' @description
#' This method does robust function minimization by doing an initial
#' optimization using the parallel tempering algorithm in the enneal package,
#' then fine tuning that optimization using the Hooke-Jeeves algorithm from the
#' dfoptim package. The first input is the objective function to be minimized
#' and the second input is the starting parameter vector. Optionally, control
#' parameters for the tempering can be input (the number of cycles, number of
#' samples per cycle, temperature vector, and proposal scaling matrix).
#' Additional inputs needed by the objective function can also be input.
#'
#' @param obj_fun The objective function to minimize
#' @param th0 The starting point for minimization
#' @param num_cyc Number of cycles for the tempering (default: 100)
#' @param samps_per_cyc Samples per cycle for the tempering (default: 20)
#' @param temp_vect Temperature vector for the tempering (default:
#'   10^(rev(seq(-1,1,by=.25)))
#' @param prop_scale_mat Proposal scale matrix for the tempering (default:
#'   ranges evenly from 0.1 for the hottest temperature to 0.001 for the
#'   coldest temperature, for all parameter values at each temperature).
#' @param verbose Whether to print out optimization information
#'   (default: FALSE)
#' @param lr The learning rate to use for fine tuning using gradient descent
#'   (default: 1e-5)
#' @param ... Additional inputs to the objective function
#'
#' @return The best fit parameter vector
#'
#' @export
temper_and_tune <- function(obj_fun,
                            th0,
                            num_cyc=100,
                            samps_per_cyc=20,
                            temp_vect = 10^(rev(seq(-1,1,by=.25))),
                            prop_scale_mat = NULL,
                            verbose=FALSE,
                            lr=1e-5,
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
  temper <- enneal::par_temper(th0,
                               obj_fun,
                               temp_vect=temp_vect,
                               prop_scale=prop_scale_mat,
                               num_cyc=num_cyc,
                               samps_per_cyc=samps_per_cyc,
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
                             fast_gradnllsiler,
                             lr=lr,
                             100,
                             1e-1,
                             rescale=TRUE,
                             verbose=verbose,
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

    if(verbose) {
      grad <- gr(th, ...)  # Update gradient
      cat(sprintf("Iteration: %d, Objective function: %.3f, Gradient: %s\r",
                  iter, f_new, toString(round(grad, 3))))
      flush.console()
    }
    
    # TODO: implement plotting
    ## Plot the progress if fn_plot is not NULL
    #if (!is.null(fn_plot) && iter %% report_period == 0) {
    #  fn_plot(...)
    #}
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


#gradient_descent <- function(th0,
                             #fn0,
                             #gr0,
                             #rescale=FALSE,
                             #lr=1e-5,
                             #func_tol=1e-6,
                             #grad_tol=1e-2,
                             #miniter=1,
                             #maxiter=10000,
                             #verbose=FALSE,
                             #fn_plot=NULL,
                             #report_period=50,
                             #...) {
  ## If is TRUE, then rescale each variable by th0 during the optimization,
  ## including for the gradient call. This requires scaling th and defining
  ## wrappers for fn0 and gr0 (hence why they have 0 added).
  #if (rescale) {
    #th <- rep(1, length(th0))
    #fn <- function(th,...) {
      #return(fn0(th*th0, ...))
    #}
    #gr <- function(th,...) {
      #grad <- gr0(th*th0,...)
      #return(grad*th0)
    #}
  #} else {
    #th <- th0
    #fn <- fn0(th,...)
    #gr <- gr0(th,...)
  #}

  #f <- fn(th, ...)
  #grad <- gr(th, ...)
  #if(verbose) {
    #cat(sprintf("Iteration: %d, Objective function: %.3f, Gradient: %s\r",
                #0, f, toString(round(grad, 3))))
    #flush.console()  # Ensure that the output is updated immediately
  #}

  #for (iter in 1:maxiter) {
    #grad <- gr(th, ...)

    ## Check the gradient norm
    #if (all(abs(grad) < grad_tol)) break

    #th_new <- th - lr * grad
    #f_new <- fn(th_new, ...)

    ## Check the change in the objective function value
    ##if (abs(fx_new - fx) < tol) break

    #th <- th_new
    #f <- f_new
    #if(verbose) {
      #cat(sprintf("Iteration: %d, Objective function: %.3f, Gradient: %s\r",
                  #iter, f, toString(round(grad, 3))))
      #flush.console()  # Ensure that the output is updated immediately
    #}
  #}
 
  #inputs <- list(th0=th0,
                 #fn0=fn0,
                 #gr0=gr0,
                 #rescale=rescale,
                 #lr=lr,
                 #func_tol=func_tol,
                 #grad_tol=grad_tol,
                 #miniter=miniter,
                 #maxiter=maxiter,
                 #verbose=verbose,
                 #fn_plot=fn_plot,
                 #report_period=report_period)

  #return(list(par = th, value = f, feval = iter, inputs=inputs))
#}