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

 # TODO: support plotting in gradient_descent
 descent <- gradient_descent(th_temper,
                             obj_fun,
                             fast_gradnllsiler,
                             lr=lr,
                             100,
                             1e-1,
                             rescale=TRUE,
                             verbose=verbose,
                             ...)

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