#' @title
#' Solve an minimization problem by doing an initial tempering than fine-tuning
#' the optimization usking the Hooke-Jeeves algorithm.
#'
#' @description
#' Many of the inputs are hard coded....
#'
#' @param obj_fun The objective function to minimize
#' @param th0 The starting point for minimization
#' @param num_cyc Number of cycles for the tempering (default: 100)
#' @param samps_per_cyc Samples percycle for the tempering (default: 20)
#' @param temp_vect Temperature vector for the tempering (default:
#'   10^(rev(seq(-1,1,by=.25)))
#' @param prop_scale_mat Proposal scale matrix for the tempering (default:
#'   ranges evenly from 0.1 for the hottest temperature to 0.001 for the
#'   coldest temperature, for all parameter values at each temperature).
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
                            ...) {
  num_param <- length(th0)

  # For the scale of the proposal distribution, use 0.1 for the highest
  # temperature and 0.001 for the coldest temperature.
  if (is.null(prop_scale_mat)) {
    prop_scale_mat <- t(replicate(num_param,
                                  rev(seq(0.001,.1,len=length(temp_vect)))))
  }
  temper <- enneal::par_temper(th0,
                               obj_fun,
                               temp_vect=temp_vect,
                               prop_scale=prop_scale_mat,
                               num_cyc=num_cyc,
                               samps_per_cyc=samps_per_cyc,
                               ...)

  n <- which.min(unlist(lapply(temper$chains,function(x){x$eta_best})))
  th_temper <- temper$chains[[n]]$theta_best

  # Refine the solution using hjk
  hjk_output <- dfoptim::hjk(th_temper,
                             obj_fun,
                             ...)

  th <- hjk_output$par

  return(list(obj_fun=obj_fun,
              th0=th0,
              temper=temper,
              th_temper=th_temper,
              hjk_output=hjk_output,
              th=th))
}