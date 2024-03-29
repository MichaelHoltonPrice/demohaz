% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/siler.R
\name{plot_siler_fit_progress}
\alias{plot_siler_fit_progress}
\title{A utility function to plot the progress of a Siler fit}
\usage{
plot_siler_fit_progress(b, x0, xvalues, xcounts)
}
\arguments{
\item{x0}{The conditional starting age [default: 0].}

\item{xvalues}{The x-values}

\item{xcounts}{The number of observations for each element in xvalues}

\item{b0}{The Siler parameter vector}
}
\description{
Plot the progress of a Siler fit. The input variables the parameter vector,
b, as well as xvalues and xcounts (not x). A histogram is made of the actual
data and the density from dsiler is plotted on top of it.
}
