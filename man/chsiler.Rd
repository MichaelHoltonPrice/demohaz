% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/siler.R
\name{chsiler}
\alias{chsiler}
\title{Calculate the Siler cumulative hazard}
\usage{
chsiler(x, b, x0 = 0)
}
\arguments{
\item{x}{The vector of ages}

\item{b}{The parameter vector}

\item{x0}{The conditional starting age [default: 0].}

\item{qvect}{The vector of quantiles}

\item{N}{The number of observations}
}
\value{
The cumulative hazard evaluated at the locations in the input vector
  x
}
\description{
Calculate the Siler cumulative hazard at the locations in the input vector x
given the parameter vector b. Optionally, a conditional starting age can be
input. See details in the documentation for hsiler for a definition of the
Siler hazard and related quantities, including the cumulative hazard.
}
