# Overview
demohaz is an R package with tools for fitting and using hazard functions.
While it is oriented towards demographic analyses, it can be used for other 
applications. Backgound, framing, and exampled can be found in the following
article and the associated github repository:

Price et al. (In review) Robust fitting of Siler models
https://github.com/MichaelHoltonPrice/price_et_al_robust_siler

# Installation
demohaz is not presently on CRAN. It can be installed from github using the
devtools package. If necessary, install devtools and testthat (testthat is
only needed if you want to run the tests as described below).

```R
install.packages(devtools)
install.packages(testthat)
```

```R
devtools::install_github('MichaelHoltonPrice/demohaz')
```

# Siler hazard (quickstart)
demohaz uses the following parameterization of the Siler hazard:

lambda(x) = b[1] * exp(-b[2]*x) + b[3] + b[5] * exp(b[5]*(x-b[4]))

This should be contrasted with the traditional parameterization:

lambda(x) = a[1] * exp(-a[2]*x) + a[3] + a[4] * exp(a[5]*x)

These are just two different ways of writing the same distribution. They are
related via

a[1] = b[1]
a[2] = b[2]
a[3] = b[3]
a[4] = b[5]*exp(-b[4]*b[5])
a[5] = b[5]

We can use the values in Gage and Dyke (1986), Table 2, Level 15 for
illustrations. The traditional parameterization is:

```R
a0 <- c(.175, 1.40, .368 * .01, .075 * .001, .917 * .1)
```

In demohaz, this can be converted to the demohaz parameterization using
trad_to_demohaz_siler_param (demohaz_to_trad_siler_param provides the reverse
transformation):

```R
b0 <- trad_to_demohaz_siler_param(a0)
```

hsiler calculates the hazard. Let's plot it:

```R
xplot <- seq(0,80,by=.1)
plot(xplot, hsiler(xplot, b0), xlab='x', ylab='Hazard', type='l', lwd=3)
```

To create random samples for the parameterization b0, use rsiler:

```R
N <- 1000
x <- rsiler(N,b0)
hist(x, freq=F, xlab='x', ylab='Density', ylim=c(0, .03))
```

Use dsiler to calculate the density, which we can add to the preceding
histogram (we redefine xplot):

```R
xplot <- seq(0,120,by=.1)
lines(xplot, dsiler(xplot, b0), lwd=3)
```

The random sample x can be fit using the fit_siler function. We jitter the
starting point using runif to show that the fitting works (i.e., without having
to start near the true parameter vector, b0).

```R
fit <- fit_siler(x, b0=b0*runif(5, min = .9, max = 1.1),verbose=TRUE, show_plot=TRUE)
```

In the preceding fit, we set two optional parameters to TRUE, verbose and
show_plot, so that (a) information about the fit is printed as it progresses
(because verbose is TRUE) and a plot of the fit is shown (because show_plot is
TRUE). Let's plot the fit:

```R
hist(x, freq=F, xlab='x', ylab='Density', ylim=c(0, .03))
lines(xplot, dsiler(xplot, fit$b), lwd=3)
```

The final piece of functionality we highlight in this quickstart guide is the
parameter x0. This allows analyses that condition on survival to the starting
age x0, which is 0 by default. Let's sample and plot the example distribution
when x0=15:

```R
x0 <- 10
x <- rsiler(N,b0,x0=x0)
hist(x, freq=F, xlab='x', ylab='Density', ylim=c(0, .03))
xplot <- seq(x0, 100,by=.1)
lines(xplot, dsiler(xplot, fit$b, x0=x0), lwd=3)
```

For additional information, including on functionality not covered in this
quick start, see the documentation for demohaz and the github repository
linked to above. In addition, by design the demohaz tests (see below) provide comprehensive coverage of the demohaz functionality. While the tests were not
designed to be user-friendly example code, they do provide comprehensive
examples of demohaz functionality.

# Tests
There are two types of tests: unit tests and functional tests. These tests
take perhaps 30 minutes to run since the functional tests require using
large numbers of observations in (a) the fitting and (b) the checks on the
random draws from the probability density. These checks are essential to show
that the package works as intended and there is no shortcut aside from using
multiple CPUs, which we want to avoid. To run all the tests, clone the package,
change directory into the package, start R, and use devtools::test().

```bash
git clone https://github.com/MichaelHoltonPrice/demohaz
cd demohaz
R
devtools::test()
```

To test set of tests in an individual file (though test_file is now deprecated
and may soon be removed):

```R
devtools::test_file('tests/testthat/test-data_io-unit.R')
```

# Local dev
If you wish to examine or modify the code, use the following sequence of steps
to rebuild the documentation, install locally, and run the tests. This assumes
you have already cloned the repository and set the working directory to the
repository root. You may need to first install the roxygen2 package.

```R
roxygen2::roxygenize()
detach("package:demohaz", unload = TRUE)
install_local('.',force=T,dep=F)
devtools::test()
```
