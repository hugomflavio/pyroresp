% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope_functions.R
\name{filter_r2}
\alias{filter_r2}
\title{Filter  slopes by threshold R2}
\usage{
filter_r2(slopes, r2 = 0.95)
}
\arguments{
\item{slopes}{a data frame obtained by using
the function \code{\link{calc_slopes}}}

\item{r2}{numeric: minimal coefficient of determination (\eqn{r^{2}})
for valid slopes. Default \eqn{r^{2}} = 0.95.}
}
\value{
A data frame with the information about extracted slopes.
}
\description{
Filter  slopes by threshold R2
}
