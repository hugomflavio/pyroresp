% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/o2_calcs_functions.R
\name{calc_delta}
\alias{calc_delta}
\title{Calculate oxygen delta for each cycle}
\usage{
calc_delta(input, zero_buffer = 3)
}
\arguments{
\item{input}{a data frame with trimmed measurements.
The output of \code{\link{trim_resp}}.}

\item{zero_buffer}{when calculating the delta, a value must be assigned to 0.
Traditionally, this is the very first value of the cycle. However, due to
natural probe noise, this can cause an upward or downard shift to the whole
delta line. This can be countered by calculating the mean of a few initial
values, instead of relying only on the very first. zero_buffer sets how many
values should be used to estimate the starting O2 concentration for the
cycle.}
}
\value{
The input data frame with calculated deltas for o2 and airsat.
}
\description{
Subtracts the initial oxygen value to the remaining,
for each probe*cycle combination.
}
