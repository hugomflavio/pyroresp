% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_mr}
\alias{plot_mr}
\title{Plot the metabolic rates}
\usage{
plot_mr(input, cycles, probes, verbose = TRUE)
}
\arguments{
\item{input}{An experiment list with calculated metabolic rates.
The output of \code{\link{process_experiment}} or
\code{\link{calc_mr}}.}

\item{cycles}{A numeric vector of which cycles to plot}

\item{probes}{A string of which probes to plot}

\item{verbose}{should argument warnings be displayed?}
}
\value{
a ggplot object
}
\description{
Plot the metabolic rates
}
