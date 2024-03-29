% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_functions.R
\name{plot_experiment}
\alias{plot_experiment}
\title{Wrapper to plot the whole experiment data for one probe}
\usage{
plot_experiment(input, cycles, probe, verbose = TRUE)
}
\arguments{
\item{input}{An experiment list with the experimental data.
The output of \code{\link{process_experiment}}.}

\item{cycles}{A numeric vector of which cycles to plot}

\item{probe}{The probe to plot (can only plot one probe at a time).}

\item{verbose}{should argument warnings be displayed?}
}
\description{
Wrapper to plot the whole experiment data for one probe
}
