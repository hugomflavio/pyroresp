% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aux_functions.R
\name{conv_w_to_ml}
\alias{conv_w_to_ml}
\title{convert weight to volume in ml}
\usage{
conv_w_to_ml(w, d = 1)
}
\arguments{
\item{w}{a units object in either grams or kilograms}

\item{d}{The density. defaults to 1 (1g = 1ml). Defaults to 1.
Higher densities lead to lower volumes, e.g. if d = 2, then 1g = 0.5ml.}
}
\value{
a units object in ml
}
\description{
convert weight to volume in ml
}
