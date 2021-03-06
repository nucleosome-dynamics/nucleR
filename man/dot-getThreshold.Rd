% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpers.R
\name{.getThreshold}
\alias{.getThreshold}
\title{Threshold getter}
\usage{
.getThreshold(threshold, data)
}
\arguments{
\item{threshold}{threshold given as an absolute value or as a string
percentage}

\item{data}{vector with values from which to derive the threshold if it's
relative}
}
\value{
a numeric vector
}
\description{
If threshold is given as a string with percentage, convert it
}
