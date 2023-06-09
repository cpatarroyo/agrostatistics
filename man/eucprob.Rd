% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/SimulatorFns.R
\name{eucprob}
\alias{eucprob}
\title{Calculates the probability of arrival by distance}
\usage{
eucprob(x, y)
}
\arguments{
\item{x}{Origin of the displacement.}

\item{y}{Arrival of the displacement.}
}
\value{
Cumulative probability of arrival from x to y given a normal probability distribution.
}
\description{
This functions returns the cumulative probability of dispersal from the quadrant x to the quadrant y. This function uses a normal probability distribution with mu = 0 and sigma = 3.
}
