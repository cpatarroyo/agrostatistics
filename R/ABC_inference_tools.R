#' Pareto distribution using \code{poweRlaw}
#'
#' @description This function calcultates the negative slope of a Pareto distribution using the \code{poweRlaw} package as described in the \code{poppr} package.
#' @importFrom poweRlaw displ
#' @importFrom poweRlaw estimate_pars
#' @importFrom stats lm
#' @importFrom stats coef
#' @returns Returns the negative slope of the Pareto distribution

power_law_beta <- function(x){
  xpow <- poweRlaw::displ(x[x > 0])           # Generate the distribution
  xpow$setPars(poweRlaw::estimate_pars(xpow)) # Estimate the parameters
  xdat <- plot(xpow, draw = FALSE)            # Extract the data
  xlm <- lm(log(y) ~ log(x), data = xdat)     # Run log-log linear model for slope
  return(-coef(xlm)[2])
}

#' Function to apply the calculation of the negative slope of the Pareto distribution as described in the \code{poppr} package.
#'
#' @description This function is a wrapper to apply the calculation of the negative slope of the Pareto distribution to a table of genotypes.
#' @returns Negative slope of the Pareto distribution from a genotype table.

Beta <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- apply(x, 1, power_law_beta)
  } else {
    res <- power_law_beta(x)
  }
  return(res)
}
