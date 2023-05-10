#' Pareto distribution using \code{poweRlaw}
#'
#' @description This function calcultates the negative slope of a Pareto distribution using the \code{poweRlaw} package as described in the \code{poppr} package.
#' @importFrom poweRlaw displ
#' @importFrom poweRlaw estimate_pars
#' @importFrom stats lm
#' @importFrom stats coef
#' @param x Value used to generate the Pareto distribution
#' @return Returns the negative slope of the Pareto distribution

power_law_beta <- function(x){
  xpow <- poweRlaw::displ(x[x > 0])           # Generate the distribution
  xpow$setPars(poweRlaw::estimate_pars(xpow)) # Estimate the parameters
  xdat <- data.frame(x=xpow$internal$values, y= xpow$internal$cum_n/max(xpow$internal$cum_n)) # Extract the data
  xlm <- lm(log(y) ~ log(x), data = xdat)     # Run log-log linear model for slope
  return(-coef(xlm)[2])
}

#' Function to apply the calculation of the negative slope of the Pareto distribution as described in the \code{poppr} package.
#'
#' @description This function is a wrapper to apply the calculation of the negative slope of the Pareto distribution to a table of genotypes.
#' @param x Abundance of a particular genotype.
#' @return Negative slope of the Pareto distribution from a genotype table.

Beta <- function(x){
  x <- drop(as.matrix(x))
  if (length(dim(x)) > 1){
    res <- apply(x, 1, power_law_beta)
  } else {
    res <- power_law_beta(x)
  }
  return(res)
}

#' Function used to generate a reference table from multiple simulations.
#'
#' @description This function is a wrapper to implement the simulation and summarizing of the populations pooled for the ABC inference of the proportion of sexual reproduction.
#' This allows to take advantage of the multi-core execution available from the \code{parallel} R package to make the production of a reference table more efficient.
#' @param x Element of a vector or list of values for the probability of sexual reproduction.
#' @param pop_number Integer. Number of individuals per population per population.
#' @param generations Integer. Number of generations to be simulated for each population.
#' @param ploidy Integer. Ploidy level of the individuals in the simulated populations.
#' @param loci Integer. Amount of microsatellite loci to be simulated for each individual.
#' @param mutation_r Double. Rate of mutation for the microsatellite loci.
#' @param spgrid Int vector. Vector of two values indicating the size of the spatial grid where the simulations take place.
#' @param sample_size Passed on to the \code{evoSim_se} function. Individuals to be sampled to calculate the summary statistics. If \code{sample_size} is more than \code{pop_number}, all individuals will be used for the calculations.
#' @importFrom poppr poppr
#' @importFrom poppr mlg.table
#' @return List of the given probability of sexual reproduction, the proportion of sexual reproduction events, and the summary statistics for each simulated population.

refTab <- function(x, pop_number, generations, ploidy, loci, mutation_r, spgrid, sample_size) {
  tempop <- evoSim_se(sexprob = x, n=pop_number, gen=generations, ploidy = ploidy, msat = loci, mutrat = mutation_r, grid=spgrid, sample_size = sample_size)
  response <- as.numeric(c(x,mean(tempop$sexprop) ,as.numeric(unlist(poppr::poppr(tempop$pop)[,6:12])), Beta(poppr::mlg.table(tempop$pop,plot = FALSE))))
  return(response)
}

#' Function to make the reference table from a pool of simulated populations.
#'
#' @description
#' This function provides an easy-to-use wrapper to generate a pool of simulated populations and calculate the summary statistics from each simulated population.
#' The function samples \code{N} values from a uniform probability distribution between 0 and 1. Each one of the sampled values will be used as a probability of sexual reproduction for each of the simulated populations.
#' The summary statistics are collected in a table that is written into a file at the working directory.
#' Because of the computationally intensive nature of the procedure this function can be executed in several different runs. If the file exists already the function adds a number to the file name to preserve previous runs.
#' @param N Integer. Number of populations to be simulated.
#' @param pop_number Integer. Limit number of individuals to reproduce per generation per population.
#' @param sample_size Integer. Passed on to the \code{evoSim_se} function. Individuals to be sampled to calculate the summary statistics. If \code{sample_size} is more than \code{pop_number}, all individuals will be used for the calculations.
#' @param generations Integer. Passed on to the \code{evoSim_se} function. Number of generations simulate for each population.
#' @param ploidy Integer. Passed on to the \code{evoSim_se} function. Ploidy level of the individuals in the simulated populations.
#' @param loci Integer. Passed on to the \code{evoSim_se} function. Number of microsatellite loci to be simulated for each individual in the population.
#' @param mutation_r Double. Passed on to the \code{evoSim_se} function. Mutation rate of the microsatellite markers.
#' @param spgrid Int vector. Passed on to the \code{evoSim_se} function. The dimensions of the spatial grid where the populations are simulated.
#' @param file_name String. Name of the folder to be created to store the parts of the reference table.
#' @param cores Integer. Amount of parallel cores to run the simulations. **Warning** Windows does not support multi-core execution, so to run this function in Windows the `cores` parameter should be set to 1.
#' @importFrom utils write.csv
#' @export
#' @returns This function does not return a value. It writes out the reference table using the name given in `file_name` to the working directory.

makeRefTable <- function(N, pop_number=100, sample_size=NULL, generations=10000, ploidy=1, loci=10, mutation_r=0.001, file_name=NULL, spgrid=c(10,10), cores=NULL) {

  #Create the vector of probability values for the simulations
  probvec <- runif(N, min=0, max=1)

  #Create the folder to store the
  tryCatch({dir.create(file_name, showWarnings = FALSE) }, error=function(cond) { stop("You need to provide a folder name!") } )

  if(is.null(cores)) {
    cores <- parallel::detectCores()
  }

  results <- parallel::mclapply(probvec, FUN=refTab, mc.cores = cores, pop_number = pop_number, generations = generations, ploidy = ploidy, loci = loci, mutation_r = mutation_r, spgrid = spgrid, sample_size = sample_size)
  tempResults <- as.data.frame(matrix(unlist(results),ncol = 10, byrow = TRUE))
  colnames(tempResults) <- c("Prob", "SexRate","H","G","lambda","E.5","Hexp", "Ia", "rbarD","Pareto")

  nameCount <- 0
  while(file.exists(paste(file_name,nameCount,".csv", sep = ""))) {
    nameCount <- nameCount + 1
  }

  write.csv(tempResults, file=paste("./",file_name,"/",file_name,nameCount,".csv", sep = ""), row.names=FALSE)
}
