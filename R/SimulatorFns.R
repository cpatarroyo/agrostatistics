#' Class individual for evolutionary simulations
#'
#' @description Class that contains a matrix of the microsatellite alleles for multiple loci and multiple ploidy levels. It also includes a position slot for spatially explicit inferences.
#' @field micro Matrix that contains the microsatellite alleles (as rows) for each locus (as columns).
#' @field pos Integer. Indicator of the position to facilitate spatially explicit models.
#' @export

setClass("indiv",slots = list(micro="matrix",pos="integer"))

#' Introduce step wise mutations to a microsatellite marker
#'
#' @description This function is used to introduce mutations to microsatellite markers following a stepwise mutation model. Following this, the mutated allele will be either a repetition smaller or a repetition bigger than the wild type allele.
#' @param y Integer. Number of repetitions of the allele.
#' @param mutrate Double. This is the mutation rate for the locus.

mutmsat <- function(y,mutrate) {
  y<-y+(rbinom(1,1,mutrate)*(2*(rbinom(1,1,0.5)-0.5)))
}

#' Converts a simulated population into a Genind object
#'
#' @description This function converts a simulated population into a Genind object. This in order to facilitate the calculation of diversity statistics using the poppr package.
#' @param x List of objects of the class "indiv" produced during the evolutionary simulation.
#' @returns Genind of the individuals of the list produced during the simulation.

sim2genind <- function(x, ploidy) {
  temp_df <- as.data.frame(do.call(rbind,lapply(x,FUN = pastealleles)))
  return(adegenet::df2genind(temp_df,sep = ";",ploidy = ploidy))
}

#' Paste the alleles of an individual separated by a semicolon
#'
#' @description This function collapses the alleles separating them with a semicolon.
#' @param x Vector of microsatellite loci to be collapsed.

pastealleles <- function(x) {
  apply(x@micro,MARGIN = 2,FUN = function(y) {paste(y,collapse = ";")})
}

#' Generate copies of the objects of class "indiv"
#'
#' @description Used to simulate asexual/clonal reproduction in a population of individuals of the class "indiv"
#' @param x Individual to be replicated

asexrep <- function(x) {
  rep(list(x),5)
}

#' Function to sample the new location of a dispersed individual
#'
#' @description This function samples the new position of an individual in the position matrix given the current position of the individual and the dispersal kernel.
#' @param x Individual of the class "indiv" with a position slot.
#' @param position Spatial matrix.
#' @param probability Dispersal kernel. Array of transition probabilities for the movement from each possible slot.
#' @returns Sampled new position within the position matrix given as a parameter.

displacement <- function(x,position,probability) {
  y <- x
  y@pos <- sample(position,1,prob = probability[,,x@pos])
  return(y)
}

#' Function to do the recombination of the alleles in a locus
#'
#' @description Selects alleles from an individual locus for a n-ploid individual.
#' @param x "indiv" individual index.
#' @param size Ploidy level of the individuals.
#' @param poplist Population list. List of the "indiv" individuals.

recombination <- function(x,size,poplist) {
  apply(rbind(poplist[[x[1]]]@micro,poplist[[x[2]]]@micro), FUN = sample,MARGIN = 2,size=size)
}

#' Create "indiv" objects from a matrix of microsatellite data.
#'
#' @description This function takes a matrix of microsatellite loci and creates an "indiv" object with this information.
#' @param x Matrix with the microsatellite repetition data.
#' @param mrows Ploidy level of the population.

addinds <- function(x,mrows) {
  new("indiv",micro=matrix(data=x[1:(length(x)-1)],nrow = mrows),pos=as.integer(x[length(x)]))
}

#' Insert mutations in microsatellite markers following a step wise mutation model
#'
#' @description Inserts a step wise mutation in a microsatellite allele using a binomial probability distribution with the probability specified as the mutationRate parameter.
#' @param x "indiv" individual where the mutations will be introduced.
#' @param mutationRate Mutation rate of the microsatellite markers.

insMutations <- function(x,mutationRate) {
  x@micro <- apply(x@micro,MARGIN = c(1,2), FUN = mutmsat, mutrate = mutationRate)
  return(x)
}

#' Calculates the probability of arrival by distance
#'
#' @description This functions returns the cumulative probability of dispersal from the quadrant x to the quadrant y. This function uses a normal probability distribution with mu = 0 and sigma = 3.
#' @param x Origin of the displacement.
#' @param y Arrival of the displacement.
#' @returns Cumulative probability of arrival from x to y given a normal probability distribution.

eucprob <- function(x,y) {
  return(dnorm((sqrt(sum((x-y)^2))),mean=0,sd=3))
}

#' Calculate the probability of dispersal to each slot
#'
#' @description This function uses the eucprob function to calculate the probability density contained in each slot depending on each origin site.
#' @param i Initial quadrant
#' @param j Arrival quadrant
#' @param mat Position matrix where the individuals are located
#' @returns The probability of dispersal from a site i to all possible quadrants of the defined position grid

distint <- function(i,j,mat) {
  return(eucprob(which(mat == i, arr.ind = TRUE),which(mat == j, arr.ind = TRUE)))
}

#' Calculate the probability of dispersal from each initial quadrant to each possible quadrant
#'
#' @description This function calculates the probabilities of displacement from an x quadrant to every other quadrant in the grid.
#' @param x The quadrant of origin.
#' @param mat The probability matrix of the displacement rates.
#' @returns It returns an array of probabilities of displacement from each of the possible origin quadrants.

distmat <- function(x, mat) {
  apply(mat,MARGIN = c(1,2),FUN = distint, j = x, mat=mat)
}

#' Obtains the position of a descendant individual from a sexual crossing
#'
#' @description After a sexual recombination, this function is used to obtain the location of the parental that is the place where the descendant individual appears.
#' @param x Parental individual.
#' @param elList Population list where all the individuals are stored.
#' @returns It returns the position of the parent x.

posRecomb <- function(x,elList) {
  elList[[x[1]]]@pos
}

#' Spatially explicit evolutionary simulation of organisms with mixed reproduction
#'
#' @description This function performs a spatially explicit evolutionary simulation of a population with mixed reproduction.
#' @details This function allows to specify the dimensions of the position matrix, the probability of sexual reproduction, the amount of microsatellite markers, the ploidy, the amount of generations and the mutation rate.
#'    This function starts by creating an individual starting in a random position of the position grid. This individual will start the first generation by reproducing clonaly. The individuals generated by clonal reproduction will disperse to other quadrants in the grid given the dispersal kernel. Once there are multiple individuals, the individuals that are located in the same quadrant are able to reproduce sexually. The sexual reproduction probability is handled as a conditional probability. The sexual reproduction then has a probability of \code{sexprob} given that the individuals of the crossing are present in the same quadrant (representing here the same site).
#'    After generating the individuals of the next generation, the mutations are introduced to the microsatellite loci with a binomial probability distribution given by p = \code{mutrat}. This will simulate the accumulation of mutations along generations of reproduction of the organisms.
#' @param n Integer. Number of individuals samples for reproduction each generation. This is the number of individuals that will be considered for reproduction. This is done to reduce the computational burden by having an infinitely exponentially growing population. This will also be a proxy of the fact that sampled populations are a small sample of all the individuals that compose an actual population.
#' @param msat Double. Amount of microsatellite loci for your individuals.
#' @param gen Integer. Amount of generations to pass before returning the resulting population.
#' @param ploidy Integer. Ploidy level of the individuals in the population.
#' @param sexprob Double. Conditional probability of sexual reproduction events.
#' @param mutrat Double. Mutation rate of the microsatellite loci.
#' @param grid Integer vector. Dimensions of the spatial grid where the populations occurs.
#' @param printAnc Boolean. If \code{TRUE} the genotype of the first individual created will be printed at the end of the simulation.
#' @returns A named list with two slots. The \code{pop} slot that contains the Genind object containing the simulated population, and the \code{sexprop} slot that contains the proportion of sexual reproduction events in each generation.
#' @export

evoSim_se<-function(n=100,msat=10,gen=100,ploidy=1,sexprob=0.5, mutrat = 0.001, grid=c(10,10), printAnc = FALSE) {

  #Create the position grid and the dispersal kernel
  posmat <- matrix(1:prod(grid), nrow = grid[1], ncol = grid[2])
  tprobmat <-as.matrix(sapply(posmat,FUN = distmat,mat=posmat))
  probarray <- array(dim = c(grid[1],grid[2],prod(grid)))
  #Create the vectors to store the amount of sexual and clonal reproduction events respectively
  sexEvents <- vector()
  clonEvents <- vector()

  for(i in 1:dim(tprobmat)[2]) {
    probarray[,,i] <- matrix(data = tprobmat[,i],nrow = grid[1], ncol = grid[2])
  }

  #Create the initial population
  ind1<-new("indiv",micro=matrix(data = as.integer(runif(msat,min=100,max=200)),nrow = ploidy,ncol = msat),pos=sample(posmat,1))
  population<-list(ind1)

  #Save the initial population
  genmemory<-ind1

  #Pass of generations
  gencnt<-1
  while(gencnt <= gen) {

    if(length(population) > n) {
      population <- sample(population,n,replace = FALSE)
    }

    #Create a temporal list to hold the new generation
    temp <- list()

    #Add the individuals product of asexual reproduction
    temp <- append(temp, unlist(lapply(sample(population, max(as.integer((1-sexprob)*length(population)),1)), FUN = asexrep)))

    #Displacement of the asexually produced offspring
    temp <- lapply(temp, FUN = displacement, position = posmat, probability = probarray)

    #Record the clonal reproduction events for this generation and put a 0 placeholder for sexual reproduction events
    clonEvents[gencnt] <- max(as.integer((1-sexprob)*length(population)),1)
    sexEvents[gencnt] <- 0

    if(sexprob>0) {

      #If the probability of sexual reproduction is bigger than 0, the individuals that are in the same quadrant are located
      tempLoc <- unlist(lapply(population, FUN = function(x) {return(x@pos)}))
      tempPairings <- lapply(unique(tempLoc), FUN = function(x,a) {which(a == x)}, a=tempLoc)

      if(max(sapply(tempPairings, length)) > 1) {

        ##If there is at least a pair of individuals in the same quadrant the recombinants are calculated
        #Only the quadrants with more than one individual are kept
        tempPairings <- tempPairings[which(sapply(tempPairings, FUN = function(x) { length(x)>1 }))]

        #All possible pairings between individuals that are in the same quadrants are arranged
        pairings <- do.call(cbind,lapply(tempPairings,combn,m=2))

        #From the possible pairings randomly select the fraction corresponding to the probability of sexual reproduction for the population size
        pairings <- pairings[, sample(1:ncol(pairings), as.integer(sexprob*length(population)), replace = TRUE)]

        if(!is.null(ncol(pairings))) {
          if(ncol(pairings)>0) {
            #Recombine the alleles for each marker in each one of the crossings that will occur
            tempRecomb <- apply(pairings,MARGIN = 2,FUN = recombination, size=ploidy,poplist = population)

            #The position attribute of the parents is stored in the last row of the matrix, since the offspring remains at the same quadrant
            tempRecomb <- rbind(tempRecomb,apply(pairings,MARGIN = 2, FUN= posRecomb, elList=population))

            #Create the individuals with the recombinant genotypes and add them to the population
            temp <- append(temp,apply(tempRecomb,MARGIN = 2,FUN = addinds, mrows = ploidy))

            #Record the amount of sexual reproduction events that occurred in the current generation
            sexEvents[gencnt] <- ncol(pairings)
          }
        }
      }
    }
    population <- temp

    #Introduce the mutations
    population <- lapply(population,FUN = insMutations,mutationRate=mutrat)

    #Step into the following generation
    gencnt<-gencnt+1
  }

  #Print the last common ancestor of the population
  if(printAnc) {
    print(genmemory)
  }

  #Return a list with a pop slot containing the created population as a genind object and the sexprop that contains the proportion of sexual reproduction events
  results <- list(pop = sim2genind(population, ploidy = ploidy), sexprop = sexEvents/(sexEvents+clonEvents))
  return(results)
}

