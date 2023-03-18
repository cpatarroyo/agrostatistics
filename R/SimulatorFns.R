#Declare the individual class
setClass("indiv",slots = list(micro="matrix",pos="integer"))

#' Introduce stepwise mutations to a microsatellite marker
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
#' @return Genind of the individuals of the list produced during the simulation.

sim2genind <- function(x, ploidy) {
  temp_df <- as.data.frame(do.call(rbind,lapply(x,FUN = pastealleles)))
  return(df2genind(temp_df,sep = ";",ploidy = ploidy))
}

#' Paste the alleles of an individual separated by a semicolon
#'
#' @description This function collapses the alleles separating them with a semicolon.
#' @param x Vector of microsatellite loci to be collapsed.

pastealleles <- function(x) {
  apply(x@micro,MARGIN = 2,FUN = function(y) {paste(y,collapse = ";")})
}


asexrep <- function(x) {
  rep(list(x),5)
}

displacement <- function(x,position,probability) {
  y <- x
  y@pos <- sample(position,1,prob = probability[,,x@pos])
  return(y)
}

recombination <- function(x,size,poplist) {
  apply(rbind(poplist[[x[1]]]@micro,poplist[[x[2]]]@micro), FUN = sample,MARGIN = 2,size=size)
}

addinds <- function(x,mrows) {
  new("indiv",micro=matrix(data=x[1:(length(x)-1)],nrow = mrows),pos=as.integer(x[length(x)]))
}

insMutations <- function(x,mutationRate) {
  x@micro <- apply(x@micro,MARGIN = c(1,2), FUN = mutmsat, mutrate = mutationRate)
  return(x)
}

eucprob <- function(x,y) {
  return(dnorm((sqrt(sum((x-y)^2))),mean=0,sd=3))
}

distint <- function(i,j,mat) {
  return(eucprob(which(mat == i, arr.ind = TRUE),which(mat == j, arr.ind = TRUE)))
}

distmat <- function(x, mat) {
  apply(mat,MARGIN = c(1,2),FUN = distint, j = x, mat=mat)
}

posRecomb <- function(x,elList) {
  elList[[x[1]]]@pos
}

simulator<-function(n=100,msat=10,gen=100,ploidy=1,sexprob=0.5, mutrat = 0.001, grid=c(10,10), printAnc = FALSE) {

  #Create the position grid and the dispersal kernel
  posmat <- matrix(1:prod(grid), nrow = grid[1], ncol = grid[2])
  tprobmat <-as.matrix(sapply(posmat,FUN = distmat,mat=posmat))
  probarray <- array(dim = c(grid[1],grid[2],prod(grid)))
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
    temp <- append(temp,unlist(lapply(population,FUN = asexrep)))

    #Displacement of the asexually produced offspring
    temp <- lapply(temp, FUN = displacement, position = posmat, probability = probarray)

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

        #From the possible pairings randomly select the fraction corresponding to the selected probability of sexual reproduction
        pairings <- pairings[, sample(1:ncol(pairings), as.integer(sexprob*ncol(pairings)))]

        if(!is.null(ncol(pairings))) {
          if(ncol(pairings)>0) {
            #Recombine the alleles for each marker in each one of the crossings that will occur
            tempRecomb <- apply(pairings,MARGIN = 2,FUN = recombination, size=ploidy,poplist = population)

            #The position attribute of the parents is stored in the last row of the matrix, since the offspring remains at the same quadrant
            tempRecomb <- rbind(tempRecomb,apply(pairings,MARGIN = 2, FUN= posRecomb, elList=population))

            #Create the individuals with the recombinant genotypes and add them to the population
            temp <- append(temp,apply(tempRecomb,MARGIN = 2,FUN = addinds, mrows = ploidy))
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

  #Return the created population as a genind object
  return(sim2genind(population, ploidy = ploidy))
}

