##### Supporting Functions for SFS Manuscript ###
## By Anikó B. Tóth
## 6/07/2017

require(stats)

# produces Diamond and Gilpin A matrix of specified size - for Type II error testing.
DGA <- function(nsp = 50, nloc = 20){ 
  m <- matrix(nrow = nsp, ncol = nloc, data = 0)
  m[1:nsp/2, 1:nloc/2] <- 1
  m[((nsp/2)+1):nsp, ((nloc/2)+1):nloc] <- 1
  return(m)
}

# SFS function: 
# input: an occupancy matrix of species (rows) by sites (columns)
# output: returns a distance matrix of pairwise associations calculated using SFS

SFS_Pairwise <- function(x){ #simpairs function, list of simpairs and fisher's test out.
  samples = ncol(x)  #S
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
  
  occs = array()
  
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #SFS Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # observed number of co-occurrences
      
      #SFS
      for (k in 0:a) # cycle through 0 through observed number of co-occurrences
          z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      
      # subtract half the probability of the observed value
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      z[j,i] = z[i,j]
    }
  }
  return(as.dist(z, diag = F, upper = F))
  
} 

#Fisher's exact test (FET)
# input: an occupancy matrix of species (rows) by sites (columns)
# output: returns a distance matrix of pairwise associations calculated using FET

FET_Pairwise <- function(x){ 
  samples = ncol(x)  #S
  f = matrix(nrow=nrow(x),ncol=nrow(x),data=0)
 
  occs = array()
  #convert to P/A. Occs = rowsums of PA matrix.
  
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # observed number of co-occurrences
      
      for (k in 0:a) # cycle through 0 through observed number of co-occurrences
          f[i,j] = phyper(a, occs[j], samples-occs[j], occs[i], lower.tail = T)
      f[j,i] <-  f[i,j]
      
    }
  }
  return(as.dist(f, diag = F, upper = F))
  
} 
