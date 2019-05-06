##### Supporting Functions for SFS Manuscript ###
## By Anikó B. Tóth
## 6/07/2017

require(stats)

# produces Diamond and Gilpin A matrix of specified size  
  # for producing the two-compartment simulated matrix.
  # also for Type II error testing. 
DGA <- function(nsp = 50, nloc = 20){ 
  m <- matrix(nrow = nsp, ncol = nloc, data = 0)
  m[1:nsp/2, 1:nloc/2] <- 1
  m[((nsp/2)+1):nsp, ((nloc/2)+1):nloc] <- 1
  return(m)
}

# FETmP function: 
# input: an occupancy matrix of species (rows) by sites (columns)
# output: returns a distance matrix of pairwise associations calculated using FETmP

FETmP_Pairwise <- function(x){ #FETmP function
  samples = ncol(x)  #number of samples
  z = matrix(nrow=nrow(x),ncol=nrow(x),data=0) #empty output matrix
  
  occs = array()
  
  #convert to P/A. Occs = rowsums of PA matrix.
  x <- x/x
  x[is.na(x)] <- 0
  occs <- rowSums(x)
  
  #FETmP Algorithm
  for (i in 2:nrow(x))  {
    for (j in 1:(i-1))
    {
      a = length(which(x[i,] > 0 & x[j,] > 0)) # observed number of co-occurrences
      
      #FETmP
      for (k in 0:a) # cycle through 0 through observed number of co-occurrences
          z[i,j] = z[i,j] + choose(occs[j] , k) * choose(samples - occs[j] , occs[i] - k) / choose(samples , occs[i])
      
      # subtract half the probability of the observed value
      z[i,j] = z[i,j] - choose(occs[j] , a) * choose(samples - occs[j] , occs[i] - a) / choose(samples , occs[i]) / 2
      z[j,i] = z[i,j] # symmetric matrix
    }
  }
  return(as.dist(z, diag = F, upper = F))
  
} 

#Fisher's exact test (FET)
# input: an occupancy matrix of species (rows) by sites (columns)
# output: returns a distance matrix of pairwise associations calculated using FET. Cutoff must be applied separately.

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

# A function for calculating Shannon's H
# input: a vector of numbers
# output: value of Shannon's H for input vector

evennessH <- function(vector){
  vector <- na.exclude(vector)
  v <- numeric()
  for(i in 1:length(vector)){
    v[i] <- vector[i]/sum(vector)*log(vector[i]/sum(vector))
  }
  return(-sum(v, na.rm = T))
}


# Function for calculating weighted connectance for postive and negative symmetric associations.
# input: the original species by site matrix, and optionally a distance matrix representing differences between species. 
    # The function calculates FETmP distances automatically if z is left as NULL.
    # The cutoff argument allows you to set an alpha level cutoff for your calculation of weighted connectance, if desired.
# output: vector of two numbers, the positive and negative weighted connectance. 

weighted_connectance_Bersier <- function(m, z= NULL, cutoff=0){
  if(is.null(z)) z <- FETmP_Pairwise(m)
  z[which(abs(z)<= cutoff)] <- NA
  
  z <- as.matrix(z, diag = F, upper = T)
  zp <- zn <- z
  zp[zp<0] <- NA
  zn[zn>0] <- NA
  
  Hpk <- numeric() 
  Hnk <- numeric()
  
  for(i in 1:nrow(z)){
    #for(j in 1:ncol(z)){
    Hpk[i] <- evennessH(zp[i,])
    Hnk[i] <- evennessH(zn[i,])}#}
  
  LDp <- numeric()
  LDn <- numeric()
  for(i in 1:nrow(z)){
    LDp[i] <- sum(zp[i,], na.rm = T)/sum(abs(z[i,]))*2^Hpk[i]
    LDn[i] <- sum(zn[i,], na.rm = T)/sum(abs(z[i,]))*2^Hnk[i]
  }
  
  return(c(sum(LDp, na.rm = T)/(2*nrow(m)), sum(LDn, na.rm = T)/(2*nrow(m))))
}

######## HELPER FUNCTIONS ##########

namerows <- function(table){
  rownames(table) <- table[,1]
  table <- table[,2:ncol(table)]
  return(table)
}

## abundance PA matrices to 0/1 format
tobinary <- function(PA.LIST){
  binary <- lapply(PA.LIST, function(x) {
    x <- x/x
    x[is.na(x)] <- 0
    return(x)
  })
  return(binary)
}

# clean empty rows and columns in 1 matrix
clean.empty <- function(x, mincol = 1, minrow = 1){
  x <- x[which(rowSums(x) > minrow-1),]
  x <- x[,which(colSums(x) > mincol-1)]
  return(x)
}

#dist object to edge list
dist2edgelist <- function(z, sppDat, linktype=NULL){  #edge list with link types attached
  k3 = as.data.frame(melt(as.matrix(z)))
  k3$X1 <- row.names(sppDat)[k3[,1]]
  k3$X2 <- row.names(sppDat)[k3[,2]] #insert species names
  k3 <- data.frame(k3$X1, k3$X2, k3$value, id = paste(k3$X1, k3$X2, sep = "-"))
  colnames(k3) <- c("Sp1", "Sp2", "Score", "id")
  return(k3)
}

# groups pairs as positive, negative, or zero.
posnegzero <- function(x){
  out <- x > 0
  out[which(x>0)] <- "Aggregation"
  out[which(x<0)] <- "Segregation"
  out[which(x==0)] <- "ZERO"
  out
}


### 
percpos <- function(x)
  length(which(x>0))/length(x)

percneg <- function(x)
  length(which(x<0))/length(x)

meanpos <- function(x)
  mean(x[which(x>0)])

meanneg <- function(x)
  mean(x[which(x<0)])
##


