## Example 2 
rm(list=ls(all=TRUE))
## needed packages
library(mvtnorm)
library(expm)

## Generate data
nm <- 200 ## sample size
pm <- 2   ## response size
km <- 10  ## explanatory size
k1m <- 5  ## true explanatory size
rhoy <- 0.8
rhox <- 0.8

Sigmay  <-0.4*((1-rhoy)*diag(pm)+rhoy*matrix(1, nrow=pm, ncol=pm))
Sigmax <- matrix(0,km-2,km-2)
for (i in 1:km-2) {
  for (j in 1:km-2) {
    Sigmax[i,j] <- rhox^(abs(i-j))
  }
}
X <- rmvnorm(nm,rep(0,km-2),Sigmax)
Xdummy <- matrix(0,nm,2)
Xdummy[1:floor(nm/3),1] <- 1 
Xdummy[(floor(nm/3)+1):floor(2*nm/3),2] <- 1 
X <- cbind(Xdummy,X)
Z <- cbind(rep(1, nm),X)
Theta1 <- rmvnorm((k1m+1),rep(0,pm),diag(pm))
Y <- (Z[,1:(k1m+1)]%*%Theta1) + rmvnorm(nm,rep(0,pm),Sigmay)

## Perform hcgcp function
group <- rep(0,km)
for (i in 1:km-2) {
  group[i+2] <- i
}
as.character(group)
group[1:2] <- 'dummy'
source("hcgcp.R")
result <- hcgcp(X,Y,group, penalty = NULL)
result