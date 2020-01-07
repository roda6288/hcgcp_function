
## Example 1 
rm(list=ls(all=TRUE))
## needed packages
library(mvtnorm)
library(expm)

## Generate data
nm <- 500  ## sample size
pm <- 100  ## response size
km <- 100  ## explanatory size
k1m <- 5   ## true explanatory size
rhoy <- 0.8
rhox <- 0.8

Sigmay  <-0.4*((1-rhoy)*diag(pm)+rhoy*matrix(1, nrow=pm, ncol=pm))
Sigmax <- matrix(0,km,km)
for (i in 1:km) {
  for (j in 1:km) {
    Sigmax[i,j] <- rhox^(abs(i-j))
  }
}
X <- rmvnorm(nm,rep(0,km),Sigmax)
Z <- cbind(rep(1, nm),X)
Theta1 <- rmvnorm((k1m+1),rep(0,pm),diag(pm))
Y <- (Z[,1:(k1m+1)]%*%Theta1) + rmvnorm(nm,rep(0,pm),Sigmay)

## Perform hcgcp function
source("hcgcp.R")
result <- hcgcp(X,Y,group = NULL, penalty = NULL)
result