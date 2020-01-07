###############################################
### hcgcp_function 2020/1/7
###
### This is the function to select explanatory variables by the EZKB selection method using the HCGCp criterion in multivariate linear regression.
### 
### Usage: 
###   hcgcp(X, Y, group = NULL, penalty = NULL)
###
### Arguments:
###   X: A n*k expalanatory matrix; n observations on the rows and k explanatory variables on the columns.
###   Y: A n*p response matrix; n observations on the rows and p response variables on the columns.
###   group: A k-dimensional vector to identify groups about explanatory variables.
###   penalty: A value of the penalty in the GCp criterion. If penalty=2, the Cp criterion is used. The default is the HCGCp criterion.
###
### Value:
###   Returns vectors of the length of the number of groups about explanatory variables, with selection results and the values of the differences of criteria.
###############################################
hcgcp <- function(X, Y, group = NULL, penalty = NULL){
  
  if (!is.matrix(X))
    stop("X has to be a matrix")
  
  if (any(is.na(X)))
    stop("Missing values in X not allowed!")
  
  if (!is.matrix(Y))
    stop("Y has to be a matrix")
  
  if (any(is.na(Y)))
    stop("Missing values in Y not allowed!")
  
  dimy <- dim(Y)
  n <- dimy[1]
  p <- dimy[2]
  k <- dim(X)[2]
  
  if (n!=dim(X)[1])
    stop("Sample sizes of Y and X have to be same")
  
  if (is.null(penalty)){
    N <- n-p-k
    if (N<5)
      stop("The default penalty cannot be used because sample size is too small")
    penalty <- ((n-k-1)/(N-2)) + ((n-k-1)*sqrt(N+p-4)*((k+1)^(1/4))*log(n))/((N-2)*sqrt(N-4)*sqrt(p))
  } else{
    if (!is.numeric(penalty)||(length(penalty)>1))
      stop("Penalty has to be a scalar")
  }
  
  ## Group label set up
  if (is.null(group)){
    group <- 1:k
  } else{
    if (!is.character(group)){
      stop("Group label has to be a character")
    }
    if (length(group)!=k)
      stop("Lengths of group label and X have to be same")
  }
  
  samegr <- intersect(group,group)
  maxgr <- length(samegr)
  
  ## Some needed calculations
  X <- cbind(X,rep(1, n))
  if (qr(X)$rank!=k+1)
    stop("X must have a full rank")
  tX <- t(X)
  invXX <- solve(tX%*%X)
  B <- X%*%invXX
  Pxc <- diag(n) - B%*%tX
  
  tY <- t(Y)
  S <- tY%*%Pxc%*%Y/(n-k-1)
  if (qr(S)$rank!=p)
    stop("Sample covariance matrix has to be non-singular")
  invS <- solve(S)
  
  ## Performe selecting variables
  boxnumber <- rep("not-select", maxgr) ## box for saving results for selection
  difgcp <- rep(0, maxgr) ## box for saving the values of differences of criteria
  elecou <- 1 ## counter
  if((p==min(p,k+1))&&(p<100)){ ## for the case that p is small
    SL <- sqrtm(invS)%*%tY%*%X
    AB <- SL%*%invXX
    for (j in samegr) {
      vec <- which(group==j)
      xc <- AB[,vec,drop=FALSE]
      R <- invXX[vec,vec]
      invR <- solve(R)
      difgcp[elecou] <- sum(diag(invR%*%t(xc)%*%xc)) - length(vec)*p*penalty ## differences of criteria
      if(difgcp[elecou]>0){
        boxnumber[elecou] <- "select"
      }
      elecou <- elecou + 1
    }
  } else{ ## for the case that p is large
    SL <- Y%*%invS%*%tY
    for (j in samegr) {
      vec <- which(group==j)
      xc <- B[,vec,drop=FALSE]
      R <- invXX[vec,vec]
      invR <- solve(R)
      difgcp[elecou] <- sum(diag(invR%*%t(xc)%*%SL%*%xc)) - length(vec)*p*penalty ## differences of criteria
      if(difgcp[elecou]>0){
        boxnumber[elecou] <- "select"
      }
      elecou <- elecou + 1
    }
  }
  
  names(boxnumber) <- samegr
  names(difgcp) <- samegr
  
  ## Out put
  return(list(selection=boxnumber,value=difgcp))
}
