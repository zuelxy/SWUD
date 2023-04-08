corr.matrix.ISO <- function(X, theta)
{
  n <- dim(X)[1]
  d <- dim(X)[2]
  Theta <- diag(theta)
  U <- matrix(apply(X^2%*%Theta, 1, sum), nrow = n, ncol = n ,byrow=F)
  V <- -2*X%*%Theta%*%t(X)
  Dist <- U + t(U) + V
  return(exp(-Dist))
}

cross.corr.matrix <- function(D.old, D.new, theta)
{
    n1 <- dim(D.new)[1]
    n2 <- dim(D.old)[1]
    d <- dim(D.new)[2]
    Theta <- diag(theta)
    
    U <- matrix(apply(D.new^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=F)
    V <- -2*D.new%*%Theta%*%t(D.old)
    W <- matrix(apply(D.old^2%*%Theta, 1, sum), nrow = n1, ncol = n2 ,byrow=T)
    Dist <- U + V + W
    
    return(exp(-Dist))
}

#**************************************************************************************
#******************* The Entropy criterion for augmented designs  ******************
#**************************************************************************************
Augmented.Entropy <- function(D.old, D.new, theta, R.old.Inv)
{
  n.new <- dim(D.new)[1]
  R.cross <- cross.corr.matrix(D.old, D.new, theta)
  R.new <- corr.matrix.ISO(D.new, theta)
  E <- try(-det(R.new - R.cross%*%R.old.Inv%*%t(R.cross), tol = 1e-16), silent = T)
  return(E)
}

#**************************************************************************************
#***************
#***************** Finding the augmented Maximum Entropy Design *****************
#**************************************************************************************
Augm.Entropy.optim <- function(D.old, n.new, d, theta, n.starts)
{
  vals <- rep(0, n.starts)
  designs <- list()
  R.old <- corr.matrix.ISO(D.old, theta)
  R.old.Inv <- solve(R.old, tol = 1e-16)
    
  Entropy.univar <- function(D.new)
  {
    D.new <- matrix(D.new, ncol = d)
    Augmented.Entropy(D.old, D.new, theta, R.old.Inv)
    }
    
  for(k in 1:n.starts)
  {
    start <- c(optimumLHS(n.new,d))
    optimum <- optim(start, Entropy.univar, method = "L-BFGS-B", lower = rep(0,n.new*d), 
                       upper = rep(1,n.new*d))
    D.optimum <- matrix(optimum$par, ncol = d)
    Entropy.optimum <- optimum$value
    designs[[k]] <- D.optimum
    vals[k] <- Entropy.optimum
    }
    
  min.val <- vals[which.min(vals)]
  best.design <- designs[[which.min(vals)]]
    
  return(list(Design = best.design, log.entropy = -min.val))
}
