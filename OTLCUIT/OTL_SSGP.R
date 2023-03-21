setwd("~/Desktop/R_script/Augmented_UD/Otlcircuit")
library(DiceKriging)

corr_matrix <- function (X, theta) 
{
  if (!is.matrix(X)) {
    X = as.matrix(X)
  }
  d = ncol(X)
  n = nrow(X)
  if (d != length(theta)) {
    stop("The dimensions of theta and X do not match. \n")
  }
  R <- matrix(0, nrow = n, ncol = n)
  rcoord <- cbind(rep(seq_len(n - 1L), times = rev(seq_len(n - 1L))), 
                  unlist(lapply(X = rev(seq_len(n - 1L)), 
                                FUN = function(nn, nm) seq_len(nn) + nm - nn, nm = n)))
  absdiff <- abs(X[rcoord[, 2L], ] - X[rcoord[, 1L], ])
  absdiff <- absdiff^2
  Theta <- matrix(theta^2, nrow = (length(absdiff)/d), ncol = d, 
                  byrow = TRUE)
  Rtemp <- Theta * absdiff
  Rtemp <- rowSums(Rtemp)
  R[rcoord] <- Rtemp
  R <- R + t(R)
  R <- exp(-R)
  R
}
otlcircuit <- function(xx)
{
  Rb1  <- 50+(150-50)*xx[1]
  Rb2  <- 25+(70-25)*xx[2]
  Rf   <- 0.5+(3-0.5)*xx[3]
  Rc1  <- 1.2+(2.5-1.2)*xx[4]
  Rc2  <- 0.25+(1.2-0.25)*xx[5]
  beta <- 50+(300-50)*xx[6]
  
  Vb1 <- 12*Rb2 / (Rb1+Rb2)
  term1a <- (Vb1+0.74) * beta * (Rc2+9)
  term1b <- beta*(Rc2+9) + Rf
  term1 <- term1a / term1b
  
  term2a <- 11.35 * Rf
  term2b <- beta*(Rc2+9) + Rf
  term2 <- term2a / term2b
  
  term3a <- 0.74 * Rf * beta * (Rc2+9)
  term3b <- (beta*(Rc2+9)+Rf) * Rc1
  term3 <- term3a / term3b
  
  Vm <- term1 + term2 + term3
  return(Vm)
}

d <- 6
X0 <- as.matrix(read.table("UD25_6.txt",header = TRUE))
y0 <- apply(X0,1,otlcircuit)
y0 <- as.matrix(y0)
##################### MLE kriging

set.seed(111)
mle_gp <- km(design = X0, response = y0, covtype ="gauss")

##################### SSVS
### initial_value
n <- nrow(X0)
p <- ncol(X0)

mu <- coef(mle_gp)$trend
sd2 <- coef(mle_gp)$sd2
theta <- 1/(coef(mle_gp)$range)
# 1.1062314 0.9881612 0.5208333 0.5208333 0.5208333 0.5208333
gm <- rep(1, p)

### hyper-parameter
B <- diag(1, p)
# diff <- apply(UDdesign, 2, function(xx) range(xx)[2]-range(xx)[1])
# names(diff) <- NULL

# delta_y <- max(resp)-min(resp) # as.vector(sd(resp))
tau <- rep(0.08, p)  #delta_y/(3*diff) 
c <- rep(50, p)
# v_gm <- 0   # sigma^2 prior parameter
# lm_gm <- 1 # can be any value
phat=rep(0.5, p)

D_gm <- function(gm, c, tau){
  a <- c^gm
  D_gm <- diag(a*tau)
  return(D_gm)
}

pi_theta <- function(theta, design, resp, mu, sigma2, D_gm, B){
  n <- nrow(design)
  R_theta <- corr_matrix(design, theta)
  R_inv <- solve(R_theta)
  
  Q1 <- -0.5*log(det(R_theta))
  Q2 <- -0.5*t(resp - mu*rep(1,n))%*%R_inv%*%(resp - mu*rep(1,n))/sigma2
  Q3 <- -0.5*t(theta)%*%solve((D_gm%*%B%*%D_gm))%*%theta
  return(Q1+Q2+Q3)
}

m <- 5000
######### keep track of the stuff
keep_mu <- rep(0, m)
keep_sigma2 <- rep(0,m)
keep_theta <- matrix(0, m, p)
keep_gm <- matrix(0, m, p)

######### Gibbs sampling
begin_time <- proc.time()
for (i in 1:m) {
  R_theta <- corr_matrix(X0, theta)
  R_inv <- solve(R_theta)
  
  ### update mu
  mu_mean <- t(rep(1,n))%*%R_inv%*%y0/sum(R_inv)
  mu_cov <- sd2/sum(R_inv)
  mu <- rnorm(1, mu_mean, sqrt(mu_cov))
  
  ### update sd2
  Q_theta <- t(y0 - mu*rep(1,n))%*%R_inv%*%(y0 - mu*rep(1,n))
  sd2 <- 1/rgamma(1, n/2,  Q_theta/2)
  
  ### update theta
  D <- D_gm(gm, c, tau)
  theta_star <- rep(NA, p)
  sig_theta <-  rep(0.05, p)
  for (k in 1:p) {
    theta_star[k] <- theta[k] + rnorm(1, 0, sig_theta[k])
  }
  rate <- exp(min(0, pi_theta(theta_star,X0,y0,mu,sd2,D,B)-
                    pi_theta(theta,X0,y0,mu,sd2,D,B)))
  d = rbinom(1,1,rate)
  theta = theta_star*d + theta*(1-d)
  #keep_rate[m,k] <- rate
  
  ### update gamma
  for (j in 1:p) {
    a <- dnorm(theta[j], mean = 0, sd = c[j]*tau[j])#*phat[j]
    b <- dnorm(theta[j], mean = 0, sd = tau[j])#*(1-phat[j])
    p_hat <- a/(a+b)
    gm[j] <- rbinom(1, 1, p_hat)
    if(is.na(gm[j]))
      stop("BUG!!!")
  }
  
  ####Store the output
  keep_mu[i] <- mu 
  keep_sigma2[i] <- sd2
  keep_theta[i,] <- theta
  keep_gm[i,] <- gm
}
end_time <- proc.time()-begin_time

gm_burn <- keep_gm[1001:5000, ]
apply(gm_burn,2,mean)
