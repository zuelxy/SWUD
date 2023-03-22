setwd("~/Desktop/R_script/Augmented_UD/Wingweight")
library(DiceKriging)

wingweight <- function(xx)
{
  Sw      <- 150+(200-150)*xx[1]
  Wfw     <- 220+(300-220)*xx[2]
  A       <- 6+(10-6)*xx[3]
  LamCaps <- ((-10)+(10-(-10))*xx[4]) * (pi/180)
  q       <- 16+(45-16)*xx[5]
  lam     <- 0.5+(1-0.5)*xx[6]
  tc      <- 0.08+(0.18-0.08)*xx[7]
  Nz      <- 2.5+(6-2.5)*xx[8]
  Wdg     <- 1700+(2500-1700)*xx[9]
  Wp      <- 0.025+(0.08-0.025)*xx[10]
  
  fact1 <- 0.036 * Sw^0.758 * Wfw^0.0035
  fact2 <- (A / ((cos(LamCaps))^2))^0.6
  fact3 <- q^0.006 * lam^0.04
  fact4 <- (100*tc / cos(LamCaps))^(-0.3)
  fact5 <- (Nz*Wdg)^0.49
  
  term1 <- Sw * Wp
  
  y <- fact1*fact2*fact3*fact4*fact5 + term1
  return(y)
}

d <- 10
X0 <- as.matrix(read.table("UD35_10.txt",header = TRUE))
y0 <- apply(X0,1,wingweight)
y0 <- as.matrix(y0)
##################### MLE kriging

set.seed(11)
mle_gp <- km(design = X0, response = y0, covtype ="gauss")

##################### SSVS
### initial_value
n <- nrow(X0)
p <- ncol(X0)

mu <- coef(mle_gp)$trend
sd2 <- coef(mle_gp)$sd2
theta <- 1/(coef(mle_gp)$range)
gm <- rep(1, p)

# [1] 0.7516920 0.5147059 0.9786522 0.5147059 0.5147059 0.5147059 0.7664742
# [8] 1.1567524 0.6091629 0.5147059
### hyper-parameter
B <- diag(1, p)
# diff <- apply(UDdesign, 2, function(xx) range(xx)[2]-range(xx)[1])
# names(diff) <- NULL

# delta_y <- max(resp)-min(resp) # as.vector(sd(resp))
tau <- rep(0.05, p)  #delta_y/(3*diff) 
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
apply(gm_burn, 2, mean)
