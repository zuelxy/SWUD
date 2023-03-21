# X0: exist runs in [0,1]^d 
# n : total number of runs including existing runs 

existed_runs <- function(X0, n, d){
  # dividing [0,1]^d to n^d new levels 
  ud_space <- matrix(nrow = n, ncol = d)
  for (k in 1:d) {
    ud_space[,k] <- seq(0,1,length.out = n)
  }
  # mapping X0 to new-level space [1,2,...,n]^d
  for (i in 1:nrow(X0)) {
    for (j in 1:ncol(X0)) {
      X0[i,j] <- which.min(abs(X0[i,j] - ud_space[, j]))
    }
  }
  # initialize D_initial = rbind(D_exist, D_random), D_exist=X0
  level <- 1:n
  exist_level <- X0[,1]
  remain_level <- level[-exist_level]
  D_random = matrix(nrow = (n-nrow(X0)), ncol = d)
  for (k in 1:d) {
    D_random[,k] = sample(remain_level)
  }
  return(list(exist_level=X0, augmt_level=D_random))
}

AUD_WD2_EP <- function(m,k,n_e){
  n <- nrow(m)
  G <- m
  i <- trunc(runif(2,n_e,n))+1   # random two elements in k-th column
  i1 = i[1]
  i2 = i[2]
  x <- G[i1,k]
  G[i1,k] <- G[i2,k]
  G[i2,k] <- x             # element exchange
  dC2 <- WD2(G)
  l = list(dC2,G)
  return(l) 
}

#####Function to extract elements from lists composed with lists#####
extract_list <- function(l){return(l[[1]])}

################ 
AugUD_WD <- function(X0, n, d, inner_it=100, it=2)
{
  initial_AUD <- existed_runs(X0, n, d)$augmt_level
  initial_Aruns <- (initial_AUD-0.5)/n
  design <- rbind(X0, initial_Aruns)
  m <- design
  
  T0 <- 0.005*WD2(m)
  J <- min(50, 0.1*nrow(m)*(nrow(m)-1))
  crit <- NULL
  temp <- NULL
  proba <- NULL
  n_e <- nrow(X0)
  
  Temperature <- T0
  Best <- m
  dC2 <- WD2(m)                   
  best <- dC2
  crit <- dC2
  
  for (q in 1:it)      # outer loop
  {
    BOLD <- Best       # BOLD = new LHS built at each iteration q
    bold <- best       # Best = new LHS built at every step         
    
    ni <- 0
    count <- 0
    na <- 0
    while(count<=inner_it)  # inner loop
    {
      count <- count+1
      
      modulo <- count%%d   # selecteed column $count mod d$
      l <- list(m)
      l <- rep(l,J)        # prepare J element exchange design
      
      g <- lapply(l, AUD_WD2_EP, k = modulo+1, n_e = n_e)
      values <- lapply(g, extract_list)   # the discrepancy of J designs
      k <- which.min(values)
      a <- values[[k]]    # the minimum discrapancy 
      
      Delta <- a-dC2
      
      if((Delta)<=(Temperature*runif(1))){  
        # higher is the temperature, higher is the probability of accepting a bad design.
        # if Delta is low, the probability is high of accepting a bad design.   
        # if Delta>Temperature, m is never accept.
        m <- g[[k]][[2]]
        dC2 <- a
        na <- na+1
        if(a<=best){
          Best <- m
          best <- a
          ni <- ni+1   #if optimization is ok, ni=ni+1
        }                       
      }
      crit <- c(crit,best)
    }
    
    v1 <- na/inner_it    # v1 <- acceptance ratio
    v2 <- ni/inner_it    # v2 <- optimization ratio
    
    temp <- c(temp,rep(Temperature,inner_it)) 
    proba <- c(proba,rep(v1,inner_it))   
    
    if (best-bold<0){
      f <- 1
      if(v1>=0.1 & v2<=v1){
        Temperature<-0.8*Temperature
      }
      else {
        if (v1>=0.1 & v2==v1){} 
        else {Temperature<-Temperature/0.8}
      }  
    }
    # if the criteria is optimized after the inner loop, then f equals 1
    else {
      f <- 0
      if (v1<=0.1){Temperature <- Temperature/0.7}
      else {Temperature <- Temperature*0.9}
    }
    # else, it is the exploratory step
  }
  List.res <- list(design,T0,inner_it,J,it,Best,crit,temp,proba)
  names(List.res) <- c("InitialDesign","TO","inner_it","J","it","design","critValues","tempValues","probaValues") 
  return(List.res)
} 
