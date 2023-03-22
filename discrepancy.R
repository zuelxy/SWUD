CD2 <- function(design){
  X <- as.matrix(design)
  dimension <- dim(X)[2] 
  n <- dim(X)[1]
  if (n < dimension) {
    stop("Warning : the number of points is lower than the dimension.")
  }
  if (min(X) < 0 || max(X) > 1) {
    warning("The design is rescaling into the unit cube [0,1]^d.")
    X <- (X-0.5)/n
  }
  s1 <- 0
  s2 <- 0
  for (i in 1:n) {
    p <- prod((1 + 0.5 * abs(X[i, ] - 0.5) - 0.5 * ((abs(X[i,] - 0.5))^2)))
    s1 <- s1 + p
    for (k in 1:n) {
      q <- prod((1 + 0.5 * abs(X[i, ] - 0.5) + 0.5 * 
                   abs(X[k, ] - 0.5) - 0.5 * abs(X[i, ] - X[k,])))
      s2 <- s2 + q
    }
  }
  DisC2 = sqrt(((13/12)^dimension) - ((2/n) * s1) + ((1/n^2) * s2))
  return(DisC2)
}

WD2 <- function(design){
  X <- as.matrix(design)
  dimension <- dim(X)[2] 
  n <- dim(X)[1]
  if (n < dimension) {
    stop("Warning : the number of points is lower than the dimension.")
  }
  if (min(X) < 0 || max(X) > 1) {
    warning("The design is rescaling into the unit cube [0,1]^d.")
    X <- (X-0.5)/n
  }
  s1 <- 0
  for (i in 1:n) {
    for (k in 1:n) {
      p <- prod((1.5 - ((abs(X[i, ] - X[k, ])) * (1 - abs(X[i, ] - X[k, ])))))
      s1 <- s1 + p
    }
  }
  DisW2 = sqrt((-((4/3)^dimension) + ((1/n^2) * s1)))
  return(DisW2)
}


WCD2 <- function(design, weight){
  X <- as.matrix(design)
  dimension <- dim(X)[2]
  n <- dim(X)[1]
  if (n < dimension) {
    stop("Warning : the number of points is lower than the dimension.")
  }
  if (min(X) < 0 || max(X) > 1) {
    warning("The design is rescaling into the unit cube [0,1]^d.")
    X <- (X-0.5)/n
  }
  s1 <- 0
  s2 <- 0
  for (i in 1:n) {
    p <- prod((1 + 0.5 * weight * (abs(X[i, ] - 0.5) - ((abs(X[i,] - 0.5))^2))))
    s1 <- s1 + p
    for (k in 1:n) {
      q <- prod((1 + 0.5 * weight * (abs(X[i, ] - 0.5) + 
                   abs(X[k, ] - 0.5) -  abs(X[i, ] - X[k,]))))
      s2 <- s2 + q
    }
  }
  DisC2 = sqrt(prod(1+weight/12) - ((2/n) * s1) + ((1/n^2) * s2))
  return(DisC2)
}

WWD2 <- function(design, weight){
  X <- as.matrix(design)
  dimension <- dim(X)[2]
  n <- dim(X)[1]
  if (n < dimension) {
    stop("Warning : the number of points is lower than the dimension.")
  }
  if (min(X) < 0 || max(X) > 1) {
    warning("The design is rescaling into the unit cube [0,1]^d.")
    X <- (X-0.5)/n
  }
  s1 <- 0
  for (i in 1:n) {
    for (k in 1:n) {
      p <- prod(1 + weight *(0.5 + ((abs(X[i, ] - X[k, ])) * (1 - abs(X[i, ] - X[k, ])))))
      s1 <- s1 + p
    }
  }
  DisW2 = sqrt(-prod(1+weight/3) + ((1/n^2) * s1))
  return(DisW2)
}


