library(CGP)
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

### prediction performance
pred_score <- function(design, output, newpoint, newout, seed){
  set.seed(seed)
  m <- CGP(design, output)
  # p <- predict(m, newpoint,"UK")
  p <- predict(m, newpoint)
  y1 <- sqrt(sum((p$Yp-newout)^2)/nrow(newpoint))
  y2 <- median(abs(p$Yp-newout))
  return(list(rmspe=y1, mar=y2))
}
### initial design 
X0 <- as.matrix(read.table("UD25_6.txt",header = T))

### test point 
newpoint <- as.matrix(read.table("testpoint.txt",header = T))
newout <- apply(newpoint,1,otlcircuit)

### max entropy
max_entropy <- as.matrix(read.table("OTL40_MaxEntropy.txt", header = F))
Aug_entrop <- rbind(X0,max_entropy)
rownames(Aug_entrop) <- 1:65
entropy_y <- apply(Aug_entrop, 1, otlcircuit)

result_entropy <- pred_score(Aug_entrop, entropy_y, newpoint, 
                             newout, seed = 213)
# [1] 0.03776645
# [1] 0.02024582

### GP_IMSE 
gp_imse <- as.matrix(read.table("OTL40_IMSPE.txt", header = F))
Aug_imse <- rbind(X0, gp_imse)
rownames(Aug_imse) <- 1:65
imse_y <- apply(Aug_imse, 1, otlcircuit)

result_imse <- pred_score(Aug_imse, imse_y, newpoint, 
                             newout, seed = 931)

# [1] 0.04232037
# [1] 0.01601855
