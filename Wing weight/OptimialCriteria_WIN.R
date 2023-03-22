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
wing_score <- function(design, output, newpoint, newout, seed){
  set.seed(seed)
  m <- km(~., design=design, response=output)
  p <- predict(m, newpoint,"UK")
  y1 <- sqrt(sum((p$mean-newout)^2)/nrow(newpoint))
  y2 <- median(abs(p$mean-newout))
  return(list(rmspe=y1, mar=y2))
}

### test point 
wingpoint <- as.matrix(read.table("testpoint10.txt",header = T))
wingout <- apply(wingpoint,1,wingweight)

### initial design 
X0 <- as.matrix(read.table("UD35_10.txt",header = T))

### max entropy
max_entropy <- as.matrix(read.table("WIN65_MaxEntropy.txt", header = F))
Aug_entrop <- rbind(X0,max_entropy)
rownames(Aug_entrop) <- 1:100
entropy_y <- apply(Aug_entrop, 1, wingweight)

result_entropy <- wing_score(Aug_entrop, entropy_y, wingpoint, 
                             wingout, seed = 172)

result_entropy$rmspe
result_entropy$mar
# [1] 2.404055
# [1] 1.632854

### GP_IMSE 
gp_imse <- as.matrix(read.table("WIN65_IMSPE.txt", header = F))
Aug_imse <- rbind(X0, gp_imse)
rownames(Aug_imse) <- 1:100
imse_y <- apply(Aug_imse, 1, wingweight)

result_imse <- wing_score(Aug_imse, imse_y, wingpoint, 
                          wingout, seed = 527)

result_imse$rmspe
result_imse$mar

# [1] 2.369156
# [1] 1.159762
