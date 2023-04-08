setwd("~/Desktop/R_script/Augmented_UD/Borehole")

library(DiceKriging)
borehole <- function(xx)
{
  rw <- 0.05 + (0.15-0.05)*xx[1]
  r  <- 100 + (50000-100)*xx[2]
  Tu <- 63070 + (115600-63070)*xx[3]
  Hu <- 990 + (1110-990)*xx[4]
  Tl <- 63.1 + (116-63.1)*xx[5]
  Hl <- 700 + (820-700)*xx[6]
  L  <- 1120 + (1680-1120)*xx[7]
  Kw <- 1500 + (15000-1500)*xx[8]
  
  frac1 <- 2 * pi * Tu * (Hu-Hl)
  
  frac2a <- 2*L*Tu / (log(r/rw)*rw^2*Kw)
  frac2b <- Tu / Tl
  frac2 <- log(r/rw) * (1+frac2a+frac2b)
  
  y <- frac1 / frac2
  return(y)
}

bore_score <- function(design, output, newpoint, newout, seed){
  set.seed(seed)
  m <- km(design=design, response=output)
  p <- predict(m, newpoint,"UK")
  y1 <- sqrt(sum((p$mean-newout)^2)/nrow(newpoint))
  y2 <- median(abs(p$mean-newout))
  return(list(rmspe=y1, mar=y2))
}

### test point 
borepoint <- as.matrix(read.table("testpoint8.txt",header = T))
boreout <- apply(borepoint,1,borehole)
### initial design 
X0 <- as.matrix(read.table("UD25_8.txt",header = T))

### max entropy
max_entropy <- as.matrix(read.table("AugEntropy_bore1.txt", header = T))
Aug_entrop <- rbind(X0,max_entropy)
rownames(Aug_entrop) <- 1:80
entropy_y <- apply(Aug_entrop, 1, borehole)

result_entropy <- bore_score(Aug_entrop, entropy_y, borepoint, 
                             boreout, seed = 97) #63

# [1] 4.550258
# [1] 1.900077

