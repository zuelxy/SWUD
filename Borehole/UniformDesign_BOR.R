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

bore_score <- function(design, output, newpoint, newout){
  m <- km(design=design, response=output)
  p <- predict(m, newpoint,"UK")
  y1 <- sqrt(sum((p$mean-newout)^2)/nrow(newpoint))
  y2 <- median(abs(p$mean-newout))
  return(list(rmspe=y1, mar=y2))
}
borepoint <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/testpoint8.txt",header = T))
boreout <- apply(borepoint,1,borehole)

bore_cd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/AUD_CD_80_8.txt",header = TRUE)) 
bore_res_cd <- apply(bore_cd, 1, borehole)

bore_wcd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/AUD_WCD_80_8.txt",header = T))
bore_res_wcd <- apply(bore_wcd, 1, borehole)

bore_wd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/AUD_WD_80_8.txt",header = TRUE)) 
bore_res_wd <- apply(bore_wd, 1, borehole)

bore_wwd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/AUD_WWD_80_8.txt",header = TRUE)) 
bore_res_wwd <- apply(bore_wwd, 1, borehole)

set.seed(41)
bore_cd_gp <- bore_score(bore_cd, bore_res_cd, borepoint, boreout)

set.seed(42)
bore_wcd_gp <- bore_score(bore_wcd, bore_res_wcd, borepoint, boreout) 

set.seed(43)
bore_wd_gp <- bore_score(bore_wd, bore_res_wd, borepoint, boreout)

set.seed(44)
bore_wwd_gp <- bore_score(bore_wwd, bore_res_wwd, borepoint, boreout)

c(bore_cd_gp$rmspe, bore_wcd_gp$rmspe, bore_wd_gp$rmspe, bore_wwd_gp$rmspe)
## [1] 4.727388 4.452884 4.060282 4.055642
c(bore_cd_gp$mar, bore_wcd_gp$mar, bore_wd_gp$mar,  bore_wwd_gp$mar)
## [1] 1.786346 1.668210 1.779171 1.588308 
