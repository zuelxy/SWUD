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
wing_score <- function(design, output, newpoint, newout){
  m <- km(~., design=design, response=output)
  p <- predict(m, newpoint,"UK")
  y1 <- sqrt(sum((p$mean-newout)^2)/nrow(newpoint))
  y2 <- median(abs(p$mean-newout))
  return(list(rmspe=y1, mar=y2))
}
wingpoint <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Wingweight/testpoint10.txt",header = T))
wingout <- apply(wingpoint,1,wingweight)

wing_cd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Wingweight/AUD_CD_100_10.txt",header = TRUE)) 
wing_res_cd <- apply(wing_cd, 1, wingweight)

wing_wcd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Wingweight/AUD_WCD_100_10_SSGP.txt",header = T))
wing_res_wcd <- apply(wing_wcd, 1, wingweight)

wing_wd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Wingweight/AUD_WD_100_10.txt",header = TRUE)) 
wing_res_wd <- apply(wing_wd, 1, wingweight)

wing_wwd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Wingweight/AUD_WWD_100_10_SSGP.txt",header = TRUE)) 
wing_res_wwd <- apply(wing_wwd, 1, wingweight)

set.seed(1)
wing_cd_gp <- wing_score(wing_cd, wing_res_cd, wingpoint, wingout)

set.seed(2)
wing_wcd_gp <- wing_score(wing_wcd, wing_res_wcd, wingpoint, wingout)

set.seed(3)
wing_wd_gp <- wing_score(wing_wd, wing_res_wd, wingpoint, wingout)

set.seed(4)
wing_wwd_gp <- wing_score(wing_wwd, wing_res_wwd, wingpoint, wingout)

c(wing_cd_gp$rmspe, wing_wcd_gp$rmspe, wing_wd_gp$rmspe, wing_wwd_gp$rmspe)
## [1] 1.572331 1.294706 1.837064 1.697304
c(wing_cd_gp$mar, wing_wcd_gp$mar, wing_wd_gp$mar, wing_wwd_gp$mar)
## [1] 0.7280769 0.7107076 0.7989006 0.7326095
