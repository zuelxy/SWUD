X0 <- as.matrix(read.table("C:/xiaoyao/SWUD/UD25_6.txt",header = TRUE))
design_cd <- AugUD_CD(X0,65,6,inner_it = 100,it=50)

weight <- c(0.9520, 0.9990, 0.1105, 0.2190, 0.0280, 0.0260)
design_wcd <- AugUD_WCD(X0,65,6,weight,inner_it = 100,it=50)



library(DiceKriging)
### computer simulator model
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

### test points
newpoint <- as.matrix(read.table("C:/xiaoyao/SWUD/testpoint.txt",header = T))
newout <- apply(newpoint,1,otlcircuit)


########
# otl_wd <- as.matrix(design_wd$design)
otl_wd <- as.matrix(read.table("C:/xiaoyao/SWUD/AUD_WD65_6.txt",header = T))
res_wd <- apply(otl_wd, 1, otlcircuit)

result_wd <- pred_score(otl_wd, res_wd, newpoint, newout,seed = 19)

########
otl_wwd <- as.matrix(read.table("C:/xiaoyao/SWUD/AUD_WWD65_6.txt",header = T))
res_wwd <- apply(otl_wwd, 1, otlcircuit)

result_wwd <- pred_score(otl_wwd, res_wwd, newpoint, newout,seed = 112)

c(result_wd$mar, result_wwd$mar)
c(result_wd$rmspe, result_wwd$rmspe)


########
otl_cd <- as.matrix(read.table("C:/xiaoyao/SWUD/AUD_CD65_6.txt",header = T))
res_cd <- apply(otl_cd, 1, otlcircuit)

result_cd <- pred_score(otl_cd, res_cd, newpoint, newout,seed = 15)

########
otl_wcd <- as.matrix(read.table("C:/xiaoyao/SWUD/AUD_WCD65_6.txt",header = T))
res_wcd <- apply(otl_wcd, 1, otlcircuit)

result_wcd <- pred_score(otl_wcd, res_wcd, newpoint, newout,seed = 18)

c(result_cd$mar, result_wcd$mar)
c(result_cd$rmspe, result_wcd$rmspe)

