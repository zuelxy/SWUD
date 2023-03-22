design_wd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/AUD_WD_80_8.txt", header = T))

design_wwd <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/AUD_WWD_80_8.txt", header = T))


par(mar=c(4, 5, 1, 1), family="STKaiti")
plot(design_wd[1:25,c(1,8)], cex=1.5, cex.lab=1.5,
     xlab = expression(r[w]%*%100) , ylab = expression(K[w]%/%100), 
     xaxt="n", yaxt="n", lwd=1.5)
points(design_wd[26:80,c(1,8)], col="red", cex=1.5, pch = 2, lwd=1.5)
axis(side = 1, at=c(0,0.25,0.5,0.75,1),
     labels = c(5,7.5,10,12.5,15), cex.axis=1.5)
axis(side = 2, at=c(0,0.2,0.4,0.6,0.8,1),
     labels = c(15,42,69,96,123,150), cex.axis=1.5)


par(mar=c(4, 5, 1, 1), family="STKaiti")
plot(design_wwd[1:25,c(1,8)], cex=1.5, cex.lab=1.5,
     xlab = expression(r[w]%*%100) , ylab = expression(K[w]%/%100), 
     xaxt="n", yaxt="n", lwd=1.5)
points(design_wwd[26:80,c(1,8)], col="red", cex=1.5, pch = 2, lwd=1.5)
axis(side = 1, at=c(0,0.25,0.5,0.75,1),
     labels = c(5,7.5,10,12.5,15), cex.axis=1.5)
axis(side = 2, at=c(0,0.2,0.4,0.6,0.8,1),
     labels = c(15,42,69,96,123,150), cex.axis=1.5)


X0 <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/UD25_8.txt", header = T))
design_entropy <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/BOR55_MaxEntropy.txt", header = F))
design_imspe <- as.matrix(read.table("~/Desktop/R_script/Augmented_UD/Borehole/BOR55_IMSPE.txt", header = F))

par(mar=c(4, 5, 1, 1), family="STKaiti")
plot(X0[,c(1,8)], cex=1.5, cex.lab=1.5,
     xlab = expression(r[w]%*%100) , ylab = expression(K[w]%/%100), 
     xaxt="n", yaxt="n", lwd=1.5)
points(design_entropy[,c(1,8)], col="red", cex=1.5, pch = 2, lwd=1.5)
axis(side = 1, at=c(0,0.25,0.5,0.75,1),
     labels = c(5,7.5,10,12.5,15), cex.axis=1.5)
axis(side = 2, at=c(0,0.2,0.4,0.6,0.8,1),
     labels = c(15,42,69,96,123,150), cex.axis=1.5)

par(mar=c(4, 5, 1, 1), family="STKaiti")
plot(X0[,c(1,8)], cex=1.5, cex.lab=1.5,
     xlab = expression(r[w]%*%100) , ylab = expression(K[w]%/%100), 
     xaxt="n", yaxt="n", lwd=1.5)
points(design_imspe[,c(1,8)], col="red", cex=1.5, pch = 2, lwd=1.5)
axis(side = 1, at=c(0,0.25,0.5,0.75,1),
     labels = c(5,7.5,10,12.5,15), cex.axis=1.5)
axis(side = 2, at=c(0,0.2,0.4,0.6,0.8,1),
     labels = c(15,42,69,96,123,150), cex.axis=1.5)

