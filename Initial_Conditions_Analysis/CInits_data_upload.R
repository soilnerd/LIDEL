##########ALL DATA UPLOAD########
#Upload all data for Initial conditions
#Nell Campbell
#Created 11/17/14
#################################

litter_nam=c("Alfalfa", "Ash", "Bluestem", "Oak", "Pine")

#load data
HWE=read.csv("HWE.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
DIFF=read.csv("DIFF.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
C6alpha=read.csv("C6alpha.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
C1_3alpha=read.csv("C1_3alpha.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
LIG=read.csv("LIG.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
C2alpha=read.csv("C2alpha.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

#create list of fraction soluble measurements
y_Fs=list(as.matrix(HWE[,2:ncol(HWE)]/100),as.matrix(DIFF[,2:ncol(DIFF)]/100))
#create list of all data
data_real_all=list(y_Fs, as.matrix(C1_3alpha[,2:ncol(C1_3alpha)]), 
                   as.matrix(C6alpha[,2:ncol(C6alpha)]),
                   as.matrix(LIG[,2:ncol(LIG)])/100,
                   as.matrix(C2alpha[,2:ncol(C2alpha)]))
names(data_real_all)<-c("y_Fs", "y_C1_3", "y_C6alph", "LIG", "C2alpha")

rownames(data_real_all$y_Fs[[1]])<-litter_nam
rownames(data_real_all$y_Fs[[2]])<-litter_nam
rownames(data_real_all$y_C1_3)<-litter_nam
rownames(data_real_all$y_C6alph)<-litter_nam
rownames(data_real_all$LIG)<-litter_nam
rownames(data_real_all$C2alpha)<-litter_nam

par(mfrow=c(1,1))
#plot of measured data for soluble fraction, HWE and by mass difference
plot(data_real_all$y_Fs[[1]][,1], col=c("orange", "lightblue", "lightgreen", "darkgrey", "pink"), type="p", 
     pch=1, main="Soluble fraction measured data", ylab="Measured soluble fraction",
     ylim=c(0,1), xaxt="n", cex=3, xlab="Litter type", lwd=3)
axis(1,1:5,label=litter_nam)
for(s in 2:ncol(data_real_all$y_Fs[[1]])){
  points(data_real_all$y_Fs[[1]][,s], col=c("orange", "lightblue", "lightgreen", "darkgrey", "pink"), 
         type="p", pch=1, cex=3, lwd=3)
}
#points for mass difference-based measurements
for(s in 1:ncol(data_real_all$y_Fs[[2]])){
  points(data_real_all$y_Fs[[2]][,s], col=c("orange", "lightblue", "lightgreen", "darkgrey", "pink"), type="p", pch=2, 
         cex=3, lwd=3)
}
legend("topright", legend=c("HWE", "Mass Diff."), pch=c(1, 2))

par(mfrow=c(2,1))
#plot of measured data for mass
plot(data_real_all$y_C1_3[,1], col=c("red", "blue", "darkgreen", "black", "purple"), type="p", 
     main="Measured Initial litter carbon", ylab="mass (g C)", lwd=3,
     xaxt="n", xlab="Litter type", cex=2, pch=2, ylim=c(min(data_real_all$y_C1_3)-10, max(data_real_all$y_C1_3)+10))
axis(1,1:5,label=litter_nam)
for(s in 2:ncol(data_real_all$y_C1_3)){
  points(data_real_all$y_C1_3[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
         type="p", cex=2, pch=2, lwd=3)
}

#plot of measured data for initial C6 DOC
plot(data_real_all$y_C6alph[,1], col=c("red", "blue", "darkgreen", "black", "purple"), type="p", 
     main="Measured Initial C6, DOC", ylab="mass (g C)", lwd=3,
     xaxt="n", xlab="Litter type", cex=2, pch=2, ylim=c(min(data_real_all$y_C6alph)-10, max(data_real_all$y_C6alph)+10))
axis(1,1:5,label=litter_nam)
for(s in 2:ncol(data_real_all$y_C6alph)){
  points(data_real_all$y_C6alph[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
         type="p", cex=2, pch=2, lwd=3)
}

#plot measured data for lignin by litter
par(mfrow=c(1,1))
plot(data_real_all$LIG[,1], col=c("red", "blue", "darkgreen", "black", "purple"), type="p", 
     main="Measured Lignin, High vs Low C", ylab="Fraction lignin C in litter C",
     xaxt="n", xlab="Litter type", cex=2, pch=2, ylim=c(0,.5), lwd=3)
axis(1,1:5,label=litter_nam)
for(s in 2:ncol(data_real_all$LIG)){
  points(data_real_all$LIG[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
         type="p", cex=2, pch=2, lwd=3)
}

par(mfrow=c(2,1))
#plot measured data for C2 (non-soluble, structural) by litter
plot(data_real_all$C2alpha[,1], col=c("red", "blue", "darkgreen", "black", "purple"), type="p", 
     main="Measured C2, Non-Soluble Structural", ylab="mass (g C)",
     xaxt="n", xlab="Litter type", cex=2, pch=2, lwd=3,
     ylim=c(min(data_real_all$C2alpha)-10, max(data_real_all$C2alpha)+10))
axis(1,1:5,label=litter_nam)
for(s in 2:ncol(data_real_all$C2alpha)){
  points(data_real_all$C2alpha[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
         type="p", cex=2, pch=2, lwd=3)
}

#plot of measured data for initial C6 DOC
plot(data_real_all$y_C6alph[,1], col=c("red", "blue", "darkgreen", "black", "purple"), type="p", 
     main="Measured Initial C6, DOC", ylab="mass (g C)", lwd=3,
     xaxt="n", xlab="Litter type", cex=2, pch=2, ylim=c(min(data_real_all$y_C6alph)-10, max(data_real_all$y_C6alph)+10))
axis(1,1:5,label=litter_nam)
for(s in 2:ncol(data_real_all$y_C6alph)){
  points(data_real_all$y_C6alph[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
         type="p", cex=2, pch=2, lwd=3)
}