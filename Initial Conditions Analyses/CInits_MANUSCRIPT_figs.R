setwd("C:/Users/eecampbe/Desktop/trial")
library(deSolve)
library(gplots)
library(pscl)
library(rootSolve)
library(coda)
library(compiler)

load("C:/Users/eecampbe/Desktop/CInit_converged.RData")

C1_params=mcmc(x_CInit[[1]][[1]][burnin:n.iter,])
C2_params=mcmc(x_CInit[[2]][[1]][burnin:n.iter,])
C3_params=mcmc(x_CInit[[3]][[1]][burnin:n.iter,])

C1_sigs=mcmc(x_CInit[[1]][[2]][burnin:n.iter,])
C2_sigs=mcmc(x_CInit[[2]][[2]][burnin:n.iter,])
C3_sigs=mcmc(x_CInit[[3]][[2]][burnin:n.iter,])

plot(C1_params)
plot(C1_sigs)

geweke.diag(C1_params)
geweke.diag(C2_params)
geweke.diag(C3_params)

geweke.diag(C1_sigs)
geweke.diag(C2_sigs)
geweke.diag(C3_sigs)

#if thinning is preferred: thin=seq(from=burnin, to=n.iter, by="thin"), then mcmc(<result>[thin,])

#create mcmc list (manual for now)
combined<-mcmc.list(C1_params, C2_params, C3_params)
combined_sigs<-mcmc.list(C1_sigs, C2_sigs, C3_sigs)
combined_var=mcmc(rbind(C1_sigs, C2_sigs, C3_sigs)^2)

print(gelman.diag(combined))
print(gelman.diag(combined_sigs))

print(summary(combined))
print(summary(combined_sigs))

#creat summary of combined chains
all_sum=summary(combined)

par(mfrow=c(2,1))
densplot(combined_sigs[,18:19], xlim=c(0.05, 0.35), ylim=c(0,30), show.obs=FALSE, main="")

#Figure 5 density plot
densplot(combined_var[,18:19],xlim=c(0.0, 0.12), ylim=c(0,200), main="", bty="n")


##########plot of Fs, measured data vs modeled estimate of soluble fraction############
###all data

#Figure 5 soluble fraction vs measured results
#pull mean and credible intervals from summary results from all chains
Mean_Fs_all=all_sum[[1]][2:6,1]
uq_all=all_sum[[2]][2:6,5]
lq_all=all_sum[[2]][2:6,1]

#Plot Fs estimates for all litter, using all data
par(mfrow=c(1,1))
plotCI(Mean_Fs_all, ui=uq_all, li=lq_all, 
       ylab="Soluble fraction", cex.lab=1.5, lwd=4,
       ylim=c(0,1), xaxt="n", cex=3, xlab="Litter type", bty="n")
#HWE data
for(s in 1:ncol(data_real_all$y_Fs[[1]])){
  points(data_real_all$y_Fs[[1]][,s], col="black", 
         type="p", pch=4, lwd=2, cex=2)
}
#Mass-difference data
for(s in 1:ncol(data_real_all$y_Fs[[2]])){
  points(data_real_all$y_Fs[[2]][,s], col="black", 
         type="p", pch=6, lwd=2, cex=2)
}
#plot mean of estimated points
points(Mean_Fs_all, type="p", pch=20, 
       cex=3)
axis(1,1:5,label=litter_nam)
legend("topright", legend=c("Estimated", "Meas-HWE", "Meas-Mass Diff"), 
       pch=c(20,4, 6), pt.cex=c(3,2,2), pt.lwd=2, cex=1.5)

