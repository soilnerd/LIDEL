#################################
#Analyses and figures for Initial Condition Analysis
#
#Nell Campbell
#10/16/15
#################################

#Setup
  rm(list=ls())
  setwd("C:/LIDEL/Initial_Conditions_Analysis")
  library(deSolve)
  library(gplots)
  library(pscl)
  library(rootSolve)
  library(coda)
  library(compiler)

#load Initial Condition Analysis workspace
  load("C:/LIDEL/Initial_Conditions_Analysis/results/CInit_converged.RData")

#create MCMC objects  
  C1_params=mcmc(x_CInit[[1]][[1]][burnin:n.iter,])
  C2_params=mcmc(x_CInit[[2]][[1]][burnin:n.iter,])
  C3_params=mcmc(x_CInit[[3]][[1]][burnin:n.iter,])
  
  C1_sigs=mcmc(x_CInit[[1]][[2]][burnin:n.iter,])
  C2_sigs=mcmc(x_CInit[[2]][[2]][burnin:n.iter,])
  C3_sigs=mcmc(x_CInit[[3]][[2]][burnin:n.iter,])

#Output geweke diagnostic results for each chain, parameter and error
  sink("C:/LIDEL/Initial_Conditions_Analysis/results/Initial_conditions_summary_CHAINS.txt")
    print(geweke.diag(C1_params))
    print(geweke.diag(C2_params))
    print(geweke.diag(C3_params))
    
    print(geweke.diag(C1_sigs))
    print(geweke.diag(C2_sigs))
    print(geweke.diag(C3_sigs))
  sink()

#create mcmc list with all chains
  #if thinning is preferred: 
  #thin=seq(from=burnin, to=n.iter, by="thin"), then mcmc(<result>[thin,])
	  combined<-mcmc.list(C1_params, C2_params, C3_params)
	  combined_sigs<-mcmc.list(C1_sigs, C2_sigs, C3_sigs)
	  combined_var=mcmc(rbind(C1_sigs, C2_sigs, C3_sigs)^2)

#print gelman diagnostics and summary of values
  sink("C:/LIDEL/Initial_Conditions_Analysis/results/Initial_conditions_summary_ALL.txt")
	  print(gelman.diag(combined))
	  print(gelman.diag(combined_sigs))
	  
	  print(summary(combined))
	  print(summary(combined_sigs))
	  print(summary(combined_var))
  sink()

#creat object with summary of combined chains
  all_sum=summary(combined)
  var_sum=summary(combined_var)

##################PLOTS#######################
#density plots of variance for HWE and mass difference measurements
  par(mfrow=c(2,1))
  densplot(combined_sigs[,18:19], xlim=c(0.05, 0.35), ylim=c(0,30), show.obs=FALSE)

#plot of Fs, measured data vs modeled estimate of soluble fraction for all chains
  #Figure 5B soluble fraction vs measured results
    #pull mean and credible intervals from summary results from all chains
      Mean_Fs_all=all_sum[[1]][2:6,1]
      uq_all=all_sum[[2]][2:6,5]
      lq_all=all_sum[[2]][2:6,1]
    
    #Print plot Fs estimates for all litter, using all measured data
      par(mfrow=c(1,1))
      plotCI(Mean_Fs_all, ui=uq_all, li=lq_all, 
             ylab="Soluble fraction", cex.lab=1.5, cex.axis=1.5, lwd=4,
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
             cex=2)
      axis(1,1:5,label=litter_nam, cex.axis=1.5)
      legend("topright", legend=c("Estimated", "Meas-HWE", "Meas-Mass Diff"), 
             pch=c(20,4, 6), pt.cex=c(2,2,2), pt.lwd=2, cex=1.5)
  
  #Export plot for Figure 5B, Fs estimates for all litter, using all measured data
    jpeg(file="C:/LIDEL/Initial_Conditions_Analysis/results/Figure5B.jpg")
    
      par(mfrow=c(1,1))
      plotCI(Mean_Fs_all, ui=uq_all, li=lq_all, 
             ylab="Soluble fraction", cex.lab=1.5, cex.axis=1.5, lwd=4,
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
             cex=2)
      axis(1,1:5,label=litter_nam, cex.axis=1.5)
      legend("topright", legend=c("Estimated", "Meas-HWE", "Meas-Mass Diff"), 
             pch=c(20,4, 6), pt.cex=c(2,2,2), pt.lwd=2, cex=1.5)
    dev.off()
    
#plot variance for figure 5A
    #pull mean and credible intervals from summary results from all chains
      Mean_var_all=var_sum[[1]][18:19,1]
      uq_var_all=var_sum[[2]][18:19,5]
      lq_var_all=var_sum[[2]][18:19,1]
  
    #Plot measurement variance estimates for soluble fraction
      par(mfrow=c(1,1))
      plotCI(Mean_var_all, ui=uq_var_all, li=lq_var_all, cex.axis=1.5,
             ylab="Variance of Soluble Fraction Measurements", cex.lab=1.5, lwd=4,
             ylim=c(0,0.1), xaxt="n", xlim=c(.8,2.2), cex=3, xlab="Measurement Type", bty="n")
      points(Mean_var_all, type="p", pch=20, 
             cex=2)
      axis(1,at=c(1,2),label=c("Hot Water Extraction", "Mass Difference"), cex.axis=1.5)
      axis(1,at=c(.8, 2.2), labels=c("",""), lwd.ticks=0)

  #Export plot for Figure 5, Fs estimates for all litter, using all measured data
    jpeg(file="C:/LIDEL/Initial_Conditions_Analysis/results/Figure5A.jpg")
    
    #Plot measurement variance estimates for soluble fraction
      par(mfrow=c(1,1))
      plotCI(Mean_var_all, ui=uq_var_all, li=lq_var_all, cex.axis=1.5,
             ylab="Variance of Soluble Fraction Measurements", cex.lab=1.5, lwd=4,
             ylim=c(0,0.1), xaxt="n", xlim=c(.8,2.2), cex=3, xlab="Measurement Type", bty="n")
      points(Mean_var_all, type="p", pch=20, 
             cex=2)
      axis(1,at=c(1,2),label=c("Hot Water Extraction", "Mass Difference"), cex.axis=1.5)
      axis(1,at=c(.8, 2.2), labels=c("",""), lwd.ticks=0)
    dev.off()