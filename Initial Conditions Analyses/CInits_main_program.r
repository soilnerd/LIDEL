#Initial DOC pool estimate
##########NOTES##########
#This script will run MCMC estimates of initial DOC parameters
#(currently only for C6)
#using all data versus HWE and mass difference separately
#and then plotting mean and CIs of the results
#Nell Campbell
#11/17/14
##########NOTES##########

  ##########SETUP FUNCTIONS AND DATA##########
  rm(list=ls())
  setwd("C:/Users/eecampbe/Dropbox/Work/getting my Phd/DOC project/DOC model")
  par(mfrow=c(3,1))
  library(deSolve)
  library(gplots)
  library(pscl)
  library(rootSolve)
  library(coda)
  
  ##########EXECUTE FUNCTION DEFINITIONS##########
  proc_C6alpha<-function(Fs,Fdoc, zC1_3a){
    C6alpha<-Fs*Fdoc*zC1_3a
    return(C6alpha)
  }
  
  proc_C2alpha<-function(Fs,Flig, zC1_3a){
    C2alpha<-(1-(Fs+Flig))*zC1_3a
    return(C2alpha)
  }
  ##########N.ITER & BURN###########
  
  source("CInits_Priors.R")
  source("CInits_Convergence.R")
  #source("Simulate_data.R")
  source("CInits_data_upload.R")
  source("CInits_functions_all.R")
  #source("CInits_functions_HWE.R")
  #source("CInits_functions_ydiff.R")
  source("CInits_likelihood.R")
  n.iter=100000
  burnin=50000
  source("CInits_MCMC_setup.R")
  ##########RUN MCMC for all data##########
  x_CInit=Run_MCMC(set=set_Cinits, num_L=length(latent_inits), num_litter=length(litter_nam), 
             data_y=data_real_all, n.iter=n.iter, params=params, param_names=param_names,
             burnin=burnin)

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

  print(gelman.diag(combined))
  print(gelman.diag(combined_sigs))

  print(summary(combined))
  print(summary(combined_sigs))

save.image(file="C:/Users/eecampbe/Desktop/CInit_converged.RData")



###HWE data
  #pull mean and credible intervals from summary results from all chains
    Mean_Fs_HWE=HWE_sum[[1]][2:6,1]
    uq_Fs_HWE=HWE_sum[[2]][2:6,5]
    lq_Fs_HWE=HWE_sum[[2]][2:6,1]
      
        #Plot Fs estimates for all litter, using HWE data
          par(mfrow=c(1,1))
          plotCI(Mean_Fs_HWE, ui=uq_Fs_HWE, li=lq_Fs_HWE, main="Soluble Litter Fraction, Measured vs Estimated, HWE Data", 
                 ylab="Soluble fraction",
                 ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
        #HWE data
          for(s in 1:ncol(data_real_all$y_Fs[[1]])){
            points(data_real_all$y_Fs[[1]][,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
                   type="p", pch=1, lwd=2)
          }
        #plot mean of estimated points
          points(Mean_Fs_HWE, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
                 cex=2)
          axis(1,1:5,label=litter_nam)
          legend("topright", legend=c("HWE", "Estimated"), pch=c(1, 2, 20))

#####ydiff
  #pull mean and credible intervals from summary results from all chains
    Mean_Fs_ydiff=ydiff_sum[[1]][2:6,1]
    uq_Fs_ydiff=ydiff_sum[[2]][2:6,5]
    lq_Fs_ydiff=ydiff_sum[[2]][2:6,1]

        #Plot Fs estimates for all litter, using ydiff data
          par(mfrow=c(1,1))
          plotCI(Mean_Fs_ydiff, ui=uq_Fs_ydiff, li=lq_Fs_ydiff, main="Soluble Litter Fraction, Measured vs Estimated, Ydiff Data", 
                 ylab="Soluble fraction",
                 ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
        #Mass-difference data
          for(s in 1:ncol(data_real_all$y_Fs[[2]])){
            points(data_real_all$y_Fs[[2]][,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
                   type="p", pch=2, lwd=2)
          }
        #plot mean of estimated points
          points(Mean_Fs_ydiff, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
                 cex=2)
          axis(1,1:5,label=litter_nam)
          legend("topright", legend=c("Mass Diff.", "Estimated"), pch=c(1, 2, 20))

##########Plot Flig, modeled vs measured##############
###all data
  #pull mean and credible interval for Flig given all chains
    Mean_Flig_all=all_sum[[1]][7:11,1]
    uq_Flig_all=all_sum[[2]][7:11,5]
    lq_Flig_all=all_sum[[2]][7:11,1]

          #Plot Flig estimates for all litter, using all data
            plotCI(Mean_Flig_all, ui=uq_Flig_all, li=lq_Flig_all, main="Lignin Litter Fraction, Measured vs Estimated", 
                   ylab="Lignin fraction",
                   ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
          #Flig data
            for(s in 1:ncol(data_real_all$LIG)){
              points(data_real_all$LIG[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
                     type="p", pch=4, lwd=2)
            }
          #plot mean of estimated points
            points(Mean_Flig_all, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
                   cex=2)
            axis(1,1:5,label=litter_nam)
            legend("topright", legend=c("Measured", "Estimated"), pch=c(4, 20))

###HWE data
  #pull mean and credible interval for Flig given all chains
    Mean_Flig_HWE=HWE_sum[[1]][7:11,1]
    uq_Flig_HWE=HWE_sum[[2]][7:11,5]
    lq_Flig_HWE=HWE_sum[[2]][7:11,1]

        #Plot Flig estimates for all litter, using hwe data
          plotCI(Mean_Flig_HWE, ui=uq_Flig_HWE, li=lq_Flig_HWE, main="HWE Lignin Litter Fraction, Measured vs HWE Estimated", 
                 ylab="Lignin fraction",
                 ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
        #Flig data
          for(s in 1:ncol(data_real_all$LIG)){
            points(data_real_all$LIG[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
                   type="p", pch=4, lwd=2)
          }
        #plot mean of estimated points
          points(Mean_Flig_HWE, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
                 cex=2)
          axis(1,1:5,label=litter_nam)
          legend("topright", legend=c("Measured", "Estimated"), pch=c(4, 20))

###ydiff data
  #pull mean and credible interval for Flig given all chains
    Mean_Flig_ydiff=ydiff_sum[[1]][7:11,1]
    uq_Flig_ydiff=ydiff_sum[[2]][7:11,5]
    lq_Flig_ydiff=ydiff_sum[[2]][7:11,1]

        #Plot Flig estimates for all litter, using ydiff data
          plotCI(Mean_Flig_ydiff, ui=uq_Flig_ydiff, li=lq_Flig_ydiff, main="Lignin Litter Fraction, Measured vs ydiff Estimated", 
                 ylab="Lignin fraction",
                 ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
        #Flig data
          for(s in 1:ncol(data_real_all$LIG)){
            points(data_real_all$LIG[,s], col=c("red", "blue", "darkgreen", "black", "purple"), 
                   type="p", pch=4, lwd=2)
          }
        #plot mean of estimated points
          points(Mean_Flig_ydiff, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
                 cex=2)
          axis(1,1:5,label=litter_nam)
          legend("topright", legend=c("Measured", "Estimated"), pch=c(4, 20))

################plot Fnss, fraction non-soluble structural##################
#calculate non-soluble structural fraction mean and credible interval given all chains
    #calculate NSS
    NSS_all=Calc_NSS(x_CInit)
    NSS_HWE=Calc_NSS(x_CInit_HWE)
    NSS_ydiff=Calc_NSS(x_CInit_ydiff)

    NSS_mc=mcmc.list(NSS_all[[1]], NSS_all[[2]], NSS_all[[3]])
    NSS_mc_HWE=mcmc.list(NSS_HWE[[1]], NSS_HWE[[2]], NSS_HWE[[3]])
    NSS_mc_ydiff=mcmc.list(NSS_ydiff[[1]], NSS_ydiff[[2]], NSS_ydiff[[3]])

    gelman.diag(NSS_mc)
    gelman.diag(NSS_mc_HWE)
    gelman.diag(NSS_mc_ydiff)

    NSS_sum=summary(NSS_mc)
    NSS_sum_HWE=summary(NSS_mc_HWE)
    NSS_sum_ydiff=summary(NSS_mc_ydiff)

###all data
    Mean_Fnss_all=NSS_sum[[1]][,1]
    uq_Fnss_all=NSS_sum[[2]][,5]
    lq_Fnss_all=NSS_sum[[2]][,1]

    #Plot Fnss estimates for all litter, using all data
      plotCI(Mean_Fnss_all, ui=uq_Fnss_all, li=lq_Fnss_all, main="Non-soluble Structural Litter Fraction, Estimated", 
             ylab="Non-soluble Structural fraction",
             ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
    #plot mean of estimated points
      points(Mean_Fnss_all, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
             cex=2)
      axis(1,1:5,label=litter_nam)
      legend("topright", legend=c("Estimated"), pch=c(20))

###HWE data
    Mean_Fnss_HWE=NSS_sum_HWE[[1]][,1]
    uq_Fnss_HWE=NSS_sum_HWE[[2]][,5]
    lq_Fnss_HWE=NSS_sum_HWE[[2]][,1]

    #Plot Fnss estimates for all litter, using HWE data
      plotCI(Mean_Fnss_HWE, ui=uq_Fnss_HWE, li=lq_Fnss_HWE, main="Non-soluble Structural Litter Fraction, HWE Estimated", 
             ylab="Non-soluble Structural fraction",
             ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
    #plot mean of estimated points
      points(Mean_Fnss_HWE, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
             cex=2)
      axis(1,1:5,label=litter_nam)
      legend("topright", legend=c("Estimated"), pch=c(20))

###ydiff data
    Mean_Fnss_HWE=NSS_sum_ydiff[[1]][,1]
    uq_Fnss_HWE=NSS_sum_ydiff[[2]][,5]
    lq_Fnss_HWE=NSS_sum_ydiff[[2]][,1]

    #Plot Fnss estimates for all litter, using ydiff data
      plotCI(Mean_Fnss_HWE, ui=uq_Fnss_HWE, li=lq_Fnss_HWE, main="Non-soluble Structural Litter Fraction, ydiff Estimated", 
             ylab="Non-soluble Structural fraction",
             ylim=c(0,1), xaxt="n", cex=2, xlab="Litter type")
    #plot mean of estimated points
      points(Mean_Fnss_HWE, col=c("red", "blue", "darkgreen", "black", "purple"), type="p", pch=20, 
             cex=2)
      axis(1,1:5,label=litter_nam)
      legend("topright", legend=c("Estimated"), pch=c(20))

################Plot Fdoc estimate, based on all 3 analyses############
#plot estimate of Fdoc, fraction of soluble material that is not microbially processed or released
    #all data, vs individual measurement-derived

      All_Fdoc=c(all_sum[[1]][1,1], HWE_sum[[1]][1,1], ydiff_sum[[1]][1,1])
      All_Fdoc_uci=c(all_sum[[2]][1,5], HWE_sum[[2]][1,5], ydiff_sum[[2]][1,5])
      All_Fdoc_lci=c(all_sum[[2]][1,1], HWE_sum[[2]][1,1], ydiff_sum[[2]][1,1])

      plotCI(All_Fdoc, ui=All_Fdoc_uci, li=All_Fdoc_lci, main="Fdoc, mean with 95% CI", 
             ylim=c(0,.5), ylab="Fdoc", xlab="Measurement data", xaxt="n", cex=2)
      points(All_Fdoc, col="blue", cex=2, pch=20)
      axis(1,1:3,label=c("All data", "HWE data", "Difference data"))

    

##########sAVE RESULTS##########
save.image("CInits_all_runs.RData")

	