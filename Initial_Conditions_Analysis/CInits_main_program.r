#Initial DOC pool estimate
##########NOTES##########
#This script will run MCMC estimates of initial LIDEL condition parameters
#using all data for HWE and mass difference measure of soluble litter fraction
#and then plotting mean and 95% BCIs of the results
#Nell Campbell
#11/17/14
##########NOTES##########

##########SETUP FUNCTIONS AND DATA##########
  rm(list=ls())
  setwd("C:/LIDEL/Initial_Conditions_Analysis")
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
	#run scripts to load priors, functions, data, and setup storage for MCMC
	  source("CInits_Priors.R")
	  source("CInits_Convergence.R")
	  source("CInits_data_upload.R")
	  source("CInits_functions_all.R")
	  source("CInits_likelihood.R")
	  n.iter=100000
	  burnin=50000
	  source("CInits_MCMC_setup.R")
  
##########RUN MCMC for all data##########
	  x_CInit=Run_MCMC(set=set_Cinits, num_L=length(latent_inits), num_litter=length(litter_nam), 
				 data_y=data_real_all, n.iter=n.iter, params=params, param_names=param_names,
				 burnin=burnin)

##########Visualize results##########		
	#create MCMC objects 
	#if thinning is preferred: thin=seq(from=burnin, to=n.iter, by="thin"), then mcmc(<result>[thin,])
		C1_params=mcmc(x_CInit[[1]][[1]][burnin:n.iter,])
		C2_params=mcmc(x_CInit[[2]][[1]][burnin:n.iter,])
		C3_params=mcmc(x_CInit[[3]][[1]][burnin:n.iter,])

		C1_sigs=mcmc(x_CInit[[1]][[2]][burnin:n.iter,])
		C2_sigs=mcmc(x_CInit[[2]][[2]][burnin:n.iter,])
		C3_sigs=mcmc(x_CInit[[3]][[2]][burnin:n.iter,])
		
	#plot results posterior distributions
		plot(C1_params)
		plot(C1_sigs)

	#run geweke diagnostics
		geweke.diag(C1_params)
		geweke.diag(C2_params)
		geweke.diag(C3_params)

		geweke.diag(C1_sigs)
		geweke.diag(C2_sigs)
		geweke.diag(C3_sigs)

		
	#create mcmc list with all chains
	  combined<-mcmc.list(C1_params, C2_params, C3_params)
	  combined_sigs<-mcmc.list(C1_sigs, C2_sigs, C3_sigs)

	  print(gelman.diag(combined))
	  print(gelman.diag(combined_sigs))

	  print(summary(combined))
	  print(summary(combined_sigs))
    

##########SAVE RESULTS##########
		save.image(file="C:/LIDEL/Initial_Conditions_Analysis/results/CInit_converged.RData")



	