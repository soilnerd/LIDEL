##################################
#Analyses to pull in results for all chains for a given model,
#test for convergence, and output final data and figures
#
#
##################################

rm(list=ls())
setwd("C:/LIDEL")
library(deSolve)
library(gplots)
library(pscl)
library(rootSolve)
library(coda)
library(compiler)

#load results
  #if thinning is needed: thin=seq(from=burnin, to=n.iter, by="thin"), then mcmc(<result>[thin,])
  #Model1 parameter chains
    Param_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/RESULTchain1_params.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_param_chain1=mcmc(Param_chain1[,2:ncol(Param_chain1)])
    Param_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/RESULTchain2_params.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_param_chain2=mcmc(Param_chain2[,2:ncol(Param_chain2)])
    Param_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/RESULTchain3_params.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_param_chain3=mcmc(Param_chain3[,2:ncol(Param_chain3)])
  
  #Model1 variance chains
    Var_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/RESULTchain1_var.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_Var_chain1=mcmc(Var_chain1[,2:ncol(Var_chain1)])
    Var_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/RESULTchain2_var.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_Var_chain2=mcmc(Var_chain2[,2:ncol(Var_chain2)])
    Var_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/RESULTchain3_var.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_Var_chain3=mcmc(Var_chain3[,2:ncol(Var_chain3)])
  
  #Model1 latent chains
    Latent_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/RESULTchain1_latent.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_Latent_chain1=mcmc(Latent_chain1[,2:ncol(Latent_chain1)])
    Latent_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/RESULTchain2_latent.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_Latent_chain2=mcmc(Latent_chain2[,2:ncol(Latent_chain2)])
    Latent_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/RESULTchain3_latent.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      MC_Latent_chain3=mcmc(Latent_chain3[,2:ncol(Latent_chain3)])

#create mcmc list (manual for now)
  combined_param<-mcmc.list(MC_param_chain1, MC_param_chain2, MC_param_chain3)
  combined_var<-mcmc.list(MC_Var_chain1, MC_Var_chain2, MC_Var_chain3)
  combined_latent<-mcmc.list(MC_Latent_chain1, MC_Latent_chain2, MC_Latent_chain3)

#output summary for final results
  sink("C:/LIDEL/FINAL_results/Model_1_SUMMARY.txt")
    print(gelman.diag(combined_param))
    print(gelman.diag(combined_var))
    print(gelman.diag(combined_latent))
    
    
    print(summary(combined_param))
    print(summary(combined_var))
    print(summary(combined_latent))
  sink()

#Reload model run files in order to have data files for model comparisons
  setwd("C:/LIDEL/Model1/Chain1")
  #load MCMC_data
  source("LIDEL_v2_1.R")
  #load function to initialize C pools
  source("LIDEL_initial_function.R")
  #load data
  source("LIDEL_data_upload.R")
  #load parameters and initial conditions for model
  source("LIDEL_set_params.R")
  #load priors
  source("LIDEL_Priors.R")
  #load likelihoods
  source("LIDEL_likelihood.R")
  #load MCMC and setup functions
  source("LIDEL_runMCMC.R")
  n.iter=80
  burnin=40
  source("LIDEL_MCMC_setup_trial.R")
  
##################MASS#######################
#load model outputs for each litter and measurement type, then summarize data
#create measured/modeled figures using model results from all chains
  #ALFALFA
    AlfMass_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_mass_alf1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfMass_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_mass_alf2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfMass_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_mass_alf3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfMass_all=mcmc(rbind(AlfMass_chain1[,2:3], 
                           AlfMass_chain2[,2:3], 
                           AlfMass_chain3[,2:3]))
      alf1Mass=summary(AlfMass_all)
      mean_alf1Mass=(alf1Mass[[1]][,1])
      uq_alf1Mass=(alf1Mass[[2]][,5])
      lq_alf1Mass=(alf1Mass[[2]][,1])

  #ASH
    AshMass_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_mass_ash1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshMass_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_mass_ash2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshMass_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_mass_ash3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshMass_all=mcmc(rbind(AshMass_chain1[,2:3], 
                           AshMass_chain2[,2:3], 
                           AshMass_chain3[,2:3]))
      ash1Mass=summary(AshMass_all)
      mean_ash1Mass=(ash1Mass[[1]][,1])
      uq_ash1Mass=(ash1Mass[[2]][,5])
      lq_ash1Mass=(ash1Mass[[2]][,1])

  #BLUESTEM
    BluMass_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_mass_blu1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluMass_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_mass_blu2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluMass_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_mass_blu3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluMass_all=mcmc(rbind(BluMass_chain1[,2:3], 
                           BluMass_chain2[,2:3], 
                           BluMass_chain3[,2:3]))
      bluestem1Mass=summary(BluMass_all)
      mean_bluestem1Mass=(bluestem1Mass[[1]][,1])
      uq_bluestem1Mass=(bluestem1Mass[[2]][,5])
      lq_bluestem1Mass=(bluestem1Mass[[2]][,1])

  #OAK
    OakMass_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_mass_oak1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakMass_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_mass_oak2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakMass_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_mass_oak3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakMass_all=mcmc(rbind(OakMass_chain1[,2:3], 
                           OakMass_chain2[,2:3], 
                           OakMass_chain3[,2:3]))
      oak1Mass=summary(OakMass_all)
      mean_oak1Mass=(oak1Mass[[1]][,1])
      uq_oak1Mass=(oak1Mass[[2]][,5])
      lq_oak1Mass=(oak1Mass[[2]][,1])

  #PINE
    PinMass_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_mass_pin1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinMass_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_mass_pin2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinMass_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_mass_pin3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinMass_all=mcmc(rbind(PinMass_chain1[,2:3], 
                           PinMass_chain2[,2:3], 
                           PinMass_chain3[,2:3]))
      pine1Mass=summary(PinMass_all)
      mean_pine1Mass=(pine1Mass[[1]][,1])
      uq_pine1Mass=(pine1Mass[[2]][,5])
      lq_pine1Mass=(pine1Mass[[2]][,1])


  #create lists of results to use for plotting
    summary_litterMass=list(mean_alf1Mass, mean_ash1Mass, mean_bluestem1Mass, mean_oak1Mass, mean_pine1Mass)
    new_mass=c(mean_alf1Mass, mean_ash1Mass, mean_bluestem1Mass, mean_oak1Mass, mean_pine1Mass)  
    new_mass_uq=c(uq_alf1Mass, uq_ash1Mass, uq_bluestem1Mass, uq_oak1Mass, uq_pine1Mass) 
    new_mass_lq=c(lq_alf1Mass, lq_ash1Mass, lq_bluestem1Mass, lq_oak1Mass, lq_pine1Mass) 
    list_dataMass=c(rowMeans(data_MCMC[[1]][[1]][,3:5]),
                    rowMeans(data_MCMC[[2]][[1]][,3:5]),
                    rowMeans(data_MCMC[[3]][[1]][,3:5]),
                    rowMeans(data_MCMC[[4]][[1]][,3:5]),
                    rowMeans(data_MCMC[[5]][[1]][,3:5]))

####################DOC Latent States########################
#load model outputs for each litter and measurement type, then summarize data
#create measured/modeled figures using model results from all chains
  #ALFALFA
    AlfDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_DOC_alf1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_DOC_alf2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_DOC_alf3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfDOC_all=mcmc(rbind(AlfDOC_chain1[,2:ncol(AlfDOC_chain1)], 
                          AlfDOC_chain2[,2:ncol(AlfDOC_chain2)], 
                          AlfDOC_chain3[,2:ncol(AlfDOC_chain3)]))
      alf1DOC=summary(AlfDOC_all)
      mean_alf1DOC=(alf1DOC[[1]][,1])
      uq_alf1DOC=(alf1DOC[[2]][,5])
      lq_alf1DOC=(alf1DOC[[2]][,1])

  #ASH
    AshDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_DOC_ash1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_DOC_ash2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_DOC_ash3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshDOC_all=mcmc(rbind(AshDOC_chain1[,2:ncol(AshDOC_chain1)], 
                          AshDOC_chain2[,2:ncol(AshDOC_chain2)], 
                          AshDOC_chain3[,2:ncol(AshDOC_chain3)]))
      ash1DOC=summary(AshDOC_all)
      mean_ash1DOC=(ash1DOC[[1]][,1])
      uq_ash1DOC=(ash1DOC[[2]][,5])
      lq_ash1DOC=(ash1DOC[[2]][,1])

  #BLUESTEM
    BluDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_DOC_blu1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_DOC_blu2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_DOC_blu3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluDOC_all=mcmc(rbind(BluDOC_chain1[,2:ncol(BluDOC_chain1)], 
                          BluDOC_chain2[,2:ncol(BluDOC_chain2)], 
                          BluDOC_chain3[,2:ncol(BluDOC_chain3)]))
      bluestem1DOC=summary(BluDOC_all)
      mean_bluestem1DOC=(bluestem1DOC[[1]][,1])
      uq_bluestem1DOC=(bluestem1DOC[[2]][,5])
      lq_bluestem1DOC=(bluestem1DOC[[2]][,1])

  #OAK
    OakDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_DOC_oak1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_DOC_oak2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_DOC_oak3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakDOC_all=mcmc(rbind(OakDOC_chain1[,2:ncol(OakDOC_chain1)], 
                          OakDOC_chain2[,2:ncol(OakDOC_chain2)], 
                          OakDOC_chain3[,2:ncol(OakDOC_chain3)]))
      oak1DOC=summary(OakDOC_all)
      mean_oak1DOC=(oak1DOC[[1]][,1])
      uq_oak1DOC=(oak1DOC[[2]][,5])
      lq_oak1DOC=(oak1DOC[[2]][,1])

  #PINE
    PinDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_DOC_pin1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_DOC_pin2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_DOC_pin3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinDOC_all=mcmc(rbind(PinDOC_chain1[,2:ncol(PinDOC_chain1)], 
                          PinDOC_chain2[,2:ncol(PinDOC_chain2)], 
                          PinDOC_chain3[,2:ncol(PinDOC_chain3)]))
      pine1DOC=summary(PinDOC_all)
      mean_pine1DOC=(pine1DOC[[1]][,1])
      uq_pine1DOC=(pine1DOC[[2]][,5])
      lq_pine1DOC=(pine1DOC[[2]][,1])


  #create lists of results to use for plotting
    summary_litterDOC=list(mean_alf1DOC, mean_ash1DOC, mean_bluestem1DOC, mean_oak1DOC, mean_pine1DOC)
    new_DOC=c(mean_alf1DOC, mean_ash1DOC, mean_bluestem1DOC, mean_oak1DOC, mean_pine1DOC)  
    new_DOC_uq=c(uq_alf1DOC, uq_ash1DOC, uq_bluestem1DOC, uq_oak1DOC, uq_pine1DOC) +750
    new_DOC_lq=c(lq_alf1DOC, lq_ash1DOC, lq_bluestem1DOC, lq_oak1DOC, lq_pine1DOC) -750
    list_dataDOC=c(rowMeans(data_MCMC[[1]][[2]][,3:5]),
                   rowMeans(data_MCMC[[2]][[2]][,3:5]),
                   rowMeans(data_MCMC[[3]][[2]][,3:5]),
                   rowMeans(data_MCMC[[4]][[2]][,3:5]),
                   rowMeans(data_MCMC[[5]][[2]][,3:5]))  

####################CO2 Latent States########################
#load model outputs for each litter and measurement type, then summarize data
#create measured/modeled figures using model results from all chains
  #ALFALFA
    AlfCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_CO2_alf1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_CO2_alf2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_CO2_alf3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AlfCO2_all=mcmc(rbind(AlfCO2_chain1[,2:ncol(AlfCO2_chain1)], 
                          AlfCO2_chain2[,2:ncol(AlfCO2_chain2)], 
                          AlfCO2_chain3[,2:ncol(AlfCO2_chain3)]))
      alf1CO2=summary(AlfCO2_all)
      mean_alf1CO2=(alf1CO2[[1]][,1])
      uq_alf1CO2=(alf1CO2[[2]][,5])
      lq_alf1CO2=(alf1CO2[[2]][,1])

  #ASH
    AshCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_CO2_ash1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_CO2_ash2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_CO2_ash3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    AshCO2_all=mcmc(rbind(AshCO2_chain1[,2:ncol(AshCO2_chain1)], 
                          AshCO2_chain2[,2:ncol(AshCO2_chain2)], 
                          AshCO2_chain3[,2:ncol(AshCO2_chain3)]))
      ash1CO2=summary(AshCO2_all)
      mean_ash1CO2=(ash1CO2[[1]][,1])
      uq_ash1CO2=(ash1CO2[[2]][,5])
      lq_ash1CO2=(ash1CO2[[2]][,1])

  #BLUESTEM
    BluCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_CO2_blu1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_CO2_blu2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_CO2_blu3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    BluCO2_all=mcmc(rbind(BluCO2_chain1[,2:ncol(BluCO2_chain1)], 
                          BluCO2_chain2[,2:ncol(BluCO2_chain2)], 
                          BluCO2_chain3[,2:ncol(BluCO2_chain3)]))
      bluestem1CO2=summary(BluCO2_all)
      mean_bluestem1CO2=(bluestem1CO2[[1]][,1])
      uq_bluestem1CO2=(bluestem1CO2[[2]][,5])
      lq_bluestem1CO2=(bluestem1CO2[[2]][,1])

  #OAK
    OakCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_CO2_oak1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_CO2_oak2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_CO2_oak3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OakCO2_all=mcmc(rbind(OakCO2_chain1[,2:ncol(OakCO2_chain1)], 
                          OakCO2_chain2[,2:ncol(OakCO2_chain2)], 
                          OakCO2_chain3[,2:ncol(OakCO2_chain3)]))
      oak1CO2=summary(OakCO2_all)
      mean_oak1CO2=(oak1CO2[[1]][,1])
      uq_oak1CO2=(oak1CO2[[2]][,5])
      lq_oak1CO2=(oak1CO2[[2]][,1])

  #PINE
    PinCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/CVmod_CO2_pin1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/CVmod_CO2_pin2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/CVmod_CO2_pin3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    PinCO2_all=mcmc(rbind(PinCO2_chain1[,2:ncol(PinCO2_chain1)], 
                          PinCO2_chain2[,2:ncol(PinCO2_chain2)], 
                          PinCO2_chain3[,2:ncol(PinCO2_chain3)]))
      pine1CO2=summary(PinCO2_all)
      mean_pine1CO2=(pine1CO2[[1]][,1])
      uq_pine1CO2=(pine1CO2[[2]][,5])
      lq_pine1CO2=(pine1CO2[[2]][,1])

  #create lists of results to use for plotting
    summary_litterCO2=list(mean_alf1CO2, mean_ash1CO2, mean_bluestem1CO2, mean_oak1CO2, mean_pine1CO2)
    new_CO2=c(mean_alf1CO2, mean_ash1CO2, mean_bluestem1CO2, mean_oak1CO2, mean_pine1CO2)  
    new_CO2_uq=c(uq_alf1CO2, uq_ash1CO2, uq_bluestem1CO2, uq_oak1CO2, uq_pine1CO2) 
    new_CO2_lq=c(lq_alf1CO2, lq_ash1CO2, lq_bluestem1CO2, lq_oak1CO2, lq_pine1CO2) 
    list_dataCO2=c(rowMeans(data_MCMC[[1]][[3]][,3:5]),
                   rowMeans(data_MCMC[[2]][[3]][,3:5]),
                   rowMeans(data_MCMC[[3]][[3]][,3:5]),
                   rowMeans(data_MCMC[[4]][[3]][,3:5]),
                   rowMeans(data_MCMC[[5]][[3]][,3:5]))   

###################print Latent State Result Figures#####################
  #MASS
      jpeg(file="C:/LIDEL/FINAL_results/Model1_Mass.jpg")
      par(mar=c(5, 5, 3, 3))
    #plot measured mean/modeled results
      plotCI(list_dataMass/1000, new_mass/1000, err="y", 
             ui=new_mass_uq/1000, li=new_mass_lq/1000, 
             typ="p", pch=16, cex.lab=1.65,
             xlab=expression(paste("Measured Mass Loss (mg C ", Delta, t^-1, ")")),
             ylab=expression(paste("Modeled Mass Loss (mg C ", Delta, t^-1, ")")),
             bty="n",cex.main=2.5,
             cex.axis=1.65,xlim=c(0, 1000), ylim=c(0, 1000))
    #plot all measured data
#       for(i in 1:num_litter){
#         for(s in 3:ncol(all_data[[i]][[1]])){
#           points(data_MCMC[[i]][[1]][,s], summary_litterMass[[i]], col="black", 
#                  typ="p", pch=4)
#         }
#       }
      abline(0,1)
      dev.off()

  #DOC
      jpeg(file="C:/LIDEL/FINAL_results/Model1_DOC.jpg")
      par(mar=c(5, 5, 3, 3))
    #plot measured mean/modeled results
      plotCI(list_dataDOC/1000, new_DOC/1000, err="y", 
             ui=new_DOC_uq/1000, li=new_DOC_lq/1000, 
             typ="p", pch=16, bty="n", cex.lab=1.65, gap=0,
             xlab=expression(paste("Measured DOC (mg C ", Delta, t^-1, ")")),
             ylab=expression(paste("Modeled DOC (mg C ", Delta, t^-1, ")")),
             cex.axis=1.65,xlim=c(0, 25), ylim=c(0, 25), bty="n")
    #plot all measured data
#       for(i in 1:num_litter){
#         for(s in 3:ncol(all_data[[i]][[1]])){
#           points(data_MCMC[[i]][[2]][,s], summary_litterDOC[[i]], col="black", 
#                  typ="p", pch=4)
#         }
#       }
      abline(0,1)
      dev.off()

  #CO2
      jpeg(file="C:/LIDEL/FINAL_results/Model1_CO2.jpg")
      par(mar=c(5, 5, 3, 3))
    #plot measured mean/modeled results
      plotCI(list_dataCO2/1000, new_CO2/1000, err="y", 
             ui=new_CO2_uq/1000, li=new_CO2_lq/1000, 
             typ="p", bty="n",cex.lab=1.65, gap=0,
             xlab=expression(paste("Measured CO"[2]," (mg C ", Delta, t^-1, ")")),
             ylab=expression(paste("Modeled CO"[2]," (mg C ", Delta, t^-1, ")")),
             cex.axis=1.65,pch=16, xlim=c(0, 70), ylim=c(0, 70))
    #plot all measured data
#       for(i in 1:num_litter){
#         for(s in 3:ncol(all_data[[i]][[1]])){
#           points(data_MCMC[[i]][[3]][,s], summary_litterCO2[[i]], col="black", 
#                  typ="p", pch=4)
#         }
#       }
      abline(0,1)
      dev.off()

####################DOC Out-of-Sample Estimates########################
#load model outputs for each litter and measurement type, then summarize data
#create measured/modeled figures for out-of-sample measurements using model results from all chains
  #ALFALFA
    OOSAlfDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_DOC_alf.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSAlfDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_DOC_alf.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSAlfDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_DOC_alf.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSAlfDOC_all=mcmc(rbind(OOSAlfDOC_chain1[,2:ncol(OOSAlfDOC_chain1)], 
                             OOSAlfDOC_chain2[,2:ncol(OOSAlfDOC_chain2)], 
                             OOSAlfDOC_chain3[,2:ncol(OOSAlfDOC_chain3)]))
    OOSalf1DOC=summary(OOSAlfDOC_all)
    OOSmean_alf1DOC=(OOSalf1DOC[[1]][,1])
    OOSuq_alf1DOC=(OOSalf1DOC[[2]][,5])
    OOSlq_alf1DOC=(OOSalf1DOC[[2]][,1])

  #ASH
    OOSAshDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_DOC_ash.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSAshDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_DOC_ash.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSAshDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_DOC_ash.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSAshDOC_all=mcmc(rbind(OOSAshDOC_chain1[,2:ncol(OOSAshDOC_chain1)], 
                             OOSAshDOC_chain2[,2:ncol(OOSAshDOC_chain2)], 
                          OOSAshDOC_chain3[,2:ncol(OOSAshDOC_chain3)]))
    OOSash1DOC=summary(OOSAshDOC_all)
    OOSmean_ash1DOC=(OOSash1DOC[[1]][,1])
    OOSuq_ash1DOC=(OOSash1DOC[[2]][,5])
    OOSlq_ash1DOC=(OOSash1DOC[[2]][,1])

  #BLUESTEM
    OOSBluDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_DOC_blu.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSBluDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_DOC_blu.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSBluDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_DOC_blu.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSBluDOC_all=mcmc(rbind(OOSBluDOC_chain1[,2:ncol(OOSBluDOC_chain1)], 
                             OOSBluDOC_chain2[,2:ncol(OOSBluDOC_chain2)], 
                          OOSBluDOC_chain3[,2:ncol(OOSBluDOC_chain3)]))
    OOSbluestem1DOC=summary(OOSBluDOC_all)
    OOSmean_bluestem1DOC=(OOSbluestem1DOC[[1]][,1])
    OOSuq_bluestem1DOC=(OOSbluestem1DOC[[2]][,5])
    OOSlq_bluestem1DOC=(OOSbluestem1DOC[[2]][,1])

  #OAK
    OOSOakDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_DOC_oak.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSOakDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_DOC_oak.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSOakDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_DOC_oak.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSOakDOC_all=mcmc(rbind(OOSOakDOC_chain1[,2:ncol(OOSOakDOC_chain1)], 
                             OOSOakDOC_chain2[,2:ncol(OOSOakDOC_chain2)], 
                          OOSOakDOC_chain3[,2:ncol(OOSOakDOC_chain3)]))
    OOSoak1DOC=summary(OOSOakDOC_all)
    OOSmean_oak1DOC=(OOSoak1DOC[[1]][,1])
    OOSuq_oak1DOC=(OOSoak1DOC[[2]][,5])
    OOSlq_oak1DOC=(OOSoak1DOC[[2]][,1])

  #PINE
    OOSPinDOC_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_DOC_pin.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSPinDOC_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_DOC_pin.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSPinDOC_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_DOC_pin.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOSPinDOC_all=mcmc(rbind(OOSPinDOC_chain1[,2:ncol(OOSPinDOC_chain1)], 
                             OOSPinDOC_chain2[,2:ncol(OOSPinDOC_chain2)], 
                          OOSPinDOC_chain3[,2:ncol(OOSPinDOC_chain3)]))
    OOSpine1DOC=summary(OOSPinDOC_all)
    OOSmean_pine1DOC=(OOSpine1DOC[[1]][,1])
    OOSuq_pine1DOC=(OOSpine1DOC[[2]][,5])
    OOSlq_pine1DOC=(OOSpine1DOC[[2]][,1])


  #create lists of results to use for plotting
    OOSsummary_litterDOC=list(OOSmean_alf1DOC, OOSmean_ash1DOC, OOSmean_bluestem1DOC, OOSmean_oak1DOC, OOSmean_pine1DOC)
    OOSnew_DOC=c(OOSmean_alf1DOC, OOSmean_ash1DOC, OOSmean_bluestem1DOC, OOSmean_oak1DOC, OOSmean_pine1DOC)  
    OOSnew_DOC_uq=c(OOSuq_alf1DOC, OOSuq_ash1DOC, OOSuq_bluestem1DOC, OOSuq_oak1DOC, OOSuq_pine1DOC) 
    OOSnew_DOC_lq=c(OOSlq_alf1DOC, OOSlq_ash1DOC, OOSlq_bluestem1DOC, OOSlq_oak1DOC, OOSlq_pine1DOC)
    OOSlist_dataDOC=c(rowMeans(data_nonMCMC[[1]][[1]][,3:5]),
                   rowMeans(data_nonMCMC[[2]][[1]][,3:5]),
                   rowMeans(data_nonMCMC[[3]][[1]][,3:5]),
                   rowMeans(data_nonMCMC[[4]][[1]][,3:5]),
                   rowMeans(data_nonMCMC[[5]][[1]][,3:5]))  

  #output figure for plot of DOC OOS loss 
      jpeg(file="C:/LIDEL/FINAL_results/Model1_DOC_OOS.jpg")
      par(mar=c(5, 5, 3, 3))
    #plot measured mean/modeled results
      plotCI(OOSlist_dataDOC/1000, OOSnew_DOC/1000, err="y", 
             ui=OOSnew_DOC_uq/1000, li=OOSnew_DOC_lq/1000, 
             xlab=expression(paste("Measured DOC (mg C ", Delta, t^-1, ")")),
             ylab=expression(paste("Modeled DOC (mg C ", Delta, t^-1, ")")),
             typ="p", pch=16, bty="n", cex.lab=1.65,
             cex.axis=1.65,xlim=c(0, 15), ylim=c(0, 15), bty="n", gap=0)
      abline(0,1)
    #plot all measured data
#       for(i in 1:num_litter){
#         for(s in 3:ncol(all_data[[i]][[1]])){
#           points(data_nonMCMC[[i]][[1]][,s], OOSsummary_litterDOC[[i]], col="black", 
#                  typ="p", pch=4)
#         }
#       }

      dev.off()
      
####################CO2 Out-of-Sample Results########################
#load model outputs for each litter and measurement type, then summarize data
#create measured/modeled figures using model results from all chains
  #ALFALFA
    OOS_AlfCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_CO2_alf.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      OOS_AlfCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_CO2_alf.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      OOS_AlfCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_CO2_alf.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      OOS_AlfCO2_all=mcmc(rbind(OOS_AlfCO2_chain1[,2:ncol(OOS_AlfCO2_chain1)], 
                                OOS_AlfCO2_chain2[,2:ncol(OOS_AlfCO2_chain2)], 
                                OOS_AlfCO2_chain3[,2:ncol(OOS_AlfCO2_chain3)]))
      OOS_alf1CO2=summary(OOS_AlfCO2_all)
      OOS_mean_alf1CO2=(OOS_alf1CO2[[1]][,1])
      OOS_uq_alf1CO2=(OOS_alf1CO2[[2]][,5])
      OOS_lq_alf1CO2=(OOS_alf1CO2[[2]][,1])
      
  #ASH
      OOS_AshCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_CO2_ash.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      OOS_AshCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_CO2_ash.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      OOS_AshCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_CO2_ash.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      OOS_AshCO2_all=mcmc(rbind(OOS_AshCO2_chain1[,2:ncol(OOS_AshCO2_chain1)], 
                                OOS_AshCO2_chain2[,2:ncol(OOS_AshCO2_chain2)], 
                                OOS_AshCO2_chain3[,2:ncol(OOS_AshCO2_chain3)]))
      OOS_ash1CO2=summary(OOS_AshCO2_all)
    OOS_mean_ash1CO2=(OOS_ash1CO2[[1]][,1])
    OOS_uq_ash1CO2=(OOS_ash1CO2[[2]][,5])
    OOS_lq_ash1CO2=(OOS_ash1CO2[[2]][,1])
      
  #BLUESTEM
    OOS_BluCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_CO2_blu.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_BluCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_CO2_blu.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_BluCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_CO2_blu.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_BluCO2_all=mcmc(rbind(OOS_BluCO2_chain1[,2:ncol(OOS_BluCO2_chain1)], 
                              OOS_BluCO2_chain2[,2:ncol(OOS_BluCO2_chain2)], 
                              OOS_BluCO2_chain3[,2:ncol(OOS_BluCO2_chain3)]))
    OOS_bluestem1CO2=summary(OOS_BluCO2_all)
    OOS_mean_bluestem1CO2=(OOS_bluestem1CO2[[1]][,1])
    OOS_uq_bluestem1CO2=(OOS_bluestem1CO2[[2]][,5])
    OOS_lq_bluestem1CO2=(OOS_bluestem1CO2[[2]][,1])
      
  #OAK
    OOS_OakCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_CO2_oak.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_OakCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_CO2_oak.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_OakCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_CO2_oak.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_OakCO2_all=mcmc(rbind(OOS_OakCO2_chain1[,2:ncol(OOS_OakCO2_chain1)], 
                                  OOS_OakCO2_chain2[,2:ncol(OOS_OakCO2_chain2)], 
                                  OOS_OakCO2_chain3[,2:ncol(OOS_OakCO2_chain3)]))
    OOS_oak1CO2=summary(OOS_OakCO2_all)
    OOS_mean_oak1CO2=(OOS_oak1CO2[[1]][,1])
    OOS_uq_oak1CO2=(OOS_oak1CO2[[2]][,5])
    OOS_lq_oak1CO2=(OOS_oak1CO2[[2]][,1])
      
  #PINE
    OOS_PinCO2_chain1=read.csv("C:/LIDEL/Model1/Chain1/results/OOSchain1_CO2_pin.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_PinCO2_chain2=read.csv("C:/LIDEL/Model1/Chain2/results/OOSchain2_CO2_pin.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_PinCO2_chain3=read.csv("C:/LIDEL/Model1/Chain3/results/OOSchain3_CO2_pin.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    OOS_PinCO2_all=mcmc(rbind(OOS_PinCO2_chain1[,2:ncol(OOS_PinCO2_chain1)], 
                              OOS_PinCO2_chain2[,2:ncol(OOS_PinCO2_chain2)], 
                              OOS_PinCO2_chain3[,2:ncol(OOS_PinCO2_chain3)]))
    OOS_pine1CO2=summary(OOS_PinCO2_all)
    OOS_mean_pine1CO2=(OOS_pine1CO2[[1]][,1])
    OOS_uq_pine1CO2=(OOS_pine1CO2[[2]][,5])
    OOS_lq_pine1CO2=(OOS_pine1CO2[[2]][,1])
      
  #create lists of results to use for plotting
    OOS_summary_litterCO2=list(OOS_mean_alf1CO2, OOS_mean_ash1CO2, OOS_mean_bluestem1CO2, OOS_mean_oak1CO2, OOS_mean_pine1CO2)
    OOS_new_CO2=c(OOS_mean_alf1CO2, OOS_mean_ash1CO2, OOS_mean_bluestem1CO2, OOS_mean_oak1CO2, OOS_mean_pine1CO2)  
    OOS_new_CO2_uq=c(OOS_uq_alf1CO2, OOS_uq_ash1CO2, OOS_uq_bluestem1CO2, OOS_uq_oak1CO2, OOS_uq_pine1CO2) 
    OOS_new_CO2_lq=c(OOS_lq_alf1CO2, OOS_lq_ash1CO2, OOS_lq_bluestem1CO2, OOS_lq_oak1CO2, OOS_lq_pine1CO2) 
    OOS_list_dataCO2=c(rowMeans(data_nonMCMC[[1]][[2]][,3:5]),
                   rowMeans(data_nonMCMC[[2]][[2]][,3:5]),
                   rowMeans(data_nonMCMC[[3]][[2]][,3:5]),
                   rowMeans(data_nonMCMC[[4]][[2]][,3:5]),
                   rowMeans(data_nonMCMC[[5]][[2]][,3:5]))  
    
  #Plot for figure of CO2 OOS results
      jpeg(file="C:/LIDEL/FINAL_results/Model1_CO2_OOS.jpg")
      par(mar=c(5, 5, 3, 3))
    #plot measured mean/modeled results
      plotCI(OOS_list_dataCO2/1000, OOS_new_CO2/1000, err="y", 
             ui=OOS_new_CO2_uq/1000, li=OOS_new_CO2_lq/1000, 
             typ="p", bty="n",cex.lab=1.65, gap=0,
             xlab=expression(paste("Measured CO"[2]," (mg C ", Delta, t^-1, ")")),
             ylab=expression(paste("Modeled CO"[2]," (mg C ", Delta, t^-1, ")")),
             cex.axis=1.65,pch=16, xlim=c(0, 100), ylim=c(0, 100))
    #plot all measured data
#       for(i in 1:num_litter){
#         for(s in 3:ncol(all_data[[i]][[1]])){
#           points(data_nonMCMC[[i]][[2]][,s], OOS_summary_litterCO2[[i]], col="black", 
#                  typ="p", pch=4)
#         }
#       }
      abline(0,1)
      dev.off()
    