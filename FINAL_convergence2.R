setwd("C:/Users/eecampbe/Desktop/trial/model_2_2")
library(deSolve)
library(gplots)
library(pscl)
library(rootSolve)
library(coda)
library(compiler)

#V2_1 parameter chains
Param_chain1=read.csv("RESULTchain1_params.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_param_chain1=mcmc(Param_chain1[,2:ncol(Param_chain1)])
Param_chain2=read.csv("RESULTchain2_params.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_param_chain2=mcmc(Param_chain2[,2:ncol(Param_chain2)])
Param_chain3=read.csv("RESULTchain3_params.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_param_chain3=mcmc(Param_chain3[,2:ncol(Param_chain3)])

#V2_1 variance chains
Var_chain1=read.csv("RESULTchain1_var.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_Var_chain1=mcmc(Var_chain1[,2:ncol(Var_chain1)])
Var_chain2=read.csv("RESULTchain2_var.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_Var_chain2=mcmc(Var_chain2[,2:ncol(Var_chain2)])
Var_chain3=read.csv("RESULTchain3_var.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_Var_chain3=mcmc(Var_chain3[,2:ncol(Var_chain3)])

#V2_1 latent chains
Latent_chain1=read.csv("RESULTchain1_latent.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_Latent_chain1=mcmc(Latent_chain1[,2:ncol(Latent_chain1)])
Latent_chain2=read.csv("RESULTchain2_latent.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_Latent_chain2=mcmc(Latent_chain2[,2:ncol(Latent_chain2)])
Latent_chain3=read.csv("RESULTchain3_latent.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
  MC_Latent_chain3=mcmc(Latent_chain3[,2:ncol(Latent_chain3)])


#if thinning is preferred: thin=seq(from=burnin, to=n.iter, by="thin"), then mcmc(<result>[thin,])

#create mcmc list (manual for now)
combined_param<-mcmc.list(MC_param_chain1, MC_param_chain2, MC_param_chain3)
combined_var<-mcmc.list(MC_Var_chain1, MC_Var_chain2, MC_Var_chain3)
combined_latent<-mcmc.list(MC_Latent_chain1, MC_Latent_chain2, MC_Latent_chain3)

sink("C:/Users/eecampbe/Desktop/trial/model_2_2_SUMMARY.txt")
  print(gelman.diag(combined_param))
  print(gelman.diag(combined_var))
  print(gelman.diag(combined_latent))
  
  
  print(summary(combined_param))
  print(summary(combined_var))
  print(summary(combined_latent))
sink()

#to do
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
#load CV model outputs
#create measured/modeled figures using all chain data for model results
##################MASS#######################
#ALFALFA
AlfMass_chain1=read.csv("CVmod_mass_alf1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfMass_chain2=read.csv("CVmod_mass_alf2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfMass_chain3=read.csv("CVmod_mass_alf3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfMass_all=mcmc(rbind(AlfMass_chain1[,2:3], AlfMass_chain2[,2:3], AlfMass_chain3[,2:3]))
  alf1Mass=summary(AlfMass_all)
  mean_alf1Mass=(alf1Mass[[1]][,1])
  uq_alf1Mass=(alf1Mass[[2]][,5])
  lq_alf1Mass=(alf1Mass[[2]][,1])

#ASH
AshMass_chain1=read.csv("CVmod_mass_ash1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshMass_chain2=read.csv("CVmod_mass_ash2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshMass_chain3=read.csv("CVmod_mass_ash3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshMass_all=mcmc(rbind(AshMass_chain1[,2:3], AshMass_chain2[,2:3], AshMass_chain3[,2:3]))
  ash1Mass=summary(AshMass_all)
  mean_ash1Mass=(ash1Mass[[1]][,1])
  uq_ash1Mass=(ash1Mass[[2]][,5])
  lq_ash1Mass=(ash1Mass[[2]][,1])

#BLUESTEM
BluMass_chain1=read.csv("CVmod_mass_blu1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluMass_chain2=read.csv("CVmod_mass_blu2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluMass_chain3=read.csv("CVmod_mass_blu3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluMass_all=mcmc(rbind(BluMass_chain1[,2:3], BluMass_chain2[,2:3], BluMass_chain3[,2:3]))
  bluestem1Mass=summary(BluMass_all)
  mean_bluestem1Mass=(bluestem1Mass[[1]][,1])
  uq_bluestem1Mass=(bluestem1Mass[[2]][,5])
  lq_bluestem1Mass=(bluestem1Mass[[2]][,1])

#OAK
OakMass_chain1=read.csv("CVmod_mass_oak1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakMass_chain2=read.csv("CVmod_mass_oak2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakMass_chain3=read.csv("CVmod_mass_oak3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakMass_all=mcmc(rbind(OakMass_chain1[,2:3], OakMass_chain2[,2:3], OakMass_chain3[,2:3]))
  oak1Mass=summary(OakMass_all)
  mean_oak1Mass=(oak1Mass[[1]][,1])
  uq_oak1Mass=(oak1Mass[[2]][,5])
  lq_oak1Mass=(oak1Mass[[2]][,1])

#PINE
PinMass_chain1=read.csv("CVmod_mass_pin1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinMass_chain2=read.csv("CVmod_mass_pin2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinMass_chain3=read.csv("CVmod_mass_pin3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinMass_all=mcmc(rbind(PinMass_chain1[,2:3], PinMass_chain2[,2:3], PinMass_chain3[,2:3]))
  pine1Mass=summary(PinMass_all)
  mean_pine1Mass=(pine1Mass[[1]][,1])
  uq_pine1Mass=(pine1Mass[[2]][,5])
  lq_pine1Mass=(pine1Mass[[2]][,1])


#for plotting
    summary_litterMass=list(mean_alf1Mass, mean_ash1Mass, mean_bluestem1Mass, mean_oak1Mass, mean_pine1Mass)
    
    new_mass=c(mean_alf1Mass, mean_ash1Mass, mean_bluestem1Mass, mean_oak1Mass, mean_pine1Mass)  
    new_mass_uq=c(uq_alf1Mass, uq_ash1Mass, uq_bluestem1Mass, uq_oak1Mass, uq_pine1Mass) 
    new_mass_lq=c(lq_alf1Mass, lq_ash1Mass, lq_bluestem1Mass, lq_oak1Mass, lq_pine1Mass) 


list_dataMass=c(rowMeans(data_MCMC[[1]][[1]][,3:5]),
                rowMeans(data_MCMC[[2]][[1]][,3:5]),
                rowMeans(data_MCMC[[3]][[1]][,3:5]),
                rowMeans(data_MCMC[[4]][[1]][,3:5]),
                rowMeans(data_MCMC[[5]][[1]][,3:5]))

par(mar=c(4,4,4,4))
plotCI(list_dataMass, new_mass, err="y", ui=new_mass_uq, li=new_mass_lq, 
       main="Mass Loss V2_2", ylab="Modeled total C loss", typ="p", pch=16,
       xlab="Measured total C loss", xlim=c(100000, 1000000), ylim=c(100000, 1000000))
#plot of litter C remaining through time
for(i in 1:num_litter){
  for(s in 3:ncol(all_data[[i]][[1]])){
    points(data_MCMC[[i]][[1]][,s], summary_litterMass[[i]], col="black", 
           typ="p", pch=4)
  }
}
abline(0,1)


####################DOC########################
#ALFALFA
AlfDOC_chain1=read.csv("CVmod_DOC_alf1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfDOC_chain2=read.csv("CVmod_DOC_alf2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfDOC_chain3=read.csv("CVmod_DOC_alf3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfDOC_all=mcmc(rbind(AlfDOC_chain1[,2:ncol(AlfDOC_chain1)], 
                      AlfDOC_chain2[,2:ncol(AlfDOC_chain2)], 
                      AlfDOC_chain3[,2:ncol(AlfDOC_chain3)]))
  alf1DOC=summary(AlfDOC_all)
  mean_alf1DOC=(alf1DOC[[1]][,1])
  uq_alf1DOC=(alf1DOC[[2]][,5])
  lq_alf1DOC=(alf1DOC[[2]][,1])

#ASH
AshDOC_chain1=read.csv("CVmod_DOC_ash1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshDOC_chain2=read.csv("CVmod_DOC_ash2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshDOC_chain3=read.csv("CVmod_DOC_ash3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshDOC_all=mcmc(rbind(AshDOC_chain1[,2:ncol(AshDOC_chain1)], 
                      AshDOC_chain2[,2:ncol(AshDOC_chain2)], 
                      AshDOC_chain3[,2:ncol(AshDOC_chain3)]))
  ash1DOC=summary(AshDOC_all)
  mean_ash1DOC=(ash1DOC[[1]][,1])
  uq_ash1DOC=(ash1DOC[[2]][,5])
  lq_ash1DOC=(ash1DOC[[2]][,1])

#BLUESTEM
BluDOC_chain1=read.csv("CVmod_DOC_blu1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluDOC_chain2=read.csv("CVmod_DOC_blu2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluDOC_chain3=read.csv("CVmod_DOC_blu3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluDOC_all=mcmc(rbind(BluDOC_chain1[,2:ncol(BluDOC_chain1)], 
                      BluDOC_chain2[,2:ncol(BluDOC_chain2)], 
                      BluDOC_chain3[,2:ncol(BluDOC_chain3)]))
  bluestem1DOC=summary(BluDOC_all)
  mean_bluestem1DOC=(bluestem1DOC[[1]][,1])
  uq_bluestem1DOC=(bluestem1DOC[[2]][,5])
  lq_bluestem1DOC=(bluestem1DOC[[2]][,1])

#OAK
OakDOC_chain1=read.csv("CVmod_DOC_oak1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakDOC_chain2=read.csv("CVmod_DOC_oak2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakDOC_chain3=read.csv("CVmod_DOC_oak3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakDOC_all=mcmc(rbind(OakDOC_chain1[,2:ncol(OakDOC_chain1)], 
                      OakDOC_chain2[,2:ncol(OakDOC_chain2)], 
                      OakDOC_chain3[,2:ncol(OakDOC_chain3)]))
  oak1DOC=summary(OakDOC_all)
  mean_oak1DOC=(oak1DOC[[1]][,1])
  uq_oak1DOC=(oak1DOC[[2]][,5])
  lq_oak1DOC=(oak1DOC[[2]][,1])

#PINE
PinDOC_chain1=read.csv("CVmod_DOC_pin1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinDOC_chain2=read.csv("CVmod_DOC_pin2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinDOC_chain3=read.csv("CVmod_DOC_pin3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinDOC_all=mcmc(rbind(PinDOC_chain1[,2:ncol(PinDOC_chain1)], 
                      PinDOC_chain2[,2:ncol(PinDOC_chain2)], 
                      PinDOC_chain3[,2:ncol(PinDOC_chain3)]))
  pine1DOC=summary(PinDOC_all)
  mean_pine1DOC=(pine1DOC[[1]][,1])
  uq_pine1DOC=(pine1DOC[[2]][,5])
  lq_pine1DOC=(pine1DOC[[2]][,1])


#for plotting
  summary_litterDOC=list(mean_alf1DOC, mean_ash1DOC, mean_bluestem1DOC, mean_oak1DOC, mean_pine1DOC)
  
  new_DOC=c(mean_alf1DOC, mean_ash1DOC, mean_bluestem1DOC, mean_oak1DOC, mean_pine1DOC)  
  new_DOC_uq=c(uq_alf1DOC, uq_ash1DOC, uq_bluestem1DOC, uq_oak1DOC, uq_pine1DOC) +750
  new_DOC_lq=c(lq_alf1DOC, lq_ash1DOC, lq_bluestem1DOC, lq_oak1DOC, lq_pine1DOC) -750

  list_dataDOC=c(rowMeans(data_MCMC[[1]][[2]][,3:5]),
                 rowMeans(data_MCMC[[2]][[2]][,3:5]),
                 rowMeans(data_MCMC[[3]][[2]][,3:5]),
                 rowMeans(data_MCMC[[4]][[2]][,3:5]),
                 rowMeans(data_MCMC[[5]][[2]][,3:5]))  


#plot DOC loss
par(mar=c(4,4,4,4))
plotCI(list_dataDOC, new_DOC, err="y", ui=new_DOC_uq, li=new_DOC_lq, 
       main="DOC V2_2", ylab="Modeled DOC", typ="p", pch=16, 
       xlab="Measured DOC", xlim=c(0, 25000), ylim=c(0, 25000), bty="n")
#plot of litter C remaining through time
for(i in 1:num_litter){
  for(s in 3:ncol(all_data[[i]][[1]])){
    points(data_MCMC[[i]][[2]][,s], summary_litterDOC[[i]], col="black", 
           typ="p", pch=4)
  }
}
abline(0,1)


####################DOC########################
#ALFALFA
AlfCO2_chain1=read.csv("CVmod_CO2_alf1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfCO2_chain2=read.csv("CVmod_CO2_alf2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfCO2_chain3=read.csv("CVmod_CO2_alf3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AlfCO2_all=mcmc(rbind(AlfCO2_chain1[,2:ncol(AlfCO2_chain1)], 
                      AlfCO2_chain2[,2:ncol(AlfCO2_chain2)], 
                      AlfCO2_chain3[,2:ncol(AlfCO2_chain3)]))
  alf1CO2=summary(AlfCO2_all)
  mean_alf1CO2=(alf1CO2[[1]][,1])
  uq_alf1CO2=(alf1CO2[[2]][,5])
  lq_alf1CO2=(alf1CO2[[2]][,1])

#ASH
AshCO2_chain1=read.csv("CVmod_CO2_ash1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshCO2_chain2=read.csv("CVmod_CO2_ash2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshCO2_chain3=read.csv("CVmod_CO2_ash3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
AshCO2_all=mcmc(rbind(AshCO2_chain1[,2:ncol(AshCO2_chain1)], 
                      AshCO2_chain2[,2:ncol(AshCO2_chain2)], 
                      AshCO2_chain3[,2:ncol(AshCO2_chain3)]))
  ash1CO2=summary(AshCO2_all)
  mean_ash1CO2=(ash1CO2[[1]][,1])
  uq_ash1CO2=(ash1CO2[[2]][,5])
  lq_ash1CO2=(ash1CO2[[2]][,1])

#BLUESTEM
BluCO2_chain1=read.csv("CVmod_CO2_blu1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluCO2_chain2=read.csv("CVmod_CO2_blu2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluCO2_chain3=read.csv("CVmod_CO2_blu3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
BluCO2_all=mcmc(rbind(BluCO2_chain1[,2:ncol(BluCO2_chain1)], 
                      BluCO2_chain2[,2:ncol(BluCO2_chain2)], 
                      BluCO2_chain3[,2:ncol(BluCO2_chain3)]))
  bluestem1CO2=summary(BluCO2_all)
  mean_bluestem1CO2=(bluestem1CO2[[1]][,1])
  uq_bluestem1CO2=(bluestem1CO2[[2]][,5])
  lq_bluestem1CO2=(bluestem1CO2[[2]][,1])

#OAK
OakCO2_chain1=read.csv("CVmod_CO2_oak1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakCO2_chain2=read.csv("CVmod_CO2_oak2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakCO2_chain3=read.csv("CVmod_CO2_oak3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
OakCO2_all=mcmc(rbind(OakCO2_chain1[,2:ncol(OakCO2_chain1)], 
                      OakCO2_chain2[,2:ncol(OakCO2_chain2)], 
                      OakCO2_chain3[,2:ncol(OakCO2_chain3)]))
  oak1CO2=summary(OakCO2_all)
  mean_oak1CO2=(oak1CO2[[1]][,1])
  uq_oak1CO2=(oak1CO2[[2]][,5])
  lq_oak1CO2=(oak1CO2[[2]][,1])

#PINE
PinCO2_chain1=read.csv("CVmod_CO2_pin1.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinCO2_chain2=read.csv("CVmod_CO2_pin2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinCO2_chain3=read.csv("CVmod_CO2_pin3.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
PinCO2_all=mcmc(rbind(PinCO2_chain1[,2:ncol(PinCO2_chain1)], 
                      PinCO2_chain2[,2:ncol(PinCO2_chain2)], 
                      PinCO2_chain3[,2:ncol(PinCO2_chain3)]))
  pine1CO2=summary(PinCO2_all)
  mean_pine1CO2=(pine1CO2[[1]][,1])
  uq_pine1CO2=(pine1CO2[[2]][,5])
  lq_pine1CO2=(pine1CO2[[2]][,1])

#for plotting
  summary_litterCO2=list(mean_alf1CO2, mean_ash1CO2, mean_bluestem1CO2, mean_oak1CO2, mean_pine1CO2)
  
  new_CO2=c(mean_alf1CO2, mean_ash1CO2, mean_bluestem1CO2, mean_oak1CO2, mean_pine1CO2)  
  new_CO2_uq=c(uq_alf1CO2, uq_ash1CO2, uq_bluestem1CO2, uq_oak1CO2, uq_pine1CO2) 
  new_CO2_lq=c(lq_alf1CO2, lq_ash1CO2, lq_bluestem1CO2, lq_oak1CO2, lq_pine1CO2) 
  
  
  list_dataCO2=c(rowMeans(data_MCMC[[1]][[3]][,3:5]),
                 rowMeans(data_MCMC[[2]][[3]][,3:5]),
                 rowMeans(data_MCMC[[3]][[3]][,3:5]),
                 rowMeans(data_MCMC[[4]][[3]][,3:5]),
                 rowMeans(data_MCMC[[5]][[3]][,3:5]))   



#plot CO2 loss
plotCI(list_dataCO2, new_CO2, err="y", ui=new_CO2_uq, li=new_CO2_lq, 
       main="CO2 V2_2", ylab="Modeled CO2", typ="p", bty="n",
       xlab="Measured CO2", pch=16, xlim=c(0, 70000), ylim=c(0, 70000))
#plot of litter C remaining through time
for(i in 1:num_litter){
  for(s in 3:ncol(all_data[[i]][[1]])){
    points(data_MCMC[[i]][[3]][,s], summary_litterCO2[[i]], col="black", 
           typ="p", pch=4)
  }
}
abline(0,1)


###################print results figures#####################
#MASS
jpeg(file="C:/Users/eecampbe/Desktop/trial/Mass_2_2.jpg")
par(mar=c(2.85, 2.85, 2.85, 2.85))
plotCI(list_dataMass, new_mass, err="y", ui=new_mass_uq, li=new_mass_lq, 
       typ="p", pch=16, bty="n",cex.main=2.5,
       cex.axis=1.85, xlim=c(100000, 1000000), ylim=c(100000, 1000000))
#plot of litter C remaining through time
for(i in 1:num_litter){
  for(s in 3:ncol(all_data[[i]][[1]])){
    points(data_MCMC[[i]][[1]][,s], summary_litterMass[[i]], col="black", 
           typ="p", pch=4)
  }
}
abline(0,1)
dev.off()

#DOC
#plot DOC loss
jpeg(file="C:/Users/eecampbe/Desktop/trial/DOC_2_2.jpg")
par(mar=c(2.85, 2.85, 2.85, 2.85))
plotCI(list_dataDOC, new_DOC, err="y", ui=new_DOC_uq, li=new_DOC_lq, 
       typ="p", pch=16, bty="n",cex.main=2.5,
       cex.axis=1.85, xlim=c(0, 25000), ylim=c(0, 25000), bty="n")
#plot of litter C remaining through time
for(i in 1:num_litter){
  for(s in 3:ncol(all_data[[i]][[1]])){
    points(data_MCMC[[i]][[2]][,s], summary_litterDOC[[i]], col="black", 
           typ="p", pch=4)
  }
}
abline(0,1)
dev.off()

#CO2
#plot CO2 loss
jpeg(file="C:/Users/eecampbe/Desktop/trial/CO2_2_2.jpg")
par(mar=c(2.85, 2.85, 2.85, 2.85))
plotCI(list_dataCO2, new_CO2, err="y", ui=new_CO2_uq, li=new_CO2_lq, 
       typ="p", bty="n",cex.main=2.5,
       cex.axis=1.85, pch=16, xlim=c(0, 65000), ylim=c(0, 65000))
#plot of litter C remaining through time
for(i in 1:num_litter){
  for(s in 3:ncol(all_data[[i]][[1]])){
    points(data_MCMC[[i]][[3]][,s], summary_litterCO2[[i]], col="black", 
           typ="p", pch=4)
  }
}
abline(0,1)
dev.off()