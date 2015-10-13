#Initial DOC pool estimate
##########NOTES##########
#This script will run the full LIDEL model
#using CIinit analysis results
#and then plotting mean and CIs of the results
#Nell Campbell
#4/12/15
##########NOTES##########

rm(list=ls())
setwd("C:/LIDEL/Model4/Chain1")
  library(deSolve)
  library(gplots)
  library(pscl)
  library(rootSolve)
  library(coda)
  library(compiler)
  

  ############START CODE##########################
  #load model
    source("LIDEL_v2_2.R")
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
  #compile functions
    enableJIT(2)

  n.iter=80000
  burnin=40000
  #use when playing with input parameters to re-run inital measured/modeled

  source("LIDEL_set_params.R")
  source("LIDEL_MCMC_setup_trial.R")

# # #for MCMC inner work
# set=set_LIDEL1
# LIDEL_INPUT=list(LIDEL_inputs)
# data_y=all_data
# chain=1
# g=1
# s=2
  ##########RUN MCMC for all data##########
      x_Cinit=Run_MCMC(set=set_LIDEL1, d_l_list=d_l_list, Init_fctn_vals=Init_fctn_vals,
                 LIDEL_INPUT=list(LIDEL_inputs), litter_nam=litter_nam,
                 data_y=all_data, data_MCMC=data_MCMC, n.iter=n.iter, params=params,  
                 param_names=param_names, Init_Cs=Init_Cs, burnin=burnin, chain=1)

save.image(file="./results/trialrun_workspace.RData")

#post-processing code to output results
	source("LIDEL_post_processing.R")