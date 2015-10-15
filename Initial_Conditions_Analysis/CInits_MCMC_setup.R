##########SETUP FOR MCMC##########
#using all data
#Nell Campbell
#moved to separate file 11/17/14
##################################

#####iterations and chains
#Set up the number of iterations. While you are getting things working, this can be small, ca 5000

n.chain=3

#####specify inputs to model
#name all inputs, including parameters and derived values
param_names=c("Fdoc","sigma_L_C6","sigma_L_C2","sigma_HWE", "sigma_DIF", "sigma_LIG",
              "Fs", "Flig", 
              "zC1_3alp","sigma_y_C1_3", "sigma_y_C6alph", "sigma_y_C2alph", 
              "LC2alph", "LC6alph")
#specify the dimension of all sites (either number of sites- 1 currently-, or sites and time, ie scalar or vector)
dim_param=c(1,1,1,1,1,1,
            length(litter_nam), length(litter_nam),
            length(litter_nam), length(litter_nam), length(litter_nam), length(litter_nam),
            length(litter_nam),length(litter_nam))
#names of parameters constant across litter types and time (KEEP ORDER THE SAME THROUGHOUT)
params=c("Fdoc", "sigma_L_C6","sigma_L_C2","sigma_HWE", "sigma_DIF", "sigma_LIG")
#specify how the parameters will be 'tuned' during the runs
tune_param=c(0.05, 1500, 2000, .05, .05, .05,
             0.05, .02, 
             1000, 500, 500, 2500,
             5000, 2000)

tune_param_HWE=c(0.08, 500, 500, .05, .05, .05,
                 0.09, .1, 
                 1000, 50, 50, 50,
                 50, 50)

tune_param_ydiff=c(0.08, 500, 500, .05, .05, .05,
                   0.09, .1, 
                   1000, 50, 50, 50,
                   50, 50)
#what sort of data are the parameters?
support_param=c("zero-one", "non-negative","non-negative","zero-one","zero-one", "zero-one",
                "zero-one", "zero-one", 
                "real", "non-negative","non-negative", "non-negative",
                "real", "real")
#create matrix for to initialize chains
parameter_inits=matrix(NA, nrow=n.chain, ncol=length(params))
#set parameter values to initialize each chain
Fdoc=c(runif(1,.1,.5), runif(1,.1,.25),runif(1,.25, .5))
sigma_L_C6=c(runif(1,100,500), runif(1,100,500), runif(1,100,500))
sigma_L_C2=c(runif(1,100,500), runif(1,100,500), runif(1,100,500))
sigma_HWE=c(runif(1,.05,.1), runif(1,.05,.2),runif(1,.05, .05))
sigma_DIF=c(runif(1,.05,.1), runif(1,.05,.2),runif(1,.05, .05))
sigma_LIG=c(runif(1,.05,.1), runif(1,.05,.2),runif(1,.05, .05))

#loop to set init values within each parameter_inits
for(j in 1:n.chain){
  parameter_inits[j,1]=Fdoc[j]
  parameter_inits[j,2]=sigma_L_C6[j]
  parameter_inits[j,3]=sigma_L_C2[j]
  parameter_inits[j,4]=sigma_HWE[j]
  parameter_inits[j,5]=sigma_DIF[j]
  parameter_inits[j,6]=sigma_LIG[j]
}	

#parameters that vary by litter type
  #generate random values for initial conditions
    Fs=matrix(runif(15, 0.1, 0.45), nrow=length(litter_nam), ncol=n.chain)
    Flig=matrix(runif(15, 0.1, 0.45), nrow=length(litter_nam), ncol=n.chain)
    zC1_3alp=matrix(runif(15, 800000, 1100000), nrow=length(litter_nam), ncol=n.chain)
    sigma_y_C1_3=matrix(runif(15, 20000, 55000), nrow=length(litter_nam), ncol=n.chain)
    sigma_y_C6alph=matrix(runif(15, 400, 6000), nrow=length(litter_nam), ncol=n.chain)
    sigma_y_C2alph=matrix(runif(15, 4000, 30000), nrow=length(litter_nam), ncol=n.chain)
  #create matrix for model initial estimate
    LC6=matrix(0, nrow=length(litter_nam), ncol=n.chain)
    LC2=matrix(0, nrow=length(litter_nam), ncol=n.chain)
  #estimate initial LC6 alpha based on other initial values
    for(b in 1:n.chain){
      for(a in 1:length(litter_nam)){
        LC6[a,b]<-proc_C6alpha(Fs=Fs[a,b], zC1_3a=zC1_3alp[a,b], Fdoc=Fdoc[b])
        LC2[a,b]<-proc_C2alpha(Fs=Fs[a,b], zC1_3a=zC1_3alp[a,b], Flig=Flig[a,b])
      }
    }

  #initial values for parameters that vary by litter type 
    latent_inits=list(Fs,Flig, zC1_3alp, sigma_y_C1_3,sigma_y_C6alph, sigma_y_C2alph, LC2, LC6)
    names(latent_inits)<-c("Fs", "Flig", "zC1_3alp","sigma_y_C1_3", "sigma_y_C6alph", "sigma_y_C2alph", "LC2alph", "LC6alph")



#run exp model based on parameters for each chain to generate initial values for first iteration

#use the above to create chain storage list (x), initialize chains, specify tuning and support
set_Cinits=setup(x1=latent_inits, num_L=length(latent_inits), n.iter=n.iter,
              n.chain=n.chain, parameter.names=param_names, dim.x=dim_param, 
              parameter_inits=parameter_inits, tune_param=tune_param, 
              support_param=support_param)

# #use the above to create chain storage list (x), initialize chains, specify tuning and support
# set_Cinits_HWE=setup(x1=latent_inits, num_L=length(latent_inits), n.iter=n.iter,
#                  n.chain=n.chain, parameter.names=param_names, dim.x=dim_param, 
#                  parameter_inits=parameter_inits, tune_param=tune_param_HWE, 
#                  support_param=support_param)
# 
# #use the above to create chain storage list (x), initialize chains, specify tuning and support
# set_Cinits_ydiff=setup(x1=latent_inits, num_L=length(latent_inits), n.iter=n.iter,
#                  n.chain=n.chain, parameter.names=param_names, dim.x=dim_param, 
#                  parameter_inits=parameter_inits, tune_param=tune_param_ydiff, 
#                  support_param=support_param)