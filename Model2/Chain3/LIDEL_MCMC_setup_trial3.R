##########SETUP FOR MCMC##########
#using all data
#Nell Campbell
#moved to separate file 11/17/14
##################################

#####iterations and chains
#Set up the number of iterations. While you are getting things working, this can be small, ca 5000
n.chain=1
num_litter=5
model_names=c("Day", "C1_5", "C6", "C7")

#####specify inputs to model
#name all inputs, including parameters and derived values
param_names=c("k1", "k2", "k4", "bet3", "lam2",
              "sigma_y_C1_5", "sigma_y_C6", "sigma_y_C7",
              "sigma_L_C1_5", "sigma_L_C6", "sigma_L_C7",
              litter_nam)
#specify the dimension of all sites (either number of sites- 1 currently-, or sites and time, ie scalar or vector)
dim_param=c(1,1,1,1,1,
            length(d_l_list[[1]][[1]]),nrow(d_l_list[[1]][[2]]),nrow(d_l_list[[1]][[3]]),
            length(d_l_list[[1]][[1]]),nrow(d_l_list[[1]][[2]]),nrow(d_l_list[[1]][[3]]))
#names of parameters constant across litter types and time (KEEP ORDER THE SAME THROUGHOUT)
params=c("k1", "k2", "k4", "bet3", "lam2",
         "sigma_y_C1_5", "sigma_y_C6", "sigma_y_C7",
         "sigma_L_C1_5", "sigma_L_C6", "sigma_L_C7")
#specify how the parameters will be 'tuned' during the runs
tune_param=c(0.1, 0.005, 0.1, 0.1, 0.05,
             5, 5, 5,
             1, 6, 6)
#what sort of data are the parameters?
support_param=c("zero-one", "zero-one","zero-one","zero-one","zero-one", 
                "non-negative", "non-negative","non-negative",
                "non-negative", "non-negative","non-negative",
                "real", "real","real", "real", "real")
#create matrix for to initialize chains
parameter_inits=matrix(NA, nrow=n.chain, ncol=length(params))
#set parameter values to initialize each chain
k1=runif(1,min=0.05, max=0.25)
k2=runif(1,min=0.002, max=0.01)
k4=runif(1,min=0.2, max=0.8)
bet3=runif(1,min=0.4,max=0.8)
lam2=runif(1,min=0.02,max=0.1)
sigma_y_C1_5=matrix(runif(1,min=10,max=1000), ncol=n.chain, nrow=length(d_l_list[[1]][[1]]))
sigma_y_C6=matrix(runif(1,min=10,max=1000), ncol=n.chain, nrow=nrow(d_l_list[[1]][[2]]))
sigma_y_C7=matrix(runif(1,min=10,max=1000), ncol=n.chain, nrow=nrow(d_l_list[[1]][[3]]))
sigma_L_C1_5=matrix(runif(1,min=10,max=1000), ncol=n.chain, nrow=length(d_l_list[[1]][[1]]))
sigma_L_C6=matrix(runif(1,min=10,max=1000), ncol=n.chain, nrow=nrow(d_l_list[[1]][[2]]))
sigma_L_C7=matrix(runif(1,min=10,max=1000), ncol=n.chain, nrow=nrow(d_l_list[[1]][[3]]))
#loop to set init values within each parameter_inits
parameter_inits=list()
for(j in 1:n.chain){
  parameter_inits[[j]]=list()
  parameter_inits[[j]][[1]]=k1[j]
  parameter_inits[[j]][[2]]=k2[j]
  parameter_inits[[j]][[3]]=k4[j]
  parameter_inits[[j]][[4]]=bet3[j]
  parameter_inits[[j]][[5]]=lam2[j]
  parameter_inits[[j]][[6]]=sigma_y_C1_5[,j]
  parameter_inits[[j]][[7]]=sigma_y_C6[,j]
  parameter_inits[[j]][[8]]=sigma_y_C7[,j]
  parameter_inits[[j]][[9]]=sigma_L_C1_5[,j]
  parameter_inits[[j]][[10]]=sigma_L_C6[,j]
  parameter_inits[[j]][[11]]=sigma_L_C7[,j]
}	

#parameters that vary by litter type
#run LIDEL model based on parameters for each chain to generate initial values for first iteration
init_mu=list()
init_proc_days=list()
for(i in 1:n.chain){
  init_mu[[i]]=list()
  init_proc_days[[i]]=list()
  for(j in 1:num_litter){
  Pars <- c(tau=LIDEL_inputs$tau[j], EM1=LIDEL_inputs$EM1, Em1=LIDEL_inputs$Em1, EM2=LIDEL_inputs$EM2, 
            Em2=LIDEL_inputs$Em2, Nmid=LIDEL_inputs$Nmid, NC=LIDEL_inputs$NC,LcM=LIDEL_inputs$LcM, 
            beta1=LIDEL_inputs$beta1, beta2=LIDEL_inputs$beta2, lambda3=LIDEL_inputs$lambda3, 
            k1=parameter_inits[[1]][[1]],
            k2=parameter_inits[[1]][[2]], 
            k4=parameter_inits[[1]][[3]], 
            beta3=parameter_inits[[1]][[4]], 
            lambda2=parameter_inits[[1]][[5]],
            Lcalph=LIDEL_inputs$Lcalph[j])
  #initial condition at day 0
  start_initial_vals=LIDEL_initials(mass=Init_fctn_vals[[j]][1:2], fdoc=Init_fctn_vals[[j]][3:4], fs=Init_fctn_vals[[j]][5:6], flig=Init_fctn_vals[[j]][7:8])
  yini <- c(C1=start_initial_vals[2], C2=start_initial_vals[3], C3=start_initial_vals[4], C4=0, 
            C5=0, C6=start_initial_vals[5], C7=0)
  print(start_initial_vals[1])
  print(sum(yini))
  out<-lsoda(func=DOCmod, y=yini, parms=Pars, times=LIDEL_inputs$days)
  #save model results for pools 1-7
  init_proc_days[[i]][[j]]=out[,1:8]
  #pull daily values for measured data
  init_mu[[i]][[j]]=cbind(as.vector(out[,1]),as.vector(apply(out[,2:6], 1, sum)), 
                          as.vector(out[,7]), as.vector(out[,8]))
  #replace 1st row so that total mass is including DOC init
  init_mu[[i]][[j]][1,]=cbind(as.vector(out[1,1]),as.vector(sum(out[1,2:8])), 
                              as.vector(out[1,7]), as.vector(out[1,8]))
  colnames(init_mu[[i]][[j]])=model_names
  }
  names(init_mu[[i]])=litter_nam
}

##############plot initial vs measured results#############
colors=c("black", "orange", "lightblue", "lightgreen", "darkgrey", "pink", "blue")

  #plot init_mu vs measured remaining litter c
   jpeg(file="C:/LIDEL/Model2/Chain3/plots/litterCinit1.jpg")

      par(mfrow=c(3,2), oma=c(0,0,2,0))
      #plot of litter C remaining through time
      for(t in 1:length(litter_nam)){
        Ctot_lim=max(rowSums(as.matrix(Init_Cs[,2:5])))+10
        plot(0, rowSums(as.matrix(Init_Cs[t,2:5])), col=colors, ylim=c(0,Ctot_lim), xlim=c(0, 365),
             main=litter_nam[t], ylab="Remaining litter C (mg)", typ="p", xlab="Time", lwd=3, pch=4)
        for(s in 3:ncol(all_data[[t]][[1]])){
          points(all_data[[t]][[1]][,2], all_data[[t]][[1]][,s], col=colors[s-1], 
                 typ="p", pch=4, lwd=3)
        }
        points(init_mu[[1]][[t]][,1], init_mu[[1]][[t]][,2],  
               typ="l", lwd=3)
      }
      mtext("chain1", outer=TRUE)
    dev.off()

  #plot of DOC accumulation through time
  jpeg(file="C:/LIDEL/Model2/Chain3/plots/DOCinit1.jpg")
    par(mfrow=c(3,2), oma=c(0,0,2,0))
    for(t in 1:num_litter){
      DOC_lim=max(cumul_DOC[[t]][,2:ncol(cumul_DOC[[t]])], init_mu[[1]][[t]][,3])+10
      plot(init_mu[[1]][[t]][,1], init_mu[[1]][[t]][,3], col=colors[1], ylim=c(0,DOC_lim),
           main=litter_nam[t], ylab="DOC (mg)", typ="l", xlab="Time", lwd=3)
      for(s in 2:ncol(cumul_DOC[[1]])){
        if(cumul_DOC[[t]][9,s]!=-99.99){
          points(cumul_DOC[[t]][,1], cumul_DOC[[t]][,s], col=colors[s], 
                 typ="l", pch=1, cex=3, lwd=3)
        }
        if(cumul_DOC[[t]][9,s]==-99.99){
          points(cumul_DOC[[t]][1:7,1], cumul_DOC[[t]][1:7,s], col=colors[s], 
                 typ="l", pch=1, cex=3, lwd=3)
        }
      }
      points(init_mu[[1]][[t]][,1], init_mu[[1]][[t]][,3], col=colors[1],
          typ="l", lwd=3)
    }
    mtext("chain1", outer=TRUE)
    #dev.off()

  #plot of CO2 accumulation through time
    jpeg(file="C:/LIDEL/Model2/Chain3/plots/CO2init1.jpg")
    par(mfrow=c(3,2), oma=c(0,0,2,0))
    for(t in 1:num_litter){
      CO2_lim=max(cumul_CO2[[t]][,2:ncol(cumul_CO2[[t]])], init_mu[[1]][[t]][,4])+10
      plot(init_mu[[1]][[t]][,1], init_mu[[1]][[t]][,4], col=colors[1], ylim=c(0,CO2_lim),
           main=litter_nam[t], ylab="CO2-C (mg)", typ="l", xlab="Time", lwd=3)
      for(s in 2:ncol(cumul_CO2[[1]])){
        if(cumul_CO2[[t]][18,s]!=-99.99){
          points(cumul_CO2[[t]][,1], cumul_CO2[[t]][,s], col=colors[s], 
                 typ="l", pch=1, cex=3, lwd=3)
        }
        if(cumul_CO2[[t]][18,s]==-99.99){
          if(t==1){
            points(cumul_CO2[[t]][1:17,1], cumul_CO2[[t]][1:17,s], col=colors[s-1], 
                   typ="l", pch=1, cex=3, lwd=3)
          }
          if(t!=1){
            points(cumul_CO2[[t]][1:12,1], cumul_CO2[[t]][1:12,s], col=colors[s-1], 
                   typ="l", pch=1, cex=3, lwd=3)
          }
        }
      }
      points(init_mu[[1]][[t]][,1], init_mu[[1]][[t]][,4], col=colors[1],
          typ="l", lwd=3)
    }
    mtext("chain1", outer=TRUE)
    dev.off()

#################Calculate Difference#########################
#calculate difference to use as initial values in setup of latent states
#Modeled difference, using measured time points
init_alfalfa=list()
init_ash=list()
init_bluestem=list()
init_oak=list()
init_pine=list()

#############for all timepoints##########
#for C1_5, measurement is in same increment as model results (mass remaining)
init_alfalfa[[1]]=init_mu[[1]][[1]][all_data[[1]][[1]][,2]+1,2]

init_ash[[1]]=init_mu[[1]][[2]][all_data[[2]][[1]][,2]+1,2]

init_bluestem[[1]]=init_mu[[1]][[3]][all_data[[3]][[1]][,2]+1,2]

init_oak[[1]]=init_mu[[1]][[4]][all_data[[4]][[1]][,2]+1,2]

init_pine[[1]]=init_mu[[1]][[5]][all_data[[5]][[1]][,2]+1,2]

#C6, DOC (calculate difference between measurement timepoints)
init_alfalfa[[2]]=subset(init_mu[[1]][[1]][,3], init_mu[[1]][[1]][,1] %in% as.numeric(all_data[[1]][[2]][,2]))-
  subset(init_mu[[1]][[1]][,3], init_mu[[1]][[1]][,1] %in% as.numeric(all_data[[1]][[2]][,1]))
  
init_ash[[2]]=subset(init_mu[[1]][[2]][,3], init_mu[[1]][[2]][,1] %in% as.numeric(all_data[[2]][[2]][,2]))-
  subset(init_mu[[1]][[2]][,3], init_mu[[1]][[2]][,1] %in% as.numeric(all_data[[2]][[2]][,1]))

init_bluestem[[2]]=subset(init_mu[[1]][[3]][,3], init_mu[[1]][[3]][,1] %in% as.numeric(all_data[[3]][[2]][,2]))-
  subset(init_mu[[1]][[3]][,3], init_mu[[1]][[3]][,1] %in% as.numeric(all_data[[3]][[2]][,1]))

init_oak[[2]]=subset(init_mu[[1]][[4]][,3], init_mu[[1]][[4]][,1] %in% as.numeric(all_data[[4]][[2]][,2]))-
  subset(init_mu[[1]][[4]][,3], init_mu[[1]][[4]][,1] %in% as.numeric(all_data[[4]][[2]][,1]))

init_pine[[2]]=subset(init_mu[[1]][[5]][,3], init_mu[[1]][[5]][,1] %in% as.numeric(all_data[[5]][[2]][,2]))-
  subset(init_mu[[1]][[5]][,3], init_mu[[1]][[5]][,1] %in% as.numeric(all_data[[5]][[2]][,1]))

#C7, CO2 (calculate difference between measurement timepoints)
init_alfalfa[[3]]=subset(init_mu[[1]][[1]][,4], init_mu[[1]][[1]][,1] %in% as.numeric(all_data[[1]][[3]][,2]))-
  subset(init_mu[[1]][[1]][,4], init_mu[[1]][[1]][,1] %in% as.numeric(all_data[[1]][[3]][,1]))

init_ash[[3]]=subset(init_mu[[1]][[2]][,4], init_mu[[1]][[2]][,1] %in% as.numeric(all_data[[2]][[3]][,2]))-
  subset(init_mu[[1]][[2]][,4], init_mu[[1]][[2]][,1] %in% as.numeric(all_data[[2]][[3]][,1]))

init_bluestem[[3]]=subset(init_mu[[1]][[3]][,4], init_mu[[1]][[3]][,1] %in% as.numeric(all_data[[3]][[3]][,2]))-
  subset(init_mu[[1]][[3]][,4], init_mu[[1]][[3]][,1] %in% as.numeric(all_data[[3]][[3]][,1]))

init_oak[[3]]=subset(init_mu[[1]][[4]][,4], init_mu[[1]][[4]][,1] %in% as.numeric(all_data[[4]][[3]][,2]))-
  subset(init_mu[[1]][[4]][,4], init_mu[[1]][[4]][,1] %in% as.numeric(all_data[[4]][[3]][,1]))

init_pine[[3]]=subset(init_mu[[1]][[5]][,4], init_mu[[1]][[5]][,1] %in% as.numeric(all_data[[5]][[3]][,2]))-
  subset(init_mu[[1]][[5]][,4], init_mu[[1]][[5]][,1] %in% as.numeric(all_data[[5]][[3]][,1]))

init_mu_all=list(init_alfalfa, init_ash, init_bluestem, init_oak, init_pine)

##############for d_l_list subset time points
init_alfalfa_d=list()
init_ash_d=list()
init_bluestem_d=list()
init_oak_d=list()
init_pine_d=list()

#for C1_5, measurement is in same increment as model results (mass remaining)
  init_alfalfa_d[[1]]=init_mu[[1]][[1]][d_l_list[[1]][[1]]+1,2]
  
  init_ash_d[[1]]=init_mu[[1]][[2]][d_l_list[[2]][[1]]+1,2]
  
  init_bluestem_d[[1]]=init_mu[[1]][[3]][d_l_list[[3]][[1]]+1,2]
  
  init_oak_d[[1]]=init_mu[[1]][[4]][d_l_list[[4]][[1]]+1,2]
  
  init_pine_d[[1]]=init_mu[[1]][[5]][d_l_list[[5]][[1]]+1,2]

#C6, DOC (calculate difference between measurement timepoints)
  init_alfalfa_d[[2]]=subset(init_mu[[1]][[1]][,3], init_mu[[1]][[1]][,1] %in% as.numeric(d_l_list[[1]][[2]][,2]))-
    subset(init_mu[[1]][[1]][,3], init_mu[[1]][[1]][,1] %in% as.numeric(d_l_list[[1]][[2]][,1]))
	
	#initial conditions tend to result in a final alfalfa DOC value that is very, very far from its converged value
    #therefore takes a long time to get there because the tuning size is small
    #set value to be more reasonable for initial conditions
    init_alfalfa_d[[2]][length(init_alfalfa_d[[2]])]=rnorm(1, 100, 10)
  
  init_ash_d[[2]]=subset(init_mu[[1]][[2]][,3], init_mu[[1]][[2]][,1] %in% as.numeric(d_l_list[[2]][[2]][,2]))-
    subset(init_mu[[1]][[2]][,3], init_mu[[1]][[2]][,1] %in% as.numeric(d_l_list[[2]][[2]][,1]))
  
  init_bluestem_d[[2]]=subset(init_mu[[1]][[3]][,3], init_mu[[1]][[3]][,1] %in% as.numeric(d_l_list[[3]][[2]][,2]))-
    subset(init_mu[[1]][[3]][,3], init_mu[[1]][[3]][,1] %in% as.numeric(d_l_list[[3]][[2]][,1]))
  
  init_oak_d[[2]]=subset(init_mu[[1]][[4]][,3], init_mu[[1]][[4]][,1] %in% as.numeric(d_l_list[[4]][[2]][,2]))-
    subset(init_mu[[1]][[4]][,3], init_mu[[1]][[4]][,1] %in% as.numeric(d_l_list[[4]][[2]][,1]))
  
  init_pine_d[[2]]=subset(init_mu[[1]][[5]][,3], init_mu[[1]][[5]][,1] %in% as.numeric(d_l_list[[5]][[2]][,2]))-
    subset(init_mu[[1]][[5]][,3], init_mu[[1]][[5]][,1] %in% as.numeric(d_l_list[[5]][[2]][,1]))

#C7, CO2 (calculate difference between measurement timepoints)
  init_alfalfa_d[[3]]=subset(init_mu[[1]][[1]][,4], init_mu[[1]][[1]][,1] %in% as.numeric(d_l_list[[1]][[3]][,2]))-
    subset(init_mu[[1]][[1]][,4], init_mu[[1]][[1]][,1] %in% as.numeric(d_l_list[[1]][[3]][,1]))
  
  init_ash_d[[3]]=subset(init_mu[[1]][[2]][,4], init_mu[[1]][[2]][,1] %in% as.numeric(d_l_list[[2]][[3]][,2]))-
    subset(init_mu[[1]][[2]][,4], init_mu[[1]][[2]][,1] %in% as.numeric(d_l_list[[2]][[3]][,1]))
  
  init_bluestem_d[[3]]=subset(init_mu[[1]][[3]][,4], init_mu[[1]][[3]][,1] %in% as.numeric(d_l_list[[3]][[3]][,2]))-
    subset(init_mu[[1]][[3]][,4], init_mu[[1]][[3]][,1] %in% as.numeric(d_l_list[[3]][[3]][,1]))
  
  init_oak_d[[3]]=subset(init_mu[[1]][[4]][,4], init_mu[[1]][[4]][,1] %in% as.numeric(d_l_list[[4]][[3]][,2]))-
    subset(init_mu[[1]][[4]][,4], init_mu[[1]][[4]][,1] %in% as.numeric(d_l_list[[4]][[3]][,1]))
  
  init_pine_d[[3]]=subset(init_mu[[1]][[5]][,4], init_mu[[1]][[5]][,1] %in% as.numeric(d_l_list[[5]][[3]][,2]))-
    subset(init_mu[[1]][[5]][,4], init_mu[[1]][[5]][,1] %in% as.numeric(d_l_list[[5]][[3]][,1]))
  
  init_mu_final=list(init_alfalfa_d, init_ash_d, init_bluestem_d, init_oak_d, init_pine_d)

#use the above to create chain storage list (x), initialize chains, specify tuning and support
set_LIDEL1=setup(n.iter=n.iter,n.chain=n.chain, parameter.names=param_names, 
          dim.x=dim_param, parameter_inits=parameter_inits, 
          init_mu_input=init_mu_final, init_proc_days=init_proc_days,
          tune_param=tune_param, Latent_tune=Latent_tune_sub, support_param=support_param)
