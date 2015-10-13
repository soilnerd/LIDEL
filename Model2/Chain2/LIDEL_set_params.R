##########LIDEL MODEL INPUTS##########
LIDEL_inputs=list()

#simulation time period
LIDEL_inputs$days=seq(0,365, by=1) 
##########set parameters
# % nitrogen (tau)
LIDEL_inputs$tau=Init_N[,2]; 
#inputs for linear equations determining liebig's law of minimum for DOC from soluble (1) and NSS (2)
LIDEL_inputs$EM1=0.1; LIDEL_inputs$Em1=0.001; LIDEL_inputs$EM2=0.5; LIDEL_inputs$Em2=0.1; 
#inputs for linear equation for internal vs external N limitation, and max Lc
LIDEL_inputs$Nmid=1.75; LIDEL_inputs$NC=3; LIDEL_inputs$LcM=0.51 
#maximum microbial growth efficiency from soluble (1) and NSS (2)
LIDEL_inputs$beta1=0.6; LIDEL_inputs$beta2=0.5; 
#DOC from lignin and microbe products
LIDEL_inputs$lambda3=0.038
# Initial LCI
LIDEL_inputs$Lcalph=Init_LCI[,2]

#########TUNING PARAMETERS FOR LATENT STATES
TUNE_DOC=read.csv("data_TUNING_DOC.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
TUNE_CO2=read.csv("data_TUNING_CO2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)

Latent_tune=list()
Latent_tune$Alfalfa=list(matrix(c(0, 95, 365, 10, 1000, 1000), ncol=2), 
                         cbind(TUNE_DOC[,1:2], TUNE_DOC[,3]),
                         cbind(TUNE_CO2[,1:2], TUNE_CO2[,3]))
names(Latent_tune$Alfalfa)<-c("Mass", "DOC", "CO2")

Latent_tune$Ash=list(matrix(c(0, 95, 365, 10, 1000, 1000), ncol=2), 
                         cbind(TUNE_DOC[,1:2], TUNE_DOC[,4]),
                         cbind(TUNE_CO2[1:18,4:5], TUNE_CO2[1:18,6]))
names(Latent_tune$Ash)<-c("Mass", "DOC", "CO2")

Latent_tune$Bluestem=list(matrix(c(0, 95, 365, 10, 1000, 1000), ncol=2), 
                     cbind(TUNE_DOC[,1:2], TUNE_DOC[,5]),
                     cbind(TUNE_CO2[1:18,4:5], TUNE_CO2[1:18,7]))
names(Latent_tune$Bluestem)<-c("Mass", "DOC", "CO2")

Latent_tune$Oak=list(matrix(c(0, 95, 365, 10, 1000, 1000), ncol=2), 
                     cbind(TUNE_DOC[,1:2], TUNE_DOC[,6]),
                     cbind(TUNE_CO2[1:18,4:5], TUNE_CO2[1:18,8]))
names(Latent_tune$Oak)<-c("Mass", "DOC", "CO2")

Latent_tune$Pine=list(matrix(c(0, 95, 365, 10, 1000, 1000), ncol=2), 
                     cbind(TUNE_DOC[,1:2], TUNE_DOC[,7]),
                     cbind(TUNE_CO2[1:18,4:5], TUNE_CO2[1:18,9]))
names(Latent_tune$Pine)<-c("Mass", "DOC", "CO2")


##########TIME SUBSET, FOR LATENT STATE AND TUNING
d_l_list=list()
Latent_tune_sub=list()
non_d_l_list=list()


#list of dates to subset
t_sub=c(7,15,28,64,95,365)
t_nonMCMC_sub=c(1,2,4,6,9,10,13,18,20,39,49,76,118,152,181,228,284)

#subset dates for each litter type, one list per litter
for(i in 1:length(litter_nam)){
  d_l_list[[i]]=list(c(95, 365), 
              subset(Latent_tune[[i]][[2]][,1:2], Latent_tune[[i]][[2]][,2] %in% t_sub),
              subset(Latent_tune[[i]][[3]][,1:2], Latent_tune[[i]][[3]][,2] %in% t_sub))
  names(d_l_list[[i]])=c("Mass", "DOC", "CO2")
  non_d_l_list[[i]]=list(subset(Latent_tune[[i]][[2]][,1:2], Latent_tune[[i]][[2]][,2] %in% t_nonMCMC_sub),
                     subset(Latent_tune[[i]][[3]][,1:2], Latent_tune[[i]][[3]][,2] %in% t_nonMCMC_sub))
  names(non_d_l_list[[i]])=c("DOC", "CO2")
  Latent_tune_sub[[i]]=list(matrix(c(95, 365, 1000, 1000), ncol=2), 
                     subset(Latent_tune[[i]][[2]][,2:3], Latent_tune[[i]][[2]][,2] %in% t_sub),
                     subset(Latent_tune[[i]][[3]][,2:3], Latent_tune[[i]][[3]][,2] %in% t_sub))
  names(Latent_tune_sub[[i]])=c("Mass", "DOC", "CO2")
}
#name list by litter
names(d_l_list)=litter_nam
names(Latent_tune_sub)=litter_nam
names(non_d_l_list)=litter_nam

##########SUBSET DATA FOR USE IN MCMC, FROM ALL MEASURED DATA AVAILABLE
data_MCMC=list()
for(i in 1:length(litter_nam)){
  data_MCMC[[i]]=list()
  for(j in 1:3){
    data_MCMC[[i]][[j]]=subset(all_data[[i]][[j]][,1:5], all_data[[i]][[j]][,2] %in% t_sub)
  }
  names(data_MCMC[[i]])=c("Mass", "DOC", "CO2")
}
names(data_MCMC)=litter_nam

#make mean of all pine DOC for day 39-64 fill the missing value for litter 1
data_MCMC[[5]][[2]][3,3]=mean(all_data[[5]][[2]][5,4:8])


data_nonMCMC=list()
  for(i in 1:length(litter_nam)){
    data_nonMCMC[[i]]=list()
    for(j in 1:2){
      f=j+1
      data_nonMCMC[[i]][[j]]=subset(all_data[[i]][[f]][,1:5], all_data[[i]][[f]][,2] %in% t_nonMCMC_sub)
    }
    names(data_nonMCMC[[i]])=c("DOC", "CO2")
  }
  names(data_nonMCMC)=litter_nam

