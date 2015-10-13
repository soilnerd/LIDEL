#code to run post-processing independently, loading saved workspace
		# setwd("C:/LIDEL/Model3/Chain3")
		# library(deSolve, lib.loc="/data/rubelscratch/eecampbe/packages")
		# library(compiler)
		# library(lattice)
		# library(coda)
		# library(MASS)
		# library(gam, lib.loc="/data/rubelscratch/eecampbe/packages")
		# library(grid)
		# library(vcd, lib.loc="/data/rubelscratch/eecampbe/packages")
		# library(pscl, lib.loc="/data/rubelscratch/eecampbe/packages")
		# library(gplots, lib.loc="/data/rubelscratch/eecampbe/packages")

		# load("trialrun_workspace.RData")

#results are organized in x_Cinit in a list as follows:
# [[1]]=parameters
# [[2]]=all variance
# [[3]]-[[7]]=latent for all litter & time
# [[8]]=list of all process results from burnin:n.iter [[1]]-[[5]]
##################Parameter results###############
#create mcmc objects
#parameter results
  RESULTchain3_params=mcmc(x_Cinit[[1]][burnin:n.iter,])
  #save param results
    write.csv(RESULTchain3_params, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_params.csv")

#variance results
  RESULTchain3_var=mcmc(x_Cinit[[2]][burnin:n.iter,])
  #save var results
    write.csv(RESULTchain3_var, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_var.csv")

#latent state results
    RESULTchain3_alfalfa=mcmc(x_Cinit[[3]][burnin:n.iter,])
      #save alfalfa results
      write.csv(RESULTchain3_alfalfa, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_alfalfa.csv")
    RESULTchain3_ash=mcmc(x_Cinit[[4]][burnin:n.iter,])
      #save ash results
      write.csv(RESULTchain3_ash, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_ash.csv")
    RESULTchain3_bluestem=mcmc(x_Cinit[[5]][burnin:n.iter,])
      #save bluestem results
      write.csv(RESULTchain3_bluestem, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_bluestem.csv")
    RESULTchain3_oak=mcmc(x_Cinit[[6]][burnin:n.iter,])
      #save oak results
      write.csv(RESULTchain3_oak, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_oak.csv")
    RESULTchain3_pine=mcmc(x_Cinit[[7]][burnin:n.iter,])
      #save pine results
      write.csv(RESULTchain3_pine, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_pine.csv")

  #all latent states
    RESULTchain3_latent=mcmc(cbind(x_Cinit[[3]][burnin:n.iter,],
                      x_Cinit[[4]][burnin:n.iter,],
                      x_Cinit[[5]][burnin:n.iter,],
                      x_Cinit[[6]][burnin:n.iter,],
                      x_Cinit[[7]][burnin:n.iter,]))
    #save all latent results
      write.csv(RESULTchain3_latent, file="C:/LIDEL/Model3/Chain3/results/RESULTchain3_latent.csv")

#print results to file
  sink("C:/LIDEL/Model3/Chain3/results/chain3_summary.txt")

    print(burnin)
    print(n.iter)
    
    print("PARAMS")
      print(summary(RESULTchain3_params))
      print(geweke.diag(RESULTchain3_params))
    print("VAR")
      print(summary(RESULTchain3_var))
      print(geweke.diag(RESULTchain3_var))
    print("LATENT")
      print(summary(RESULTchain3_latent))
      print(geweke.diag(RESULTchain3_latent))

  sink()


################plot MCMC results################
#parameters
  jpeg(file="C:/LIDEL/Model3/Chain3/plots/RESULTchain3_PARAMS_A.jpg")
    plot(RESULTchain3_params[,1:3])
  dev.off()
  jpeg(file="C:/LIDEL/Model3/Chain3/plots/RESULTchain3_PARAMS_B.jpg")
    plot(RESULTchain3_params[,4:5])
  dev.off()

#variance
  for(p in 1:26){
    a=(p*3)-2
    c=(a+2)
    fileNam=paste("C:/LIDEL/Model3/Chain3/plots/RESULTchain3_var", p, ".jpg", sep="")
    jpeg(file=fileNam)
      plot(RESULTchain3_var[,a:c])
    dev.off()
  }

#latent, cols 1:63
for(p in 1:21){
  a=(p*3)-2
  c=(a+2)
  fileNam=paste("C:/LIDEL/Model3/Chain3/plots/RESULTchain3_latent", p, ".jpg", sep="")
  jpeg(file=fileNam)
    plot(RESULTchain3_latent[,a:c])
  dev.off()
}
#latent, cols 64:65
fileNam=paste("C:/LIDEL/Model3/Chain3/plots/RESULTchain3_latent", 22, ".jpg", sep="")
jpeg(file=fileNam)
  plot(RESULTchain3_latent[,64:65])
dev.off()


####################RMSE################################
all_rmse=vector("numeric", length=(length(burnin:n.iter)-1))

rmse<-function(error){
  sqrt(mean(error^2))
}

est_fin_alf_DOC=matrix(0, nrow=length(d_l_list[[1]][[2]][,1]), ncol=length(burnin:n.iter))
est_fin_ash_DOC=matrix(0, nrow=length(d_l_list[[2]][[2]][,1]), ncol=length(burnin:n.iter))
est_fin_blu_DOC=matrix(0, nrow=length(d_l_list[[3]][[2]][,1]), ncol=length(burnin:n.iter))
est_fin_oak_DOC=matrix(0, nrow=length(d_l_list[[4]][[2]][,1]), ncol=length(burnin:n.iter))
est_fin_pin_DOC=matrix(0, nrow=length(d_l_list[[5]][[2]][,1]), ncol=length(burnin:n.iter))

est_fin_alf_CO2=matrix(0, nrow=length(d_l_list[[1]][[3]][,1]), ncol=length(burnin:n.iter))
est_fin_ash_CO2=matrix(0, nrow=length(d_l_list[[2]][[3]][,1]), ncol=length(burnin:n.iter))
est_fin_blu_CO2=matrix(0, nrow=length(d_l_list[[3]][[3]][,1]), ncol=length(burnin:n.iter))
est_fin_oak_CO2=matrix(0, nrow=length(d_l_list[[4]][[3]][,1]), ncol=length(burnin:n.iter))
est_fin_pin_CO2=matrix(0, nrow=length(d_l_list[[5]][[3]][,1]), ncol=length(burnin:n.iter))


for(i in 1:length(all_rmse)){
#C6, DOC (calculate difference between measurement timepoints)
  #alfalfa non-estimated time points
    final_alfalfa_DOC=subset(x_Cinit[[8]][[1]][,7,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(non_d_l_list[[1]][[1]][,2]))-
      subset(x_Cinit[[8]][[1]][,7,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(non_d_l_list[[1]][[1]][,1]))
  #alfalfa error DOC
    error_alfalfa_DOC=data_nonMCMC[[1]][[1]][,3:5]-final_alfalfa_DOC[[1]]
  #alfalfa, estimated time points
    est_fin_alf_DOC[,i]=subset(x_Cinit[[8]][[1]][,7,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(d_l_list[[1]][[2]][,2]))-
      subset(x_Cinit[[8]][[1]][,7,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(d_l_list[[1]][[2]][,1]))

  #ash non-estimated time points
    final_ash_DOC=subset(x_Cinit[[8]][[2]][,7,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(non_d_l_list[[2]][[1]][,2]))-
      subset(x_Cinit[[8]][[2]][,7,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(non_d_l_list[[2]][[1]][,1]))
  #ash error DOC
    error_ash_DOC=data_nonMCMC[[2]][[1]][,3:5]-final_ash_DOC[[1]]
  #ash, estimated time points
    est_fin_ash_DOC[,i]=subset(x_Cinit[[8]][[2]][,7,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(d_l_list[[2]][[2]][,2]))-
      subset(x_Cinit[[8]][[2]][,7,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(d_l_list[[2]][[2]][,1]))
  
  #bluestem non-estimated time points
    final_bluestem_DOC=subset(x_Cinit[[8]][[3]][,7,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(non_d_l_list[[3]][[1]][,2]))-
      subset(x_Cinit[[8]][[3]][,7,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(non_d_l_list[[3]][[1]][,1]))
  #bluestem error DOC
    error_bluestem_DOC=data_nonMCMC[[3]][[1]][,3:5]-final_bluestem_DOC[[1]]
  #bluestem, estimated time points
    est_fin_blu_DOC[,i]=subset(x_Cinit[[8]][[3]][,7,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(d_l_list[[3]][[2]][,2]))-
      subset(x_Cinit[[8]][[3]][,7,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(d_l_list[[3]][[2]][,1]))
  
  #oak non-estimated time points
    final_oak_DOC=subset(x_Cinit[[8]][[4]][,7,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(non_d_l_list[[4]][[1]][,2]))-
      subset(x_Cinit[[8]][[4]][,7,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(non_d_l_list[[4]][[1]][,1]))
  #oak error DOC
    error_oak_DOC=data_nonMCMC[[4]][[1]][,3:5]-final_oak_DOC[[1]]
  #oak, estimated time points
    est_fin_oak_DOC[,i]=subset(x_Cinit[[8]][[4]][,7,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(d_l_list[[4]][[2]][,2]))-
      subset(x_Cinit[[8]][[4]][,7,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(d_l_list[[4]][[2]][,1]))
  
  #pine non-estimated time points  
    final_pine_DOC=subset(x_Cinit[[8]][[5]][,7,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(non_d_l_list[[5]][[1]][,2]))-
      subset(x_Cinit[[8]][[5]][,7,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(non_d_l_list[[5]][[1]][,1]))
  #pine error DOC
    error_pine_DOC=data_nonMCMC[[5]][[1]][,3:5]-final_pine_DOC[[1]]
  #pine, estimated time points
    est_fin_pin_DOC[,i]=subset(x_Cinit[[8]][[5]][,7,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(d_l_list[[5]][[2]][,2]))-
      subset(x_Cinit[[8]][[5]][,7,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(d_l_list[[5]][[2]][,1]))

  #
  #
  #C7, CO2 (calculate difference between measurement timepoints)
  #
  #
  #alfalfa non-estimated time points 
    final_alfalfa_CO2=subset(x_Cinit[[8]][[1]][,8,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(non_d_l_list[[1]][[2]][,2]))-
      subset(x_Cinit[[8]][[1]][,8,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(non_d_l_list[[1]][[2]][,1]))
  #alfalfa error CO2
    error_alfalfa_CO2=data_nonMCMC[[1]][[2]][,3:5]-final_alfalfa_CO2[[1]]
  #alfalfa, estimated time points
    est_fin_alf_CO2[,i]=subset(x_Cinit[[8]][[1]][,8,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(d_l_list[[1]][[3]][,2]))-
      subset(x_Cinit[[8]][[1]][,8,i], x_Cinit[[8]][[1]][,1,i] %in% as.numeric(d_l_list[[1]][[3]][,1]))
    
  #ash non-estimated time points 
    final_ash_CO2=subset(x_Cinit[[8]][[2]][,8,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(non_d_l_list[[2]][[2]][,2]))-
      subset(x_Cinit[[8]][[2]][,8,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(non_d_l_list[[2]][[2]][,1]))
  #ash error CO2
    error_ash_CO2=data_nonMCMC[[2]][[2]][,3:5]-final_ash_CO2[[1]]
  #ash, estimated time points
    est_fin_ash_CO2[,i]=subset(x_Cinit[[8]][[2]][,8,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(d_l_list[[2]][[3]][,2]))-
      subset(x_Cinit[[8]][[2]][,8,i], x_Cinit[[8]][[2]][,1,i] %in% as.numeric(d_l_list[[2]][[3]][,1]))
    
  #bluestem non-estimated time points 
    final_bluestem_CO2=subset(x_Cinit[[8]][[3]][,8,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(non_d_l_list[[3]][[2]][,2]))-
      subset(x_Cinit[[8]][[3]][,8,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(non_d_l_list[[3]][[2]][,1]))
  #bluestem error CO2 
    error_bluestem_CO2=data_nonMCMC[[3]][[2]][,3:5]-final_bluestem_CO2[[1]]
  #bluestem, estimated time points
    est_fin_blu_CO2[,i]=subset(x_Cinit[[8]][[3]][,8,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(d_l_list[[3]][[3]][,2]))-
      subset(x_Cinit[[8]][[3]][,8,i], x_Cinit[[8]][[3]][,1,i] %in% as.numeric(d_l_list[[3]][[3]][,1]))
    
  #oak non-estimated time points 
    final_oak_CO2=subset(x_Cinit[[8]][[4]][,8,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(non_d_l_list[[4]][[2]][,2]))-
      subset(x_Cinit[[8]][[4]][,8,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(non_d_l_list[[4]][[2]][,1]))
  #oak error CO2
    error_oak_CO2=data_nonMCMC[[4]][[2]][,3:5]-final_oak_CO2[[1]]
  #oak, estimated time points
    est_fin_oak_CO2[,i]=subset(x_Cinit[[8]][[4]][,8,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(d_l_list[[4]][[3]][,2]))-
      subset(x_Cinit[[8]][[4]][,8,i], x_Cinit[[8]][[4]][,1,i] %in% as.numeric(d_l_list[[4]][[3]][,1]))
    
  #pine non-estimated time points 
    final_pine_CO2=subset(x_Cinit[[8]][[5]][,8,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(non_d_l_list[[5]][[2]][,2]))-
      subset(x_Cinit[[8]][[5]][,8,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(non_d_l_list[[5]][[2]][,1]))
  #pine error CO2
    error_pine_CO2=data_nonMCMC[[5]][[2]][,3:5]-final_pine_CO2[[1]]
  #alfalfa, estimated time points
    est_fin_pin_CO2[,i]=subset(x_Cinit[[8]][[5]][,8,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(d_l_list[[5]][[3]][,2]))-
      subset(x_Cinit[[8]][[5]][,8,i], x_Cinit[[8]][[5]][,1,i] %in% as.numeric(d_l_list[[5]][[3]][,1]))

all_error=rbind(error_alfalfa_DOC, error_ash_DOC, error_bluestem_DOC, error_oak_DOC, error_pine_DOC,
                error_alfalfa_CO2, error_ash_CO2, error_bluestem_CO2, error_oak_CO2, error_pine_CO2)

all_rmse[i]=rmse(all_error)

}

	#save as MCMC result
		all_rmse_MCMC=mcmc(all_rmse)
	
	#write MCMC results
		write.csv(all_rmse, file="C:/LIDEL/Model3/Chain3/results/RMSE_data_chain3.csv")
		
	#write summary results to file
	  sink("C:/LIDEL/Model3/Chain3/results/RMSE_summary3.txt")
		  print(summary(all_rmse_MCMC))
	  sink()
	  
	#write plot of RMSE for MCMC to file
		jpeg(file="C:/LIDEL/Model3/Chain3/results/RMSE_results3.jpg")
		  plot(all_rmse_MCMC)
		dev.off()

############plot meas/mod comparison Mass####################
    #MASS
    mod_mass_alf=mcmc(cbind(apply(t(x_Cinit[[8]][[1]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[1]][366,2:6,]), 1, sum)))
    
    mod_mass_ash=mcmc(cbind(apply(t(x_Cinit[[8]][[2]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[2]][366,2:6,]), 1, sum)))
    
    mod_mass_blu=mcmc(cbind(apply(t(x_Cinit[[8]][[3]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[3]][366,2:6,]), 1, sum)))
    
    mod_mass_oak=mcmc(cbind(apply(t(x_Cinit[[8]][[4]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[4]][366,2:6,]), 1, sum)))
    
    mod_mass_pin=mcmc(cbind(apply(t(x_Cinit[[8]][[5]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[5]][366,2:6,]), 1, sum)))
							
	#code to print .csv files of modeled results to compare to measured estimates						
	CVmod_mass_alf3=cbind(apply(t(x_Cinit[[8]][[1]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[1]][366,2:6,]), 1, sum))
	CVmod_mass_ash3=cbind(apply(t(x_Cinit[[8]][[2]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[2]][366,2:6,]), 1, sum))
	CVmod_mass_blu3=cbind(apply(t(x_Cinit[[8]][[3]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[3]][366,2:6,]), 1, sum))
	CVmod_mass_oak3=cbind(apply(t(x_Cinit[[8]][[4]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[4]][366,2:6,]), 1, sum))
	CVmod_mass_pin3=cbind(apply(t(x_Cinit[[8]][[5]][96,2:6,]), 1, sum),
                            apply(t(x_Cinit[[8]][[5]][366,2:6,]), 1, sum))
		#write output for modeled mass measurements
		write.csv(CVmod_mass_alf3, file="C:/LIDEL/Model3/Chain3/results/CVmod_mass_alf3.csv")
		write.csv(CVmod_mass_ash3, file="C:/LIDEL/Model3/Chain3/results/CVmod_mass_ash3.csv")
		write.csv(CVmod_mass_blu3, file="C:/LIDEL/Model3/Chain3/results/CVmod_mass_blu3.csv")
		write.csv(CVmod_mass_oak3, file="C:/LIDEL/Model3/Chain3/results/CVmod_mass_oak3.csv")
		write.csv(CVmod_mass_pin3, file="C:/LIDEL/Model3/Chain3/results/CVmod_mass_pin3.csv")
		
		#write output for modeled DOC measurements
		write.csv(t(est_fin_alf_DOC), file="C:/LIDEL/Model3/Chain3/results/CVmod_DOC_alf3.csv")
		write.csv(t(est_fin_ash_DOC), file="C:/LIDEL/Model3/Chain3/results/CVmod_DOC_ash3.csv")
		write.csv(t(est_fin_blu_DOC), file="C:/LIDEL/Model3/Chain3/results/CVmod_DOC_blu3.csv")
		write.csv(t(est_fin_oak_DOC), file="C:/LIDEL/Model3/Chain3/results/CVmod_DOC_oak3.csv")
		write.csv(t(est_fin_pin_DOC), file="C:/LIDEL/Model3/Chain3/results/CVmod_DOC_pin3.csv")
		
		#write output for modeled CO2 measurements
		write.csv(t(est_fin_alf_CO2), file="C:/LIDEL/Model3/Chain3/results/CVmod_CO2_alf3.csv")
		write.csv(t(est_fin_ash_CO2), file="C:/LIDEL/Model3/Chain3/results/CVmod_CO2_ash3.csv")
		write.csv(t(est_fin_blu_CO2), file="C:/LIDEL/Model3/Chain3/results/CVmod_CO2_blu3.csv")
		write.csv(t(est_fin_oak_CO2), file="C:/LIDEL/Model3/Chain3/results/CVmod_CO2_oak3.csv")
		write.csv(t(est_fin_pin_CO2), file="C:/LIDEL/Model3/Chain3/results/CVmod_CO2_pin3.csv")
    
    alf1Mass=summary(mod_mass_alf)
    mean_alf1Mass=(alf1Mass[[1]][,1])
    uq_alf1Mass=(alf1Mass[[2]][,5])
    lq_alf1Mass=(alf1Mass[[2]][,1])
    
    #ash
    ash1Mass=summary(mod_mass_ash)
    mean_ash1Mass=(ash1Mass[[1]][,1])
    uq_ash1Mass=(ash1Mass[[2]][,5])
    lq_ash1Mass=(ash1Mass[[2]][,1])
    
    #bluestem
    bluestem1Mass=summary(mod_mass_blu)
    mean_bluestem1Mass=(bluestem1Mass[[1]][,1])
    uq_bluestem1Mass=(bluestem1Mass[[2]][,5])
    lq_bluestem1Mass=(bluestem1Mass[[2]][,1])
    
    #oak
    oak1Mass=summary(mod_mass_oak)
    mean_oak1Mass=(oak1Mass[[1]][,1])
    uq_oak1Mass=(oak1Mass[[2]][,5])
    lq_oak1Mass=(oak1Mass[[2]][,1])
    
    #pine
    pine1Mass=summary(mod_mass_pin)
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
    
  jpeg(file="C:/LIDEL/Model3/Chain3/results/Meas_Mod_MASS3.jpg")
    colors=c("orange", "lightblue", "lightgreen", "darkgrey", "pink", "blue")
    #plot mass loss
    plotCI(list_dataMass, new_mass, err="y", ui=new_mass_uq, li=new_mass_lq, 
           main="Mass Loss V1_4", ylab="Modeled total C loss", typ="p", 
           xlab="Measured total C loss")
    #plot of litter C remaining through time
    for(i in 1:num_litter){
      for(s in 3:ncol(all_data[[i]][[1]])){
        points(data_MCMC[[i]][[1]][,s], summary_litterMass[[i]], col=colors[s], 
               typ="p", pch=1)
      }
    }
    abline(0,1)
  dev.off()

############plot meas/mod comparison DOC####################
  fin_fin_alf_DOC=mcmc(t(est_fin_alf_DOC))
  fin_fin_ash_DOC=mcmc(t(est_fin_ash_DOC))
  fin_fin_blu_DOC=mcmc(t(est_fin_blu_DOC))
  fin_fin_oak_DOC=mcmc(t(est_fin_oak_DOC))
  fin_fin_pin_DOC=mcmc(t(est_fin_pin_DOC))
  
      #DOC
      alf1DOC=summary(fin_fin_alf_DOC)
      mean_alf1DOC=(alf1DOC[[1]][,1])
      uq_alf1DOC=(alf1DOC[[2]][,5])
      lq_alf1DOC=(alf1DOC[[2]][,1])
      
      #ash
      ash1DOC=summary(fin_fin_ash_DOC)
      mean_ash1DOC=(ash1DOC[[1]][,1])
      uq_ash1DOC=(ash1DOC[[2]][,5])
      lq_ash1DOC=(ash1DOC[[2]][,1])
      
      #bluestem
      bluestem1DOC=summary(fin_fin_blu_DOC)
      mean_bluestem1DOC=(bluestem1DOC[[1]][,1])
      uq_bluestem1DOC=(bluestem1DOC[[2]][,5])
      lq_bluestem1DOC=(bluestem1DOC[[2]][,1])
      
      #oak
      oak1DOC=summary(fin_fin_oak_DOC)
      mean_oak1DOC=(oak1DOC[[1]][,1])
      uq_oak1DOC=(oak1DOC[[2]][,5])
      lq_oak1DOC=(oak1DOC[[2]][,1])
      
      #pine
      pine1DOC=summary(fin_fin_pin_DOC)
      mean_pine1DOC=(pine1DOC[[1]][,1])
      uq_pine1DOC=(pine1DOC[[2]][,5])
      lq_pine1DOC=(pine1DOC[[2]][,1])
      
      #for plotting
      summary_litterDOC=list(mean_alf1DOC, mean_ash1DOC, mean_bluestem1DOC, mean_oak1DOC, mean_pine1DOC)
      
      new_DOC=c(mean_alf1DOC, mean_ash1DOC, mean_bluestem1DOC, mean_oak1DOC, mean_pine1DOC)  
      new_DOC_uq=c(uq_alf1DOC, uq_ash1DOC, uq_bluestem1DOC, uq_oak1DOC, uq_pine1DOC) 
      new_DOC_lq=c(lq_alf1DOC, lq_ash1DOC, lq_bluestem1DOC, lq_oak1DOC, lq_pine1DOC) 
      
      
      list_dataDOC=c(rowMeans(data_MCMC[[1]][[2]][,3:5]),
                      rowMeans(data_MCMC[[2]][[2]][,3:5]),
                      rowMeans(data_MCMC[[3]][[2]][,3:5]),
                      rowMeans(data_MCMC[[4]][[2]][,3:5]),
                      rowMeans(data_MCMC[[5]][[2]][,3:5]))  

    jpeg(file="C:/LIDEL/Model3/Chain3/results/Meas_Mod_DOC3.jpg")
      colors=c("orange", "lightblue", "lightgreen", "darkgrey", "pink", "blue")
      #plot mass loss
      plotCI(list_dataDOC, new_DOC, err="y", ui=new_DOC_uq, li=new_DOC_lq, 
             main="DOC V1_4", ylab="Modeled DOC", typ="p", 
             xlab="Measured DOC")
      #plot of litter C remaining through time
      for(i in 1:num_litter){
        for(s in 3:ncol(all_data[[i]][[1]])){
          points(data_MCMC[[i]][[2]][,s], summary_litterDOC[[i]], col=colors[s], 
                 typ="p", pch=1)
        }
      }
      abline(0,1)
    dev.off()

############plot meas/mod comparison CO2####################
  fin_fin_alf_CO2=mcmc(t(est_fin_alf_CO2))
  fin_fin_ash_CO2=mcmc(t(est_fin_ash_CO2))
  fin_fin_blu_CO2=mcmc(t(est_fin_blu_CO2))
  fin_fin_oak_CO2=mcmc(t(est_fin_oak_CO2))
  fin_fin_pin_CO2=mcmc(t(est_fin_pin_CO2))
      
      #DOC
      alf1CO2=summary(fin_fin_alf_CO2)
      mean_alf1CO2=(alf1CO2[[1]][,1])
      uq_alf1CO2=(alf1CO2[[2]][,5])
      lq_alf1CO2=(alf1CO2[[2]][,1])
      
      #ash
      ash1CO2=summary(fin_fin_ash_CO2)
      mean_ash1CO2=(ash1CO2[[1]][,1])
      uq_ash1CO2=(ash1CO2[[2]][,5])
      lq_ash1CO2=(ash1CO2[[2]][,1])
      
      #bluestem
      bluestem1CO2=summary(fin_fin_blu_CO2)
      mean_bluestem1CO2=(bluestem1CO2[[1]][,1])
      uq_bluestem1CO2=(bluestem1CO2[[2]][,5])
      lq_bluestem1CO2=(bluestem1CO2[[2]][,1])
      
      #oak
      oak1CO2=summary(fin_fin_oak_CO2)
      mean_oak1CO2=(oak1CO2[[1]][,1])
      uq_oak1CO2=(oak1CO2[[2]][,5])
      lq_oak1CO2=(oak1CO2[[2]][,1])
      
      #pine
      pine1CO2=summary(fin_fin_pin_CO2)
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
      
      
    jpeg(file="C:/LIDEL/Model3/Chain3/results/Meas_Mod_CO23.jpg")
      colors=c("orange", "lightblue", "lightgreen", "darkgrey", "pink", "blue")
      #plot mass loss
      plotCI(list_dataCO2, new_CO2, err="y", ui=new_CO2_uq, li=new_CO2_lq, 
             main="CO2 V1_4", ylab="Modeled CO2", typ="p", 
             xlab="Measured CO2")
      #plot of litter C remaining through time
      for(i in 1:num_litter){
        for(s in 3:ncol(all_data[[i]][[1]])){
          points(data_MCMC[[i]][[3]][,s], summary_litterCO2[[i]], col=colors[s], 
                 typ="p", pch=1)
        }
      }
      abline(0,1)
    dev.off()
