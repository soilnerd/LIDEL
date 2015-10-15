##########TEST RESULTS FOR CONVERGENCE##########
Calc_NSS=function(x){
  results=list()
  for(s in 1:n.chain){
  #calculate non-soluble structural pool
  Fnss_alf=1-(x[[s]][burnin:n.iter,2]+x[[s]][burnin:n.iter,7])
  Fnss_ash=1-(x[[s]][burnin:n.iter,3]+x[[s]][burnin:n.iter,8])
  Fnss_blue=1-(x[[s]][burnin:n.iter,4]+x[[s]][burnin:n.iter,9])
  Fnss_oak=1-(x[[s]][burnin:n.iter,5]+x[[s]][burnin:n.iter,10])
  Fnss_pine=1-(x[[s]][burnin:n.iter,6]+x[[s]][burnin:n.iter,11])
  results[[s]]=mcmc(cbind(Fnss_alf, Fnss_ash, Fnss_blue, Fnss_oak, Fnss_pine))
  }
  return(results)
}

all_results_converge_test=function(x){
    #to check key parameters
    burn_results=list()
    #convert to mcmc and summarize results for each chain
    for(s in 1:n.chain){
      #moment match for mean of lignin fraction, based on alpha and beta shape parameters
      	est_Flig_alf=x[[s]][burnin:n.iter,42]/(x[[s]][burnin:n.iter,42]+
      	                                               x[[s]][burnin:n.iter,47])
      	est_Flig_ash=x[[s]][burnin:n.iter,43]/(x[[s]][burnin:n.iter,43]+
      	                                               x[[s]][burnin:n.iter,48])
      	est_Flig_blue=x[[s]][burnin:n.iter,44]/(x[[s]][burnin:n.iter,44]+
      	                                                x[[s]][burnin:n.iter,49])
      	est_Flig_oak=x[[s]][burnin:n.iter,45]/(x[[s]][burnin:n.iter,45]+
      	                                               x[[s]][burnin:n.iter,50])
      	est_Flig_pine=x[[s]][burnin:n.iter,46]/(x[[s]][burnin:n.iter,46]+
      	                                               x[[s]][burnin:n.iter,51])
      	#moment match for mean of soluble fraction, based on alpha and beta shape parameters
      	est_Fs_alf=x[[s]][burnin:n.iter,32]/(x[[s]][burnin:n.iter,32]+
      	                                               x[[s]][burnin:n.iter,37])
      	est_Fs_ash=x[[s]][burnin:n.iter,33]/(x[[s]][burnin:n.iter,33]+
      	                                               x[[s]][burnin:n.iter,38])
      	est_Fs_blue=x[[s]][burnin:n.iter,34]/(x[[s]][burnin:n.iter,34]+
      	                                                x[[s]][burnin:n.iter,39])
      	est_Fs_oak=x[[s]][burnin:n.iter,35]/(x[[s]][burnin:n.iter,35]+
      	                                               x[[s]][burnin:n.iter,40])
      	est_Fs_pine=x[[s]][burnin:n.iter,36]/(x[[s]][burnin:n.iter,36]+
      	                                                x[[s]][burnin:n.iter,41])
        
      	x_new=cbind(x[[s]][burnin:n.iter,1:31],
      	            est_Fs_alf,est_Fs_ash, est_Fs_blue, est_Fs_oak,est_Fs_pine, 
                    est_Flig_alf,est_Flig_ash,est_Flig_blue,est_Flig_oak, est_Flig_pine)
      burn_results[[s]]=mcmc(x_new)
      print(s)
      print(summary(burn_results[[s]]))
      print(geweke.diag(burn_results[[s]]))
    }
    return(burn_results)
    } 



MCMClist_converge_test=function(x){
    #create mcmc list (manual for now)
    combined<-mcmc.list(x[[1]], x[[2]], x[[3]])
    print(gelman.diag(combined))
    print(summary(combined))
    return(combined)
}
