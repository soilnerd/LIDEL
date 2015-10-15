##########LIKELIHOOD FUNCTIONS##########  
#for use with all data
#Nell Campbell
#moved to separate file 11/17/14
##################################
#likelihood functions for Latent states	
#x is vector or matrix of measured values

#likelihood function for Fraction Soluble in litter
LikeFs=function(L_C6, proc_mu, L_C2, proc_mu_C2, param_update){

  #Log likelihood of latent C6alpha given model result
  LogFs_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  
  #Log likelihood of latent C2alpha given model result
  LogFs_proc_C2<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  
  #total
  return(LogFs_proc+LogFs_proc_C2)
}

#likelihood function for alpha & beta parameters
  LikeAB=function(x_dif, alp, bet, Fs){
        #--->LikeAB=function(x_Hwe, x_dif, alp, bet, Fs){
  
  #Log likelihood of both datasets given latent soluble pool
        #--->Log_HWE<-dbeta(x_Hwe, alp, bet, log=TRUE)
  Log_ydiff<-dbeta(x_dif, alp, bet, log=TRUE)
  #Log likelihood of latent C6alpha given model result
  LogFs<-dbeta(Fs, alp, bet, log=TRUE)
  #total
      #--->return(sum(Log_HWE[is.finite(Log_HWE)])+sum(Log_ydiff[is.finite(Log_ydiff)])+LogFs)
  return(sum(Log_ydiff[is.finite(Log_ydiff)])+LogFs)
}

#likelihood function for alpha & beta lignin parameters
LikeABlig=function(x_lig, alplig, betlig, Flig){
  
  #Log likelihood of dataset given lignin measurement
  Log_ylig<-dbeta(x_lig, alplig, betlig, log=TRUE)
  #Log likelihood of latent C6alpha given model result
  LogFlig<-dbeta(Flig, alplig, betlig, log=TRUE)
  #total
  return(sum(Log_ylig[is.finite(Log_ylig)])+LogFlig)
}

#likelihood function for Latent intial mass
LikeZ1_3<-function(x, Z1_3, L_C6, proc_mu, L_C2, proc_mu_C2, param_update){
  
  #Log likelihood of dataset given latent initial litter mass
  LogZ1_3_y<-dnorm(x, mean=Z1_3, sd=param_update$sigma_y_C1_3, log=TRUE)
  #Log likelihood of latent C6alpha given model result
  LogZ1_3_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  #Log likelihood of latent C2alpha given model result
  LogZ1_3_proc_C2<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  #total
  return(sum(LogZ1_3_y[is.finite(LogZ1_3_y)])+LogZ1_3_proc+LogZ1_3_proc_C2)
}

#likelihood function for Latent DOC
LikeLc6<-function(x, L_C6, param_update){
  #moment match for alpha and beta parameters for gamma distributions
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_L=dnorm(x, L_C6, param_update$sigma_y_C6alph, log=TRUE)
  #total
  return(sum(LogL_L[is.finite(LogL_L)]))
}

#likelihood function for Latent C2, structural non-soluble
LikeLc2<-function(x_C2, L_C2, param_update){
  #moment match for alpha and beta parameters for gamma distributions
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_C2=dnorm(x_C2, L_C2, param_update$sigma_y_C2alph, log=TRUE)
  #total
  return(sum(LogL_C2[is.finite(LogL_C2)]))
}

#likelihood function for DOC fraction parameter
LikeFdOC<-function(L_C6, proc_mu, param_update){
  #moment match for gamma distribution
  
  LogLc6_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  return(sum(LogLc6_proc[is.finite(LogLc6_proc)]))
}

#likelihood function for lignin fraction parameter
LikeFlig<-function(L_C2, proc_mu_C2, param_update){
  #moment match for gamma distribution
  
  #Log likelihood of latent C2alpha given model result
  LogFlig_C2<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  return(sum(LogFlig_C2[is.finite(LogFlig_C2)]))
}

#likelihood function for sigma_yC1_3
Likesig_yC1_3<-function(x, zC1_3, param_update){
  LogL_sig_yC1_3<-dnorm(x, mean=zC1_3, sd=param_update$sigma_y_C1_3, log=TRUE)
  return(sum(LogL_sig_yC1_3[is.finite(LogL_sig_yC1_3)]))
}

#likelihood function for sigma_yC6alph
Likesig_yC6alp<-function(x, L_C6, param_update){
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_sigma_yC6=dnorm(x, mean=L_C6, sd=param_update$sigma_y_C6alph, log=TRUE)
  return(sum(LogL_sigma_yC6[is.finite(LogL_sigma_yC6)]))
}

#likelihood function for sigma_yC2alph
Likesig_yC2alp<-function(x, L_C2, param_update){
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_sigma_yC2=dnorm(x, mean=L_C2, sd=param_update$sigma_y_C2alph, log=TRUE)
  return(sum(LogL_sigma_yC2[is.finite(LogL_sigma_yC2)]))
}

#likelihood function for sigma_LC6alph
Likesig_LC6alp<-function(L_C6, proc_mu, param_update){
  
  LogsigLC6_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  return(sum(LogsigLC6_proc[is.finite(LogsigLC6_proc)])) 
}

#likelihood function for sigma_LC2alph
Likesig_LC2alp<-function(L_C2, proc_mu_C2, param_update){
  
  LogsigLC2_proc<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  return(sum(LogsigLC2_proc[is.finite(LogsigLC2_proc)])) 
}
