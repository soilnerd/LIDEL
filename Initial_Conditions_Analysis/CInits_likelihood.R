##########LIKELIHOOD FUNCTIONS##########  
#for use with all data
#Nell Campbell
#moved to separate file 11/17/14
##################################
#likelihood functions for Latent states	
#x is vector or matrix of measured values

#likelihood function for Fraction Soluble in litter
LikeFs=function(fs, x_Hwe, x_dif, shapes, L_C6, proc_mu, L_C2, proc_mu_C2, param_update){

  #log likelihood of data given moment matched shape parameters
    Log_HWE<-dbeta(x_Hwe, shapes[1], shapes[2], log=TRUE)
    Log_ydiff<-dbeta(x_dif, shapes[3], shapes[4], log=TRUE)
    
  #Log likelihood of latent C6alpha given model result
    LogFs_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  
  #Log likelihood of latent C2alpha given model result
    LogFs_proc_C2<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  
  #total
  return(sum(Log_HWE[is.finite(Log_HWE)])+sum(Log_ydiff[is.finite(Log_ydiff)])+LogFs_proc+LogFs_proc_C2)
}

#likelihood function for lignin fraction parameter
LikeFlig<-function(x_lig, a_LIG, b_LIG, L_C2, proc_mu_C2, param_update){
  
  #log likelihood of data given moment matched shape parameters
  Log_LIG<-dbeta(x_lig, a_LIG, b_LIG, log=TRUE)
  
  #Log likelihood of latent C2alpha given model result
  LogFlig_C2<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  
  return(sum(Log_LIG[is.finite(Log_LIG)])+LogFlig_C2)
}

#likelihood function for Latent intial mass
LikeZ1_3<-function(x, Z1_3, L_C6, proc_mu, L_C2, proc_mu_C2, param_update, sig_y_1_3){
  
  #Log likelihood of dataset given latent initial litter mass
  LogZ1_3_y<-dnorm(x, mean=Z1_3, sd=sig_y_1_3, log=TRUE)
  #Log likelihood of latent C6alpha given model result
  LogZ1_3_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  #Log likelihood of latent C2alpha given model result
  LogZ1_3_proc_C2<-dnorm(L_C2, mean=proc_mu_C2, sd=param_update$sigma_L_C2, log=TRUE)
  #total
  return(sum(LogZ1_3_y[is.finite(LogZ1_3_y)])+LogZ1_3_proc+LogZ1_3_proc_C2)
}

Likesig_y_1_3<-function(x_1_3, z_1_3, sig_y_1_3){
  
  #Log likelihood of dataset given latent initial litter mass
  LogZ1_3_y_sig<-dnorm(x_1_3, mean=z_1_3, sd=sig_y_1_3, log=TRUE)

  return(sum(LogZ1_3_y_sig[is.finite(LogZ1_3_y_sig)]))
}

#likelihood function for sigma_yC6alph
Likesig_yC6alp<-function(xC6, L_C6, sigma_y_C6alph){
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_sigma_yC6=dnorm(xC6, mean=L_C6, sd=sigma_y_C6alph, log=TRUE)
  
  return(sum(LogL_sigma_yC6[is.finite(LogL_sigma_yC6)]))
}

#likelihood function for sigma_yC2alph
Likesig_yC2alp<-function(xC2, L_C2, sigma_yC2alph){
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_sigma_yC2=dnorm(xC2, mean=L_C2, sd=sigma_yC2alph, log=TRUE)
  
  return(sum(LogL_sigma_yC2[is.finite(LogL_sigma_yC2)]))
}

#likelihood function for Latent C2, structural non-soluble
LikeLc2<-function(x_C2, L_C2, sigma_yC2alph){
  #moment match for alpha and beta parameters for gamma distributions
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_C2=dnorm(x_C2, L_C2, sigma_yC2alph, log=TRUE)
  #total
  return(sum(LogL_C2[is.finite(LogL_C2)]))
}

#likelihood function for Latent DOC
LikeLc6<-function(x_C6, L_C6, sigma_yC6alph){
  #moment match for alpha and beta parameters for gamma distributions
  
  #Log likelihood of DOC data given latent C6alpha
  LogL_L=dnorm(x_C6, L_C6, sigma_yC6alph, log=TRUE)
  #total
  return(sum(LogL_L[is.finite(LogL_L)]))
}

#likelihood function for DOC fraction parameter
LikeFdOC<-function(L_C6, proc_mu, param_update){
  #moment match for gamma distribution
  
  LogFdoc_proc<-dnorm(L_C6, mean=proc_mu, sd=param_update$sigma_L_C6, log=TRUE)
  
  return(sum(LogFdoc_proc[is.finite(LogFdoc_proc)]))
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

Like_sig_HWE<-function(x_Hwe, a_all_HWE, b_all_HWE){
  
  Logsig_HWE<-dbeta(x_Hwe, a_all_HWE, b_all_HWE, log=TRUE)
  
  return(sum(Logsig_HWE[is.finite(Logsig_HWE)]))
}

Like_sig_DIF<-function(x_dif, a_all_DIF, b_all_DIF){
  
  Logsig_DIF<-dbeta(x_dif, a_all_DIF, b_all_DIF, log=TRUE)
  
  return(sum(Logsig_DIF[is.finite(Logsig_DIF)]))
}

Like_sig_LIG<-function(x_lig, a_all_LIG, b_all_LIG){
  
  Logsig_LIG<-dbeta(x_lig, a_all_LIG, b_all_LIG, log=TRUE)
  
  return(sum(Logsig_LIG[is.finite(Logsig_LIG)]))
}
