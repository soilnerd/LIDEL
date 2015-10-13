DOCmod <- function(t, y, parms){

  #set initial state values
  C1=y[1]; C2=y[2]; C3=y[3]; C4=y[4];
  C5=y[5]; C6=y[6]; C7=y[7]
  
  #set parameter values
  tau=parms[1]
  EM1=parms[2]
  Em1=parms[3]
  EM2=parms[4]
  Em2=parms[5]
  Nmid=parms[6] ###option, V1.1 & V1.2, logistic N control
  NC=parms[7] ###option, V1.1 & V1.2, logistic N control
  ###Nm=parms[6] ###option, V1.3 & V1.4, linear N control
  ###NM=parms[7] ###option, V1.3 & V1.4, linear N control
  LcM=parms[8]
  beta1=parms[9]
  beta2=parms[10]
  lambda3=parms[11]
  #(fitted parameters)
  k1=parms[12]
  k2=parms[13]
  k4=parms[14]
  beta3=parms[15]
  lambda2max=parms[16]
  Lcalp=parms[17]
  
  
  #calculated values
    Lc <- C3/(C2+C3)
    #N limitation
        if(tau>=3){
          gamma=1
        }  
        if(tau<3){
            #
            ###gamma <- Nm+((1-Nm)/NM)*tau ###option, V1.3 & V1.4, linear N control
             gamma <- 1/(1+(exp(-NC*(tau-Nmid)))) ###option, V1.1 & V1.2, logistic N control
           }
    #lignin limitation
      #limitation on C1 & C2 decay
        eps_k <- if(Lc>.7){
                exp(-3.0*.7)
              }else{
                exp(-3.0*Lc)
              }
      #limitation on microbial uptake
        eps_beta <- if(Lc>.7){
                0
              }else {
                (1-exp(-0.7*(abs(Lc-0.7)*10)))
              }
    #minimum N and lignin limitation
      mu_k <- min(gamma, eps_k)
      mu_beta <- min(gamma, eps_beta)

    #decay rates
      k3 <- if(Lc>.7){
                k2*exp(-2.1)
              }else { if(Lc<.4){
                k2*(0.2/(1+200/exp(8.15*Lc)))
                #0 ###option to more closely match moorehead
              }else{
                k2*(0.2/(1+200/exp(8.15*Lc)))
              }
              }
  
      k5 <- k2*exp(-3.0*.7)
      #k5 = k3 ###option to have microbial products decay exactly match lignin decay

  
    #DOC transfers
      #non-soluble structural (C2) transfer to DOC (C6)
        lambda1 <- if(tau<3){
                min((EM1-((EM1-Em1)/LcM)*Lc),(EM1-((EM1-Em1)/3)*tau))
                }else {
                min((EM1-((EM1-Em1)/LcM)*Lc),(EM1-((EM1-Em1)/3)*3))
                }
        #set restriction to lambda1 cannot go below 0
          if(lambda1>0){
            lambda1=lambda1
          }else{
            lambda1=0
          }
  
      #microbe biomass (C4) transfer to DOC (C6)
        lambda2<-lambda2max ###option, V1.1 & V1.3 making transfer of microbe biomass to DOC not controlled by nitrogen
        #lambda2 <- lambda2max*gamma ###option, V1.2 & V1.4, making transfer of microbe biomass to DOC controlled by nitrogen

      #soluble (C1) transfer to DOC (C6)
        lambda4 <- if(tau<3){
                  min((EM2-((EM2-Em2)/LcM)*Lc),(EM2-((EM2-Em2)/3)*tau))
                }else {
                  min((EM2-((EM2-Em2)/LcM)*Lc),(EM2-((EM2-Em2)/3)*3))
                }
        #set restriction to lambda4 cannot go below 0
          if(lambda4>0){
            lambda4=lambda4
          }else{
            lambda4=0
          }
    
    #Set increment to 0 if next step is sufficiently small, to prevent negative values
      solubleC <- if(C1-mu_k*k1*C1<=0.0000000001){0}
        else{mu_k*k1*C1}
    #Set increment to 0 if next step is sufficiently small, to prevent negative values
      nonSolStrC <- if(C2-mu_k*k2*C2<=0.0000000001){0}
        else{mu_k*k2*C2}
      ligC <- k3*C3
    #Set increment to 0 if next step is sufficiently small, to prevent negative values
      micC <- if(C4-k4*C4<=0.0000000001){0}
        else{k4*C4}
      prodC <- k5*C5

      dC1 <- -solubleC 
      dC2 <- -nonSolStrC
      dC3 <- -ligC
      dC4 <- -micC+mu_beta*beta1*(1-lambda4)*solubleC+mu_beta*beta2*(1-lambda1)*nonSolStrC
      dC5 <- -prodC + beta3*(1-lambda2)*micC
      dC6 <- lambda1*nonSolStrC+lambda3*ligC+lambda2*micC+lambda3*prodC+lambda4*solubleC
      dC7 <- (1-mu_beta*beta1)*(1-lambda4)*solubleC+(1-mu_beta*beta2)*(1-lambda1)*nonSolStrC+
              (1-lambda3)*ligC+(1-beta3)*(1-lambda2)*micC+(1-lambda3)*prodC
       
       return(list(c(dC1, dC2, dC3, dC4, dC5, dC6, dC7)))
  
}
