##########SETUP FOR PRIORS##########
#priors for all initial values
#Nell Campbell
#updated 11/19/14
###################################

prior=function(param,theta){  
  if(param=="Fs") { #print("Fs")
    #print(theta)
    return(dbeta(theta, 1, 1, log=TRUE))
  }
  if(param=="Flig") { #print("Flig")
    #print(theta)
    return(dbeta(theta, 1, 1, log=TRUE))
  }
  if(param=="zC1_3alp") { #print("zC1_3alp")
    #print(theta)
    return(dunif(theta,min=1, max=1500000, log=TRUE))
  }
  if(param == "Fdoc") {
    return(dbeta(theta, 1, 1, log=TRUE))
  }
  if(param == "sigma_L_C6") { #print("sigma_L_C6")
    #print(theta)
    return(dunif(theta, min=1, max=50000, log=TRUE))
  }
  if(param == "sigma_L_C2") { #print("sigma_L_C2")
    #print(theta)
    return(dunif(theta, min=1, max=100000, log=TRUE))
  }
  if(param == "sigma_HWE") { #print("sigma_L_C2")
    #print(theta)
    return(dbeta(theta, 1, 1, log=TRUE))
  }
  if(param == "sigma_DIF") { #print("sigma_L_C2")
    #print(theta)
    return(dbeta(theta, 1, 1, log=TRUE))
  }
  if(param == "sigma_LIG") { #print("sigma_L_C2")
    #print(theta)
    return(dbeta(theta, 1, 1, log=TRUE))
  }
}

y_c1_3_priors=read.csv("C1_3alpha_priors.csv", header=FALSE, sep=",", stringsAsFactors=FALSE)
y_c2_priors=read.csv("C2alpha_priors.csv", header=FALSE, sep=",", stringsAsFactors=FALSE)
y_c6_priors=read.csv("C6alpha_priors.csv", header=FALSE, sep=",", stringsAsFactors=FALSE)

priori=function(param,theta, i){ 
#   if(param=="sigma_y_C1_3"){ #print("sigma_y_C1_3")
#     #print(theta)
#     x=densigamma(theta,y_c1_3_priors[1,i], y_c1_3_priors[2,i])
#     if(x==0){
#       return(prior=-100000)
#       print("sigma_y_C1_3 is 0")
#     } else {
#     return(log(x))
#     }
  if(param=="sigma_y_C1_3"){ #print("sigma_y_C1_3")
    #print(theta)
    return(dunif(theta,min=1, max=55000, log=TRUE))
    }
  if(param=="sigma_y_C6alph"){
    x2=densigamma(theta,y_c6_priors[1,i], y_c6_priors[2,i])
    if(x2==0){
      return(prior=-100000)
      print("sigma_y_c6alph is 0")
    } else {
      return(log(x2))
  } 
  }
  if(param=="sigma_y_C2alph"){ #print("sigma_y_C2alpha")
    #print(theta)
    x3=densigamma(theta,y_c2_priors[1,i], y_c2_priors[2,i])
    if(x3==0){
      return(prior=-100000)
      print("sigma_y_C2alph is 0")
    } else {
      return(log(x3))
  }
}
}
