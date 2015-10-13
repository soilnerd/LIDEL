##########SETUP FOR PRIORS##########
#priors for LIDEL model parameters
#Nell Campbell
#updated 12/8/14
###################################
#read in prior shape and rate parameters
prior_params=read.csv("data_PRIORS.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
    CM_prior_alf=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="CM" & prior_params[,2]=="Alfalfa"))
    CM_prior_ash=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="CM" & prior_params[,2]=="Ash"))
    CM_prior_blue=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="CM" & prior_params[,2]=="Bluestem"))
    CM_prior_oak=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="CM" & prior_params[,2]=="Oak"))
    CM_prior_pine=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="CM" & prior_params[,2]=="Pine"))

    C6_prior_alf=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C6" & prior_params[,2]=="Alfalfa"))
    C6_prior_ash=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C6" & prior_params[,2]=="Ash"))
    C6_prior_blue=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C6" & prior_params[,2]=="Bluestem"))
    C6_prior_oak=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C6" & prior_params[,2]=="Oak"))
    C6_prior_pine=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C6" & prior_params[,2]=="Pine"))
    
    C7_prior_alf=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C7" & prior_params[,2]=="Alfalfa"))
    C7_prior_ash=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C7" & prior_params[,2]=="Ash"))
    C7_prior_blue=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C7" & prior_params[,2]=="Bluestem"))
    C7_prior_oak=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C7" & prior_params[,2]=="Oak"))
    C7_prior_pine=as.matrix(subset(prior_params[,3:5], prior_params[,1]=="C7" & prior_params[,2]=="Pine"))

conj_prior=function(param, litter){
  if(param=="sigma_y_C1_5" & litter=="Alfalfa"){ #print("sigma_yC1_5")
    #print(theta)
    return(CM_prior_alf[,2:3])
  }
  if(param=="sigma_y_C1_5" & litter=="Ash"){ #print("sigma_yC1_5")
    #print(theta)
    return(CM_prior_ash[,2:3])
  }
  if(param=="sigma_y_C1_5" & litter=="Bluestem"){ #print("sigma_yC1_5")
    #print(theta)
    return(CM_prior_blue[,2:3])
  }
  if(param=="sigma_y_C1_5" & litter=="Oak"){ #print("sigma_yC1_5")
    #print(theta)
    return(CM_prior_oak[,2:3])
  }
  if(param=="sigma_y_C1_5" & litter=="Pine"){ #print("sigma_yC1_5")
    #print(theta)
    return(CM_prior_pine[,2:3])
  }
  
  if(param=="sigma_y_C6" & litter=="Alfalfa"){ #print("sigma_yC1_5")
    #print(theta)
    return(C6_prior_alf[,2:3])
  }
  if(param=="sigma_y_C6" & litter=="Ash"){ #print("sigma_yC1_5")
    #print(theta)
    return(C6_prior_ash[,2:3])
  }
  if(param=="sigma_y_C6" & litter=="Bluestem"){ #print("sigma_yC1_5")
    #print(theta)
    return(C6_prior_blue[,2:3])
  }
  if(param=="sigma_y_C6" & litter=="Oak"){ #print("sigma_yC1_5")
    #print(theta)
    return(C6_prior_oak[,2:3])
  }
  if(param=="sigma_y_C6" & litter=="Pine"){ #print("sigma_yC1_5")
    #print(theta)
    return(C6_prior_pine[,2:3])
  }
  
  if(param=="sigma_y_C7" & litter=="Alfalfa"){ #print("sigma_yC1_5")
    #print(theta)
    return(C7_prior_alf[,2:3])
  }
  if(param=="sigma_y_C7" & litter=="Ash"){ #print("sigma_yC1_5")
    #print(theta)
    return(C7_prior_ash[,2:3])
  }
  if(param=="sigma_y_C7" & litter=="Bluestem"){ #print("sigma_yC1_5")
    #print(theta)
    return(C7_prior_blue[,2:3])
  }
  if(param=="sigma_y_C7" & litter=="Oak"){ #print("sigma_yC1_5")
    #print(theta)
    return(C7_prior_oak[,2:3])
  }
  if(param=="sigma_y_C7" & litter=="Pine"){ #print("sigma_yC1_5")
    #print(theta)
    return(C7_prior_pine[,2:3])
  }
}

conj_prior_uninf=function(param){
  if(param=="sigma_L_C1_5"){
    return(matrix(0.001, ncol=2, nrow=length(d_l_list[[1]][[1]])))
  }
  if(param=="sigma_L_C6"){
    return(matrix(0.001, ncol=2, nrow=nrow(d_l_list[[1]][[2]])))
  }
  if(param=="sigma_L_C7"){
    return(matrix(0.001, ncol=2, nrow=nrow(d_l_list[[1]][[3]])))
  }
}

prior=function(theta){  
    return(dbeta(theta,1, 1, log=TRUE))
}
