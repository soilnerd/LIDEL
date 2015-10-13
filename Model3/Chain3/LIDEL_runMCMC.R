############################SETUP###################
setup=function(n.iter, init_mu_input, init_proc_days, n.chain, parameter.names,dim.x,parameter_inits,
               tune_param,support_param, Latent_tune){
  #set up storage for chains
  x=list()
  for(i in 1:length(parameter.names)) {
    #this sets up the storage list x. Note the double brackets, which are the way we access 
    #elements in a list.
    if(i<=length(params)){
      x[[i]]=array(NA,dim=c(dim.x[i],n.iter,n.chain)) 
      x[[i]][,1,n.chain]=parameter_inits[[n.chain]][[i]]
    }
    if(i>5 & i<=8){
      alfalfa=array(NA, dim=c(dim.x[i], n.iter, n.chain))
      alfalfa[,1,n.chain]=parameter_inits[[n.chain]][[i]]
      ash=array(NA, dim=c(dim.x[i], n.iter, n.chain))
      ash[,1,n.chain]=parameter_inits[[n.chain]][[i]]
      bluestem=array(NA, dim=c(dim.x[i], n.iter, n.chain))
      bluestem[,1,n.chain]=parameter_inits[[n.chain]][[i]]
      oak=array(NA, dim=c(dim.x[i], n.iter, n.chain))
      oak[,1,n.chain]=parameter_inits[[n.chain]][[i]]
      pine=array(NA, dim=c(dim.x[i], n.iter, n.chain))
      pine[,1,n.chain]=parameter_inits[[n.chain]][[i]]
      x[[i]]=list(alfalfa, ash, bluestem, oak, pine)
      names(x[[i]])=litter_nam
    }
    if(i>8 & i<12){
      x[[i]]=array(NA, dim=c(dim.x[i], n.iter, n.chain))
      x[[i]][,1,n.chain]=parameter_inits[[n.chain]][[i]]
    }
    if(i==12){
      #storage for mass loss
      Alfalfa_mass=matrix(0,nrow=length(init_mu_input[[1]][[1]]), ncol=n.iter)
      Alfalfa_mass[,1]=init_mu_input[[1]][[1]]
      #storage for DOC
      Alfalfa_DOC=matrix(0,nrow=length(init_mu_input[[1]][[2]]), ncol=n.iter)
      Alfalfa_DOC[,1]=init_mu_input[[1]][[2]]
      #storage for CO2
      Alfalfa_CO2=matrix(0,nrow=length(init_mu_input[[1]][[3]]), ncol=n.iter)
      Alfalfa_CO2[,1]=init_mu_input[[1]][[3]]
      #storage for process model results, only for burnin:n.iter
      Alfalfa_proc=array(0, dim=c(dim(init_proc_days[[1]][[1]]), (n.iter-burnin)))
      #initialize with current result just to start proc_curr before s loop in MCMC
      #will be over-written at burnin = s
      Alfalfa_proc[,,1]=init_proc_days[[1]][[1]]
      x[[i]]=list(Alfalfa_mass,Alfalfa_DOC,Alfalfa_CO2, Alfalfa_proc)
      names(x[[i]])=c("Mass", "DOC", "CO2", "Process_model")
  
    }
    if(i==13){
      #storage for mass loss
      Ash_mass=matrix(0,nrow=length(init_mu_input[[2]][[1]]), ncol=n.iter)
      Ash_mass[,1]=init_mu_input[[2]][[1]]
      #storage for DOC
      Ash_DOC=matrix(0,nrow=length(init_mu_input[[2]][[2]]), ncol=n.iter)
      Ash_DOC[,1]=init_mu_input[[2]][[2]]
      #storage for CO2
      Ash_CO2=matrix(0,nrow=length(init_mu_input[[2]][[3]]), ncol=n.iter)
      Ash_CO2[,1]=init_mu_input[[2]][[3]]
      #storage for process model results, only for burnin:n.iter
      Ash_proc=array(0, dim=c(dim(init_proc_days[[1]][[2]]), (n.iter-burnin)))
      #initialize with current result just to start proc_curr before s loop in MCMC
      #will be over-written at burnin = s
      Ash_proc[,,1]=init_proc_days[[1]][[2]]
      x[[i]]=list(Ash_mass,Ash_DOC,Ash_CO2, Ash_proc)
      names(x[[i]])=c("Mass", "DOC", "CO2", "Process_model")

    }
    if(i==14){
      #storage for mass loss
      Bluestem_mass=matrix(0,nrow=length(init_mu_input[[3]][[1]]), ncol=n.iter)
      Bluestem_mass[,1]=init_mu_input[[3]][[1]]
      #storage for DOC
      Bluestem_DOC=matrix(0,nrow=length(init_mu_input[[3]][[2]]), ncol=n.iter)
      Bluestem_DOC[,1]=init_mu_input[[3]][[2]]
      #storage for CO2
      Bluestem_CO2=matrix(0,nrow=length(init_mu_input[[3]][[3]]), ncol=n.iter)
      Bluestem_CO2[,1]=init_mu_input[[3]][[3]]
      #storage for process model results, only for burnin:n.iter
      Bluestem_proc=array(0, dim=c(dim(init_proc_days[[1]][[3]]), (n.iter-burnin)))
      #initialize with current result just to start proc_curr before s loop in MCMC
      #will be over-written at burnin = s
      Bluestem_proc[,,1]=init_proc_days[[1]][[3]]
      x[[i]]=list(Bluestem_mass,Bluestem_DOC,Bluestem_CO2, Bluestem_proc)
      names(x[[i]])=c("Mass", "DOC", "CO2", "Process_model")
      
    }
    if(i==15){
      #storage for mass loss
      Oak_mass=matrix(0,nrow=length(init_mu_input[[4]][[1]]), ncol=n.iter)
      Oak_mass[,1]=init_mu_input[[4]][[1]]
      #storage for DOC
      Oak_DOC=matrix(0,nrow=length(init_mu_input[[4]][[2]]), ncol=n.iter)
      Oak_DOC[,1]=init_mu_input[[4]][[2]]
      #storage for CO2
      Oak_CO2=matrix(0,nrow=length(init_mu_input[[4]][[3]]), ncol=n.iter)
      Oak_CO2[,1]=init_mu_input[[4]][[3]]
      #storage for process model results, only for burnin:n.iter
      Oak_proc=array(0, dim=c(dim(init_proc_days[[1]][[4]]), (n.iter-burnin)))
      #initialize with current result just to start proc_curr before s loop in MCMC
      #will be over-written at burnin = s
      Oak_proc[,,1]=init_proc_days[[1]][[4]]
      x[[i]]=list(Oak_mass,Oak_DOC,Oak_CO2, Oak_proc)
      names(x[[i]])=c("Mass", "DOC", "CO2", "Process_model")
      
    }
    if(i==16){
      #storage for mass loss
      Pine_mass=matrix(0,nrow=length(init_mu_input[[5]][[1]]), ncol=n.iter)
      Pine_mass[,1]=init_mu_input[[5]][[1]]
      #storage for DOC
      Pine_DOC=matrix(0,nrow=length(init_mu_input[[5]][[2]]), ncol=n.iter)
      Pine_DOC[,1]=init_mu_input[[5]][[2]]
      #storage for CO2
      Pine_CO2=matrix(0,nrow=length(init_mu_input[[5]][[3]]), ncol=n.iter)
      Pine_CO2[,1]=init_mu_input[[5]][[3]]
      #storage for process model results, only for burnin:n.iter
      Pine_proc=array(0, dim=c(dim(init_proc_days[[1]][[5]]), (n.iter-burnin)))
      #initialize with current result just to start proc_curr before s loop in MCMC
      #will be over-written at burnin = s
      Pine_proc[,,1]=init_proc_days[[1]][[5]]
      x[[i]]=list(Pine_mass,Pine_DOC,Pine_CO2, Pine_proc)
      names(x[[i]])=c("Mass", "DOC", "CO2", "Process_model")
      
    }
  }
  #assign parameter names to elements of the list
  names(x)=parameter.names
  
  #assign tuning parameters
  for(l in 1:length(tune_param)){
    x$tune[parameter.names[l]]=list(tune_param[l])
  }
    x$tune["Latent"]=list(Latent_tune)
  #assign support parameters
  for(m in 1:length(support_param)){
    x$support[parameter.names[m]]=list(support_param[m])
  }
  return(x)
} #end of setup function

############################CHOOSE###################
#Function to choose current or proposed value
#Choose function is generic, referring to prior list in main program
choose=function(x, z, Like_z, Like_x, param, tune, support){
  if(param=="Latent"){  

    #no prior on Latent
    #these are both logs so we add rather than multiply
    numerator = Like_z  
    denominator = Like_x
    
    q.ratio = q_sel(theta=x,mu=z,tune=tune, type="density", support=support)/
      q_sel(theta=z,mu=x,tune=tune, type="density", support=support)
    
    #because these are logs, we exponetiate the difference to ge the ratio.

	  e_num_dom_L = exp(numerator - denominator)
	  
		if(e_num_dom_L==0 & q.ratio==Inf){
			b= paste(param, "curr:", x, "new:", z, "like param x", Like_x, "like param z", Like_z, "q", q.ratio, "exp", e_num_dom_L)
			write(b, "./results/R_L_Inf_fail.txt")
			R=0
		  } else {
			R = e_num_dom_L  * q.ratio 
		  }
    
	if(R=="NaN"){
			a= paste(param, "curr:", x, "new:", z, "like param x", Like_x, "like param z", Like_z, "q", q.ratio, "exp", e_num_dom_L)
			write(a, "./results/R_L_NAN_fail.txt")
			R=0
	}
	
    if(R > runif(1,min=0,max=1)){
      new_draw_L =z
    } else {
      new_draw_L = x
    }
    return(new_draw_L)
  } else	{
    # these are both logs so we add rather than multiply
    numerator = Like_z + prior(theta=z)	
    denominator = Like_x + prior(theta=x)	
  
  q.ratio = q_sel(theta=x,mu=z,tune=tune, type="density", support=support)/
    q_sel(theta=z,mu=x, tune=tune, type="density", support=support)
  
  #because these are logs, we exponetiate the difference to ge the ratio.
  e_num_dom = exp(numerator - denominator)
  if(e_num_dom==0 & q.ratio==Inf){
    a= paste(param, "curr:", x, "new:", z, "like param x", Like_x, "like param z", Like_z, "q", q.ratio, "exp", e_num_dom)
    write(a, "./results/R_proc_fail.txt")
    R_proc=0
  } else {
    R_proc = e_num_dom  * q.ratio 
  }
  
  	if(R_proc=="NaN"){
			a= paste(param, "curr:", x, "new:", z, "like param x", Like_x, "like param z", Like_z, "q", q.ratio, "exp", e_num_dom)
			write(a, "./results/R_proc_NAN_fail.txt")
			R_proc=0
	}
  
  if(R_proc > runif(1,min=0,max=1)) {
    new_draw =z
  } else {
    new_draw = x
  }
  return(new_draw)
  }
}

############################Q###################
#Function to generate a proposed parameter value
q_sel=function(theta,mu,tune, type, support){
  sigma=tune
  if(support == "non-negative"){
    if (type == "density") return (dgamma(theta,mu^2/sigma^2, mu/sigma^2))
    if (type == "draw"){
      new_g=rgamma(1,mu^2/sigma^2, mu/sigma^2)
      #if sigma is much larger than mu, new_g will always return as 0
      #redraw to center mu on uniform range between minimum - tuning parameter
      if(new_g<0.01){
        while(new_g<0.01){
          mu=runif(1, min=0.01, max=sigma)
          new_g=rgamma(1,mu^2/sigma^2, mu/sigma^2)
        }
      }
      return (new_g)
    }
  } #end of non-negative support block
  if(support == "real"){
    if (type == "density") return (dnorm(theta,mu,sigma))
    if (type == "draw") return (rnorm(1,mu,sigma))
  } #end of real support block
  if(support == "zero-one"){
    #do moment matching for beta distribution
    a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
    b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2 
    #a & b can become negative if sigma is too large relative to mu
    #notify mu and sigma values if this occurs
    if(a<0 | b<0){
      # print("beta a b fault")
      # print("mu")
      # print(mu)
      # print("tune")
      # print(sigma)
      while(a<0 | b<0){
        sigma=runif(1,min=0, max=sigma)
        a <-(mu^2-mu^3-mu*sigma^2)/sigma^2
        b <- (mu-2*mu^2+mu^3-sigma^2+mu*sigma^2)/sigma^2
      }
      # print("new tune")
      #f=paste("new sigma", sigma)
      #print(f)
    }
    #return density
    if (type == "density") return (dbeta(theta,a,b))		
    #select possible new value, constrained by boundary
    new_b=rbeta(1,a,b)
    if(new_b>0.9999 | new_b<0.0001){
      while(new_b>0.9999 | new_b<0.0001){
        new_b=rbeta(1,a+1,b+1)
      }
    } 
    #return draw for new value
    if (type == "draw") return (new_b)
  } #end of zero-one support block
}

############################RUN_MCMC###################
Run_MCMC = function(set, d_l_list, data_y, data_MCMC, n.iter, params, param_names, Init_fctn_vals,
                    LIDEL_INPUT, litter_nam, burnin, Init_Cs, chain){
  
  #set up list to count number of time proposal is accepted.
  #set up list to track parameters in each iteration  
  #For most efficient results, acceptance be about 50% for few parameters, about 25% 
  #for models with many.  These are just rules of thumb.
  
  #set storage for all chains of results
    final=list()
  #assign non-fitted model inputs to single list, for use in mapply
    Latent_support=list(set$support$Alfalfa, set$support$Ash, set$support$Bluestem, set$support$Oak, set$support$Pine)
  
    #set acceptance indices to 0
    #acceptance list
      accept=list()
      accept_latent=list()
    #acceptance for parameters
      for(i in 1:5){
        accept[[i]]=0
      }
      names(accept)=params[1:5]
    #acceptance for latent states of each litter type
      accept_latent$alfalfa=list(vector(mode="numeric", length=length(d_l_list[[1]][[1]])),
                          vector(mode="numeric", length=nrow(d_l_list[[1]][[2]])),
                          vector(mode="numeric", length=nrow(d_l_list[[1]][[3]])))
      accept_latent$ash=list(vector(mode="numeric", length=length(d_l_list[[2]][[1]])),
                          vector(mode="numeric", length=nrow(d_l_list[[2]][[2]])),
                          vector(mode="numeric", length=nrow(d_l_list[[2]][[3]])))
      accept_latent$bluestem=list(vector(mode="numeric", length=length(d_l_list[[3]][[1]])),
                          vector(mode="numeric", length=nrow(d_l_list[[3]][[2]])),
                          vector(mode="numeric", length=nrow(d_l_list[[3]][[3]])))
      accept_latent$oak=list(vector(mode="numeric", length=length(d_l_list[[4]][[1]])),
                          vector(mode="numeric", length=nrow(d_l_list[[4]][[2]])),
                          vector(mode="numeric", length=nrow(d_l_list[[4]][[3]])))
      accept_latent$pine=list(vector(mode="numeric", length=length(d_l_list[[5]][[1]])),
                          vector(mode="numeric", length=nrow(d_l_list[[5]][[2]])),
                          vector(mode="numeric", length=nrow(d_l_list[[5]][[3]])))

    #function to divide acceptance by iterations
      accept_print=function(x, iter){
        y=x/iter
        return(y)
      }
    
  #initialize full model results on daily timestep
  #organized by litter, then full results for day & all 7 carbon pools (8 columns)
  #in set function, full model results stored in array of dim=c(days, 8, burnin:n.iter)
  #initial value [,,1] will be overwritten once s loop reaches burnin and starts storing full model results until n.iter
      proc_curr=list()
    #Alfalfa
      proc_curr[[1]]=set[[length(params)+1]][[4]][,,1]
    #Ash
      proc_curr[[2]]=set[[length(params)+2]][[4]][,,1]
    #Bluestem
      proc_curr[[3]]=set[[length(params)+3]][[4]][,,1]
    #Oak
      proc_curr[[4]]=set[[length(params)+4]][[4]][,,1]
    #Pine
      proc_curr[[5]]=set[[length(params)+5]][[4]][,,1]
  
    for(s in 2:n.iter){                  

      #iteration counter, prints every 1000th i
#        if(s%%100==0) cat(s,"\n");flush.console() 
#         #ALTERNATE-RUBEL write 1000 iteration counter to text document
             if(s%%100==0){
				prog=paste("chain", chain,s, sep=" ")
				write(prog,"./results/progress3.txt")
		}
      
      #print acceptance up to end of burnin
        if(s==burnin){
          #calculate % acceptance for parameters
            a=list()
            for(i in 1:num_litter){
              a[[i]]=as.numeric(accept[[i]])/(burnin-2)
            }
            names(a)=params[1:5]
    #######not functioning need to fix
          #calculate % acceptance for latent states
            a_latent=list(lapply(accept_latent[[1]], accept_print, iter=burnin-2), 
                          lapply(accept_latent[[2]], accept_print, iter=burnin-2),
                          lapply(accept_latent[[3]], accept_print, iter=burnin-2),
                          lapply(accept_latent[[4]], accept_print, iter=burnin-2),
                          lapply(accept_latent[[5]], accept_print, iter=burnin-2))
          
#            print("burnin acceptance rate")
#            print(a) 
 #           print(a_latent)

#           #ALTERNATE-RUBEL write acceptance to text file
#             write(a, "curr_burn_accept.txt")

          #reset accept indices to 0
            accept=list()
            for(i in 1:5){
              accept[[i]]=0
            }
            names(accept)=params[1:5]
          #reset latent vectors to 0
            accept_latent$alfalfa=list(vector(mode="numeric", length=length(d_l_list[[1]][[1]])),
                                       vector(mode="numeric", length=nrow(d_l_list[[1]][[2]])),
                                       vector(mode="numeric", length=nrow(d_l_list[[1]][[3]])))
            accept_latent$ash=list(vector(mode="numeric", length=length(d_l_list[[2]][[1]])),
                                   vector(mode="numeric", length=nrow(d_l_list[[2]][[2]])),
                                   vector(mode="numeric", length=nrow(d_l_list[[2]][[3]])))
            accept_latent$bluestem=list(vector(mode="numeric", length=length(d_l_list[[3]][[1]])),
                                        vector(mode="numeric", length=nrow(d_l_list[[3]][[2]])),
                                        vector(mode="numeric", length=nrow(d_l_list[[3]][[3]])))
            accept_latent$oak=list(vector(mode="numeric", length=length(d_l_list[[4]][[1]])),
                                   vector(mode="numeric", length=nrow(d_l_list[[4]][[2]])),
                                   vector(mode="numeric", length=nrow(d_l_list[[4]][[3]])))
            accept_latent$pine=list(vector(mode="numeric", length=length(d_l_list[[5]][[1]])),
                                    vector(mode="numeric", length=nrow(d_l_list[[5]][[2]])),
                                    vector(mode="numeric", length=nrow(d_l_list[[5]][[3]])))
        } #end print of % acceptance start to burnin iterations

##########CREATE OBJECTS FOR UPDATING###########
      #Initialize estimate of true conditions specific to litter and data type from previous iteration accepted values
      #outer list organized by litter, inner lists organized by measurements of mass(1), DOC (2), and CO2 (3)
        L_update=list()
          #previous estimate of Alfalfa
          L_update[[1]]=list(set[[length(params)+1]][[1]][,s-1], set[[length(params)+1]][[2]][,s-1], 
                             set[[length(params)+1]][[3]][,s-1])
          #previous estimate of Ash
          L_update[[2]]=list(set[[length(params)+2]][[1]][,s-1], set[[length(params)+2]][[2]][,s-1], 
                             set[[length(params)+2]][[3]][,s-1])
          #previous estimate of Bluestem
          L_update[[3]]=list(set[[length(params)+3]][[1]][,s-1], set[[length(params)+3]][[2]][,s-1], 
                             set[[length(params)+3]][[3]][,s-1])
          #previous estimate of Oak
          L_update[[4]]=list(set[[length(params)+4]][[1]][,s-1], set[[length(params)+4]][[2]][,s-1], 
                             set[[length(params)+4]][[3]][,s-1])
          #previous estimate of Pine
          L_update[[5]]=list(set[[length(params)+5]][[1]][,s-1], set[[length(params)+5]][[2]][,s-1], 
                             set[[length(params)+5]][[3]][,s-1])
        
      #Initialize parameters update vector from previous iteration accepted values
        param_update=list(c(set[[1]][1,s-1,1], set[[2]][1,s-1,1], set[[3]][1,s-1,1],
                            set[[4]][1,s-1,1], set[[5]][1,s-1,1]))
      #Initialize error terms update list from previous iteration conjugate prior results
      #list within list, in order to be used properly in mapply function for latent update
        err_y_C1_5_update=list(set[[6]][[1]][,s-1,1], set[[6]][[2]][,s-1,1], set[[6]][[3]][,s-1,1], 
                                    set[[6]][[4]][,s-1,1], set[[6]][[5]][,s-1,1])
        err_y_C6_update=list(set[[7]][[1]][,s-1,1], set[[7]][[2]][,s-1,1], set[[7]][[3]][,s-1,1], 
                                  set[[7]][[4]][,s-1,1], set[[7]][[5]][,s-1,1])
        err_y_C7_update=list(set[[8]][[1]][,s-1,1], set[[8]][[2]][,s-1,1], set[[8]][[3]][,s-1,1], 
                                  set[[8]][[4]][,s-1,1], set[[8]][[5]][,s-1,1])
        err_L_C1_5_update=list(set[[9]][,s-1,1])
        err_L_C6_update=list(set[[10]][,s-1,1])
        err_L_C7_update=list(set[[11]][,s-1,1])
        err_update=list(err_y_C1_5_update,err_y_C6_update, err_y_C7_update, err_L_C1_5_update, err_L_C6_update, err_L_C7_update)
      #name components of update objects
        names(param_update[[1]])=params[1:5]

      #convert input terms to list format to work in mapply function for latent update
        LIDEL_tau=list(LIDEL_INPUT[[1]]$tau[1], LIDEL_INPUT[[1]]$tau[2], LIDEL_INPUT[[1]]$tau[3],
                       LIDEL_INPUT[[1]]$tau[4], LIDEL_INPUT[[1]]$tau[5])
        LIDEL_lcalph=list(LIDEL_INPUT[[1]]$Lcalph[1], LIDEL_INPUT[[1]]$Lcalph[2], LIDEL_INPUT[[1]]$Lcalph[3],
                       LIDEL_INPUT[[1]]$Lcalph[4], LIDEL_INPUT[[1]]$Lcalph[5])

###############UPDATE LATENT STATES FOR ALL LITTERS AND ALL TIMES##########
        new_Latent_apply=mapply(Latents_likelihood, #<--FUNCTION
                              #parameters by litter type
                                #time, %N in litter, initial LCI
                                d_l=d_l_list, 
                                #current Latent states, current full model results
                                L_curr=L_update, proc_current=proc_curr,
                                #data, tuning value for measurement times, support term
                                all_data=data_MCMC, L_tuning=set$tune$Latent, L_support=Latent_support,
                                #error terms
                                err_y_C1_5=err_y_C1_5_update,
                                err_y_C6=err_y_C6_update,
                                err_y_C7=err_y_C7_update,
                                err_L_C1_5=err_L_C1_5_update,
                                err_L_C6=err_L_C6_update,
                                err_L_C7=err_L_C7_update,
                              #universal parameters
                                #non-estimated model inputs, estimated model inputs, estimated error values
                                SIMPLIFY=FALSE)                        

# ################Latent_Likelihood function testing#############
         #TESTING: mapply function running Latents_likelihood function for all litter types
#                  ###to speed things up testing/changing internal Latents_likelihood
#                 b=1
#                 d_l=d_l_list[[b]]
#                 L_curr=L_update[[b]]
#                 proc_current=proc_curr[[b]]
#                 all_data=data_MCMC[[b]]
#                 L_tuning=set$tune$Latent[[b]]
#                 L_support=Latent_support[[b]]
#                 err_y_C1_5=err_y_C1_5_update[[b]]
#                 err_y_C6=err_y_C6_update[[b]]
#                 err_y_C7=err_y_C7_update[[b]]
#                 err_L_C1_5=err_L_C1_5_update[[b]]
#                 err_L_C6=err_L_C6_update[[b]]
#                 err_L_C7=err_L_C7_update[[b]]
# #  
#         #TESTING: to run Latent_likelihood for individual litter, to check results of mapply
#                 new_Latent=list()
#                 d=1
#                 new_Latent[[d]]=Latents_likelihood( 
#                         #parameters by litter type
#                           #time, %N in litter, initial LCI
#                           d_l=d_l_list[[d]],
#                           #current Latent states, current full model results
#                           L_curr=L_update[[d]], proc_current=proc_curr[[d]],
#                           #data, tuning value for measurement times, support term
#                           all_data=data_MCMC[[d]], L_tuning=set$tune$Latent[[d]], L_support=Latent_support[[d]],
#                         #universal parameters
#                           #non-estimated model inputs, estimated model inputs, estimated error values
#                           err_y_C1_5=err_y_C1_5_update[[d]],
#                           err_y_C6=err_y_C6_update[[d]],
#                           err_y_C7=err_y_C7_update[[d]],
#                           err_L_C1_5=err_L_C1_5_update[[d]],
#                           err_L_C6=err_L_C6_update[[d]],
#                           err_L_C7=err_L_C7_update[[d]])

      #pull new latent values from results
        L_update=list(new_Latent_apply[[1]][[1]], new_Latent_apply[[2]][[1]], new_Latent_apply[[3]][[1]], 
                          new_Latent_apply[[4]][[1]], new_Latent_apply[[5]][[1]])

      #update accept_latent for mass
        accept_latent[[1]][[1]]=accept_latent[[1]][[1]]+new_Latent_apply[[1]][[3]][[1]]
        accept_latent[[2]][[1]]=accept_latent[[2]][[1]]+new_Latent_apply[[2]][[3]][[1]]
        accept_latent[[3]][[1]]=accept_latent[[3]][[1]]+new_Latent_apply[[3]][[3]][[1]]
        accept_latent[[4]][[1]]=accept_latent[[4]][[1]]+new_Latent_apply[[4]][[3]][[1]]
        accept_latent[[5]][[1]]=accept_latent[[5]][[1]]+new_Latent_apply[[5]][[3]][[1]]
      #update accept_latent for DOC
        accept_latent[[1]][[2]]=accept_latent[[1]][[2]]+new_Latent_apply[[1]][[3]][[2]]
        accept_latent[[2]][[2]]=accept_latent[[2]][[2]]+new_Latent_apply[[2]][[3]][[2]]
        accept_latent[[3]][[2]]=accept_latent[[3]][[2]]+new_Latent_apply[[3]][[3]][[2]]
        accept_latent[[4]][[2]]=accept_latent[[4]][[2]]+new_Latent_apply[[4]][[3]][[2]]
        accept_latent[[5]][[2]]=accept_latent[[5]][[2]]+new_Latent_apply[[5]][[3]][[2]]
      #update accept_latent for CO2
        accept_latent[[1]][[3]]=accept_latent[[1]][[3]]+new_Latent_apply[[1]][[3]][[3]]
        accept_latent[[2]][[3]]=accept_latent[[2]][[3]]+new_Latent_apply[[2]][[3]][[3]]
        accept_latent[[3]][[3]]=accept_latent[[3]][[3]]+new_Latent_apply[[3]][[3]][[3]]
        accept_latent[[4]][[3]]=accept_latent[[4]][[3]]+new_Latent_apply[[4]][[3]][[3]]
        accept_latent[[5]][[3]]=accept_latent[[5]][[3]]+new_Latent_apply[[5]][[3]][[3]]

##############UPDATE MODEL PARAMETERS ACROSS ALL LITTER TYPES###########
      #process model results need to be re-run given full suite of updated latent states
        #if, for testing purposes, this function is run using s=2, initial latent states prior to updating with the latents_likelihood function
        #results will be very similar but not exactly the same as proc_curr, since numerical solver will give slightly altered results
        #when run in Parameter_new_mod
#         proc_curr_apply=mapply(Parameter_new_mod,
#                               d_l=d_l_list, LIDEL_inpt=LIDEL_INPUT, Init_vals=Init_fctn_vals, tau=LIDEL_tau, Lcalph=LIDEL_lcalph, 
#                               param_ups=param_update, L_curr=L_update, SIMPLIFY=FALSE)
#for running deterministically from initial condition
proc_curr_apply=mapply(Parameter_new_mod,
                       d_l=d_l_list, LIDEL_inpt=LIDEL_INPUT, Init_vals=Init_fctn_vals, tau=LIDEL_tau, Lcalph=LIDEL_lcalph, 
                       param_ups=param_update, SIMPLIFY=FALSE)

# ################Parameter_new_mod function testing#############
#         #TESTING: mapply function running Parameter_new_mod for all litter types
#                 ###to speed things up testing/changing internal Latents_likelihood
#                 c=1
#                 d_l=d_l_list[[c]]
#                 LIDEL_inpt=LIDEL_INPUT[[c]]
#                 LIDEL_init=LIDEL_inits[[c]]
#                 tau=LIDEL_tau[[c]]
#                 Lcalph=LIDEL_lcalph[[c]]
#                 param_ups=param_update[[c]]
#                 L_curr=L_updaten[[c]]
#                 Init_vals=Init_fctn_vals[[c]]
#               
# 
#        #TESTING: to run Parameter_new_mod for individual litter, to check results of mapply
#           proc_curr_apply=list()
#           proc_curr_apply[[1]]=Parameter_new_mod(d_l=d_l_list[[1]], LIDEL_inpt=LIDEL_INPUT[[1]], Init_vals=Init_fctn_vals[[1]],
#                                    tau=LIDEL_tau[[1]], Lcalph=LIDEL_lcalph[[1]], 
#                                    param_ups=param_update[[1]], L_curr=L_update[[1]])


    #save current full model results 
      proc_curr=list(proc_curr_apply[[1]][[1]], proc_curr_apply[[2]][[1]], proc_curr_apply[[3]][[1]],
                          proc_curr_apply[[4]][[1]], proc_curr_apply[[5]][[1]])

    #save model results that correspond to measurement intervals for latent states
      proc_Latent=list(proc_curr_apply[[1]][[2]], proc_curr_apply[[2]][[2]], proc_curr_apply[[3]][[2]],
                       proc_curr_apply[[4]][[2]], proc_curr_apply[[5]][[2]])
      #for updating parameters
      for(p in 1:length(params)){
        if(p<=5){
          #current parameter value
            curr=param_update[[1]][p]
          #generate likelihood based on current parameter values
            Like_curr = Like_proc_param(proc_Lat=proc_Latent, Latent=L_update, errs_L_C1_5=err_L_C1_5_update, 
                                        errs_L_C6=err_L_C6_update, errs_L_C7=err_L_C7_update, num_litter=num_litter)
          #select potential new parameter value
            z=q_sel(mu=curr,tune=as.numeric(set$tune[p]), support=set$support[p], type="draw")
          #check that beta 3 + lambda 2 are not greater than 1
            if(p==4){
            if((z+param_update[[1]][5])>=1){
              ztrialbet3=0
              while((z+param_update[[1]][5])>=1){
                z=q_sel(mu=curr,tune=as.numeric(set$tune[p]), support=set$support[p], type="draw")
                ztrialbet3=ztrialbet3+1
              }
              # print("ztrialbet3")
              # print(ztrialbet3)
            }
            }
            if(p==5){
              if((z+param_update[[1]][4])>=1){
                ztriallam2=0
                while((z+param_update[[1]][4])>=1){
                  z=q_sel(mu=curr,tune=as.numeric(set$tune[p]), support=set$support[p], type="draw")
                  ztriallam2=ztriallam2+1
                }
                # print("ztriallam2")
                # print(ztriallam2)
              }
            }
          #update parameters to reflect potential new parameter value
            param_update[[1]][p]=z
          #run model for current estimates
#             proc_new_apply=mapply(Parameter_new_mod,
#                                    d_l=d_l_list, LIDEL_inpt=LIDEL_INPUT, Init_vals=Init_fctn_vals, tau=LIDEL_tau, Lcalph=LIDEL_lcalph, 
#                                    param_ups=param_update, L_curr=L_update, SIMPLIFY=FALSE)
#for simulating model deterministically
          proc_new_apply=mapply(Parameter_new_mod,
                                d_l=d_l_list, LIDEL_inpt=LIDEL_INPUT, Init_vals=Init_fctn_vals, tau=LIDEL_tau, Lcalph=LIDEL_lcalph, 
                                param_ups=param_update, SIMPLIFY=FALSE)

          #save model results that correspond to measurement intervals for latent states
            proc_Latent_new=list(proc_new_apply[[1]][[2]], proc_new_apply[[2]][[2]], proc_new_apply[[3]][[2]],
                                 proc_new_apply[[4]][[2]], proc_new_apply[[5]][[2]])
          #generate likelihood based on potential new parameter
            Like_z = Like_proc_param(proc_Lat=proc_Latent_new, Latent=L_update, errs_L_C1_5=err_L_C1_5_update, 
                                     errs_L_C6=err_L_C6_update, errs_L_C7=err_L_C7_update, num_litter=num_litter)
          #select parameter using choose function
            param_update[[1]][p] =choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, 
                                      tune=as.numeric(set$tune[p]), support=set$support[p], param=params[p])
          #update number of accepted proposals and replace current init_mu and init_mu_diff with results using new parameter
            if(param_update[[1]][p] != set[[p]][1,s-1,1]) {
              accept[p]=as.numeric(accept[p])+1  
              proc_curr=list(proc_new_apply[[1]][[1]], proc_new_apply[[2]][[1]], proc_new_apply[[3]][[1]],
                                proc_new_apply[[4]][[1]], proc_new_apply[[5]][[1]])
              proc_Latent=proc_Latent_new
            }
      } else { if(p>5 & p<9){ #for sigma_y estimates, using conjugate prior
          #call matrix of alpha and beta prior values for variances
              inform_priorsAlf=conj_prior(param=params[p], litter="Alfalfa")
              inform_priorsAsh=conj_prior(param=params[p], litter="Ash")
              inform_priorsBlu=conj_prior(param=params[p], litter="Bluestem")
              inform_priorsOak=conj_prior(param=params[p], litter="Oak")
              inform_priorsPin=conj_prior(param=params[p], litter="Pine")
          #posterior alpha calculations, vector of length nrow(d_l[[s]])
            #n = 5 litter types * 3 measurements per litter for each individual s (state) and d_l[[s]] (time of measurement)
            #therefore, n=15
              alpha_sigxAlf=inform_priorsAlf[,1]+(3/2)
              alpha_sigxAsh=inform_priorsAsh[,1]+(3/2)
              alpha_sigxBlu=inform_priorsBlu[,1]+(3/2)
              alpha_sigxOak=inform_priorsOak[,1]+(3/2)
              alpha_sigxPin=inform_priorsPin[,1]+(3/2)
            
          #determine sum component of beta, difference between measurement vs true state for vector of length nrow(d_l[[p-5]])
            #apply function sums values across each row (for each d_l[[s]])
              sumBeta_sigxAlf=apply((data_MCMC[[1]][[p-5]][,3:ncol(data_MCMC[[1]][[p-5]])]-L_update[[1]][[p-5]])^2, 1, sum)
              sumBeta_sigxAsh=apply((data_MCMC[[2]][[p-5]][,3:ncol(data_MCMC[[2]][[p-5]])]-L_update[[2]][[p-5]])^2, 1, sum)
              sumBeta_sigxBlu=apply((data_MCMC[[3]][[p-5]][,3:ncol(data_MCMC[[3]][[p-5]])]-L_update[[3]][[p-5]])^2, 1, sum)
              sumBeta_sigxOak=apply((data_MCMC[[4]][[p-5]][,3:ncol(data_MCMC[[4]][[p-5]])]-L_update[[4]][[p-5]])^2, 1, sum)
              sumBeta_sigxPin=apply((data_MCMC[[5]][[p-5]][,3:ncol(data_MCMC[[5]][[p-5]])]-L_update[[5]][[p-5]])^2, 1, sum) 


        #posterior beta calculations, vector of length nrow(d_l[[s]])
            beta_sigxAlf=inform_priorsAlf[,2]+(sumBeta_sigxAlf/2)
            beta_sigxAsh=inform_priorsAsh[,2]+(sumBeta_sigxAsh/2)
            beta_sigxBlu=inform_priorsBlu[,2]+(sumBeta_sigxBlu/2)
            beta_sigxOak=inform_priorsOak[,2]+(sumBeta_sigxOak/2)
            beta_sigxPin=inform_priorsPin[,2]+(sumBeta_sigxPin/2)
        #create matrix of alpha [,1] and beta [,2] for posterior
            final_paramAlf=as.matrix(cbind(alpha_sigxAlf, beta_sigxAlf))
            final_paramAsh=as.matrix(cbind(alpha_sigxAsh, beta_sigxAsh))
            final_paramBlu=as.matrix(cbind(alpha_sigxBlu, beta_sigxBlu))
            final_paramOak=as.matrix(cbind(alpha_sigxOak, beta_sigxOak))
            final_paramPin=as.matrix(cbind(alpha_sigxPin, beta_sigxPin))
        #random draws from posterior distribution to be new sigma_y SD
            err_update[[p-5]][[1]]=apply(final_paramAlf, 1, function(x) sqrt(rigamma(1, x[1], x[2]))) 
            err_update[[p-5]][[2]]=apply(final_paramAsh, 1, function(x) sqrt(rigamma(1, x[1], x[2])))
            err_update[[p-5]][[3]]=apply(final_paramBlu, 1, function(x) sqrt(rigamma(1, x[1], x[2])))
            err_update[[p-5]][[4]]=apply(final_paramOak, 1, function(x) sqrt(rigamma(1, x[1], x[2])))
            err_update[[p-5]][[5]]=apply(final_paramPin, 1, function(x) sqrt(rigamma(1, x[1], x[2])))
          
          
        } else { if(p>=9){ #for sigma_L estimates, using conjugate prior
          #call alpha and beta uninformative prior values
            uninform_priors=conj_prior_uninf(param=params[p])
          #posterior alpha calculation
            alpha_sigL=uninform_priors[,1]+(5/2)
            
          #loop to determine sum component of beta, difference between true state vs model simulation
          #for each row of d_l[[s]]
            sumBeta_sigL1=(L_update[[1]][[p-8]]-proc_Latent[[1]][[p-8]])^2
            sumBeta_sigL2=(L_update[[2]][[p-8]]-proc_Latent[[2]][[p-8]])^2
            sumBeta_sigL3=(L_update[[3]][[p-8]]-proc_Latent[[3]][[p-8]])^2
            sumBeta_sigL4=(L_update[[4]][[p-8]]-proc_Latent[[4]][[p-8]])^2
            sumBeta_sigL5=(L_update[[5]][[p-8]]-proc_Latent[[5]][[p-8]])^2
          
            sumBeta_sigL=sumBeta_sigL1+sumBeta_sigL2+sumBeta_sigL3+sumBeta_sigL4+sumBeta_sigL5
          
          #posterior beta calculation
            beta_sigL=uninform_priors[,2]+(sumBeta_sigL/2)
          
          #create matrix of alpha [,1] and beta [,2] for posterior
            final_param_unin=as.matrix(cbind(alpha_sigL, beta_sigL))

          #random draw from posterior distribution to be new sigma_y SD
            err_update[[p-5]][[1]]=apply(final_param_unin, 1, function(x) sqrt(rigamma(1, x[1], x[2])))
          
        } #end last if 
        } #end else (l>=9)
        } #end else (5<l<9) 
      } #end of 'p' loop
      
        
        #update parameter values in chain with the results of the MCMC selection & Gibbs
          set[[1]][1,s,1]=param_update[[1]][[1]]
          set[[2]][1,s,1]=param_update[[1]][[2]]
          set[[3]][1,s,1]=param_update[[1]][[3]]
          set[[4]][1,s,1]=param_update[[1]][[4]]
          set[[5]][1,s,1]=param_update[[1]][[5]]
        #update variance values
          set[[6]][[1]][,s,1]=err_update[[1]][[1]]
          set[[6]][[2]][,s,1]=err_update[[1]][[2]]
          set[[6]][[3]][,s,1]=err_update[[1]][[3]]
          set[[6]][[4]][,s,1]=err_update[[1]][[4]]
          set[[6]][[5]][,s,1]=err_update[[1]][[5]]

          set[[7]][[1]][,s,1]=err_update[[2]][[1]]
          set[[7]][[2]][,s,1]=err_update[[2]][[2]]
          set[[7]][[3]][,s,1]=err_update[[2]][[3]]
          set[[7]][[4]][,s,1]=err_update[[2]][[4]]
          set[[7]][[5]][,s,1]=err_update[[2]][[5]]

          set[[8]][[1]][,s,1]=err_update[[3]][[1]]
          set[[8]][[2]][,s,1]=err_update[[3]][[2]]
          set[[8]][[3]][,s,1]=err_update[[3]][[3]]
          set[[8]][[4]][,s,1]=err_update[[3]][[4]]
          set[[8]][[5]][,s,1]=err_update[[3]][[5]]

          set[[9]][,s,1]=err_update[[4]][[1]]

          set[[10]][,s,1]=err_update[[5]][[1]]

          set[[11]][,s,1]=err_update[[6]][[1]]

        #update mass
          set[[length(params)+1]][[1]][,s]=L_update[[1]][[1]]
          set[[length(params)+2]][[1]][,s]=L_update[[2]][[1]]
          set[[length(params)+3]][[1]][,s]=L_update[[3]][[1]]
          set[[length(params)+4]][[1]][,s]=L_update[[4]][[1]]
          set[[length(params)+5]][[1]][,s]=L_update[[5]][[1]]
        #update DOC
          set[[length(params)+1]][[2]][,s]=L_update[[1]][[2]]
          set[[length(params)+2]][[2]][,s]=L_update[[2]][[2]]
          set[[length(params)+3]][[2]][,s]=L_update[[3]][[2]]
          set[[length(params)+4]][[2]][,s]=L_update[[4]][[2]]
          set[[length(params)+5]][[2]][,s]=L_update[[5]][[2]]
        #update CO2
          set[[length(params)+1]][[3]][,s]=L_update[[1]][[3]]
          set[[length(params)+2]][[3]][,s]=L_update[[2]][[3]]
          set[[length(params)+3]][[3]][,s]=L_update[[3]][[3]]
          set[[length(params)+4]][[3]][,s]=L_update[[4]][[3]]
          set[[length(params)+5]][[3]][,s]=L_update[[5]][[3]]
        #update proc if s>=burnin
        if(s>=burnin){
          burnin_adj=burnin+1
          set[[length(params)+1]][[4]][,,s-burnin_adj]=proc_curr[[1]]
          set[[length(params)+2]][[4]][,,s-burnin_adj]=proc_curr[[2]]
          set[[length(params)+3]][[4]][,,s-burnin_adj]=proc_curr[[3]]
          set[[length(params)+4]][[4]][,,s-burnin_adj]=proc_curr[[4]]
          set[[length(params)+5]][[4]][,,s-burnin_adj]=proc_curr[[5]]
        }

    }#end of s loop of n.iter
    
    #print acceptance rate
        a=list()
        for(i in 1:num_litter){
          a[[i]]=as.numeric(accept[[i]])/(n.iter-(burnin-1))
        }
        names(a)=params[1:5]
		a1=c(a[[1]], a[[2]], a[[3]], a[[4]], a[[5]])
        #calculate % acceptance for latent states
          a_latent2=list(lapply(accept_latent[[1]], accept_print, iter=(n.iter-(burnin-1))), 
                      lapply(accept_latent[[2]], accept_print, iter=(n.iter-(burnin-1))),
                      lapply(accept_latent[[3]], accept_print, iter=(n.iter-(burnin-1))),
                      lapply(accept_latent[[4]], accept_print, iter=(n.iter-(burnin-1))),
                      lapply(accept_latent[[5]], accept_print, iter=(n.iter-(burnin-1))))
					  
          alfalfa_accept=c("alfalfa", a_latent2[[1]][[1]], a_latent2[[1]][[2]],a_latent2[[1]][[3]])
		  ash_accept=c("ash", a_latent2[[2]][[1]], a_latent2[[2]][[2]],a_latent2[[2]][[3]])
		  bluestem_accept=c("blue", a_latent2[[3]][[1]], a_latent2[[3]][[2]],a_latent2[[3]][[3]])
		  oak_accept=c("oak", a_latent2[[4]][[1]], a_latent2[[4]][[2]], a_latent2[[4]][[3]])
		  pine_accept=c("pine", a_latent2[[5]][[1]], a_latent2[[5]][[2]], a_latent2[[5]][[3]])
		 
		  
			write(alfalfa_accept,"./results/Latent_alf_acceptance.txt")
			write(ash_accept,"./results/Latent_ash_acceptance.txt")
			write(bluestem_accept,"./results/Latent_blue_acceptance.txt")
			write(oak_accept,"./results/Latent_oak_acceptance.txt")
			write(pine_accept,"./results/Latent_pine_acceptance.txt")
		  
#        #ALTERNATE: RUBEL print to txt file
          accept_nam=paste("./results/acceptance", chain, ".txt", sep="")
          write(a1, accept_nam)
		  

    
    #pull together list of final results, organized
    final<-list()
    
    #list of all parameters estimated across litter types
      final$parameters=cbind(as.vector(set[[1]]), as.vector(set[[2]]),as.vector(set[[3]]), as.vector(set[[4]]), as.vector(set[[5]]))
        colnames(final[[1]])=params[1:5]
    #matrix of all variance terms, labeled by litter type and time point
        var_ALF=cbind(t(set[[6]][[1]][,,1]), t(set[[7]][[1]][,,1]), t(set[[8]][[1]][,,1]))
          colnames(var_ALF)=paste("Alf", c(paste("meas_m", d_l_list[[1]][[1]]), 
                                           paste("meas_d", d_l_list[[1]][[2]][,2]), 
                                           paste("meas_c", d_l_list[[1]][[3]][,2])))
                                  
        var_ASH=cbind(t(set[[6]][[2]][,,1]), t(set[[7]][[2]][,,1]), t(set[[8]][[2]][,,1]))
        
          colnames(var_ASH)=paste("Ash", c(paste("meas_m", d_l_list[[2]][[1]]), 
                                           paste("meas_d", d_l_list[[2]][[2]][,2]), 
                                           paste("meas_c", d_l_list[[2]][[3]][,2])))
        
        var_BLU=cbind(t(set[[6]][[3]][,,1]), t(set[[7]][[3]][,,1]), t(set[[8]][[3]][,,1]))
        
          colnames(var_BLU)=paste("Blue", c(paste("meas_m", d_l_list[[3]][[1]]), 
                                           paste("meas_d", d_l_list[[3]][[2]][,2]), 
                                           paste("meas_c", d_l_list[[3]][[3]][,2])))
        
        var_OAK=cbind(t(set[[6]][[4]][,,1]), t(set[[7]][[4]][,,1]), t(set[[8]][[4]][,,1]))
        
          colnames(var_OAK)=paste("Oak", c(paste("meas_m", d_l_list[[4]][[1]]), 
                                           paste("meas_d", d_l_list[[4]][[2]][,2]), 
                                           paste("meas_c", d_l_list[[4]][[3]][,2])))
        
        var_PIN=cbind(t(set[[6]][[5]][,,1]), t(set[[7]][[5]][,,1]), t(set[[8]][[5]][,,1]))
        
          colnames(var_PIN)=paste("Pine", c(paste("meas_m", d_l_list[[5]][[1]]), 
                                           paste("meas_d", d_l_list[[5]][[2]][,2]), 
                                           paste("meas_c", d_l_list[[5]][[3]][,2])))

        var_MOD=cbind(t(set[[9]][,,1]), t(set[[10]][,,1]), t(set[[11]][,,1]))

          colnames(var_MOD)=c(paste("mod_mass", d_l_list[[1]][[1]]), 
                              paste("mod_doc", d_l_list[[1]][[2]][,2]), 
                              paste("mod_co2", d_l_list[[1]][[3]][,2]))

      final$var=cbind(var_ALF, var_ASH, var_BLU, var_OAK, var_PIN, var_MOD)

    #list of measurements within litter and time
      final$Alfalfa=cbind(t(set[[12]][[1]]), t(set[[12]][[2]]), t(set[[12]][[3]]))
        colnames(final$Alfalfa)=c(paste("L Alf Ma", d_l_list[[1]][[1]]), paste("L Alf D", d_l_list[[1]][[2]][,2]), paste("L Alf CO2", d_l_list[[1]][[3]][,2]))

      final$Ash=cbind(t(set[[13]][[1]]), t(set[[13]][[2]]), t(set[[13]][[3]]))
        colnames(final$Ash)=c(paste("L Ash Ma", d_l_list[[2]][[1]]), paste("L Ash D", d_l_list[[2]][[2]][,2]), paste("L Ash CO2", d_l_list[[2]][[3]][,2]))

      final$Bluestem=cbind(t(set[[14]][[1]]), t(set[[14]][[2]]), t(set[[14]][[3]]))
        colnames(final$Bluestem)=c(paste("L Blue Ma", d_l_list[[3]][[1]]), paste("L Blue D", d_l_list[[3]][[2]][,2]), paste("L Blue CO2", d_l_list[[3]][[3]][,2]))

      final$Oak=cbind(t(set[[15]][[1]]), t(set[[15]][[2]]), t(set[[15]][[3]]))
        colnames(final$Oak)=c(paste("L Oak Ma", d_l_list[[4]][[1]]), paste("L Oak D", d_l_list[[4]][[2]][,2]), paste("L Oak CO2", d_l_list[[4]][[3]][,2]))

      final$Pine=cbind(t(set[[16]][[1]]), t(set[[16]][[2]]), t(set[[16]][[3]])) 
        colnames(final$Pine)=c(paste("L Pine Ma", d_l_list[[5]][[1]]), paste("L Pine D", d_l_list[[5]][[2]][,2]), paste("L Pine CO2", d_l_list[[5]][[3]][,2]))
        
      final$procs=list(set[[12]][[4]], set[[13]][[4]], set[[14]][[4]], set[[15]][[4]], set[[16]][[4]])
        names(final$procs)=litter_nam

  return(final)
  
} #end of Run_MCMC function