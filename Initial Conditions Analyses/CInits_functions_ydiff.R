############################SETUP###################
#function to set up storage for chains and name them--entries must be made in the function body
setup=function(x1,num_L, n.iter, n.chain,parameter.names,dim.x,parameter_inits,
               tune_param,support_param){
	#set up storage for chains
	x=list()
	for(i in 1:length(parameter.names)) {
		#this sets up the storage list x. Note the double brackets, which are the way we access 
      #elements in a list.
		x[[i]]=array(NA,dim=c(dim.x[i],n.iter,n.chain)) 
	}
	#assign parameter names to elements of the list
	names(x)=parameter.names
  
	#initialize parameter values
	for(j in 1:(ncol(parameter_inits))){
		for(k in 1:n.chain){
			x[[j]][1,1,k]=parameter_inits[k,j]
		}
	}

	#assign initial values for latent states
  for(p in (1:num_L)){
  for(g in (1:n.chain)){
		x[[p+ncol(parameter_inits)]][,1,g]=x1[[p]][,g]
    #####----->for across time points this would be a list of lists
		  #####----->as: for(j in (1:time)){ 
		  #####----->x[[p+ncol(parameter_inits)]][g][litter,init(n.iter),j]=x1[[L_inits]][g][litter,p,j]}
  }
  }
	#assign tuning parameters
	for(l in 1:length(tune_param)){
		x$tune[parameter.names[l]]=list(tune_param[l])
	}
	#assign support parameters
	for(m in 1:length(support_param)){
		x$support[parameter.names[m]]=list(support_param[m])
	}
return(x)
} #end of setup function

############################CHOOSE###################
#Function to choose current or proposed value

#ChooseL function has a prior that depends on the current L
  #sent to the function
chooseL=function(x, z, Like_z, Like_x, proc_mu_26, sigma_L, tune, support){
#   print("x")
#   print(x)
#   print("z")
#   print(z)
#   print("proc")
#   print(proc_mu_26)
#   print("sig")
#   print(sigma_L)
  u_gam=proc_mu_26
    sig_gam=sigma_L
    
    alph_gam=u_gam^2/sig_gam^2
    beta_gam=u_gam/sig_gam^2
    
  #likelihood + prior, for original x and new draw z
  #these are both logs so we add rather than multiply
    numerator = Like_z + dgamma(z, alph_gam, beta_gam, log=TRUE)
    denominator = Like_x + dgamma(x, alph_gam, beta_gam, log=TRUE)  
  
  #because these are logs, we exponetiate the difference to get the ratio
    q.ratio = q(theta=x,mu=z,tune=tune, type="density", support=support)/
      q(theta=z,mu=x,tune=tune, type="density", support=support)
    R = exp(numerator - denominator) * q.ratio
        ####---->should check that this is correct...seems to be
  
  #accept/reject function
    if(is.finite(R)){
      if(R > runif(1,min=0,max=1)) {
        new =z
      }
      else new = x
      return(new)
    }
    else return(x)
  }

#Choose function is generic, referring to prior list in main program
choose=function(x, z, Like_z, Like_x, param, tune, support){
  
  #likelihood + prior, for original x and new draw z
    #these are both logs so we add rather than multiply
    numerator = Like_z + prior(param,theta=z)	
    denominator = Like_x + prior(param,theta=x)	
  
  #because these are logs, we exponetiate the difference to get the ratio. 
    q.ratio = q(theta=x,mu=z,tune=tune, type="density", support=support)/
      q(theta=z,mu=x,tune=tune, type="density", support=support)
    R = exp(numerator - denominator) * q.ratio

  #accept/reject function
  	if(is.finite(R)){
      	if(R > runif(1,min=0,max=1)) {new =z
  		}
     	    else new = x
      	return(new)
      	}
      else return(x)
  }

#Choose function for soluble fraction
chooseFs_Flig=function(x, z, Like_z, Like_x, param, alp, bet, tune, support){
  
  #likelihood + prior, for original x and new draw z
  #these are both logs so we add rather than multiply
  numerator = Like_z + dbeta(z, alp, bet)  
  denominator = Like_x + dbeta(x, alp, bet)	
  
  #because these are logs, we exponetiate the difference to get the ratio. 
  q.ratio = q(theta=x,mu=z,tune=tune, type="density", support=support)/
    q(theta=z,mu=x,tune=tune, type="density", support=support)
  R = exp(numerator - denominator) * q.ratio
  
  #accept/reject function
  if(is.finite(R)){
    if(R > runif(1,min=0,max=1)) {new =z
    }
    else new = x
    return(new)
  }
  else return(x)
}

############################Q###################
#Function to generate a proposed parameter value
q=function(theta,mu,tune, type, support){
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
    		  print("beta a b fault")
    		  print("mu")
    		  print(mu)
    		  print("tune")
    		  print(sigma)
    		}
    #return density
      if (type == "density") return (dbeta(theta,a,b))		
    #select possible new value, constrained by boundary
      new_b=rbeta(1,a,b)
  		if(new_b>0.99 | new_b<0.01){
  		  while(new_b>0.99 | new_b<0.01){
  		    new_b=rbeta(1,a+1,b+1)
         }
       } 
    #return draw for new value
		  if (type == "draw") return (new_b)
	} #end of zero-one support block
}

############################RUN_MCMC###################
Run_MCMC_ydiff = function(set, num_L, num_litter, data_y, n.iter, params, param_names, burnin){
#set up list to count number of time proposal is accepted.
#set up list to track parameters in each iteration  
#For most efficient results, acceptance be about 50% for few parameters, about 25% 
  #for models with many.  These are just rules of thumb.

  #set storage for all chains of results
	final=list()

for(g in 1:n.chain){
  print("Chain:")
  print(g)
  #reset acceptance list and parameter update list
  accept=vector(mode="numeric", length=length(param_names))
  names(accept)=(param_names)
  
  param_update=list()
  for(b in 1:length(params)){
    param_update[[b]]=0
  }
  names(param_update)=(params)
  
  for(s in 2:n.iter){
                     
            #print acceptance up to end of burnin
            if(s==burnin){
              a=vector(mode="numeric", length=length(param_names))
              names(a)<-param_names
              for(o in 1:length(params)){
                a[o]=accept[o]/burnin
              }
              for(p in 1:num_L){
                a[p+length(params)]=accept[p+length(params)]/(burnin*num_litter)
              }
              print("burnin acceptance rate")
              print(a) 
              #reset accept vector to 0
              accept=vector(mode="numeric", length=length(param_names))
              names(accept)=(param_names)
            }
    #for(s in 2:2){
  	if(s%%1000==0) cat(s,"\n");flush.console() #iteration counter, prints every 1000th i
  
    #Initialize parameters specific to litter type from previous iteration accepted values
          ####---->for DOC model this would need to be 2 matrices with
            #2 loops each, stepping through initial L's (each column) to run process model
            #as well as each L where measured data exist (e.g. L[j,i], j=sample period, i=litter)
  		L_update=matrix(NA, ncol=num_L, nrow=num_litter)
  		L_update[,1]=set$Fs[,s-1,g] 
    	L_update[,2]=set$alpha[,s-1,g] 
  	  L_update[,3]=set$beta[,s-1,g] 
  	  L_update[,4]=set$zC1_3alp[,s-1,g] 
  	  L_update[,5]=set$LC6alph[,s-1,g] 
      L_update[,6]=set$Flig[,s-1,g]
      L_update[,7]=set$alphalig[,s-1,g]
      L_update[,8]=set$betalig[,s-1,g]
      L_update[,9]=set$LC2alph[,s-1,g]
        #####----->for full model, L_update=array(NA, dim=(num_litter, num_L, time)) might not need to define
        	#####----->L_update[,L1,]=set[[L_inits]][g][,L1,]
        	#####----->or L_update=set[[L_inits]][g]
  
    #Initialize parameters update vector from previous iteration accepted values
      for(k in 1:length(params)){
    		param_update[k]=set[[k]][1,s-1,g]
    	}
  	      #}###---->this bracket for testing only
    
    #Loop to update parameters specific to litter type values across litter types (i loop) and parameter (u loop)
      for(u in 1:num_L){
        for(i in 1:num_litter){
          #run model for current estimate of initial DOC
          proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], Fdoc=as.numeric(param_update[1]))
          proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
          if(u==1){ #this loop for updating Fs parameter for each litter
                  curr=L_update[i,u]
      			    #generate likelihood based on current Latent state value
      			      Like_curr = LikeFs(L_C6=L_update[i,5], proc_mu=proc_mu, param_update=param_update,
      			                         L_C2=L_update[i,9], proc_mu_C2=proc_mu_C2)
      			    #select potential new Latent state value
      				    z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                  #check that Flig+Fs!=1 or greater
                  if((z+L_update[i,6])>=1){
                    while((z+L_update[i,6])>=1){
            		      z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
            		    }
            		    }
                #update Latent state to reflect potential new parameter value
          				L_update[i,u]=z
      			    #re-run model for current estimate of initial DOC
      			      proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], Fdoc=as.numeric(param_update[1]))
      			      proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
          			#generate likelihood based on potential new parameter
          				Like_z = LikeFs(L_C6=L_update[i,5], proc_mu=proc_mu, param_update=param_update,
          				                L_C2=L_update[i,9], proc_mu_C2=proc_mu_C2)
    			      #select parameter using choose function
# ###for troubleshooting
#         print("Fs")
#         print(curr)
# 		    print(Like_curr)
#         print(z)
# 		    print(Like_z)
#         print("alpha")
#         print(L_update[i,2])
#         print("beta")
#         print(L_update[i,3])

      				    L_update[i,u]=chooseFs_Flig(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="Fs", alp=L_update[i,2], bet=L_update[i,3],
                              tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
                #update number of accepted proposals for a given parameter
          				if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
                
            } else{ if(u==2){ #this loop is for estimating alpha values
                  curr=L_update[i,u]
                #generate likelihood based on current Latent state value
                      #--->Like_curr = LikeAB(x_Hwe=data_y$y_Fs[[1]][i,], x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])
                  Like_curr = LikeAB(x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])
                #select potential new Latent state value
                  z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                #update Latent state to reflect potential new parameter value
                  L_update[i,u]=z
                #generate likelihood based on potential new parameter
                      #--->Like_z = LikeAB(x_Hwe=data_y$y_Fs[[1]][i,], x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])  
                  Like_z = LikeAB(x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])  
                #select parameter using choose function
                  L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="alpha",
                                       tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
                #update number of accepted proposals for a given parameter
                  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
                
            } else{ if(u==3){ #this loop is for estimating beta values
                  curr=L_update[i,u]
                #generate likelihood based on current Latent state value
                      #--->Like_curr = LikeAB(x_Hwe=data_y$y_Fs[[1]][i,], x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])
                  Like_curr = LikeAB(x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])
                #select potential new Latent state value
                  z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                #update Latent state to reflect potential new parameter value
                  L_update[i,u]=z
                #generate likelihood based on potential new parameter
                      #--->Like_z = LikeAB(x_Hwe=data_y$y_Fs[[1]][i,], x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])
                  Like_z = LikeAB(x_dif=data_y$y_Fs[[2]][i,], alp=L_update[i,2], bet=L_update[i,3], Fs=L_update[i,1])
                #select parameter using choose function
                  L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="beta",
                                     tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
                #update number of accepted proposals for a given parameter
                  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
              
            } else{ if(u==4){ #this loop is for initial litter mass estimate (zC1_3alp)
                  curr=L_update[i,u]
            #print("curr zC1_3")
            #print(curr)
                #re-run model for current estimate of initial DOC
                  proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], Fdoc=as.numeric(param_update[1]))
                  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
                #generate likelihood based on current Latent state value
                  Like_curr = LikeZ1_3(x=data_y[[2]][i,], Z1_3=curr, L_C6=L_update[i,5], proc_mu=proc_mu, 
                                       L_C2=L_update[i,9], proc_mu_C2=proc_mu_C2,
                                       param_update=param_update)
            #print(Like_curr)
                #select potential new Latent state value
                  z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                  #check that new draw is not <=0
                  if(z<=0){
                    while(z<=0){
                      z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                    }
                  }
            # print("z zC1_3")
            # print(z)
                #update Latent state to reflect potential new parameter value
                  L_update[i,u]=z
                #re-run model for current estimate of initial DOC
                  proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], Fdoc=as.numeric(param_update[1]))
                  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
                #generate likelihood based on potential new parameter
                  Like_z = LikeZ1_3(x=data_y[[2]][i,], Z1_3=L_update[i,4], L_C6=L_update[i,5], proc_mu=proc_mu, 
                                    L_C2=L_update[i,9], proc_mu_C2=proc_mu_C2,
                                    param_update=param_update)
            #print(Like_z)   
            #select parameter using choose function
                  L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="zC1_3alp",
                                     tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
                #update number of accepted proposals for a given parameter
                  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
              
            } else{ if(u==5){#loop for initial DOC latent state (LC6alph)
                  curr=L_update[i,u]
                #generate likelihood based on current Latent state value
                  Like_curr = LikeLc6(x=data_y[[3]][i,], L_C6=L_update[i,5], 
                                         param_update=param_update)
                #select potential new Latent state value
                  z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                #update Latent state to reflect potential new parameter value
                  L_update[i,u]=z
                #generate likelihood based on potential new parameter
                  Like_z = LikeLc6(x=data_y[[3]][i,], L_C6=L_update[i,5], 
                                   param_update=param_update)
                #re-run model for current estimate of initial DOC
                  proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], Fdoc=as.numeric(param_update[1]))
                #select parameter using choose function
#           print("LC6alph")
#           if(proc_mu<0){
#             print(L_update)
#             print(param_update)
#           }
                  L_update[i,u]=chooseL(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, proc_mu_26=proc_mu, sigma_L=param_update$sigma_L_C6,
                                       tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
                #update number of accepted proposals for a given parameter
                  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
                
          } else{ if(u==6){ #this loop is for estimating Flig parameter for each litter
                  curr=L_update[i,u]
                #run model to estimate current proc_mu_C2
                  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
                #generate likelihood based on current Latent state value
                  Like_curr = LikeFlig(L_C2=L_update[i,9], proc_mu_C2=proc_mu_C2, param_update=param_update)
                #select potential new Latent state value
                  z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                      #check that Flig+Fs!=1 or greater
                      if((z+L_update[i,1])>=1){
                        while((z+L_update[i,1])>=1){
                          z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
                        }
                      }
                #update Latent state to reflect potential new parameter value
                  L_update[i,u]=z
                #re-run model for current estimate of initial DOC
                  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
                #generate likelihood based on potential new parameter
                  Like_z = LikeFlig(L_C2=L_update[i,9], proc_mu_C2=proc_mu_C2, param_update=param_update)
                #select parameter using choose function
            # ###for troubleshooting
            #         print("Fs")
            #         print(curr)
            # 		    print(Like_curr)
            #         print(z)
            # 		    print(Like_z)
            #         print("alpha")
            #         print(L_update[i,2])
            #         print("beta")
            #         print(L_update[i,3])
                
                L_update[i,u]=chooseFs_Flig(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="Flig", alp=L_update[i,7], bet=L_update[i,8],
                                            tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
                #update number of accepted proposals for a given parameter
                if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
                
          } else{ if(u==7){ #this loop is for estimating alphalig values
              curr=L_update[i,u]
            #generate likelihood based on current Latent state value
              Like_curr = LikeABlig(x_lig=data_y[[4]][i,], alplig=L_update[i,7], betlig=L_update[i,8], Flig=L_update[i,6])
            #select potential new Latent state value
              z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
            #update Latent state to reflect potential new parameter value
             L_update[i,u]=z
            #generate likelihood based on potential new parameter
              Like_z =  LikeABlig(x_lig=data_y[[4]][i,], alplig=L_update[i,7], betlig=L_update[i,8], Flig=L_update[i,6])  
            #select parameter using choose function
              L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="alphalig",
                                   tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
            #update number of accepted proposals for a given parameter
            if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
            
          } else{ if(u==8){ #this loop is for estimating betalig values
              curr=L_update[i,u]
            #generate likelihood based on current Latent state value
              Like_curr = LikeABlig(x_lig=data_y[[4]][i,], alplig=L_update[i,7], betlig=L_update[i,8], Flig=L_update[i,6])
            #select potential new Latent state value
              z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
            #update Latent state to reflect potential new parameter value
              L_update[i,u]=z
            #generate likelihood based on potential new parameter
              Like_z =  LikeABlig(x_lig=data_y[[4]][i,], alplig=L_update[i,7], betlig=L_update[i,8], Flig=L_update[i,6])  
            #select parameter using choose function
              L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="betalig",
                                   tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
            #update number of accepted proposals for a given parameter
              if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
            
          } else{ #loop for estimating initial condition of C2, non-soluble structural
              curr=L_update[i,u]
            #generate likelihood based on current Latent state value
              Like_curr = LikeLc2(x_C2=data_y[[5]][i,], L_C2=L_update[i,9], 
                                  param_update=param_update)
            #select potential new Latent state value
              z=q(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
            #update Latent state to reflect potential new parameter value
              L_update[i,u]=z
            #generate likelihood based on potential new parameter
              Like_z = LikeLc2(x_C2=data_y[[5]][i,], L_C2=L_update[i,9], 
                               param_update=param_update)
            #re-run model for current estimate of initial DOC
              proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], zC1_3a=L_update[i,4])
            #select parameter using choose function
#       print("LC2alph")
#       if(proc_mu_C2<0){
#         print(L_update)
#       }
              L_update[i,u]=chooseL(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, proc_mu_26=proc_mu_C2, sigma_L=param_update$sigma_L_C2,
                                    tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
            #update number of accepted proposals for a given parameter
            if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
          } #end final else 
          } #end u==8 else
          } #end u==7 else
          } #end u==6 else
          } #end u==5 else
          } #end u==4 else
          } #end u==3 else
  		}#end of all else (u==2)
      }#end of i loop
      }#end of u loop
    
  	#Loop to update parameter values not specific to litter types
  		for(l in 1:length(params)){
        
    			#set initial parameter value
    			  curr=param_update[[l]]
          #run model for all litter types
            proc_all=vector(mode="numeric", length=num_litter)
    			  proc_all_C2=vector(mode="numeric", length=num_litter)
            for(i in 1:num_litter){
              proc_all[i]=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], 
                                       Fdoc=as.numeric(param_update[1]))
              proc_all_C2[i]=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,6], 
                                          zC1_3a=L_update[i,4])
            }

        if(l==1){ #for Fdoc
      			#generate likelihood based on current parameter values
      			  Like_curr = LikeFdOC(L_C6=L_update[,5], proc_mu=proc_all, param_update=param_update)
      			#select potential new parameter value
      			  z=q(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
      			#update parameters to reflect potential new parameter value
      			  param_update[[l]]=z
            #run model for new Fdoc draw
        			for(i in 1:num_litter){
          			  proc_all[i]=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,4], 
                                           Fdoc=as.numeric(param_update[1]))
          			}
      			#generate likelihood based on potential new parameter
      			  Like_z = LikeFdOC(L_C6=L_update[,5], proc_mu=proc_all, param_update=param_update)
  
          } else { if(l==2){ #for sigma_y_C1_3
              #generate likelihood based on current parameter values
                Like_curr = Likesig_yC1_3(x=data_y[[2]], zC1_3=L_update[,4], param_update=param_update)
         # print("sigma_y_C1_3")
        #  print(curr)
         # print(Like_curr)
              #select potential new parameter value
                z=q(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
        #  print("z, sigma_y_C1_3")
        #  print(z)
            #update parameters to reflect potential new parameter value
                param_update[[l]]=z
              #generate likelihood based on potential new parameter
                Like_z = Likesig_yC1_3(x=data_y[[2]], zC1_3=L_update[,4], param_update=param_update)
         #  print(Like_z) 
          } else { if(l==3){ #for sigma_y_C6alph
              #generate likelihood based on current parameter values
                Like_curr = Likesig_yC6alp(x=data_y[[3]], L_C6=L_update[,5], param_update=param_update)
              #select potential new parameter value
                z=q(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
              #update parameters to reflect potential new parameter value
                param_update[[l]]=z
              #generate likelihood based on potential new parameter
                Like_z = Likesig_yC6alp(x=data_y[[3]], L_C6=L_update[,5], param_update=param_update)
              
          } else { if(l==4){ #for sigma_L_C6alph
            #generate likelihood based on current parameter values
              Like_curr = Likesig_LC6alp(L_C6=L_update[,5], proc_mu=proc_all, param_update=param_update)
                    #print("sigma_L_C6")
                    #print(curr)
                    #print(Like_curr)
            #select potential new parameter value
              z=q(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
            #update parameters to reflect potential new parameter value
              param_update[[l]]=z
            #generate likelihood based on potential new parameter
              Like_z = Likesig_LC6alp(L_C6=L_update[,5], proc_mu=proc_all, param_update=param_update)
                    #print(z)
                    #print(Like_z)
            
          } else { if(l==5){ #for sigma_y_C2alph
            #generate likelihood based on current parameter values
              Like_curr = Likesig_yC2alp(x=data_y[[5]], L_C2=L_update[,9], param_update=param_update)
            #select potential new parameter value
              z=q(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
            #update parameters to reflect potential new parameter value
              param_update[[l]]=z
            #generate likelihood based on potential new parameter
              Like_z = Likesig_yC2alp(x=data_y[[5]], L_C2=L_update[,9], param_update=param_update)
            
          } else{ #for sigma_L_C2
            #generate likelihood based on current parameter values
              Like_curr = Likesig_LC2alp(L_C2=L_update[,9], proc_mu_C2=proc_all_C2, param_update=param_update)
                  #print("sigma_L_C2")
                  #print(curr)
                  #print(Like_curr)
            #select potential new parameter value
              z=q(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
            #update parameters to reflect potential new parameter value
              param_update[[l]]=z
            #generate likelihood based on potential new parameter
              Like_z = Likesig_LC2alp(L_C2=L_update[,9], proc_mu_C2=proc_all_C2, param_update=param_update)
                  #print(z)
                  #print(Like_z)
          } #end else
          } #end else (l==5) 
          } #end else (l==4)
          } #end else (l==3)
          }#end of if/else
        
  			#select parameter using choose function
  			  param_update[[l]] =choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param=params[l], 
  			                            tune=as.numeric(set$tune[l]), support=set$support[l])
  			#update number of accepted proposals for a given parameter
  			  if(param_update[[l]] != set[[l]][1,s-1,g]) accept[l]=accept[l]+1	
  		}#end of 'l' loop
    
    
  	#update parameter values in chain with the results of the MCMC selection
  	###################need to add update in chain for XY_update final parameter values
  	for(n in 1:length(params)){
  		set[[n]][1,s,g]=param_update[[n]]
  	}	
  	set$Fs[,s,g]=L_update[,1]
  	set$alpha[,s,g]=L_update[,2]
  	set$beta[,s,g]=L_update[,3]
    set$zC1_3alp[,s,g]=L_update[,4]
  	set$LC6alph[,s,g]=L_update[,5]
    set$Flig[,s,g]=L_update[,6]
    set$alphalig[,s,g]=L_update[,7]
    set$betalig[,s,g]=L_update[,8]
    set$LC2alph[,s,g]=L_update[,9]
  }#end of s loop of n.iter
  
  	#print acceptance rate
        print(accept)
        a=vector(mode="numeric", length=length(param_names))
        names(a)<-param_names
        for(o in 1:length(params)){
          a[o]=accept[o]/(n.iter-burnin)
        }
        for(p in 1:num_L){
          a[p+length(params)]=accept[p+length(params)]/((n.iter-burnin)*num_litter)
        }
        print("acceptance rate")
        print(a) 
 
  	#print filled in chains
    #create matrix (temporary, should be array) of all parameter results
      #####----->this should be set up in the structure of the code above...
  
  final[[g]]<-cbind(as.vector(set$Fdoc[,,g]), t(as.matrix(set$Fs[,,g])), t(as.matrix(set$Flig[,,g])),
                    t(as.matrix(set$zC1_3alp[,,g])), t(as.matrix(set$LC6alph[,,g])),
                    t(as.matrix(set$LC2alph[,,g])),
                    as.vector(set$sigma_y_C1_3[,,g]), as.vector(set$sigma_y_C6alph[,,g]), 
                    as.vector(set$sigma_L_C6[,,g]), as.vector(set$sigma_y_C2alph[,,g]),
                    as.vector(set$sigma_L_C2[,,g]), 
                    t(as.matrix(set$alpha[,,g])), t(as.matrix(set$beta[,,g])),
                    t(as.matrix(set$alphalig[,,g])), t(as.matrix(set$betalig[,,g])))    

  colnames(final[[g]])<-c(params[1], "Fs_alfalfa", "Fs_ash", "Fs_blue", "Fs_oak", "Fs_pine", 
                          "Flig_alfalfa", "Flig_ash", "Flig_blue", "Flig_oak", "Flig_pine",
                          "Init_1_3_alfalfa", "Init_1_3_ash", "Init_1_3_blue", "Init_1_3_oak", "Init_1_3_pine",
                          "Init_C6_alfalfa", "Init_C6_ash", "Init_C6_blue", "Init_C6_oak", "Init_C6_pine",
                          "Init_C2_alfalfa", "Init_C2_ash", "Init_C2_blue", "Init_C2_oak", "Init_C2_pine",
                          params[2:length(params)],
                          "alp_alfalfa", "alp_ash", "alp_blue", "alp_oak", "alp_pine",
                          "bet_alfalfa", "bet_ash", "bet_blue", "bet_oak", "bet_pine",
                          "alpL_alfalfa", "alpL_ash", "alpL_blue", "alpL_oak", "alpL_pine",
                          "betL_alfalfa", "betL_ash", "betL_blue", "betL_oak", "betL_pine")
  }
return(final)
#return(set)
}



