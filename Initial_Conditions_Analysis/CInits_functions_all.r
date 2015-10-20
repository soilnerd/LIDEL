###########################
#Setup function and MCMC functions for LIDEL initial condition analysis
#Nell Campbell
#
###########################

############################SETUP###################
	#function to set up storage for chains and name them--entries must be made in the function body
		setup=function(x1,num_L, n.iter, n.chain,parameter.names,dim.x,parameter_inits,
					   tune_param,support_param){
			#set up storage for chains
				x=list()
				for(i in 1:length(parameter.names)) {
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
#Functions to choose current or proposed value

	#ChooseL function has a prior that depends on the current L
	#sent to the function
		chooseL=function(x, z, Like_z, Like_x, proc_mu_26, sigma_L, tune, support){
			#likelihood + prior, for original x and new draw z
			#these are both logs so we add rather than multiply
				numerator = Like_z + dnorm(z, proc_mu_26, sigma_L, log=TRUE)
				denominator = Like_x + dnorm(x, proc_mu_26, sigma_L, log=TRUE)  
			#because these are logs, we exponetiate the difference to get the ratio
				q.ratio = q_sel(theta=x,mu=z,tune=tune, type="density", support=support)/
				  q_sel(theta=z,mu=x,tune=tune, type="density", support=support)
				
				R = exp(numerator - denominator) * q.ratio
		  
			#accept/reject function
				if(R > runif(1,min=0,max=1)) {
					  new_draw =z
				} else {
					  new_draw = x
				}
			return(new_draw)
		} #end of chooseL function

	#Choose function is generic, referring to prior list in main program
		choose=function(x, z, Like_z, Like_x, param, tune, support){
			#likelihood + prior, for original x and new draw z
			#these are both logs so we add rather than multiply
				numerator = Like_z + prior(param,theta=z)	
				denominator = Like_x + prior(param,theta=x)	
			#because these are logs, we exponetiate the difference to get the ratio. 
				q.ratio = q_sel(theta=x,mu=z,tune=tune, type="density", support=support)/
				  q_sel(theta=z,mu=x, tune=tune, type="density", support=support)
				  
				R = exp(numerator - denominator) * q.ratio

			#accept/reject function
				if(R > runif(1,min=0,max=1)) {
					  new_draw =z
					} else {
					new_draw = x
					}
			return(new_draw)
		} #end of choose function

	#Choosei for informative priors by litter type (i) 
		choosei=function(x, z, Like_z, Like_x, param, i, tune, support){
		  
		    #likelihood + prior, for original x and new draw z
		    #these are both logs so we add rather than multiply
				numerator = Like_z + priori(param,theta=z, i)  
				denominator = Like_x + priori(param,theta=x, i)	
		  
		    #because these are logs, we exponetiate the difference to get the ratio. 
				q.ratio = q_sel(theta=x,mu=z,tune=tune, type="density", support=support)/
				  q_sel(theta=z,mu=x, tune=tune, type="density", support=support)
				  
				R = exp(numerator - denominator) * q.ratio
		  
			#accept/reject function
				if(R > runif(1,min=0,max=1)) {
						new_draw =z
					} else {
						new_draw = x
					}
			return(new_draw)
		} #end of choosei function

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
		} #end of q_sel function

############################RUN_MCMC###################
	Run_MCMC = function(set, num_L, num_litter, data_y, n.iter, params, param_names, burnin){
	###notes
		#Set up list to count number of time proposal is accepted.
		#Set up list to track parameters in each iteration  
		#For most efficient results, acceptance be about 50% for few parameters, about 25% 
		# for models with many.  

		#set storage for all chains of results
			final=list()
		#outer most loop, for chain (g)
			for(g in 1:n.chain){
			#print status
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
			#inner loop, for MCMC iteration (s)
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
						} #end of if loop for s==burnin
						
				###FOR TESTING- allows full run of single MCMC iteration
				###for(s in 2:2){
				
				#iteration counter, prints every 1000th i
					if(s%%1000==0) cat(s,"\n");flush.console() 
		  
				#Initialize parameters specific to litter type from previous iteration accepted values
					L_update=matrix(NA, ncol=num_L, nrow=num_litter)
					L_update[,1]=set$Fs[,s-1,g] 
					L_update[,2]=set$Flig[,s-1,g] 
					L_update[,3]=set$zC1_3alp[,s-1,g] 
					L_update[,4]=set$sigma_y_C1_3[,s-1,g] 
					L_update[,5]=set$sigma_y_C6alph[,s-1,g] 
					L_update[,6]=set$sigma_y_C2alph[,s-1,g]
					L_update[,7]=set$LC2alph[,s-1,g]
					L_update[,8]=set$LC6alph[,s-1,g]

				#Initialize parameters update vector from previous iteration accepted values
				  for(k in 1:length(params)){
						param_update[k]=set[[k]][1,s-1,g]
					}
				###FOR TESTING- bracket to test initialization code
				###}###
			
				#Loop to update parameters specific to litter type values across litter types (i loop) and parameter (u loop)
					  for(u in 1:num_L){
						for(i in 1:num_litter){
						  #run model for current estimate of initial DOC
							proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
							proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
						  if(u==1){ #this loop for updating Fs parameter for each litter
								  curr=L_update[i,u]
								  #moment matching for beta distribution parameters
									a_HWE <-(curr^2-curr^3-curr*param_update$sigma_HWE^2)/param_update$sigma_HWE^2
									b_HWE <- (curr-2*curr^2+curr^3-param_update$sigma_HWE^2+curr*param_update$sigma_HWE^2)/param_update$sigma_HWE^2 
								  
									a_DIF <-(curr^2-curr^3-curr*param_update$sigma_DIF^2)/param_update$sigma_DIF^2
									b_DIF <- (curr-2*curr^2+curr^3-param_update$sigma_DIF^2+curr*param_update$sigma_DIF^2)/param_update$sigma_DIF^2 
								  
									#generate likelihood based on current Latent state value
									  Like_curr = LikeFs(fs=curr, x_Hwe=data_y$y_Fs[[1]][i,], x_dif=data_y$y_Fs[[2]][i,], 
													 shapes=c(a_HWE,b_HWE,a_DIF,b_DIF),
													 L_C6=L_update[i,8], proc_mu=proc_mu, param_update=param_update,
														 L_C2=L_update[i,7], proc_mu_C2=proc_mu_C2)
									#select potential new draw for Fs
										z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
							  
								  #check that Flig+Fs!=1 or greater
									if((z+L_update[i,2])>=1){
									  while((z+L_update[i,2])>=1){
										  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
										}
										}
									a_HWE_new <-(z^2-z^3-z*param_update$sigma_HWE^2)/param_update$sigma_HWE^2
									b_HWE_new <- (z-2*z^2+z^3-param_update$sigma_HWE^2+z*param_update$sigma_HWE^2)/param_update$sigma_HWE^2 
									
									a_DIF_new <-(z^2-z^3-z*param_update$sigma_DIF^2)/param_update$sigma_DIF^2
									b_DIF_new <- (z-2*z^2+z^3-param_update$sigma_DIF^2+z*param_update$sigma_DIF^2)/param_update$sigma_DIF^2 
									
								#make sure variance term and new draw don't result in negative shape values
									if(a_HWE_new<0 | b_HWE_new<0 | a_DIF_new<0 | b_DIF_new<0){
									  while(a_HWE_new<0 | b_HWE_new<0 | a_DIF_new<0 | b_DIF_new<0){
										z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
										a_HWE_new <-(z^2-z^3-z*param_update$sigma_HWE^2)/param_update$sigma_HWE^2
										b_HWE_new <- (z-2*z^2+z^3-param_update$sigma_HWE^2+z*param_update$sigma_HWE^2)/param_update$sigma_HWE^2 
										
										a_DIF_new <-(z^2-z^3-z*param_update$sigma_DIF^2)/param_update$sigma_DIF^2
										b_DIF_new <- (z-2*z^2+z^3-param_update$sigma_DIF^2+z*param_update$sigma_DIF^2)/param_update$sigma_DIF^2
									#check that Flig+Fs!=1 or greater
										if((z+L_update[i,2])>=1){
										  while((z+L_update[i,2])>=1){
											z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
										
											a_HWE_new <-(z^2-z^3-z*param_update$sigma_HWE^2)/param_update$sigma_HWE^2
											b_HWE_new <- (z-2*z^2+z^3-param_update$sigma_HWE^2+z*param_update$sigma_HWE^2)/param_update$sigma_HWE^2 
											
											a_DIF_new <-(z^2-z^3-z*param_update$sigma_DIF^2)/param_update$sigma_DIF^2
											b_DIF_new <- (z-2*z^2+z^3-param_update$sigma_DIF^2+z*param_update$sigma_DIF^2)/param_update$sigma_DIF^2
											print("hello")
										  }
										}
									  }
									}
								#update Latent state to reflect potential new parameter value
										L_update[i,u]=z
									#re-run model for current estimate of initial DOC
										proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
										proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
									#generate likelihood based on potential new parameter
										Like_z = LikeFs(fs=z, x_Hwe=data_y$y_Fs[[1]][i,], x_dif=data_y$y_Fs[[2]][i,], 
														shapes=c(a_HWE_new,b_HWE_new,a_DIF_new,b_DIF_new),
														L_C6=L_update[i,8], proc_mu=proc_mu, param_update=param_update,
														L_C2=L_update[i,7], proc_mu_C2=proc_mu_C2)
									  #select parameter using choose function
										L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="Fs", 
											  tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
										if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
								
							} else{ if(u==2){ #this loop is for estimating Flig
								  curr=L_update[i,u]
								  #moment matching for beta distribution parameters
									a_LIG <-(curr^2-curr^3-curr*param_update$sigma_LIG^2)/param_update$sigma_LIG^2
									b_LIG <- (curr-2*curr^2+curr^3-param_update$sigma_LIG^2+curr*param_update$sigma_LIG^2)/param_update$sigma_LIG^2 
								  #run model to estimate current proc_mu_C2
									proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
								  #generate likelihood based on current Latent state value
									Like_curr = LikeFlig(x_lig=data_y[[4]][i,], a_LIG=a_LIG, b_LIG=b_LIG, L_C2=L_update[i,7], 
														 proc_mu_C2=proc_mu_C2, param_update=param_update)
								  #select potential new Latent state value
									z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								  #check that Flig+Fs!=1 or greater
									if((z+L_update[i,1])>=1){
									  while((z+L_update[i,1])>=1){
										z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
									  }
									}
								  #moment matching for beta distribution parameters
									a_LIG_new <-(z^2-z^3-z*param_update$sigma_LIG^2)/param_update$sigma_LIG^2
									b_LIG_new <- (z-2*z^2+z^3-param_update$sigma_LIG^2+z*param_update$sigma_LIG^2)/param_update$sigma_LIG^2 
								  #make sure variance term and new draw don't result in negative shape values
								  if(a_LIG_new<0 | b_LIG_new<0){
									while(a_LIG_new<0 | b_LIG_new<0){
									  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
										a_LIG_new <-(z^2-z^3-z*param_update$sigma_LIG^2)/param_update$sigma_LIG^2
										b_LIG_new <- (z-2*z^2+z^3-param_update$sigma_LIG^2+z*param_update$sigma_LIG^2)/param_update$sigma_LIG^2
									  #check that Flig+Fs!=1 or greater
									  if((z+L_update[i,1])>=1){
										while((z+L_update[i,1])>=1){
										  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
										  a_LIG_new <-(z^2-z^3-z*param_update$sigma_LIG^2)/param_update$sigma_LIG^2
										  b_LIG_new <- (z-2*z^2+z^3-param_update$sigma_LIG^2+z*param_update$sigma_LIG^2)/param_update$sigma_LIG^2
										}
									  }
									}
								  }
								  #update Latent state to reflect potential new parameter value
									L_update[i,u]=z
								  #re-run model for current estimate of initial DOC
									proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
								  #generate likelihood based on potential new parameter
									Like_z = LikeFlig(x_lig=data_y[[4]][i,], a_LIG=a_LIG_new, b_LIG=b_LIG_new, L_C2=L_update[i,7], 
													proc_mu_C2=proc_mu_C2, param_update=param_update)

									L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="Flig",
															  tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								  #update number of accepted proposals for a given parameter
									if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
								  
							} else{ if(u==3){ #this loop is for z1_3
								  curr=L_update[i,u]
								#re-run model for current estimate of initial DOC
								  proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
								  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
								#generate likelihood based on current Latent state value
								  Like_curr = LikeZ1_3(x=data_y[[2]][i,], Z1_3=curr, L_C6=L_update[i,8], proc_mu=proc_mu, 
													   L_C2=L_update[i,7], proc_mu_C2=proc_mu_C2,
													   param_update=param_update, sig_y_1_3=L_update[i,4])
								#select potential new Latent state value
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								#check that new draw is not <=0 (truncate)
								  if(z<=0){
									while(z<=0){
									  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
									}
								  }
								#update Latent state to reflect potential new parameter value
								  L_update[i,u]=z
								#re-run model for current estimate of initial DOC
								  proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
								  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
								#generate likelihood based on potential new parameter
								  Like_z = LikeZ1_3(x=data_y[[2]][i,], Z1_3=z, L_C6=L_update[i,8], proc_mu=proc_mu, 
													L_C2=L_update[i,7], proc_mu_C2=proc_mu_C2,
													param_update=param_update, sig_y_1_3=L_update[i,4])
								#select parameter using choose function
								  L_update[i,u]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="zC1_3alp",
													 tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
								  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
							  
							} else{ if(u==4){ #this loop is for sigma_y_1_3
								  curr=L_update[i,u]
								#generate likelihood based on current Latent state value
								  Like_curr = Likesig_y_1_3(x_1_3=data_y[[2]][i,], z_1_3=L_update[i,3], 
															sig_y_1_3=curr)
								#select potential new Latent state value
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								#update Latent state to reflect potential new parameter value
								  L_update[i,u]=z
								#generate likelihood based on new value
								  Like_z = Likesig_y_1_3(x_1_3=data_y[[2]][i,], z_1_3=L_update[i,3], 
															sig_y_1_3=z)
								#select parameter using choose function
								  L_update[i,u]=choosei(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="sigma_y_C1_3", i=i,
													   tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
								  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
								
							} else{ if(u==5){ #this loop is for sigma_y_C6alph
								  curr=L_update[i,u]
								#generate likelihood based on current Latent state value
								  Like_curr = Likesig_yC6alp(xC6=data_y[[3]][i,], L_C6=L_update[i,8], 
															sigma_y_C6alph=curr)
								#select potential new Latent state value
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								#update Latent state to reflect potential new parameter value
								  L_update[i,u]=z
								#generate likelihood based on new value
								  Like_z = Likesig_yC6alp(xC6=data_y[[3]][i,], L_C6=L_update[i,8], 
														  sigma_y_C6alph=z)
								#select parameter using choose function
								  L_update[i,u]=choosei(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="sigma_y_C6alph", i=i,
														tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
								  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
								  
						  } else{ if(u==6){ #this loop is for sigma_y_C2alph
								  curr=L_update[i,u]
								#generate likelihood based on current Latent state value
								  Like_curr = Likesig_yC2alp(xC2=data_y[[5]][i,], L_C2=L_update[i,7], 
															 sigma_yC2alph=curr)
								#select potential new Latent state value
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								#update Latent state to reflect potential new parameter value
								  L_update[i,u]=z
								#generate likelihood based on new value
								  Like_z = Likesig_yC2alp(xC2=data_y[[5]][i,], L_C2=L_update[i,7], 
														  sigma_yC2alph=z)
								#select parameter using choose function
								  L_update[i,u]=choosei(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, param="sigma_y_C2alph", i=i,
														tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
								  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
								
						  } else{ if(u==7){ #this loop is for estimating initial condition of C2, non-soluble structural
								  curr=L_update[i,u]
								#generate likelihood based on current Latent state value
								  Like_curr = LikeLc2(x_C2=data_y[[5]][i,], L_C2=curr, 
													  sigma_yC2alph=L_update[i,6])
								#select potential new Latent state value
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								#update Latent state to reflect potential new parameter value
								  L_update[i,u]=z
								#generate likelihood based on potential new parameter
								  Like_z = LikeLc2(x_C2=data_y[[5]][i,], L_C2=z, 
												   sigma_yC2alph=L_update[i,6])
								#re-run model for current estimate of initial DOC
								  proc_mu_C2=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
								#select parameter using choose function
								  L_update[i,u]=chooseL(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, proc_mu_26=proc_mu_C2, sigma_L=param_update$sigma_L_C2,
														tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
								  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1
							
						  } else{ if(u==8){ #this loop is for estimating initial condition of C6, DOC
								  curr=L_update[i,u]
								#generate likelihood based on current Latent state value
								  Like_curr = LikeLc6(x_C6=data_y[[3]][i,], L_C6=curr, 
													  sigma_yC6alph=L_update[i,5])
								#select potential new Latent state value
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u], type="draw")
								#update Latent state to reflect potential new parameter value
								  L_update[i,u]=z
								#generate likelihood based on potential new parameter
								  Like_z = LikeLc6(x_C6=data_y[[3]][i,], L_C6=z, 
												   sigma_yC6alph=L_update[i,5])
								#re-run model for current estimate of initial DOC
								  proc_mu=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
								#select parameter using choose function
								  L_update[i,u]=chooseL(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, proc_mu_26=proc_mu, sigma_L=param_update$sigma_L_C6,
														tune=as.numeric(set$tune[length(params)+u]), support=set$support[length(params)+u])
								#update number of accepted proposals for a given parameter
								  if(L_update[i,u] != set[[length(params)+u]][i,s-1,g]) accept[length(params)+u]=accept[length(params)+u]+1

						  } 
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
								  proc_all[i]=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
								  proc_all_C2[i]=proc_C2alpha(Fs=L_update[i,1],Flig=L_update[i,2], zC1_3a=L_update[i,3])
								}

						if(l==1){ #for Fdoc
							  #generate likelihood based on current parameter values
								Like_curr = LikeFdOC(L_C6=L_update[,8], proc_mu=proc_all, param_update=param_update)
							  #select potential new parameter value
								z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
							  #update parameters to reflect potential new parameter value
								param_update[[l]]=z
							  #run model for new Fdoc draw
									for(i in 1:num_litter){
									  proc_all[i]=proc_C6alpha(Fs=L_update[i,1], zC1_3a=L_update[i,3], Fdoc=as.numeric(param_update[1]))
									}
							  #generate likelihood based on potential new parameter
								Like_z = LikeFdOC(L_C6=L_update[,8], proc_mu=proc_all, param_update=param_update)
				  
						  } else { if(l==2){ #for sigma_L_C6alph
							  #generate likelihood based on current parameter values
								Like_curr = Likesig_LC6alp(L_C6=L_update[,8], proc_mu=proc_all, param_update=param_update)
							  #select potential new parameter value
								z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
							  #update parameters to reflect potential new parameter value
								param_update[[l]]=z
							  #generate likelihood based on potential new parameter
								Like_z = Likesig_LC6alp(L_C6=L_update[,8], proc_mu=proc_all, param_update=param_update)

						  } else { if(l==3){ #for sigma_L_C2
							  #generate likelihood based on current parameter values
								Like_curr = Likesig_LC2alp(L_C2=L_update[,7], proc_mu_C2=proc_all_C2, param_update=param_update)
							  #select potential new parameter value
								z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
							  #update parameters to reflect potential new parameter value
								param_update[[l]]=z
							  #generate likelihood based on potential new parameter
								Like_z = Likesig_LC2alp(L_C2=L_update[,7], proc_mu_C2=proc_all_C2, param_update=param_update)

						  } else { if(l==4){ #for sigma_HWE
							  #moment match shape parameters for beta distribution
								a_all_HWE <-(L_update[,1]^2-L_update[,1]^3-L_update[,1]*curr^2)/curr^2
								b_all_HWE <- (L_update[,1]-2*L_update[,1]^2+L_update[,1]^3-curr^2+
											L_update[,1]*curr^2)/curr^2 
							  #generate likelihood based on current parameter values
								Like_curr = Like_sig_HWE(x_Hwe=data_y$y_Fs[[1]], a_all_HWE=a_all_HWE, b_all_HWE=b_all_HWE)
							  #select potential new parameter value
								z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
							  #moment match shape parameters for beta distribution
								a_all_HWE_new <-(L_update[,1]^2-L_update[,1]^3-L_update[,1]*z^2)/z^2
								b_all_HWE_new <- (L_update[,1]-2*L_update[,1]^2+L_update[,1]^3-z^2+
												L_update[,1]*z^2)/z^2 
							  #check that there isn't a problem with moment matching
								check=sum(ifelse(a_all_HWE_new<0 | b_all_HWE_new<0, 1, 0))
								  while(check!=0){
									z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
									a_all_HWE_new <-(L_update[,1]^2-L_update[,1]^3-L_update[,1]*z^2)/z^2
									b_all_HWE_new <- (L_update[,1]-2*L_update[,1]^2+L_update[,1]^3-z^2+
														L_update[,1]*z^2)/z^2 
									check=sum(ifelse(a_all_HWE_new<0 | b_all_HWE_new<0, 1, 0))
								  }
							  #update parameters to reflect potential new parameter value
								param_update[[l]]=z
							  #generate likelihood based on potential new parameter
								Like_z = Like_sig_HWE(x_Hwe=data_y$y_Fs[[1]], a_all_HWE=a_all_HWE_new, b_all_HWE=b_all_HWE_new)
							
						  } else { if(l==5){ #for sigma_Dif
							  #moment match shape parameters for beta distribution
								a_all_DIF <-(L_update[,1]^2-L_update[,1]^3-L_update[,1]*curr^2)/curr^2
								b_all_DIF <- (L_update[,1]-2*L_update[,1]^2+L_update[,1]^3-curr^2+
												L_update[,1]*curr^2)/curr^2 
							  #generate likelihood based on current parameter values
								Like_curr = Like_sig_DIF(x_dif=data_y$y_Fs[[2]], a_all_DIF=a_all_DIF, b_all_DIF=b_all_DIF)
							  #select potential new parameter value
								z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
							  #moment match shape parameters for beta distribution
								a_all_DIF_new <-(L_update[,1]^2-L_update[,1]^3-L_update[,1]*z^2)/z^2
								b_all_DIF_new <- (L_update[,1]-2*L_update[,1]^2+L_update[,1]^3-z^2+
												L_update[,1]*z^2)/z^2 
							  #check that there isn't a problem with moment matching
							  check_dif=sum(ifelse(a_all_DIF_new<0 | b_all_DIF_new<0, 1, 0))
								while(check_dif!=0){
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
								  a_all_DIF_new <-(L_update[,1]^2-L_update[,1]^3-L_update[,1]*z^2)/z^2
								  b_all_DIF_new <- (L_update[,1]-2*L_update[,1]^2+L_update[,1]^3-z^2+
													  L_update[,1]*z^2)/z^2 
								  check_dif=sum(ifelse(a_all_DIF_new<0 | b_all_DIF_new<0, 1, 0))
								}
							  #update parameters to reflect potential new parameter value
								param_update[[l]]=z
							  #generate likelihood based on potential new parameter
								Like_z = Like_sig_DIF(x_dif=data_y$y_Fs[[2]], a_all_DIF=a_all_DIF_new, b_all_DIF=b_all_DIF_new)

						  } else{ #for sigma_LIG
							  #moment match shape parameters for beta distribution
								a_all_LIG <-(L_update[,2]^2-L_update[,2]^3-L_update[,2]*curr^2)/curr^2
								b_all_LIG <- (L_update[,2]-2*L_update[,2]^2+L_update[,2]^3-curr^2+
												L_update[,2]*curr^2)/curr^2 
							  #generate likelihood based on current parameter values
								Like_curr = Like_sig_LIG(x_lig=data_y[[4]], a_all_LIG=a_all_LIG, b_all_LIG=b_all_LIG)
							  #select potential new parameter value
								z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
							  #moment match shape parameters for beta distribution
								a_all_LIG_new <-(L_update[,2]^2-L_update[,2]^3-L_update[,2]*z^2)/z^2
								b_all_LIG_new <- (L_update[,2]-2*L_update[,2]^2+L_update[,2]^3-z^2+
												L_update[,2]*z^2)/z^2 
							  #check that there isn't a problem with moment matching
								check_lig=sum(ifelse(a_all_LIG_new<0 | b_all_LIG_new<0, 1, 0))
								while(check_lig!=0){
								  z=q_sel(mu=curr,tune=as.numeric(set$tune[l]), support=set$support[l], type="draw")
								  a_all_LIG_new <-(L_update[,2]^2-L_update[,2]^3-L_update[,2]*z^2)/z^2
								  b_all_LIG_new <- (L_update[,2]-2*L_update[,2]^2+L_update[,2]^3-z^2+
													  L_update[,2]*z^2)/z^2 
								  check_lig=sum(ifelse(a_all_LIG_new<0 | b_all_LIG_new<0, 1, 0))
								}
							  #update parameters to reflect potential new parameter value
								param_update[[l]]=z
							  #generate likelihood based on potential new parameter
								Like_z = Like_sig_LIG(x_lig=data_y[[4]], a_all_LIG=a_all_LIG_new, b_all_LIG=b_all_LIG_new)
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
					for(n in 1:length(params)){
						set[[n]][1,s,g]=param_update[[n]]
					}	
					
					set$Fs[,s,g]=L_update[,1]
					set$Flig[,s,g]=L_update[,2]
					set$zC1_3alp[,s,g]=L_update[,3]
					set$sigma_y_C1_3[,s,g]=L_update[,4]
					set$sigma_y_C6alph[,s,g]=L_update[,5]
					set$sigma_y_C2alph[,s,g]=L_update[,6]
					set$LC2alph[,s,g]=L_update[,7]
					set$LC6alph[,s,g]=L_update[,8]

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
			  final[[g]]=list()
			  final[[g]][[1]]<-cbind(as.vector(set$Fdoc[,,g]), t(as.matrix(set$Fs[,,g])), t(as.matrix(set$Flig[,,g])),
								t(as.matrix(set$zC1_3alp[,,g])), t(as.matrix(set$LC6alph[,,g])),
								t(as.matrix(set$LC2alph[,,g])))
								
			  colnames(final[[g]][[1]])<-c(params[1], "Fs_alfalfa", "Fs_ash", "Fs_blue", "Fs_oak", "Fs_pine", 
									  "Flig_alfalfa", "Flig_ash", "Flig_blue", "Flig_oak", "Flig_pine",
									  "Init_1_3_alfalfa", "Init_1_3_ash", "Init_1_3_blue", "Init_1_3_oak", "Init_1_3_pine",
									  "Init_C6_alfalfa", "Init_C6_ash", "Init_C6_blue", "Init_C6_oak", "Init_C6_pine",
									  "Init_C2_alfalfa", "Init_C2_ash", "Init_C2_blue", "Init_C2_oak", "Init_C2_pine")
			  final[[g]][[2]]<-cbind(t(as.matrix(set$sigma_y_C1_3[,,g])),t(as.matrix(set$sigma_y_C2alph[,,g])),
									 t(as.matrix(set$sigma_y_C6alph[,,g])), 
									 as.vector(set$sigma_L_C6[,,g]), 
									 as.vector(set$sigma_L_C2[,,g]), 
									 as.vector(set$sigma_HWE[,,g]), 
									 as.vector(set$sigma_DIF[,,g]),
									 as.vector(set$sigma_LIG[,,g]))
			  colnames(final[[g]][[2]])=c("sigy1_3_alfalfa", "sigy1_3_ash", "sigy1_3_blue", "sigy1_3_oak", "sigy1_3_pine",
									  "sigy2_alfalfa", "sigy2_ash", "sigy2_blue", "sigy2_oak", "sigy2_pine",
									  "sigy6_alfalfa", "sigy6_ash", "sigy6_blue", "sigy6_oak", "sigy6_pine",
									  params[2:length(params)])
		  } # end of outer chain loop (g)
	#return results
		return(final)
	} #end of MCMC function




