##########LIKELIHOOD FUNCTIONS##########  
#for use with all data
#Nell Campbell
#moved to separate file 11/17/14
#updated to consider time steps independently 2/28/15
##################################
#likelihood function for latent states
#function to step through time points, update latent states 2:L, and model results
#created 2/12/15
#Nell Campbell

Latents_likelihood <- function(d_l, L_curr, err_y_C1_5, err_y_C6, err_y_C7,
                               err_L_C1_5, err_L_C6, err_L_C7, proc_current, all_data, L_tuning, L_support){
  
  #set final time index, L, for each measurement type (1=mass, 2=DOC, 3=cO2)
    L1=length(d_l[[1]])
    L2=nrow(d_l[[2]])
    L3=nrow(d_l[[3]])
  
  #set error terms, as vectors of length(d_l[[s]])
    Err_yC15=err_y_C1_5
    Err_yC6=err_y_C6
    Err_yC7=err_y_C7
    Err_LC15=err_L_C1_5
    Err_LC6=err_L_C6
    Err_LC7=err_L_C7

  #list to track acceptance rate of latent states
    acceptL=list(vector(mode="numeric", length=length(d_l[[1]])),
                          vector(mode="numeric", length=nrow(d_l[[2]])),
                          vector(mode="numeric", length=nrow(d_l[[3]])))
  
  #GIBBS STEP: conjugate prior to determine posterior for mass measurements, normal/normal (see pg 180 in Bayesian models for Ecology)
    #Mass time indexed in d_l[[1]]. Only date of measurement because it is in same increment as model result (mass remaining)
    #all d_l[[1]] values +1 to accomodate row shift in proc_current with day 0
      #for l=1 (day 95), calculate mean and sd of conditional posterior distribution 
        mu1=(
          (sum(proc_current[d_l[[1]][1]+1,2:6])/(Err_LC15[1]^2))+
             (sum(all_data[[1]][1,3:5])/(Err_yC15[1]^2))
          )/(
            (1/(Err_LC15[1]^2))+
              (ncol(all_data[[1]][,3:5])/(Err_yC15[1]^2))
             )
        sigZ1=sqrt(1/(
          (1/Err_LC15[1]^2)+
            (ncol(all_data[[1]][,3:5])/(Err_yC15[1]^2))
                      )
                   )
      #update L_curr, d_l[[1]][1]
        L_curr[[1]][1]=rnorm(1,mu1,sigZ1)
      #for l=2 (day 365), calculate mean and sd of conditional posterior distribution
        mu2=(
          (sum(proc_current[d_l[[1]][2]+1,2:6])/((Err_LC15[2])^2))+
            (sum(all_data[[1]][2,3:5])/((Err_yC15[2])^2))
          )/(
            (1/(Err_LC15[2]^2))+
              (ncol(all_data[[1]][,3:5])/(Err_yC15[2]^2))
            )
        sigZ2=sqrt(1/(
          (1/Err_LC15[2]^2)+
            (ncol(all_data[[1]][,3:5])/(Err_yC15[2]^2))
                      )
                  )
      #update L_curr, d_l[[1]][2]
        L_curr[[1]][2]=rnorm(1,mu2,sigZ2)
      #100% acceptance
        acceptL[[1]]=c(1,1)
  
  #METROPOLIS HASTINGS: loop to update DOC latent states
    for(l in 1:L2){
      if(l<L2){
        #current value of latent DOC at day index 'l'
          curr=L_curr[[2]][l]
        #generate likelihood based on current parameter values
                      #subtract Final_day (d_l[[2]][l,2]) DOC from Initial_day DOC (d_[[2]][l,1]) to calculate modeled DOC change for that time interval
                      #row reference for proc_current is + 1 to accomodate day 0
          Like_curr = dnorm(curr, mean=proc_current[d_l[[2]][l,2]+1,7]-proc_current[d_l[[2]][l,1]+1,7], sd=Err_LC6[l], log=TRUE)+
                      #l+1 latent state, given subsequent model prediction, Final_day (d_l[[2]][l+1,2]) minus Initial_day (d_l[[2]][l+1,1]))
                      #sd at l+1 Err_Lc6
                      dnorm(L_curr[[2]][l+1], mean=proc_current[d_l[[2]][l+1,2]+1,7]-proc_current[d_l[[2]][l+1,1]+1,7], sd=Err_LC6[l+1], log=TRUE)+
                      #all measurements at l (all_data[[2]][l,3:5])
                      sum(dnorm(all_data[[2]][l,3:ncol(all_data[[2]])], mean=curr, sd=Err_yC6[l], log=TRUE))
        #select potential new parameter value
          z=q_sel(mu=curr,tune=L_tuning[[2]][l,2], support=L_support, type="draw")
        #generate likelihood based on potential new parameter
                    #subtract Final_day (d_l[[2]][l,2]) DOC from Initial_day DOC (d_[[2]][l,1]) to calculate modeled DOC change for that time interval
                    #row reference for proc_current is + 1 to accomodate day 0
          Like_z =  dnorm(z, mean=proc_current[d_l[[2]][l,2]+1,7]-proc_current[d_l[[2]][l,1]+1,7], sd=Err_LC6[l], log=TRUE)+
                    #l+1 latent state, given subsequent model prediction, Final_day (d_l[[2]][l+1,2]) minus Initial_day (d_l[[2]][l+1,1]))
                    #sd at l+1 Err_Lc6
                    dnorm(L_curr[[2]][l+1], mean=proc_current[d_l[[2]][l+1,2]+1,7]-proc_current[d_l[[2]][l+1,1]+1,7], sd=Err_LC6[l+1], log=TRUE)+
                    #all measurements at l (all_data[[2]][l,3:5])
                    sum(dnorm(all_data[[2]][l,3:ncol(all_data[[2]])], mean=z, sd=Err_yC6[l], log=TRUE))
        #select parameter using choose function
          L_curr[[2]][l]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, 
                                param="Latent", tune=L_tuning[[2]][l,2], support=L_support)
        #update number of accepted proposals for a given parameter
          if(L_curr[[2]][l]==z) acceptL[[2]][l]=1
        
      }  else {
          curr=L_curr[[2]][l]
        #generate likelihood based on current parameter values
          Like_curr = dnorm(curr, mean=proc_current[d_l[[2]][l,2]+1,7]-proc_current[d_l[[2]][l,1]+1,7], sd=Err_LC6[l], log=TRUE)+
                      sum(dnorm(all_data[[2]][l,3:ncol(all_data[[2]])], mean=curr, sd=Err_yC6[l], log=TRUE))
        #select potential new parameter value
          z=q_sel(mu=curr,tune=L_tuning[[2]][l,2], support=L_support, type="draw")
        #generate likelihood based on potential new parameter
          Like_z =  dnorm(z, mean=proc_current[d_l[[2]][l,2]+1,7]-proc_current[d_l[[2]][l,1]+1,7], sd=Err_LC6[l], log=TRUE)+
                    sum(dnorm(all_data[[2]][l,3:ncol(all_data[[2]])], mean=z, sd=Err_yC6[l], log=TRUE))
        #select parameter using choose function
          L_curr[[2]][l]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, 
                                param="Latent", tune=L_tuning[[2]][l,2], support=L_support)
        #update number of accepted proposals for a given parameter
          if(L_curr[[2]][l]==z) acceptL[[2]][l]=1
      } #end else statements
    } #end for loop
  
  #function to update latent states for CO2 measurements
  for(l in 1:L3){
    if(l==1){
      #current value of latent CO2 at day index 'l'
        curr=L_curr[[3]][l]
      #generate likelihood based on current parameter values
                    #subtract Final_day (d_l[[3]][l,2]) CO2 from Initial_day CO2 (d_[[3]][l,1]) to calculate modeled CO2 change for that time interval
                    #row reference for proc_current is + 1 to accomodate day 0
        Like_curr = dnorm(curr, mean=proc_current[d_l[[3]][l,2]+1,8]-proc_current[d_l[[3]][l,1]+1,8], sd=Err_LC7[l], log=TRUE)+
                    #l+1 latent state, given subsequent model prediction, Final_day (d_l[[3]][l+1,2]) minus Initial_day (d_l[[3]][l+1,1]))
                    #sd at l+1 Err_Lc7
                    dnorm(L_curr[[3]][l+1], mean=proc_current[d_l[[3]][l+1,2]+1,8]-proc_current[d_l[[3]][l+1,1]+1,8], sd=Err_LC7[l+1], log=TRUE)+
                    #all measurements at l (all_data[[3]][l,3:5])
                    sum(dnorm(all_data[[3]][l,3:ncol(all_data[[3]])], mean=curr, sd=Err_yC7[l], log=TRUE))
      #select potential new parameter value
        z=q_sel(mu=curr,tune=L_tuning[[3]][l,2], support=L_support, type="draw")
      #generate likelihood based on potential new parameter
                  #subtract Final_day (d_l[[3]][l,2]) CO2 from Initial_day CO2 (d_[[3]][l,1]) to calculate modeled CO2 change for that time interval
                  #row reference for proc_current is + 1 to accomodate day 0
        Like_z =  dnorm(z, mean=proc_current[d_l[[3]][l,2]+1,8]-proc_current[d_l[[3]][l,1]+1,8], sd=Err_LC7[l], log=TRUE)+
                  #l+1 latent state, given subsequent model prediction, Final_day (d_l[[3]][l+1,2]) minus Initial_day (d_l[[3]][l+1,1]))
                  #sd at l+1 Err_Lc7
                  dnorm(L_curr[[3]][l+1], mean=proc_current[d_l[[3]][l+1,2]+1,8]-proc_current[d_l[[3]][l+1,1]+1,8], sd=Err_LC7[l+1], log=TRUE)+
                  #all measurements at l (all_data[[3]][l,3:5])
                  sum(dnorm(all_data[[3]][l,3:ncol(all_data[[3]])], mean=z, sd=Err_yC7[l], log=TRUE))

      #select parameter using choose function
        L_curr[[3]][l]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, 
                              param="Latent", tune=L_tuning[[3]][l,2], support=L_support)
      #update number of accepted proposals for a given parameter
        if(L_curr[[3]][l]==z) acceptL[[3]][l]=1
      
    }  else { if(l>1 & l<L3){
      #current value of latent CO2 at day index 'l'
        curr=L_curr[[3]][l]
      #generate likelihood based on current parameter values
                    #subtract Final_day (d_l[[3]][l,2]) CO2 from Initial_day CO2 (d_[[3]][l,1]) to calculate modeled CO2 change for that time interval
                    #row reference for proc_current is + 1 to accomodate day 0
        Like_curr = dnorm(curr, mean=proc_current[d_l[[3]][l,2]+1,8]-proc_current[d_l[[3]][l,1]+1,8], sd=Err_LC7[l], log=TRUE)+
                    #l+1 latent state, given subsequent model prediction, Final_day (d_l[[3]][l+1,2]) minus Initial_day (d_l[[3]][l+1,1]))
                    #sd at l+1 Err_Lc7
                    dnorm(L_curr[[3]][l+1], mean=proc_current[d_l[[3]][l+1,2]+1,8]-proc_current[d_l[[3]][l+1,1]+1,8], sd=Err_LC7[l+1], log=TRUE)+
                    #all measurements at l (all_data[[3]][l,3:5])
                    sum(dnorm(all_data[[3]][l,3:ncol(all_data[[3]])], mean=curr, sd=Err_yC7[l], log=TRUE))
      #select potential new parameter value
        z=q_sel(mu=curr,tune=L_tuning[[3]][l,2], support=L_support, type="draw")
      #generate likelihood based on potential new parameter
                  #subtract Final_day (d_l[[3]][l,2]) CO2 from Initial_day CO2 (d_[[3]][l,1]) to calculate modeled CO2 change for that time interval
                  #row reference for proc_current is + 1 to accomodate day 0
        Like_z =  dnorm(z, mean=proc_current[d_l[[3]][l,2]+1,8]-proc_current[d_l[[3]][l,1]+1,8], sd=Err_LC7[l], log=TRUE)+
                  #l+1 latent state, given subsequent model prediction, Final_day (d_l[[3]][l+1,2]) minus Initial_day (d_l[[3]][l+1,1]))
                  #sd at l+1 Err_Lc7
                  dnorm(L_curr[[3]][l+1], mean=proc_current[d_l[[3]][l+1,2]+1,8]-proc_current[d_l[[3]][l+1,1]+1,8], sd=Err_LC7[l+1], log=TRUE)+
                  #all measurements at l (all_data[[3]][l,3:5])
                  sum(dnorm(all_data[[3]][l,3:ncol(all_data[[3]])], mean=z, sd=Err_yC7[l], log=TRUE))

      #select parameter using choose function
        L_curr[[3]][l]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, 
                              param="Latent", tune=L_tuning[[3]][l,2], support=L_support)
      #update number of accepted proposals for a given parameter
        if(L_curr[[3]][l]==z) acceptL[[3]][l]=1
      
    } else {
        curr=L_curr[[3]][l]
      #generate likelihood based on current parameter values
        Like_curr = dnorm(curr, mean=proc_current[d_l[[3]][l,2]+1,8]-proc_current[d_l[[3]][l,1]+1,8], sd=Err_LC7[l], log=TRUE)+
                    sum(dnorm(all_data[[3]][l,3:ncol(all_data[[3]])], mean=curr, sd=Err_yC7[l], log=TRUE))
      #select potential new parameter value
        z=q_sel(mu=curr,tune=L_tuning[[3]][l,2], support=L_support, type="draw")
      #generate likelihood based on potential new parameter
        Like_z =  dnorm(z, mean=proc_current[d_l[[3]][l,2]+1,8]-proc_current[d_l[[3]][l,1]+1,8], sd=Err_LC7[l], log=TRUE)+
                  sum(dnorm(all_data[[3]][l,3:ncol(all_data[[3]])], mean=z, sd=Err_yC7[l], log=TRUE))
      #select parameter using choose function
        L_curr[[3]][l]=choose(x=curr, z=z, Like_z=Like_z, Like_x=Like_curr, 
                              param="Latent", tune=L_tuning[[3]][l,2], support=L_support)
      #update number of accepted proposals for a given parameter
        if(L_curr[[3]][l]==z) acceptL[[3]][l]=1
    } #end last else 
    } #end all else
  } #end for loop
  
      #######to do: calculation of mass balance & add that as a vector returned to the analysis
  return(list(L_curr, proc_current, acceptL))
} #end likelihood function


#likelihood function for model parameters
Like_proc_param=function(proc_Lat, Latent, errs_L_C1_5, errs_L_C6, errs_L_C7, num_litter){

  LP_paramSum=0
  #generate across litter types
    for(i in 1:num_litter){
      LP_paramMass<-dnorm(Latent[[i]][[1]], mean=proc_Lat[[i]][[1]], sd=errs_L_C1_5[[1]], log=TRUE)
      LP_paramDOC<-dnorm(Latent[[i]][[2]], mean=proc_Lat[[i]][[2]], sd=errs_L_C6[[1]], log=TRUE)
      LP_paramCO2<-dnorm(Latent[[i]][[3]], mean=proc_Lat[[i]][[3]], sd=errs_L_C7[[1]], log=TRUE)
      
      LP_paramSum=LP_paramSum+sum(LP_paramMass, LP_paramDOC, LP_paramCO2)
      #print(LP_paramSum)
    }
  #total
  return(LP_paramSum)
}

#Function to generate proc_current model results based on latent states and new draw for a given parameter value
#Parameter_new_mod <- function(d_l, LIDEL_inpt, Init_vals, tau, Lcalph, param_ups, L_curr){
Parameter_new_mod <- function(d_l, LIDEL_inpt, Init_vals, tau, Lcalph, param_ups){
  
  #Setup for analysis
    #Currently, CO2 has more measured days than DOC, therefore sets length of L
      #L=nrow(d_l[[3]])

    #set parameters for LIDEL model
      Pars <- c(tau=tau, EM1=LIDEL_inpt$EM1, Em1=LIDEL_inpt$Em1, EM2=LIDEL_inpt$EM2, 
                Em2=LIDEL_inpt$Em2, Nmid=LIDEL_inpt$Nmid, NC=LIDEL_inpt$NC,LcM=LIDEL_inpt$LcM, 
                beta1=LIDEL_inpt$beta1, beta2=LIDEL_inpt$beta2, lambda3=LIDEL_inpt$lambda3, 
                #new draw for parameter in question will be read with full parameter_ups list
                k1=as.numeric(param_ups[1]),
                k2=as.numeric(param_ups[2]), 
                k4=as.numeric(param_ups[3]), 
                beta3=as.numeric(param_ups[4]), 
                lambda2=as.numeric(param_ups[5]),
                Lcalph=Lcalph)
  
    #simulate full initial model results
    #initial condition at day 0
      call_initial_vals=LIDEL_initials(mass=Init_vals[1:2], fdoc=Init_vals[3:4], fs=Init_vals[5:6], flig=Init_vals[7:8])
      yini <- c(C1=call_initial_vals[2], C2=call_initial_vals[3], C3=call_initial_vals[4], C4=0, 
                C5=0, C6=call_initial_vals[5], C7=0)
    #run model from initial values to final day
    Time <- LIDEL_inpt$days 
    #save model results to proc_new
    proc_new=lsoda(func=DOCmod, y=yini, parms=Pars, times=Time)
    
#   for(l in 1:L){
#     #currently, only 1 CO2 measurement before DOC measurement. Therefore run out from initial (day 0) to day of 1st CO2 measurement.
#     #then initialize based on latent CO2 estimate, and run until next measured day.
#     if(l==1){
#       #simulate full initial model results
#         #initial condition at day 0
#             call_initial_vals=LIDEL_initials(mass=Init_vals[1:2], fdoc=Init_vals[3:4], fs=Init_vals[5:6], flig=Init_vals[7:8])
#             yini <- c(C1=call_initial_vals[2], C2=call_initial_vals[3], C3=call_initial_vals[4], C4=0, 
#                       C5=0, C6=call_initial_vals[5], C7=0)
#         #run model from initial values to final day
#             Time <- LIDEL_inpt$days 
#         #save model results to proc_new
#             proc_new=lsoda(func=DOCmod, y=yini, parms=Pars, times=Time)
#       
#       #simulate from 1st latent CO2 to day of next latent CO2/DOC measurement
#           #pull non-latent initial values from proc_new for d_l[[3]][l,2] (Final_day for l)
#           #each row called by +1 to associated d_l day, to accomodate shift in rows caused by day 0 in model simulation results
#             yini <- c(C1=proc_new[d_l[[3]][l,2]+1,2], C2=proc_new[d_l[[3]][l,2]+1,3], 
#                       C3=proc_new[d_l[[3]][l,2]+1,4], C4=proc_new[d_l[[3]][l,2]+1,5], 
#                       C5=proc_new[d_l[[3]][l,2]+1,6], C6=proc_new[d_l[[3]][l,2]+1,7], 
#                     #add L_curr CO2 gain at Final_day (L_curr[[3]][l]), to proc_new CO2 at Initial_day of d_l[[3]][l,1] (proc_new row ref: d_l[[3]][l,1]+1)
#                       C7=proc_new[d_l[[3]][l,1]+1,8]+L_curr[[3]][l])
#         #time points from Final_day of d_l[[3]] row l to Final_day of d_l[[3]] row L
#             Time <- LIDEL_inpt$days[(d_l[[3]][l,2]+1):(d_l[[3]][L,2]+1)] 
#         #update from one time step beyond d_l[[3]][l,2], to final day at l=L3
#         #current value at d_l[[3]][l,2] should remain the same
#             mod_out=lsoda(func=DOCmod, y=yini, parms=Pars, times=Time)
#             proc_new[Time[2:length(Time)]+1,]=mod_out[2:nrow(mod_out),]
#       
#     } else { if(l>1 & l<L){
#       
#       #simulate from day at d_l[[3]][l,2] to Final_day of next latent CO2/DOC measurement
#           #pull non-latent initial values from proc_new for d_l[[3]][l,2] (Final_day for l)
#           #each row called by +1 to associated d_l day, to accomodate shift in rows caused by day 0 in model simulation results
#             yini <- c(C1=proc_new[d_l[[3]][l,2]+1,2], C2=proc_new[d_l[[3]][l,2]+1,3], 
#                       C3=proc_new[d_l[[3]][l,2]+1,4], C4=proc_new[d_l[[3]][l,2]+1,5], 
#                       C5=proc_new[d_l[[3]][l,2]+1,6], 
#                     #add L_curr DOC gain at Final_day (L_curr[[2]][l-1]), to proc_new DOC at Initial_day of d_l[[2]][l-1,1] (proc_new row ref: d_l[[2]][l-1,1]+1)
#                     #l reference for DOC d_l and L_curr is (l-1) because DOC has 1 fewer time points measured for this analysis than DOC
#                       C6=proc_new[d_l[[2]][l-1,1]+1,7]+L_curr[[2]][l-1], 
#                     #add L_curr CO2 gain at Final_day (L_curr[[3]][l]), to proc_new CO2 at Initial_day of d_l[[3]][l,1] (proc_new row ref: d_l[[3]][l,1]+1)
#                       C7=proc_new[d_l[[3]][l,1]+1,8]+L_curr[[3]][l])
#         #time points from Final_day of d_l[[3]] row l to Final_day of d_l[[3]] row L
#              Time <- LIDEL_inpt$days[(d_l[[3]][l,2]+1):(d_l[[3]][L,2]+1)] 
#         #update from one time step beyond d_l[[3]][l,2], to final day at l=L3
#         #current value at d_l[[3]][l,2] should remain the same
#             mod_out=lsoda(func=DOCmod, y=yini, parms=Pars, times=Time)
#             proc_new[Time[2:length(Time)]+1,]=mod_out[2:nrow(mod_out),]
#       
#     } } #end else statements
# 
#   } #end for loop
  
  #create list to store subset of model results that correspond to d_l measurement intervals
    mu_diff=list()
      #mass remaining equates to sum of model pools 1:5
        mu_diff[[1]]=apply(proc_new[(d_l[[1]]+1),2:6],1, sum)
      #calculate difference between measured times for DOC
        mu_diff[[2]]=subset(proc_new[,7], proc_new[,1] %in% as.numeric(d_l[[2]][,2]))-
          subset(proc_new[,7], proc_new[,1] %in% as.numeric(d_l[[2]][,1]))
      #calculate difference between measured times for CO2
        mu_diff[[3]]=subset(proc_new[,8], proc_new[,1] %in% as.numeric(d_l[[3]][,2]))-
          subset(proc_new[,8], proc_new[,1] %in% as.numeric(d_l[[3]][,1]))
  
  
  return(list(proc_new, mu_diff))
} #end likelihood function

