##########ALL DATA UPLOAD########
#Upload data for full LIDEL model
#Nell Campbell
#Created 12/2/14
#revised 2/26/15
#################################

###main factors for data structure
#all litter types
  litter_nam=c("Alfalfa", "Ash", "Bluestem", "Oak", "Pine")
#data names
  data_name=c("y_C1_3", "y_C6","y_C7")

###load data
#load initial Nitrogen (tau)
  Init_N=read.csv("data_Init_N.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
#load initial LCI (Lcalph)
  Init_LCI=read.csv("data_Init_LCI.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
#Load mean and sd of initial mass and fraction soluble, doc, and lignin from Cinit_bayes analysis
  Init_Cs=read.csv("data_Bayes_initial_C1_2_3_6.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
      #mean and sd of mass
      Init_mass_distrib=as.matrix(Init_Cs[,2:3])
      #alpha and beta values for initial fdoc, fs, and flig
      Init_fdoc_distrib= cbind((Init_Cs[,4]^2-Init_Cs[,4]^3-(Init_Cs[,4]*Init_Cs[,5]^2))/(Init_Cs[,5]^2),
                                (Init_Cs[,4]-2*Init_Cs[,4]^2+Init_Cs[,4]^3-Init_Cs[,5]^2+(Init_Cs[,4]*Init_Cs[,5]^2))/(Init_Cs[,5]^2))
      Init_fs_distrib=cbind((Init_Cs[,6]^2-Init_Cs[,6]^3-(Init_Cs[,6]*Init_Cs[,7]^2))/(Init_Cs[,7]^2),
                            (Init_Cs[,6]-2*Init_Cs[,6]^2+Init_Cs[,6]^3-Init_Cs[,7]^2+(Init_Cs[,6]*Init_Cs[,7]^2))/(Init_Cs[,7]^2))
      Init_flig_distrib=cbind((Init_Cs[,8]^2-Init_Cs[,8]^3-(Init_Cs[,8]*Init_Cs[,9]^2))/(Init_Cs[,9]^2),
                              (Init_Cs[,8]-2*Init_Cs[,8]^2+Init_Cs[,8]^3-Init_Cs[,9]^2+(Init_Cs[,8]*Init_Cs[,9]^2))/(Init_Cs[,9]^2))

      #call set of initial values
      Initial_c=list()
      for(i in 1:5){
        Initial_c[[i]]=LIDEL_initials(mass=Init_mass_distrib[i,], fs=Init_fs_distrib[i,], 
                                    fdoc=Init_fdoc_distrib[i,], flig=Init_flig_distrib[i,])
        names(Initial_c[[i]])=c("Mass", "C1", "C2", "C3", "C6")
      }
      #list of values for LIDEL_initials function by litter type, for use in mapply
      #mass mean & sd, fdoc alpha & beta, fs alpha & beta, flig alpha & beta
        Init_fctn_vals=list(c(Init_mass_distrib[1,], Init_fdoc_distrib[1,], Init_fs_distrib[1,], Init_flig_distrib[1,]),
                            c(Init_mass_distrib[2,], Init_fdoc_distrib[2,], Init_fs_distrib[2,], Init_flig_distrib[2,]),
                            c(Init_mass_distrib[3,], Init_fdoc_distrib[3,], Init_fs_distrib[3,], Init_flig_distrib[3,]),
                            c(Init_mass_distrib[4,], Init_fdoc_distrib[4,], Init_fs_distrib[4,], Init_flig_distrib[4,]),
                            c(Init_mass_distrib[5,], Init_fdoc_distrib[5,], Init_fs_distrib[5,], Init_flig_distrib[5,]))
        for(i in 1:5){
          names(Init_fctn_vals[[i]])=c("Init_mass", "Init_mass_sd", "fdoc_alpha", "fdoc_beta", "fs_alpha", "fs_beta", "flig_alpha", "flig_beta")
        }
        names(Init_fctn_vals)=litter_nam
#loss of carbon
  data_Ctot=read.csv("data_Ctot.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
#gain of DOc
  data_DOC=read.csv("data_DOC.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)
#gain of CO2
  data_CO2=read.csv("data_CO2.csv", header=TRUE, sep=",", stringsAsFactors=FALSE)


#################ALFALFA#################
#create list of measured data for alfalfa
  Alfalfa_dat=list(as.matrix(data_Ctot[data_Ctot[,1]=="Alfalfa",2:ncol(data_Ctot)]), 
               as.matrix(data_DOC[data_DOC[,1]=="Alfalfa",2:ncol(data_DOC)]),
               as.matrix(data_CO2[data_CO2[,1]=="Alfalfa",2:ncol(data_CO2)]))
  names(Alfalfa_dat)=data_name

#create matrix of cumulative DOC
  Alf_DOC=matrix(0, nrow=nrow(Alfalfa_dat$y_C6)+1, ncol=ncol(Alfalfa_dat$y_C6)-1)
  #fill initial values
    #initial value of day 10, proxy for initial DOC
      Alf_DOC[1,1]=10
    #fill initial DOC value from initial values for C pools
      Alf_DOC[1,2:(ncol(Alfalfa_dat$y_C6)-1)]=Initial_c[[1]][5]
  #fill all remaining time points for measurements
    Alf_DOC[2:nrow(Alf_DOC),1]=Alfalfa_dat$y_C6[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Alfalfa_dat$y_C6)){
      for(i in 1:nrow(Alfalfa_dat$y_C6)){
        if(Alfalfa_dat$y_C6[i,j]!=-99.99){
        Alf_DOC[i+1,j-1]=Alf_DOC[i,j-1]+Alfalfa_dat$y_C6[i,j]
      } 
      if(Alfalfa_dat$y_C6[i,j]==-99.99 & i<8){
        Alf_DOC[i+1,j-1]=Alf_DOC[i,j-1]
      }
      if(Alfalfa_dat$y_C6[i,j]==-99.99 & i>=8){
        Alf_DOC[i+1,j-1]=-99.99
      }
      }
    }

#create matrix of cumulative CO2
  Alf_CO2=matrix(0, nrow=nrow(Alfalfa_dat$y_C7), ncol=ncol(Alfalfa_dat$y_C7)-1)
  #fill initial values
    Alf_CO2[1,2:(ncol(Alfalfa_dat$y_C7)-1)]=Alfalfa_dat$y_C7[1,3:ncol(Alfalfa_dat$y_C7)]
  #fill time point
    Alf_CO2[,1]=Alfalfa_dat$y_C7[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Alfalfa_dat$y_C7)){
      for(i in 2:nrow(Alfalfa_dat$y_C7)){
        if(Alfalfa_dat$y_C7[i,j]!=-99.99){
          Alf_CO2[i,j-1]=Alf_CO2[i-1,j-1]+Alfalfa_dat$y_C7[i,j]
        } 
        if(Alfalfa_dat$y_C7[i,j]==-99.99 & i<8){
          Alf_CO2[i,j-1]=Alf_CO2[i-1,j-1]
        }
        if(Alfalfa_dat$y_C7[i,j]==-99.99 & i>=8){
          Alf_CO2[i,j-1]=-99.99
        }
      }
    }

#################ASH#################
#create list of measured data for ash
Ash_dat=list(as.matrix(data_Ctot[data_Ctot[,1]=="Ash",2:ncol(data_Ctot)]), 
             as.matrix(data_DOC[data_DOC[,1]=="Ash",2:ncol(data_DOC)]),
             as.matrix(data_CO2[data_CO2[,1]=="Ash",2:ncol(data_CO2)]))
names(Ash_dat)=data_name

#create matrix of cumulative DOC
  Ash_DOC=matrix(0, nrow=nrow(Ash_dat$y_C6)+1, ncol=ncol(Ash_dat$y_C6)-1)
  #fill initial values
    #initial value of day 10, proxy for initial DOC
      Ash_DOC[1,1]=10
    #fill initial DOC value from initial values for C pools
      Ash_DOC[1,2:(ncol(Ash_dat$y_C6)-1)]=Initial_c[[2]][5]
  #fill time points
    Ash_DOC[2:nrow(Ash_DOC),1]=Ash_dat$y_C6[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Ash_dat$y_C6)){
      for(i in 1:nrow(Ash_dat$y_C6)){
        if(Ash_dat$y_C6[i,j]!=-99.99){
          Ash_DOC[i+1,j-1]=Ash_DOC[i,j-1]+Ash_dat$y_C6[i,j]
        } 
        if(Ash_dat$y_C6[i,j]==-99.99 & i<8){
          Ash_DOC[i+1,j-1]=Ash_DOC[i,j-1]
        }
        if(Ash_dat$y_C6[i,j]==-99.99 & i>=8){
          Ash_DOC[i+1,j-1]=-99.99
        }
      }
    }

#create matrix of cumulative CO2
  Ash_CO2=matrix(0, nrow=nrow(Ash_dat$y_C7), ncol=ncol(Ash_dat$y_C7)-1)
  #fill initial values
    Ash_CO2[1,2:(ncol(Ash_dat$y_C7)-1)]=Ash_dat$y_C7[1,3:ncol(Ash_dat$y_C7)]
  #fill time point
    Ash_CO2[,1]=Ash_dat$y_C7[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Ash_dat$y_C7)){
      for(i in 2:nrow(Ash_dat$y_C7)){
        if(Ash_dat$y_C7[i,j]!=-99.99){
          Ash_CO2[i,j-1]=Ash_CO2[i-1,j-1]+Ash_dat$y_C7[i,j]
        } 
        if(Ash_dat$y_C7[i,j]==-99.99 & i<8){
          Ash_CO2[i,j-1]=Ash_CO2[i-1,j-1]
        }
        if(Ash_dat$y_C7[i,j]==-99.99 & i>=8){
          Ash_CO2[i,j-1]=-99.99
        }
      }
    }

#################BLUESTEM#################
#create list of measured data for bluestem
Bluestem_dat=list(as.matrix(data_Ctot[data_Ctot[,1]=="Bluestem",2:ncol(data_Ctot)]), 
         as.matrix(data_DOC[data_DOC[,1]=="Bluestem",2:ncol(data_DOC)]),
         as.matrix(data_CO2[data_CO2[,1]=="Bluestem",2:ncol(data_CO2)]))
names(Bluestem_dat)=data_name

#create matrix of cumulative DOC
  Blue_DOC=matrix(0, nrow=nrow(Bluestem_dat$y_C6)+1, ncol=ncol(Bluestem_dat$y_C6)-1)
  #fill initial values
    #initial value of day 10, proxy for initial DOC
      Blue_DOC[1,1]=10
    #fill initial DOC value from initial values for C pools
      Blue_DOC[1,2:(ncol(Bluestem_dat$y_C6)-1)]=Initial_c[[3]][5]
  #fill time points
    Blue_DOC[2:nrow(Blue_DOC),1]=Bluestem_dat$y_C6[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Bluestem_dat$y_C6)){
      for(i in 1:nrow(Bluestem_dat$y_C6)){
        if(Bluestem_dat$y_C6[i,j]!=-99.99){
          Blue_DOC[i+1,j-1]=Blue_DOC[i,j-1]+Bluestem_dat$y_C6[i,j]
        } 
        if(Bluestem_dat$y_C6[i,j]==-99.99 & i<8){
          Blue_DOC[i+1,j-1]=Blue_DOC[i,j-1]
        }
        if(Bluestem_dat$y_C6[i,j]==-99.99 & i>=8){
          Blue_DOC[i+1,j-1]=-99.99
        }
      }
    }

#create matrix of cumulative CO2
  Blue_CO2=matrix(0, nrow=nrow(Bluestem_dat$y_C7), ncol=ncol(Bluestem_dat$y_C7)-1)
  #fill initial values
    Blue_CO2[1,2:(ncol(Bluestem_dat$y_C7)-1)]=Bluestem_dat$y_C7[1,3:ncol(Bluestem_dat$y_C7)]
  #fill time point
    Blue_CO2[,1]=Bluestem_dat$y_C7[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Bluestem_dat$y_C7)){
      for(i in 2:nrow(Bluestem_dat$y_C7)){
        if(Bluestem_dat$y_C7[i,j]!=-99.99){
          Blue_CO2[i,j-1]=Blue_CO2[i-1,j-1]+Bluestem_dat$y_C7[i,j]
        } 
        if(Bluestem_dat$y_C7[i,j]==-99.99 & i<8){
          Blue_CO2[i,j-1]=Blue_CO2[i-1,j-1]
        }
        if(Bluestem_dat$y_C7[i,j]==-99.99 & i>=8){
          Blue_CO2[i,j-1]=-99.99
        }
      }
    }

#################OAK#################
#create list of measured data for oak
Oak_dat=list(as.matrix(data_Ctot[data_Ctot[,1]=="Oak",2:ncol(data_Ctot)]), 
         as.matrix(data_DOC[data_DOC[,1]=="Oak",2:ncol(data_DOC)]),
         as.matrix(data_CO2[data_CO2[,1]=="Oak",2:ncol(data_CO2)]))
names(Oak_dat)=data_name

#create matrix of cumulative DOC
  Oak_DOC=matrix(0, nrow=nrow(Oak_dat$y_C6)+1, ncol=ncol(Oak_dat$y_C6)-1)
  #fill initial values
    #initial value of day 10, proxy for initial DOC
      Oak_DOC[1,1]=10
    #fill initial DOC value from initial values for C pools
      Oak_DOC[1,2:(ncol(Oak_dat$y_C6)-1)]=Initial_c[[4]][5]
  #fill time points
    Oak_DOC[2:nrow(Oak_DOC),1]=Oak_dat$y_C6[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Oak_dat$y_C6)){
      for(i in 1:nrow(Oak_dat$y_C6)){
        if(Oak_dat$y_C6[i,j]!=-99.99){
          Oak_DOC[i+1,j-1]=Oak_DOC[i,j-1]+Oak_dat$y_C6[i,j]
        } 
        if(Oak_dat$y_C6[i,j]==-99.99 & i<8){
          Oak_DOC[i+1,j-1]=Oak_DOC[i,j-1]
        }
        if(Oak_dat$y_C6[i,j]==-99.99 & i>=8){
          Oak_DOC[i+1,j-1]=-99.99
        }
      }
    }

#create matrix of cumulative CO2
  Oak_CO2=matrix(0, nrow=nrow(Oak_dat$y_C7), ncol=ncol(Oak_dat$y_C7)-1)
  #fill initial values
    Oak_CO2[1,2:(ncol(Oak_dat$y_C7)-1)]=Oak_dat$y_C7[1,3:ncol(Oak_dat$y_C7)]
  #fill time point
    Oak_CO2[,1]=Oak_dat$y_C7[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Oak_dat$y_C7)){
      for(i in 2:nrow(Oak_dat$y_C7)){
        if(Oak_dat$y_C7[i,j]!=-99.99){
          Oak_CO2[i,j-1]=Oak_CO2[i-1,j-1]+Oak_dat$y_C7[i,j]
        } 
        if(Oak_dat$y_C7[i,j]==-99.99 & i<8){
          Oak_CO2[i,j-1]=Oak_CO2[i-1,j-1]
        }
        if(Oak_dat$y_C7[i,j]==-99.99 & i>=8){
          Oak_CO2[i,j-1]=-99.99
        }
      }
    }

#################PINE#################
#create list of measured data for pine
Pine_dat=list(as.matrix(data_Ctot[data_Ctot[,1]=="Pine",2:ncol(data_Ctot)]), 
         as.matrix(data_DOC[data_DOC[,1]=="Pine",2:ncol(data_DOC)]),
         as.matrix(data_CO2[data_CO2[,1]=="Pine",2:ncol(data_CO2)]))
names(Pine_dat)=data_name

#create matrix of cumulative DOC
  Pine_DOC=matrix(0, nrow=nrow(Pine_dat$y_C6)+1, ncol=ncol(Pine_dat$y_C6)-1)
  #fill initial values
    #initial value of day 10, proxy for initial DOC
      Pine_DOC[1,1]=10
    #fill initial DOC value from initial values for C pools
      Pine_DOC[1,2:(ncol(Pine_dat$y_C6)-1)]=Initial_c[[5]][5]
  #fill time points
    Pine_DOC[2:nrow(Pine_DOC),1]=Pine_dat$y_C6[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Pine_dat$y_C6)){
      for(i in 1:nrow(Pine_dat$y_C6)){
        if(Pine_dat$y_C6[i,j]!=-99.99){
          Pine_DOC[i+1,j-1]=Pine_DOC[i,j-1]+Pine_dat$y_C6[i,j]
        } 
        if(Pine_dat$y_C6[i,j]==-99.99 & i<8){
          Pine_DOC[i+1,j-1]=Pine_DOC[i,j-1]
        }
        if(Pine_dat$y_C6[i,j]==-99.99 & i>=8){
          Pine_DOC[i+1,j-1]=-99.99
        }
      }
    }

#create matrix of cumulative CO2
  Pine_CO2=matrix(0, nrow=nrow(Pine_dat$y_C7), ncol=ncol(Pine_dat$y_C7)-1)
  #fill initial values
    Pine_CO2[1,2:(ncol(Pine_dat$y_C7)-1)]=Pine_dat$y_C7[1,3:ncol(Pine_dat$y_C7)]
  #fill time point
    Pine_CO2[,1]=Pine_dat$y_C7[,2]
  #calculate cumulative DOC for each replicate
    for(j in 3:ncol(Pine_dat$y_C7)){
      for(i in 2:nrow(Pine_dat$y_C7)){
        if(Pine_dat$y_C7[i,j]!=-99.99){
          Pine_CO2[i,j-1]=Pine_CO2[i-1,j-1]+Pine_dat$y_C7[i,j]
        } 
        if(Pine_dat$y_C7[i,j]==-99.99 & i<8){
          Pine_CO2[i,j-1]=Pine_CO2[i-1,j-1]
        }
        if(Pine_dat$y_C7[i,j]==-99.99 & i>=8){
          Pine_CO2[i,j-1]=-99.99
        }
      }
    }

all_data=list(Alfalfa_dat,Ash_dat,Bluestem_dat,Oak_dat,Pine_dat)
names(all_data)=litter_nam

cumul_DOC=list(Alf_DOC, Ash_DOC, Blue_DOC, Oak_DOC, Pine_DOC)
names(cumul_DOC)=litter_nam

cumul_CO2=list(Alf_CO2, Ash_CO2, Blue_CO2, Oak_CO2, Pine_CO2)
names(cumul_CO2)=litter_nam

# #################PLOTS#################
# par(mfrow=c(3,2))
# colors=c("black", "orange", "lightblue", "lightgreen", "darkgrey", "pink")


# par(mfrow=c(3,2))
# #plot of litter C remaining through time
  # for(t in 1:length(litter_nam)){
    # Ctot_lim=Initial_c[[5]][1]
    # plot(c(0), c(Initial_c[[t]][1]),
              # col=colors, ylim=c(0,Ctot_lim), xlim=c(0, 365),
         # main=litter_nam[t], ylab="Remaining litter C (mg)", typ="p", xlab="Time", lwd=3, pch=4)
    # for(s in 3:ncol(all_data[[t]][[1]])){
      # points(all_data[[t]][[1]][,2], all_data[[t]][[1]][,s], col=colors[s-1], 
             # typ="p", pch=4, lwd=3)
    # }
  # }

# par(mfrow=c(3,2))
# #plot of DOC accumulation through time
# for(t in 1:length(litter_nam)){
  # DOC_lim=max(cumul_DOC[[t]][,2:ncol(cumul_DOC[[t]])])+10
  # plot(cumul_DOC[[t]][,1], cumul_DOC[[t]][,2], col=colors[1], ylim=c(0,DOC_lim),
       # main=litter_nam[t], ylab="DOC (mg)", typ="p", xlab="Time", lwd=3, pch=4)
  # for(s in 2:ncol(cumul_DOC[[1]])){
    # if(cumul_DOC[[t]][9,s]!=-99.99){
      # points(cumul_DOC[[t]][,1], cumul_DOC[[t]][,s], col=colors[s-1], 
             # typ="p", pch=4, lwd=3)
    # }
    # if(cumul_DOC[[t]][9,s]==-99.99){
    # points(cumul_DOC[[t]][1:7,1], cumul_DOC[[t]][1:7,s], col=colors[s-1], 
           # typ="p", pch=4, lwd=3)
    # }
  # }
# }

# par(mfrow=c(3,2))
# #plot of CO2 accumulation through time
# for(t in 1:length(litter_nam)){
  # CO2_lim=max(cumul_CO2[[t]][,2:ncol(cumul_CO2[[t]])])+10
  # plot(cumul_CO2[[t]][,1], cumul_CO2[[t]][,2], col=colors[1], ylim=c(0,CO2_lim),
       # main=litter_nam[t], ylab="CO2-C (mg)", typ="p", xlab="Time", lwd=3, pch=4)
  # for(s in 2:ncol(cumul_CO2[[1]])){
    # if(cumul_CO2[[t]][18,s]!=-99.99){
      # points(cumul_CO2[[t]][,1], cumul_CO2[[t]][,s], col=colors[s-1], 
             # typ="p", pch=4, lwd=3)
    # }
    # if(cumul_CO2[[t]][18,s]==-99.99){
      # if(t==1){
      # points(cumul_CO2[[t]][1:17,1], cumul_CO2[[t]][1:17,s], col=colors[s-1], 
             # typ="p", pch=4, lwd=3)
      # }
      # if(t!=1){
        # points(cumul_CO2[[t]][1:12,1], cumul_CO2[[t]][1:12,s], col=colors[s-1], 
               # typ="p", pch=4, lwd=3)
      # }
    # }
  # }
# }