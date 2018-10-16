############################################################
#Project:       Master Research project
#Description :  Multispecies Competition Lotka-Voterra model Euler simulation functions

# By:     Gilbert K Langat

setwd("E:/MSC/Code Repo/MastersCode")
############################################################
## FUNCTION DEFINITIONS

# lvm(t,pop,parms)
# Use:	Function to calculate derivative of multispecies Lotka-Volterra equations
# Input:
# 	t: time (not used here, because there is no explicit time dependence)
# 	pop: vector containing the current abundance of all species
# 	parms: 	dummy variable,(used to pass on parameter values), Here, int_mat,str_mat,carry_cap,int_growth are parameter values
#           with: int_mat     - interaction matrix,
#                 str_mat     - competition strength matrix
#                 carry_cap   -  carrying capacity
#                 int_growth  - intrinsic growth rate
# Output:
#	dN: derivative of the Modified Lotka-Volterra equations
#############################################################

LVM <- function(t,pop,int_mat,str_mat,carry_cap,int_growth){ 
  
  dN=int_growth*pop*((carry_cap-(int_mat*str_mat)%*%pop)/carry_cap)
  
  return(dN)
  
}


## own sorting function

top <- function(x, n=S){
  result <- numeric()
  for(i in 1:n){
    j <- which.max(x)
    result[i] <- x[j]
    x[j] <- -Inf
  }
  result
}

###########################################################
## This function implement the LVM functions in Euler numerical simulation 

# Non-switching community

###########################################################

NonSwitch_Euler = function( t_int, y_int, stepsize, t_end,int_mat,str_mat,carry_cap,int_growth){     
  m = length(y_int)+1
  
  # Number of time steps 
  nsteps = ceiling((t_end-t_int)/stepsize) 
  
  # dimensions of output matrix
  Y_out = matrix(ncol =m, nrow = nsteps+1)
  Y_out[1,] = c(t_int, y_int) # initial output matrix
  
  # loop for implementing the Euler function over nsteps
  for (i in 1:nsteps) {
    Y_out[i+1,]= Y_out[i,]+ stepsize*c(1, LVM(Y_out[i,1],Y_out[i,2:m],int_mat,str_mat,carry_cap,int_growth))
    
  }
  
  return(data.frame(Y_out))
  
}

#######################################################################################
# Elimination switching euler simulation

#######################################################################################

Elimination_switch_Euler = function(t_int, y_int, stepsize, t_end, int_mat, str_mat,NestA,carry_cap,int_growth){
  m = length(y_int)+1
  
  # Number of steps 
  
  nsteps = ceiling((t_end-t_int)/stepsize)
  
  #Output matrix dimensions
  
  Y_out = matrix(ncol =m, nrow = nsteps+1)
  Y_out[1,] = c(t_int, y_int)
  
  #Updating the population after each timestep
  current_pop=y_int
  
  # loop for implementing the function over nsteps
  for (i in 1:nsteps) {
    
    ###############################################
    # Elimination switching rule implemention
    #############################################
    
    #Community matrix
    
    Elim_Intstr=(int_mat*str_mat)*current_pop 
    
    #Ensuring the selected species interacts with more than one species 
    repeat{
      IntRowsums=rowSums(int_mat) 
      introwsums_greater1= which(IntRowsums>1,arr.ind = T) # select rows with more than one species interactions
      j_elim=sample_g(introwsums_greater1) # From the species selected, choose one 
      if (sum(int_mat[j_elim,]) < S){
        break
      }
    }
    
    #the vector selected based on overal commnity matrix
    Vec_elim= Elim_Intstr[j_elim,] 
    
    #Select a maximum value but not on  the main diagonal
    Max_Vec_elim =((int_mat*OffDiag_Strmat)*current_pop)[j_elim,]
    j= which(Max_Vec_elim == max(Max_Vec_elim))
    
    #Check the indices of all zeros the sample one element
    NI=which(int_mat[j_elim,]!=1,arr.ind = T)
    f=sample_g(NI)
    
    # swapping between the zero selected the original value
    
    int_mat[j_elim,c(j,f)]= int_mat[j_elim,c(f,j)]
    
    # Calculate the nestedness at the beginning and end of time steps
    
    if ( i== 0 || i == nsteps){
      Nesta=nested(int_mat,method="NODF",rescale = FALSE, normalised = TRUE)
    } else {
      Nesta = 0
    }
    NestA <- c(NestA,Nesta)
    
    # Euler iteration
    
    
    Y_out[i+1,]= Y_out[i,]+ stepsize*c(1, LVM(Y_out[i,1],Y_out[i,2:m],int_mat,str_mat,carry_cap,int_growth))
    
    current_pop= Y_out[i+1,2:m]
    
  }
  
  return(data.frame(Y_out,NestA))
  
  
} # End of Elimination switching euler function

###############################################################################
# Optimization switching euler simulation

##############################################################################

Optimization_switch_Euler = function(t_int, y_int, stepsize, t_end, int_mat, str_mat,Nest_opt,carry_cap,int_growth){
  m = length(y_int)+1
  # Number of time steps 
  
  nsteps = ceiling((t_end-t_int)/stepsize)
  
  # Output matrix dimensions
  
  Y_out = matrix(ncol =m, nrow = nsteps+1)
  Y_out[1,] = c(t_int, y_int)
  
  #Updating species population at end of each time step
  new_pop=y_int
  
  # Iteration loop over  nsteps
  
  for (i in 1:nsteps) {
    
    ############################################
    # Optimization switching rule implementation at each time step
    ############################################
    
    #Updating Community matrix
    Opt_Intstr=(int_mat*str_mat)*new_pop               
    # Select species interaction with more than one partner
    
    repeat{
      IntRowsums=rowSums(int_mat)
      introwsums_greater=which(IntRowsums>2,arr.ind = T) # choose species with more than two interactions
      if (length(introwsums_greater)>1){
        j_opt= sample_g(introwsums_greater)
      } else{
        j_opt=introwsums_greater
      }
      if (sum(int_mat[j_opt,])< S){
        break
      }
      
    }
    
    Vec_opt= Opt_Intstr[j_opt,] # the choosen vector in the overall community matrix
    
    # Choose the  2 non-interacting partners (zeros in int_mat) & sample_g one
    
    
    j_k=which(Vec_opt!=Vec_opt[j_opt] & Vec_opt!=0,arr.ind = T) ; k_opt=sample_g(j_k)
    j_m=which(Vec_opt!=Vec_opt[j_opt]& Vec_opt!=Vec_opt[k_opt],arr.ind = T)
    
    new_removed_ind=c()
    removed_ind=c(j_opt,k_opt)
    
    dN_switch_Original=int_growth*new_pop*((carry_cap-(int_mat*str_mat)%*%new_pop)/carry_cap)
    
    # Checking if switching increase the species growth
    
    enter=TRUE
    while(enter==TRUE){
      m_opt=sample_g(j_m)
      
      #Updating the interaction and strength matrices
      switch_int_mat=int_mat
      switch_str_mat=str_mat
      #Swapping the selected partners 
      
      int_mat[j_opt,c(k_opt,m_opt)]=int_mat[j_opt,c(m_opt,k_opt)]
      str_mat[j_opt,c(k_opt,m_opt)]=str_mat[j_opt,c(m_opt,k_opt)]
      
      dN_switch=int_growth*new_pop*((carry_cap-(int_mat*str_mat)%*%new_pop)/carry_cap)
      
      if ( (dN_switch[j_opt]-dN_switch_Original[j_opt]) > 0 
           && length(dN_switch[j_opt]-dN_switch_Original[j_opt])>0 
           &&!is.null(dN_switch[j_opt]-dN_switch_Original[j_opt])
           && !is.na(dN_switch[j_opt]-dN_switch_Original[j_opt])){
        
        enter=FALSE
        
      }
      
      else{
        
        int_mat=switch_int_mat
        
        str_mat=switch_str_mat
        
        new_removed_ind=c(new_removed_ind,m_opt)
        j_m=j_m[!(j_m %in% new_removed_ind)]
        
        if (length(j_m)==0){
          
          enter=FALSE
          
        }
        
      } #End of else
      
    }# End of while loop
    
    ## Calcualting the nestedness of matrix  at the beginning and end of the time steps 
    if ( i== 0 || i == nsteps){
      Nestopt=nested(int_mat,method="NODF",rescale = FALSE, normalised = TRUE)
    } else {
      Nestopt = 0
    }
    Nest_opt <- c(Nest_opt,Nestopt)   
    
    #Euler iteration over nsteps 
    
    Y_out[i+1,]= Y_out[i,]+ stepsize*c(1, LVM(Y_out[i,1],Y_out[i,2:m],int_mat,str_mat,carry_cap,int_growth))
    new_pop= Y_out[i+1,2:m]
    
  }
  return(data.frame(Y_out,Nest_opt))
  
} # End of optimization switching euler fuction

#################################################
#Stability computation function 


LVM_stability <- function(pop){
  
  dN1=int_growth*pop*((carry_cap-(IM*SM)%*%pop)/carry_cap)
  return(dN1)
}
#######################

#Sampling function
sample_g <- function(x) {
  if (length(x) <= 1) {
    return(x)
  } 
  else {
    return(sample(x,1))
  }
}

######################################################################################################
# SAVING THE FUNCTIONS 
save(LVM,NonSwitch_Euler,Elimination_switch_Euler,Optimization_switch_Euler,
     LVM_stability,sample_g,file="SingleCommunityFunctions.RData")

######################################################################################################

###END
