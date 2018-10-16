#########################################################################################
#CODE FOR RESEARCH MASTERS THESIS
#Project: Implementattion of Modified Lotka-Volterra model on Euler method
#By: Gilbert K Langat



##############################################
#Clearing memory and Saving code to folder
##############################################

rm(list=ls(all=TRUE)) 


################################################################################
## IMPORTING FUNCTIONs 
################################################################################

source("E:/MSC/Code Repo/MastersCode/MetacommunityFunctions.R")

setwd("E:/MSC/Code Repo/MastersCode/Meta/Results 5")

#set.seed(123)
##########################
# LOAD LIBRARIES

library(tidyverse)
library(BiodiversityR)
library(pracma)
library(vegan)
library(bipartite)
library(rnetcarto)
library(ggthemes)
library(reshape2)
library(magic)



################################################################################
#Appending containers
Total_Abundance_Nonswitch_datalist=list(); Total_Abundance_Elim_datalist=list()
Total_Abundance_Opt_datalist=list()

Non_switchingdatalist=list()
Elimination_switchingdatalist=list()
Optimization_switchingdatalist=list()
Local_Stability_datalist=list();  Abundance_datalist= list()
GlobalStability_datalist = list();Global_Productivity_datalist=list()

combinedabundancelist = list()

###############################################################
#initial value for parameters

SPECIES=c(20,50,100)
for (S in SPECIES ){
  p=S+1
  sites=3
  
  int_growth=runif(S,min =0.5,max =0.5)
  carry_cap= runif(S, min=1, max=1)
  
  pop=runif(S,min =0.1,max =0.9)
  Original_pop=replicate(sites,pop)
  
  d=0.2
  
  Emigration_Prop=matrix(runif(sites*sites,min =d,max =d), nrow = sites, ncol=sites)
  
  Mortality= matrix(runif(sites*sites,min =0.5,max =0.5), nrow = sites)
  for (emmirate in 1:sites) {
    Mortality[emmirate,emmirate]=-1
  }
  steps=10001
  ww=replicate(101,0)
  cc=c(0:100)
  sd=0.1
  
  ##D for str_mat (competition strength matrix)
  D = matrix(rnorm(S*S,0,sd),nrow=S);diag(D)=0
  OffDiag_Strmat = abs(D) ; diag(D)= 1
  str_mat = abs(D)
  
  Meta_Strength_Matrix= array(replicate(sites,str_mat), dim=c(S,S,sites))
  
  len =length(seq(0.1,0.9,0.05))
  con_abund_NonSw = c();
  con_abund_elimination = c();con_abund_optimization = c()
  Meta_max_eigen_NonSw=c();Global_maxeigen_Nonswitch=c();Lead_Eigen_Nonswitch=c()
  Meta_max_eigen_elim=c();Global_maxeigen_Elim=c();Lead_Eigen_Elim=c()
  Meta_max_eigen_opt=c();Global_maxeigen_Opt=c() ;Lead_Eigen_Opt=c()
  ###################
  
  ####
  #Appending containers
  ##
  connectance_values=c() ; connect_abundance=c()
  connect_abund_elim_switch=c();connect_abund_opt_switch=c()
  connect_nest_elim_switch= c();connect_nest_opt_switch= c()
  connect_stab_elim_switch= c();connect_stab_opt_switch= c()
  connect_nested= c() ;connect_stability= c() 
  
  connectA=c()
  
  #### Loop for geneting species dynamics over time for each connectance value
  # z - connectance
  #con_seq - connectance sequence, form o.1 to 0.9 with step of 0.05
  Total_abundance=list()
  con_seq=seq(0.1,0.9,0.1)
  for (z in 1:length(con_seq)){
    s = S*S; nc = round(con_seq[z]*S*S, digits = 0)
    int_mat = matrix(1,S,S); diag(int_mat)=0                                       
    while (s>nc){                                           
      a <- sample(S,1)                                        
      b <- sample(S,1)                                        
      if (int_mat[a,b]>0 && int_mat[b,a]>0 ){      
        int_mat[a,b] <- 0
        int_mat[b,a] <- 0
      }
      s <- s-1
    }
    Meta_Interaction_Matrix_Offdiag= array(replicate(sites,int_mat), dim=c(S,S,sites))
    
    diag(int_mat)=1
    
    Meta_Interaction_Matrix= array(replicate(sites,int_mat), dim=c(S,S,sites))
    
    #######################################################################################################
    Meta_Non_Switching=Spatial_euler_meth(0,pop, 0.01, 100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,site)
    
    ##Non-Switching total abundance per time step for each patch
    
    Meta_Nonswitch_TAbundnance <- colSums(aperm(Meta_Non_Switching[,2:(S+1),], c(2,1,3)))
    connect_abundance=data.frame(rbind(connect_abundance,colMeans(Meta_Nonswitch_TAbundnance)))
    
    Total_Abundance_Nonswitch =transform(connect_abundance,TAbundance_Nonswitch=rowSums(connect_abundance[,1:sites]))
    
    ########################################################################################################
    
    Meta_elimination_switch=Meta_euler_meth_elimination_switch(0, pop, 0.01,100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,site)
    
    ##Elimination total abundance per time step for each patch
    
    Meta_Elimination_TAbundnance <- colSums(aperm(Meta_elimination_switch[,2:(S+1),], c(2,1,3)))
    connect_abund_elim_switch=data.frame(rbind(connect_abund_elim_switch,colMeans(Meta_Elimination_TAbundnance)))
    
    Total_Abundance_Elim =transform(connect_abund_elim_switch,TAbundance_Elim=rowSums(connect_abund_elim_switch[,1:sites]))
    
    ######################################################################################################
    
    ## Optimization switching population dynamics plots over the time step for every connectance value.
    ## Plots saved in a folder.
    
    
    Meta_optimization_switch=Meta_euler_meth_optimization_switch(0, pop, 0.01,100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,site)
    
    ##Optimization total abundance per time step for each patch
    
    Meta_Optimization_TAbundnance <- colSums(aperm(Meta_optimization_switch[,2:(S+1),], c(2,1,3)))
    connect_abund_opt_switch=data.frame(rbind(connect_abund_opt_switch,colMeans(Meta_Optimization_TAbundnance)))
    
    Total_Abundance_Opt =transform(connect_abund_opt_switch,TAbundance_Opt=rowSums(connect_abund_opt_switch[,1:sites]))
    
    #########################################################################################
    
    ## STABILITY COMPUTATION
    
    ########################################################################################
    
    max_eigen_NonSw=c();max_eigen_elim=c();max_eigen_NonSpatialSw=c()
    max_eigen_opt=c();Max_eigen_NonSW=c();Max_eigen_NonSpatialSw=c()
    Max_eigen_elim=c();Max_eigen_opt=c()
    Jacobian_Nonswitch = list();Jacobian_Nonspatial=list()
    Jacobian_Elim = list();Jacobian_Opt=list()
    
    
    ## STABLITY PER SITE
    
    for (site in 1:sites){
      
      # IM   : Interaction matrix
      # SM  : Strength matrix
      
      IM=Meta_Interaction_Matrix[,,site]
      
      SM =Meta_Strength_Matrix[,,site]
      
      ## species poplation densities for each of the switching criterion
      
      VSpatialNonSW=Meta_Non_Switching[steps,2:(S+1),site]
      V_elim_sw=Meta_elimination_switch[steps,2:(S+1),site]
      V_opt_sw=Meta_optimization_switch[steps,2:(S+1),site]
      
      #Jacobian per site
      Jac_NonSW=jacobian( Modified_Spatial_LVM_stab,VSpatialNonSW)
      Jacobian_Nonswitch[[site]] = jacobian( Modified_Spatial_LVM_stab,VSpatialNonSW)
      
      Jac_Meta_elimination_switch=jacobian(Modified_Spatial_LVM_stab,V_elim_sw)
      Jacobian_Elim[[site]]=jacobian(Modified_Spatial_LVM_stab,V_elim_sw)
      
      Jac_optimization_Switch=jacobian(Modified_Spatial_LVM_stab,V_opt_sw)
      Jacobian_Opt[[site]]=jacobian(Modified_Spatial_LVM_stab,V_opt_sw)
      
      #Eigenvalues for each of the switching criterion
      eigen_NonSW=eigen(Jac_NonSW)$values
      eigen_elim=eigen(Jac_Meta_elimination_switch)$values
      eigen_opt=eigen(Jac_optimization_Switch)$values
      
      # Selecting the lead eigenvalues.
      Max_eigen_NonSW=data.frame(max(Re(eigen_NonSW)))
      Max_eigen_elim =data.frame(max(Re(eigen_elim)))
      Max_eigen_opt =data.frame(max(Re(eigen_opt)))
      
      max_eigen_NonSw=append(max_eigen_NonSw,Max_eigen_NonSW)
      max_eigen_elim=append(max_eigen_elim,Max_eigen_elim)
      max_eigen_opt=append(max_eigen_opt,Max_eigen_opt)
      
    }
    
    Meta_max_eigen_NonSw = data.frame(cbind(Meta_max_eigen_NonSw,unlist(max_eigen_NonSw)))
    Meta_max_eigen_elim = data.frame(cbind(Meta_max_eigen_elim,unlist(max_eigen_elim)))
    Meta_max_eigen_opt = data.frame(cbind(Meta_max_eigen_opt,unlist(max_eigen_opt)))
    
    
    
    ######################################
    
    #Global stability computation
    
    ##Non switch
    
    site=1
    
    JacobibianNonswitchPatch=list()
    
    while(site< sites){
      
      JacobibianNonswitchPatch[site]=list(Jacobian_Nonswitch[[site]],Jacobian_Nonswitch[[(site+1)]])
      site=site+1
      JacobibianNonswitchPatch[[sites]]<-Jacobian_Nonswitch[[sites]]
      
    }
    
    
    Jacobians_Nonswitch=do.call(adiag,JacobibianNonswitchPatch)
    
    
    ##Elimination
    
    site=1
    
    JacobibianElimPatch=list()
    
    while(site< sites){
      
      JacobibianElimPatch[site]=list(Jacobian_Elim[[site]],Jacobian_Elim[[(site+1)]])
      
      site=site+1
      
      JacobibianElimPatch[[sites]]<-Jacobian_Elim[[sites]]
      
    }
    
    
    Jacobians_Elim=do.call(adiag,JacobibianElimPatch)
    
    
    ##Optimization
    
    site=1
    
    JacobibianOptPatch=list()
    
    while(site< sites){
      JacobibianOptPatch[site]=list(Jacobian_Opt[[site]],Jacobian_Opt[[(site+1)]])
      site=site+1
      JacobibianOptPatch[[sites]]<-Jacobian_Opt[[sites]]
      
    }
    
    
    Jacobians_Opt=do.call(adiag,JacobibianOptPatch)
    
    Intra_specific_Matrix =diag(nrow=(S*sites))
    
    Interpatch_dispersal = `diag<-`(matrix(0,S,S),d/(sites-1))
    Inter_patchdispersal=matrix(rep(t(Interpatch_dispersal),sites),nrow = S*sites,byrow=T)
    Dispersal_matrix = matrix(Inter_patchdispersal,nrow = S*sites,ncol = S*sites)
    
    diag(Dispersal_matrix)<--d
    
    
    Global_Jacobian_Nonswitch=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Nonswitch
    Global_eigen_Nonswitch = eigen(Global_Jacobian_Nonswitch)$values
    Lead_Eigen_Nonswitch = data.frame(max(Re(Global_eigen_Nonswitch)))
    Global_maxeigen_Nonswitch = append(Global_maxeigen_Nonswitch,Lead_Eigen_Nonswitch)
    
    Global_Jacobian_Elim=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Elim
    Global_eigen_Elim = eigen(Global_Jacobian_Elim)$values
    Lead_Eigen_Elim = data.frame(max(Re(Global_eigen_Elim)))
    Global_maxeigen_Elim = append(Global_maxeigen_Elim,Lead_Eigen_Elim)
    
    Global_Jacobian_Opt=-Intra_specific_Matrix+Dispersal_matrix +Jacobians_Opt
    Global_eigen_Opt = eigen(Global_Jacobian_Opt)$values
    Lead_Eigen_Opt = data.frame(max(Re(Global_eigen_Opt)))
    Global_maxeigen_Opt = append(Global_maxeigen_Opt,Lead_Eigen_Opt)
    
    
  } ## End of Connectance loop
  
  
  Global_Productivity = data.frame(cbind(Total_Abundance_Nonswitch,Total_Abundance_Elim,Total_Abundance_Opt))
  Total_Abundance_Nonswitch$S=S; Total_Abundance_Elim$S = S; Total_Abundance_Opt$S=S
  Total_Abundance_Nonswitch_datalist[[S]] =Total_Abundance_Nonswitch
  Total_Abundance_Elim_datalist[[S]] = Total_Abundance_Elim
  Total_Abundance_Opt_datalist[[S]] = Total_Abundance_Opt
  Total_Abundance_Nonswitch_data= bind_rows(Total_Abundance_Nonswitch_datalist)
  Total_Abundance_Elim_data= bind_rows(Total_Abundance_Elim_datalist)
  Total_Abundance_Opt_data = bind_rows(Total_Abundance_Opt_datalist)
  
  
  Global_Productivity$S = S
  Global_Productivity_datalist[[S]]=Global_Productivity
  Global_Productivity_data=bind_rows(Global_Productivity_datalist)
  
  
  Local_Stability = data.frame(cbind(stability_Nonswitchsite=t(Meta_max_eigen_NonSw),
                                     stability_Elimsite=t(Meta_max_eigen_elim),stability_Optsite=t(Meta_max_eigen_opt)))
  
  Local_Stability$S = S
  Local_Stability_datalist[[S]]=Local_Stability
  Local_Stability_data=bind_rows(Local_Stability_datalist)
  
  Global_Stability=data.frame(rbind(unlist(Global_maxeigen_Nonswitch),unlist(Global_maxeigen_Elim),unlist(Global_maxeigen_Opt)))
  
  Meta_Stability=t(Global_Stability)
  
  row.names(Meta_Stability)<-NULL
  GLOBAL_META_STABILITY<- data.frame(cbind(Connectance=con_seq,Meta_Stability))
  
  GLOBAL_META_STABILITY$S = S
  GlobalStability_datalist[[S]]=GLOBAL_META_STABILITY
  GlobalStability_data=bind_rows(GlobalStability_datalist)
  
  Abundance_Connectance = data.frame(cbind(Connectance=con_seq ,Global_Productivity,Local_Stability, t(Global_Stability)))
  row.names(Abundance_Connectance)<- NULL
  
  Abundance_Connectance$S = S
  Abundance_datalist[[S]]=Abundance_Connectance
  AbundanceCommunity_data=bind_rows(Abundance_datalist)
  
  
  ##Productivity per site
  
  for (Patch in 1:sites){
    
    Productivity_Site = cbind(Abundance_Connectance[1], Total_Abundance_Nonswitch[,Patch],Total_Abundance_Elim[,Patch],Total_Abundance_Opt[,Patch])
    ProductivitySite_melted=melt(Productivity_Site,id.vars = "Connectance")
    ggplot(data=ProductivitySite_melted,aes(x=Connectance,y=value, group=variable))+
      geom_point(aes(shape=variable,color=variable),size=5)+
      geom_smooth(method = "loess",se=F,aes(color=variable))+
      theme_few()+
      labs(y="Productivity",x="Connectance")+
      theme(legend.position = "None",aspect.ratio = 1)+
      theme(axis.title=element_text(size=20,face="plain"),
            axis.text =element_text(face = "plain",size=15))
    
    ggsave(paste0("Productivitysite_",toString(c(Patch,S)),".pdf"),width = 5,height = 5)
    ggsave(paste0("Productivitysite_",toString(c(Patch,S)),".tiff"),width = 5,height = 5)
    
    LocalStability_Site = cbind(Abundance_Connectance[1], Local_Stability[c((Patch),(Patch+sites),(Patch+(sites*2)))])
    LocalStabilitySite_melted=melt(LocalStability_Site,id.vars = "Connectance")
    ggplot(data=LocalStabilitySite_melted,aes(x=Connectance,y=value, group=variable))+
      geom_hline(aes(yintercept=0),size=1,colour='red',linetype='dotted')+
      geom_point(aes(shape=variable,color=variable),size=5)+
      geom_smooth(method = "loess",se=F,aes(color=variable))+
      theme_few()+
      labs(y=expression("Re"~(lambda["Max"])),x="Connectance")+
      theme(legend.position = "None",aspect.ratio = 1)+
      theme(axis.title=element_text(size=20,face="plain"),
            axis.text =element_text(face = "plain",size=15))
    
    ggsave(paste0("localstabilitysite_",toString(c(Patch,S)),".pdf"),width = 5,height = 5)
    ggsave(paste0("localstabilitysite_",toString(c(Patch,S)),".tiff"),width = 5,height = 5)
    
  }
  
  
  ##Global Productivity 
  
  Productivity= cbind(Abundance_Connectance[1],Total_Abundance_Nonswitch[,(sites+1)],Total_Abundance_Elim[,(sites+1)],Total_Abundance_Opt[,(sites+1)])
  Productivity_melted=melt(Productivity,id.vars = "Connectance")
  ggplot(data=Productivity_melted,aes(x=Connectance,y=value, group=variable))+
    geom_point(aes(shape=variable,color=variable),size=5)+
    geom_smooth(method = "loess",se=F,aes(color=variable))+
    theme_few()+
    labs(y="Productivity",x="Connectance")+
    theme(legend.position = "None",aspect.ratio = 1)+
    #       legend.key.width = unit(1.5,"cm"))+
    theme(axis.title=element_text(size=20,face="plain"),
          axis.text =element_text(face = "plain",size=15))
  ggsave(paste0("Productivity_",toString(S),".pdf"),width = 5,height = 5)
  ggsave(paste0("Productivity_",toString(S),".tiff"),width = 5,height = 5)
  
  
  ##Global Stability
  GlobalStability_melted=melt(GLOBAL_META_STABILITY,id.vars = c("Connectance","S"))
  ggplot(data=GlobalStability_melted,aes(x=Connectance,y=value,group=variable))+
    #geom_hline(aes(yintercept=0),size=1,colour='red',linetype='dotted')+
    geom_point(aes(shape=variable,color=variable),size=5)+
    geom_smooth(method = "loess",se=F,aes(color=variable))+
    theme_few()+
    labs(y=expression("Re"~(lambda["Max"])),x="Connectance")+
    theme(legend.position = "None",aspect.ratio = 1)+
    theme(axis.title=element_text(size=20,face="plain"),
          axis.text =element_text(face = "plain",size=15))
  ggsave(paste0("GlobalStability_",toString(S),".pdf"),width = 5,height = 5)
  ggsave(paste0("GlobalStability_",toString(S),".tiff"),width = 5,height = 5)
  
}# End of species number loop

##Binding the all species data together
AbundanceCommunity_data=bind_rows(Abundance_datalist)
GlobalStability_data=bind_rows(GlobalStability_datalist)
Global_Productivity_data=bind_rows(Global_Productivity_datalist)
Local_Stability_data=bind_rows(Local_Stability_datalist)
Total_Abundance_Nonswitch_data= bind_rows(Total_Abundance_Nonswitch_datalist)
Total_Abundance_Elim_data= bind_rows(Total_Abundance_Elim_datalist)
Total_Abundance_Opt_data = bind_rows(Total_Abundance_Opt_datalist)




##Saving the data to a folder
save(Meta_Non_Switching,Meta_elimination_switch,Meta_optimization_switch,Abundance_Connectance,
     GlobalStability_data,AbundanceCommunity_data,Local_Stability_data,Global_Productivity_data,
     Total_Abundance_Elim_data,Total_Abundance_Opt_data,Total_Abundance_Nonswitch_data, file="FullMetacommunity.RData")

##
#load("FullMetacommunity.RData")
