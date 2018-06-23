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

source("E:/MSC/Code Repo/MastersCode/Spatial model functions.R")

setwd("E:/MSC/Code Repo/MastersCode/Meta/Results 2")

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

SPECIES=10 #c(50,100)
for (S in SPECIES ){
p=S+1
sites=3

int_growth=runif(S,min =0.5,max =0.5)
carry_cap= runif(S, min=1, max=1)
#Original_pop=matrix(runif(S*sites,min =0,max =1), nrow=S, ncol=sites)
#pop=Original_pop[,site]
  # if(site==1){
  # pop=Original_pop[,1] 
  # } else if (site==2){
  # pop=Original_pop[,2]
  # } else {pop=Original_pop[,3]}


#pop=Original_pop[1];pop=Original_pop[2];pop=Original_pop[3]
pop=runif(S,min =0,max =1)

Original_pop=replicate(sites,pop)

d=0.2

Emigration_Prop=matrix(runif(sites*sites,min =d,max =d), nrow = sites, ncol=sites)
#  for (emmirate in 1:sites) {
#    Emigration_Prop[emmirate,emmirate]=0
# }


Mortality= matrix(runif(sites*sites,min =0.75,max =0.75), nrow = sites)
 for (emmirate in 1:sites) {
   Mortality[emmirate,emmirate]=-1
 }

#(Emigration_Prop*Mortality)%*%t(Original_pop)


steps=10001
#steps=nsteps+1
ww=replicate(101,0)
cc=c(0:100)
sd=0.1
##D for str_mat (competition strength matrix)
str_mat = abs(matrix(rnorm(S*S,0,sd),nrow=S));diag(str_mat)=1
Meta_Strength_Matrix= array(replicate(sites,str_mat), dim=c(S,S,sites))

# Meta_Strength_Matrix= array(rep(abs(rnorm(S*S*sites,0,0.2)),sites),dim = c(S,S,sites))
# 
# for (j in 1:S){
#   for(patch1 in 1:sites){
#     Meta_Strength_Matrix[j,j,patch1]=1
#   }
# }

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
  #######################
  # Meta_Interaction_Matrix=array(rep(rbinom(S*S*sites,1,con_seq[z]),sites), dim=c(S,S,sites))
  # 
  # for (locality in 1:sites){
  #   Local_Meta_Interaction_Matrix=Meta_Interaction_Matrix[,,locality]
  #   ind<-lower.tri(Local_Meta_Interaction_Matrix)
  #   Local_Meta_Interaction_Matrix[ind]<-t(Local_Meta_Interaction_Matrix)[ind]
  #   Meta_Interaction_Matrix[,,locality]=Local_Meta_Interaction_Matrix
  # }
  # MM=Meta_Interaction_Matrix
  # 
  # for (meta1 in 1:S){for(patch_int in 1:sites){MM[meta1,meta1,patch_int]=0 }}
  # 
  # for (meta in 1:S){
  #   for(patch_int in 1:sites){
  #     Meta_Interaction_Matrix[meta,meta,patch_int]=1
  #   }
  # }
  ########################## 
  #######################################################################################################
  Meta_Non_Switching=Spatial_euler_meth(0,pop, 0.01, 100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site)
  
  #######################################################################################################
  
  Meta_Nonswitch_TAbundnance = c()
  
  
  for(patch in 1: sites){
    Total_AbundNonswitch = c()
    for(i in 1:length(Meta_Non_Switching[,1,1])){
      
      meta_abund_nonswitch=sum(Meta_Non_Switching[i,2:(S+1), patch])
      Total_AbundNonswitch=c(Total_AbundNonswitch,meta_abund_nonswitch)
    }
    
    Meta_Nonswitch_TAbundnance=data.frame(cbind(Meta_Nonswitch_TAbundnance,Total_AbundNonswitch))
    
  }
  
  connect_abundance=data.frame(rbind(connect_abundance,Meta_Nonswitch_TAbundnance[steps,]))
  Total_Abundance_Nonswitch =transform(connect_abundance,TAbundance_Nonswitch=rowSums(connect_abundance[,1:sites]))
  ########################################################################################################
  Meta_elimination_switch=Meta_euler_meth_elimination_switch(0, pop, 0.01,100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site)
  ########################################################################################################
  
  Meta_Elimination_TAbundnance = c()
  
  
  for(patch in 1: sites){
    Total_AbundElim = c()
    for(i in 1:length(Meta_elimination_switch[,1,1])){
      
      meta_abund_elimation=sum(Meta_elimination_switch[i,2:(S+1), patch])
      Total_AbundElim=c(Total_AbundElim,meta_abund_elimation)
    }
    
    Meta_Elimination_TAbundnance=data.frame(cbind(Meta_Elimination_TAbundnance,Total_AbundElim))
    
  }
  
  connect_abund_elim_switch=data.frame(rbind(connect_abund_elim_switch, Meta_Elimination_TAbundnance[steps,]))
  Total_Abundance_Elim =transform(connect_abund_elim_switch,TAbundance_Elim=rowSums(connect_abund_elim_switch[,1:sites]))
  ######################################################################################################
  # ## Optimization switching population dynamics plots over the time step for every connectance value.
  # ## Plots saved in a folder.
  
  Meta_optimization_switch=Meta_euler_meth_optimization_switch(0, pop, 0.01,100,Meta_Interaction_Matrix,Meta_Strength_Matrix,carry_cap,int_growth,Original_pop,Emigration_Prop,site)
  ######################################################################################################
  
  Meta_Optimization_TAbundnance = c()
  
  
  for(patch in 1: sites){
    Total_AbundOpt = c()
    for(i in 1:length(Meta_optimization_switch[,1,1])){
      
      meta_abund_Opt=sum(Meta_optimization_switch[i,2:(S+1), patch])
      Total_AbundOpt=c(Total_AbundOpt,meta_abund_Opt)
    }
    
    Meta_Optimization_TAbundnance=data.frame(cbind(Meta_Optimization_TAbundnance,Total_AbundOpt))
    
  }
  
  connect_abund_opt_switch=data.frame(rbind(connect_abund_opt_switch,Meta_Optimization_TAbundnance[steps,]))
  Total_Abundance_Opt =transform(connect_abund_opt_switch,TAbundance_Opt=rowSums(connect_abund_opt_switch[,1:sites]))
  #Total abunance over time step
  
  Total_abund= data.frame(cbind(Meta_Non_Switching[,1,1], Meta_Nonswitch_TAbundnance,Meta_Elimination_TAbundnance, Meta_Optimization_TAbundnance))
  Total_Abund=plyr::rename(Total_abund,c(Meta_Non_Switching...1..1. ="time", Total_AbundNonswitch="TotalAbundance_NonSwitchS1",Total_AbundNonswitch.1="TotalAbundance_NonSwitchS2",Total_AbundNonswitch.2="TotalAbundance_NonSwitchS3",
                                         Total_AbundElim ="Total_AbundanceElim_S1", Total_AbundElim.1 ="Total_AbundanceElim_S2",Total_AbundElim.2 ="Total_AbundanceElim_S3",
                                         Total_AbundOpt = "Total_AbundanceOpt_S1",Total_AbundOpt.1 = "Total_AbundanceOpt_S2",Total_AbundOpt.2 = "Total_AbundanceOpt_S3"))
  
  for (site in 1:sites){
    Total_Abund_site = cbind(Total_Abund[1],Meta_Nonswitch_TAbundnance[,site],Meta_Elimination_TAbundnance[,site],Meta_Optimization_TAbundnance[,site])
    ##################################
    Total_Abund_site_melted=melt(Total_Abund_site,id.vars = "time")
    # ggplot(data=Total_Abund_site_melted,aes(x=time,y=value, group=variable))+
    #   #geom_hline(aes(yintercept=0),size=1,colour='red',linetype="dotted")+
    #   geom_line(aes(color=variable),size=1)+
    #   scale_color_discrete(breaks=c("TotalAbundance_NonSwitchS1","Total_AbundanceElim_S1","Total_AbundanceOpt_S1"),
    #                        labels=c("Non-switch","Elimination", "Optimization"))+
    #   theme_few()+
    # 
    #   labs(y="Total Abundance",x="Time",aspect.ratio = 2)+
    #   theme(legend.position = c(0.85,0.15),legend.title = element_blank(), aspect.ratio = 1)+
    #   theme(axis.title=element_text(size=12,face="bold"),
    #         axis.text =element_text(face = "bold",size=12))
    # ggsave(paste0("TotalAbundance_",toString(c(site,z)),".pdf"),width = 5,height = 5)
    #######################3
  }
  
  
  ##Stability computatation
  
  max_eigen_NonSw=c();max_eigen_elim=c();max_eigen_NonSpatialSw=c()
  max_eigen_opt=c();Max_eigen_NonSW=c();Max_eigen_NonSpatialSw=c()
  Max_eigen_elim=c();Max_eigen_opt=c()
  Jacobian_Nonswitch = list();Jacobian_Nonspatial=list()
  Jacobian_Elim = list();Jacobian_Opt=list()
  
  for (site in 1:sites){
    
    IM=Meta_Interaction_Matrix[,,site]
    SM =Meta_Strength_Matrix[,,site]
    ######################################################################################
    VSpatialNonSW=Meta_Non_Switching[steps,2:(S+1),site]
    
    V_elim_sw=Meta_elimination_switch[steps,2:(S+1),site]
    V_opt_sw=Meta_optimization_switch[steps,2:(S+1),site]
    
    Jac_NonSW=jacobian( Modified_Spatial_LVM_stab,VSpatialNonSW)
    Jacobian_Nonswitch[[site]] = jacobian( Modified_Spatial_LVM_stab,VSpatialNonSW)
    
    Jac_Meta_elimination_switch=jacobian(Modified_Spatial_LVM_stab,V_elim_sw)
    Jacobian_Elim[[site]]=jacobian(Modified_Spatial_LVM_stab,V_elim_sw)
    
    Jac_optimization_Switch=jacobian(Modified_Spatial_LVM_stab,V_opt_sw)
    Jacobian_Opt[[site]]=jacobian(Modified_Spatial_LVM_stab,V_opt_sw)
    
    eigen_NonSW=eigen(Jac_NonSW)$values
    eigen_elim=eigen(Jac_Meta_elimination_switch)$values
    eigen_opt=eigen(Jac_optimization_Switch)$values
    
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
  
  
  
  ###########################################################################
  ### Commputing stability over time for different values of connectance
  Meta_elim_max_eigen=c()
  Meta_opt_max_eigen=c()
  Meta_max_eigen_NonSWT=c()
  
  for (site in 1:sites){
    
    Max_eigen_NonSWT=c()
    max_eigen_NonSwT=c()
    
    for(nsw in 1:length(Meta_Non_Switching[,1,1])){
      
      VSpatialNonSWT=Meta_Non_Switching[nsw,2:(S+1),site]
      
      Jac_NonSWT=jacobian(Modified_Spatial_LVM_stab,VSpatialNonSWT)
      eigen_NonSWT=eigen(Jac_NonSWT)$values
      Max_eigen_NonSWT=data.frame(max(Re(eigen_NonSWT)))
      max_eigen_NonSwT=append(max_eigen_NonSwT,Max_eigen_NonSWT)
    }
    Meta_max_eigen_NonSWT = data.frame(cbind(Meta_max_eigen_NonSWT,unlist(max_eigen_NonSwT)))
  }
  
  Meta_max_eigen_NonSWT_Stab=plyr::rename(Meta_max_eigen_NonSWT,c(cbind.Meta_max_eigen_NonSWT..unlist.max_eigen_NonSwT..="StabNonSwitchsite1",
                                                                  unlist.max_eigen_NonSwT.="StabNonSwitchsite2",unlist.max_eigen_NonSwT..1="StabNonSwitchsite3"))
  
  for (site in 1:sites){
    
    elim_max_eigen=c()
    elim_Max_eigen=c()
    
    
    for(elim_sw in 1:length(Meta_elimination_switch[,1,site])){
      
      V_elim=Meta_elimination_switch[elim_sw,2:(S+1),site]
      ##########################
      Jac_elim=jacobian(Modified_Spatial_LVM_stab,V_elim)
      elim_eigen=eigen(Jac_elim)$values
      elim_Max_eigen=data.frame(max(Re(elim_eigen)))
      elim_max_eigen=append(elim_max_eigen,elim_Max_eigen)
      
    }
    Meta_elim_max_eigen = data.frame(cbind(Meta_elim_max_eigen,unlist(elim_max_eigen)))
  }
  
  Meta_elim_max_eigen_Stab=plyr::rename(Meta_elim_max_eigen,c(cbind.Meta_elim_max_eigen..unlist.elim_max_eigen..="StabElimsite1",
                                                              unlist.elim_max_eigen.="StabElimsite2",unlist.elim_max_eigen..1="StabElimsite3"))
  
  for (site in 1:sites){
    
    opt_max_eigen =c()
    opt_Max_eigen = c()
    
    for(opt_sw in 1:length(Meta_optimization_switch[,1,site])){
      
      V_opt=Meta_optimization_switch[opt_sw,2:(S+1),site]
      ##########################
      Jac_opt=jacobian(Modified_Spatial_LVM_stab,V_opt)
      opt_eigen=eigen(Jac_opt)$values
      opt_Max_eigen=data.frame(max(Re(opt_eigen)))
      opt_max_eigen=append(opt_max_eigen,opt_Max_eigen)
    }
    
    Meta_opt_max_eigen = data.frame(cbind(Meta_opt_max_eigen, unlist(opt_max_eigen)))
    
  }
  
  
  Meta_opt_max_eigen_Stab=plyr::rename(Meta_opt_max_eigen,c(cbind.Meta_opt_max_eigen..unlist.opt_max_eigen..="StabOptsite1",
                                                            unlist.opt_max_eigen.="StabOptsite2",unlist.opt_max_eigen..1="StabOptsite3"))
  
  Stabb_time= data.frame(cbind(Meta_Non_Switching[,1,1],Meta_max_eigen_NonSWT_Stab,Meta_elim_max_eigen_Stab,Meta_opt_max_eigen_Stab))
  Stab_time=plyr::rename(Stabb_time,c(Meta_Non_Switching...1..1.="time"))

  for (site in 1:sites){
    Stab_time_site = cbind(Stab_time[1],Meta_max_eigen_NonSWT_Stab[,site],Meta_elim_max_eigen_Stab[,site],Meta_opt_max_eigen_Stab[,site])
    ############################################
    Stab_time_site1_melted=melt(Stab_time_site,id.vars = "time")
    #  ggplot(data=Stab_time_site1_melted,aes(x=time,y=value, group=variable))+
    #    geom_hline(aes(yintercept=0),size=1,colour='red',linetype="dotted")+
    #    geom_line(aes(color=variable),size=1)+
    #    scale_color_discrete(breaks=c("StabNonSwitchsite1","StabElimsite1","StabOptsite1"),
    #                         labels=c("Non-switch","Elimination", "Optimization"))+
    #    theme_few()+
    # 
    #    labs(y=expression("Re"~(lambda)),x="Time",aspect.ratio = 2)+
    #    theme(legend.position = c(0.85,0.85),legend.title = element_blank(), aspect.ratio = 1)+
    #    theme(axis.title=element_text(size=12,face="bold"),
    #         axis.text =element_text(face = "bold",size=12))
    # ggsave(paste0("Stabtime_",toString(c(site,z)),".pdf"),width = 5,height = 5)
################################
  }
  
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
  
  #Jacobians_Nonswitch = adiag(Jacobian_Nonswitch[[1]],Jacobian_Nonswitch[[2]],Jacobian_Nonswitch[[3]])
  
  ##Elimination
  site=1
  JacobibianElimPatch=list()
  while(site< sites){
    JacobibianElimPatch[site]=list(Jacobian_Elim[[site]],Jacobian_Elim[[(site+1)]])
    site=site+1
    JacobibianElimPatch[[sites]]<-Jacobian_Elim[[sites]]
  }
  
  Jacobians_Elim=do.call(adiag,JacobibianElimPatch)
  #Jacobians_Elim = adiag(Jacobian_Elim[[1]],Jacobian_Elim[[2]],Jacobian_Elim[[3]])
  
  ##Optimization
  site=1
  JacobibianOptPatch=list()
  while(site< sites){
    JacobibianOptPatch[site]=list(Jacobian_Opt[[site]],Jacobian_Opt[[(site+1)]])
    site=site+1
    JacobibianOptPatch[[sites]]<-Jacobian_Opt[[sites]]
  }
  
  Jacobians_Opt=do.call(adiag,JacobibianOptPatch)
  
  #Jacobians_Opt = adiag(Jacobian_Opt[[1]],Jacobian_Opt[[2]],Jacobian_Opt[[3]])
  
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

#Global Productivity For each connectance value
#Global_TotalAbundance_Nonswitch = transform(connect_abundance,GlobalAbundance_Noswitch=rowSums(connect_abundance[1:sites]))
#Global_TotalAbundance_Elimination =transform(connect_abund_elim_switch,GlobalAbundance_Elim= rowSums(connect_abund_elim_switch[1:sites]))
#Global_TotalAbundance_Optimization = transform(connect_abund_opt_switch,GlobalAbundance_Opt=rowSums(connect_abund_opt_switch[1:sites]))

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
#global_StabNospatial=t(unlist(Global_maxeigen_Nonspatial)),global_StabNonswitch=t(unlist(Global_maxeigen_Nonswitch)),
#global_StabEli=t(unlist(Global_maxeigen_Elim)),Global_StabOpt=t(unlist(Global_maxeigen_Opt))))

row.names(Abundance_Connectance)<- NULL

Abundance_Connectance$S = S
Abundance_datalist[[S]]=Abundance_Connectance
AbundanceCommunity_data=bind_rows(Abundance_datalist)

##Productivity per site
for (Patch in 1:sites){
  
  Productivity_Site = cbind(Abundance_Connectance[1], Total_Abundance_Nonswitch[,Patch],Total_Abundance_Elim[,Patch],Total_Abundance_Opt[,Patch])
  ProductivitySite_melted=melt(Productivity_Site,id.vars = "Connectance")
  ggplot(data=ProductivitySite_melted,aes(x=Connectance,y=value, group=variable))+
    geom_point(aes(shape=variable,color=variable),size=2)+
    geom_smooth(method = "loess",se=F,aes(color=variable))+
    theme_few()+
    labs(y="Productivity",x="Connectance")+
    theme(legend.position = "None",aspect.ratio = 1)+
    theme(axis.title=element_text(size=12,face="bold"),
          axis.text =element_text(face = "bold",size=12))
  
  ggsave(paste0("Productivitysite_",toString(c(Patch,S)),".pdf"),width = 5,height = 5)
  ggsave(paste0("Productivitysite_",toString(c(Patch,S)),".tiff"),width = 5,height = 5)
  
  LocalStability_Site = cbind(Abundance_Connectance[1], Local_Stability[c((Patch),(Patch+sites),(Patch+(sites*2)))])
  LocalStabilitySite_melted=melt(LocalStability_Site,id.vars = "Connectance")
  ggplot(data=LocalStabilitySite_melted,aes(x=Connectance,y=value, group=variable))+
    geom_hline(aes(yintercept=0),size=1,colour='red',linetype='dotted')+
    geom_point(aes(shape=variable,color=variable),size=2)+
    geom_smooth(method = "loess",se=F,aes(color=variable))+
    theme_few()+
    labs(y=expression("Re"~(lambda)),x="Connectance")+
    theme(legend.position = "None",aspect.ratio = 1)+
    theme(axis.title=element_text(size=12,face="bold"),
          axis.text =element_text(face = "bold",size=12))
  
  ggsave(paste0("localstabilitysite_",toString(c(Patch,S)),".pdf"),width = 5,height = 5)
  ggsave(paste0("localstabilitysite_",toString(c(Patch,S)),".tiff"),width = 5,height = 5)
  
}

##Global Productivity 

Productivity= cbind(Abundance_Connectance[1],Total_Abundance_Nonswitch[,(sites+1)],Total_Abundance_Elim[,(sites+1)],Total_Abundance_Opt[,(sites+1)])
Productivity_melted=melt(Productivity,id.vars = "Connectance")
ggplot(data=Productivity_melted,aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y="Productivity",x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))
ggsave(paste0("Productivity_",toString(S),".pdf"),width = 5,height = 5)
ggsave(paste0("Productivity_",toString(S),".pdf"),width = 5,height = 5)

##Global Stability
#GlobalStability= cbind(Abundance_Connectance[1],Abundance_Connectance[23:25])
GlobalStability_melted=melt(GLOBAL_META_STABILITY,id.vars = c("Connectance","S"))
ggplot(data=GlobalStability_melted,aes(x=Connectance,y=value,group=variable))+
  #geom_hline(aes(yintercept=0),size=1,colour='red',linetype='dotted')+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y=expression("Re"~(lambda)),x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))
ggsave(paste0("GlobalStability_",toString(S),".pdf"),width = 5,height = 5)
ggsave(paste0("GlobalStability_",toString(S),".pdf"),width = 5,height = 5)


}# End of species number loop


##Binding the all species data together
AbundanceCommunity_data=bind_rows(Abundance_datalist)
GlobalStability_data=bind_rows(GlobalStability_datalist)
Global_Productivity_data=bind_rows(Global_Productivity_datalist)
Local_Stability_data=bind_rows(Local_Stability_datalist)
Total_Abundance_Nonswitch_data= bind_rows(Total_Abundance_Nonswitch_datalist)
Total_Abundance_Elim_data= bind_rows(Total_Abundance_Elim_datalist)
Total_Abundance_Opt_data = bind_rows(Total_Abundance_Opt_datalist)


GlobalStability_data$MayRule = with(GlobalStability_data, sd* sqrt(Connectance*S))
GlobalStability_data$GravelRule = with(GlobalStability_data,(sd* sqrt(Connectance*(S-1)))-d)

GlobalStabilityMelted_data = melt(GlobalStability_data, id.vars = c("Connectance","S","MayRule","GravelRule"))
ggplot(GlobalStabilityMelted_data, aes(x=MayRule,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y=expression("Re"~(lambda)),x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))
ggsave("MaysGlobalStability.pdf",width = 5,height = 5)
ggsave("MaysGlobalStability.pdf",width = 5,height = 5)

ggplot(GlobalStabilityMelted_data, aes(x=GravelRule,y=value, group=variable))+
  geom_point(aes(shape=variable,color=variable),size=2)+
  geom_smooth(method = "loess",se=F,aes(color=variable))+
  theme_few()+
  labs(y=expression("Re"~(lambda)),x="Connectance")+
  theme(legend.position = "None",aspect.ratio = 1)+
  theme(axis.title=element_text(size=12,face="bold"),
        axis.text =element_text(face = "bold",size=12))
ggsave("GravelGlobalStability.pdf",width = 5,height = 5)
ggsave("GravelGlobalStability.pdf",width = 5,height = 5)

##Saving the data to a folder
save(Meta_Non_Switching,Meta_elimination_switch,Meta_optimization_switch,Abundance_Connectance,
     GlobalStability_data,AbundanceCommunity_data,Local_Stability_data,Global_Productivity_data,
     Total_Abundance_Elim_data,Total_Abundance_Opt_data,Total_Abundance_Nonswitch_data, file="Sites.RData")
#save(Modified_NoSpatial_LVM,Modified_Spatial_LVM_stab,NoSpatial_LVM, Non_euler_meth, Spatial_LVM,Spatial_euler_meth,Meta_euler_meth_elimination_switch,Meta_euler_meth_optimization_switch,file="Metafunctions.RData")

load(file="Sites.RData")
