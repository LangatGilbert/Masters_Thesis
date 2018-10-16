#########################################################################################
#CODE FOR RESEARCH MASTERS THESIS
#PROJECT: Implementattion of Modified Lotka-Volterra model on Euler method for SINGLE COMMUNITY
#By: Gilbert K Langat



##############################################
#Clearing R memory and Saving code to folder
##############################################

rm(list=ls(all=TRUE)) 


################################################################################
## IMPORTING FUNCTIONs 
################################################################################

source("E:/MSC/Code Repo/MastersCode/Single_Community_Functions.R")

##SET WORKING DIRECTORY, choose the path to where you have saved your work

setwd("E:/MSC/Code Repo/MastersCode/Final Results Part 1/Pubs/sample/Final")

#set.seed(123)
##########################
# LOAD R LIBRARIES

library(BiodiversityR)  # For nestedness computation.
library(pracma)         # Numerical computations 
library(vegan)          # computing nestedness
library(bipartite)      # Community structuring 
library(tidyverse)      # Data manipulation and visualization
library(ggthemes)       # Data visualization 
library(reshape2)       # for organizing data
library("plotrix")      # for plotting 
library("Matrix")       # Matrix manipulation

#####################################################################################
#*ComNetworks - No. of Networks
#####################################################################################

CommNetworks = 50

##Creating Dummy Network matrices

Net_Community <- function(S,CommNetworks){
  Netmat <- array(0, dim=c(S,S,CommNetworks))
  return(Netmat)
}

##Appending containers

Network_Abundsigma_datalist=list()    ;Network_Compet_sigma_datalist1=list()
Network_Compet_sigma_datalist2=list() ;Network_Species_rankdatalist=list()

######################################################################################
##LOOP THROUGH ALL THE NETWORKS
# Ntwork - each network
######################################################################################
for(Ntwork in 1:CommNetworks){
  
  ##print(Newtk) - for tracking networks steps.
  print(Ntwork)
  ################################################################################
  ## empty list for storing data at for each network.
  
  CombinedAbund_datalist=list()   ;CombinedCompet_datalist1=list()
  CombinedCompet_datalist2=list() ;CompetCombineddatalist=list()
  Overall_Speciesranking=list()   ;con_abund_NonSw = c()
  con_abund_elimination = c()     ;con_abund_optimization = c()
  
  len =length(seq(0.1,0.9,0.05))
  
  ################################################################################
  #Loop through the  different values of S (number of species/community sizes)
  #S        :   number of species (community size)
  #species  :   reps the different community size groups
  ###############################################################################
  
  species= c(20,50,100)
  
  for (S in c(species)){
    
    # Net_int_mat reps the community interaction matrix for all networks
    # Net_str_mat reps the community strength matrix for all networks
    
    Net_int_mat = Net_Community(S,CommNetworks)
    Net_str_mat = Net_Community(S,CommNetworks)
    
    # Initializing the parameters
    p=S+1
    int_growth=runif(S,min =0.5,max =0.5)
    carry_cap= runif(S, min=1, max=1)
    pop=runif(S,min =0,max =1)
    
    steps=2001
    
    
    ## Empty list for storing data for each community size
    
    Abundsigma_datalist=list()    ;Compet_sigma_datalist1=list()
    Compet_sigma_datalist2=list() ;Species_rankdatalist=list()
    con_abund_elimination=c();con_abund_optimization=c();con_abund_NonSw=c();connect_abundance=c()
    connect_nest_elim_switch= c();connect_nest_opt_switch= c()
    connect_stab_elim_switch= c();connect_stab_opt_switch= c()
    
    connect_nested= c(); max_eigen_NonSw=c();max_eigen_elim=c()
    max_eigen_opt=c();Max_eigen_NonSW=c(); Max_eigen_elim=c();Max_eigen_opt=c()
    connect_compet_nonswitch_competitiveness=c(); connect_compet_elim_competitiveness=c()
    connect_compet_opt_competitiveness=c();connect_compet_nonswitch_endurance=c()
    connect_compet_elim_endurance =c(); connect_compet_opt_endurance =c()
    
    ###############################################################################
    # Loop over different values of sigma(Standard deviation)
    # SD        : standard deviation
    #*std_dvtn  : SD sequence, from 0.1 to 0.3
    ###############################################################################
    
    std_dvtn=c(0.1,0.2,0.3)
    
    for (SD in 1:length(std_dvtn)){
      
      Net_str_mat[,,Ntwork]= matrix(rnorm(S*S,0,std_dvtn[SD]),nrow=S)
      diag(Net_str_mat[,,Ntwork])=0 ##All diagonal values=zero
      OffDiag_Strmat = abs(Net_str_mat[,,Ntwork]) 
      
      diag(Net_str_mat[,,Ntwork])=1
      str_mat=abs(Net_str_mat[,,Ntwork])
      
      
      #For appending data for each connectance value
      con_abund_elimination=c();con_abund_optimization=c();con_abund_NonSw=c(); connect_abundance=c()
      connect_nest_elim_switch= c()
      connect_nest_opt_switch= c();connect_stab_elim_switch= c();connect_stab_opt_switch= c()
      
      connect_nested= c() ;max_eigen_NonSw=c();max_eigen_elim=c()
      max_eigen_opt=c();Max_eigen_NonSW=c();Max_eigen_elim=c();Max_eigen_opt=c()
      connect_compet_nonswitch_competitiveness=c(); connect_compet_elim_competitiveness=c()
      connect_compet_opt_competitiveness=c(); connect_compet_nonswitch_endurance=c()
      connect_compet_elim_endurance =c(); connect_compet_opt_endurance =c()
      
      #############################################################################
      # Loop for geneting species dynamics over time for each connectance value
      # connectz            :  connectance
      #connectance_sequence :  connectance sequence, from 0.1 to 0.9 with step of 0.05
      ############################################################################
      
      connectance_sequence=seq(0.1,0.9,0.05)
      
      for (connectz in connectance_sequence){
        
        Net_int_mat[,,Ntwork]=matrix(rbinom(S*S,1,connectz),S,S)
        
        diag(Net_int_mat[,,Ntwork]) = 1
        int_mat <- Net_int_mat[,,Ntwork]
        ind<-lower.tri(int_mat)
        int_mat[ind]<-t(int_mat)[ind]  ## Ensuring upper triangular matrix == lower triangular matrix
        
        ########################################################################
        # Figures saved in a folder.
        # Euler method implemented over time 0-100 with 0.01 time step.
        # plotting the species dynamics over the time interval.
        # Non-swithing population dynamics plots for all connectance values.
        #######################################################################
        
        Non_Switching_Dynamics=NonSwitch_Euler(0,pop, 0.1, 200,int_mat,str_mat,carry_cap,int_growth)
        
        ##NON_SWITCHIN SPECIES COMPETITIVENESS, ENDURANCE, TOTAL BIOMASS & NESTEDNESS
        
        Non_Switching_row_wise = c() 
        Non_Switching_column_wise = c()
        
        Non_Switching_row_wise=rowSums(int_mat*str_mat*unlist(Non_Switching_Dynamics[steps,2:(S+1)]))
        Non_Switching_column_wise=colSums(int_mat*str_mat*unlist(Non_Switching_Dynamics[steps,2:(S+1)]))
        
        connect_compet_nonswitch_competitiveness=rbind(connect_compet_nonswitch_competitiveness,Non_Switching_row_wise)
        connect_compet_nonswitch_endurance=rbind(connect_compet_nonswitch_endurance,Non_Switching_column_wise)
        
        
        #NonSwitching Biomass per timestep
        NonSwitching_TotalAbundance =transform(Non_Switching_Dynamics, T_abundanceNonsw=rowSums(Non_Switching_Dynamics[2:(S+1)]))
        
        ## NonSwitching Total Biomass for each connectance
        con_abund_NonSw1 = c()
        con_abund_NonSw1=unlist(Non_Switching_Dynamics[steps,2:(S+1)])
        con_abund_NonSw =rbind(con_abund_NonSw,con_abund_NonSw1)
        
        #NonSwitching Nestedness for each connectance value
        initial_nest=c()
        initial_nest = nested(int_mat,method='NODF',rescale=FALSE,normalised = TRUE)
        connect_nested = append(connect_nested,initial_nest)
        
        # Elimination switching Nestedness
        
        NestedElim = c()
        Nested_Elim= nested(int_mat,method = "NODF",rescale = FALSE,normalised = TRUE)
        NestedElim= c(NestedElim,Nested_Elim)
        
        ##################################################################
        # Elimination switching population dynamics plots over the time step for every connectance value.
        # Plots saved in a folder.
        ##################################################################
        
        Elimination_Switching_Dynamics=Elimination_switch_Euler(0, pop, 0.1,200,int_mat, str_mat,NestedElim,carry_cap,int_growth)
        
        ## ElIMINATION SWITCHING SPECIES COMPETITIVENESS, ENDURANCE, TOTAL BIOMASS 
        
        elimination_row_wise=c()
        elimination_column_wise = c()
        
        elimination_row_wise=rowSums(int_mat*str_mat*unlist(Elimination_Switching_Dynamics[steps,2:(S+1)]))
        elimination_column_wise=colSums(int_mat*str_mat*unlist(Elimination_Switching_Dynamics[steps,2:(S+1)]))
        
        connect_compet_elim_competitiveness=rbind(connect_compet_elim_competitiveness,elimination_row_wise)
        connect_compet_elim_endurance=rbind(connect_compet_elim_endurance,elimination_column_wise)
        
        ##Biomass at end of timesteps
        con_abund_elim = c()
        con_abund_elim=unlist(Elimination_Switching_Dynamics[steps,2:(S+1)]) 
        
        ## Total Biomass for each connectance value
        con_abund_elimination =rbind(con_abund_elimination,con_abund_elim)
        
        #Species Total Biomass appended to the species dynamics dataframe.
        TBiomass_Elimination=transform(Elimination_Switching_Dynamics, T_abundanceElim=rowSums(Elimination_Switching_Dynamics[2:(S+1)]))
        
        #Elimination switching Nestedness for each connectance value
        connect_nest_elim_switch= append(connect_nest_elim_switch,(sum(TBiomass_Elimination[,(S+2)]))/2)
        
        #Optimization Switching Nestedness
        Nest_opt = c()
        Nestopt= nested(int_mat,method = "NODF",rescale = FALSE,normalised = TRUE)
        Nest_opt= c(Nest_opt,Nestopt)
        
        
        ### Optimization switching population dynamics plots over the time step for every connectance value.
        
        Optimization_Switching_Dynamics=Optimization_switch_Euler(0, pop, 0.1,200,int_mat, str_mat,Nest_opt,carry_cap,int_growth) 
        
        ######---OPTIMIZATION SWITCHING SPECIES COMPETITIVESS, ENDURANCE, TOTAL BIOMASS & NESTEDNESS --###############
        optimization_row_wise =c()
        optimization_column_wise =c()
        
        optimization_row_wise=rowSums(int_mat*str_mat*unlist(Optimization_Switching_Dynamics[steps,2:(S+1)]))
        optimization_column_wise=colSums(int_mat*str_mat*unlist(Optimization_Switching_Dynamics[steps,2:(S+1)]))
        
        connect_compet_opt_competitiveness=rbind(connect_compet_opt_competitiveness,optimization_row_wise)
        connect_compet_opt_endurance=rbind(connect_compet_opt_endurance,optimization_column_wise)
        
        # Optimization switching biomass for each connectance value
        con_abund_opt = c()
        con_abund_opt=unlist(Optimization_Switching_Dynamics[steps,2:(S+1)])
        con_abund_optimization =rbind(con_abund_optimization,con_abund_opt)
        
        # Species total biomass appended to population dynamics data frame.
        
        TBiomass_Optimization=transform(Optimization_Switching_Dynamics, T_abundanceOpt=rowSums(Optimization_Switching_Dynamics[2:(S+1)]))
        
        # Optimization Switching nestedness for each connectance value
        connect_nest_opt_switch= append(connect_nest_opt_switch,(sum(TBiomass_Optimization[,(S+2)]))/2)
        
        ########--COMMUNITY STABILITY COMPUTATION--##########
        
        IM=int_mat
        SM = str_mat
        
        # Species biomass at end of timesteps for each switching types
        NonSwitchingBiomass=c(do.call("cbind",Non_Switching_Dynamics[steps,2:(S+1)]))
        ElimSwitchingBiomass=c(do.call("cbind",Elimination_Switching_Dynamics[steps,2:(S+1)]))
        OptSwitchingBiomass=c(do.call("cbind",Optimization_Switching_Dynamics[steps,2:(S+1)]))
        
        # THE JACOBIAN 
        Jac_NonSW=jacobian(LVM_stability,NonSwitchingBiomass)
        Jac_elimination_Switch=jacobian(LVM_stability,ElimSwitchingBiomass)
        Jac_optimization_Switch=jacobian(LVM_stability,OptSwitchingBiomass)
        
        # THE EIGENVALUES
        eigen_NonSW=eigen(Jac_NonSW)$values
        eigen_elim=eigen(Jac_elimination_Switch)$values
        eigen_opt=eigen(Jac_optimization_Switch)$values
        
        # MAXIMUM REAL EIGENVALUE
        Max_eigen_NonSW=data.frame(max(Re(eigen_NonSW)))
        Max_eigen_elim =data.frame(max(Re(eigen_elim)))
        Max_eigen_opt =data.frame(max(Re(eigen_opt)))
        
        ##APPEND THE MAXIMUM EIGENVALUE FOR EACH CONNECTANCE VALUE
        max_eigen_NonSw=append(max_eigen_NonSw,Max_eigen_NonSW)
        max_eigen_elim=append(max_eigen_elim,Max_eigen_elim)
        max_eigen_opt=append(max_eigen_opt,Max_eigen_opt)
        
        
      }# End of Connectance Loop
      
      
      #DATAFRAME STORING TOTAL BIOMASS,NESTEDNESS, STABILITY FOR EACH CONNECTANCE VALUE.
      
      BIONESTAB_CONNECTANCE=do.call(rbind, Map(data.frame,Connectance=connectance_sequence, Abundance_Nonswitch=rowSums(con_abund_NonSw),
                                               Abund_elim_switch=rowSums(con_abund_elimination),Abund_opt_switch=rowSums(con_abund_optimization),
                                               Nested_Nonswitch=connect_nested,Nest_elim_switch=connect_nest_elim_switch,Nest_opt_switch=connect_nest_opt_switch,
                                               Stability_NonSW=max_eigen_NonSw,Stability_Elim_sw=max_eigen_elim,Stability_Opt_sw=max_eigen_opt))
      
      BIOM_NEST_STAB_CONNECTANCE=sort.data.frame(BIONESTAB_CONNECTANCE)
      
      BIOM_NEST_STAB_CONNECTANCE$SD = SD
      Abundsigma_datalist[[SD]]=BIOM_NEST_STAB_CONNECTANCE
      Abundsigma_data=bind_rows(Abundsigma_datalist)
      
      #############################################################################################
      
      #TOTAL BIOMASS |NESTEDNESS|STABILITY PLOTS VS CONNECTANCE
      ##TOTAL BIOMASS 
      ############################################################################################# 
      
      Productivity= cbind(BIOM_NEST_STAB_CONNECTANCE[1],BIOM_NEST_STAB_CONNECTANCE[2:4])
      Productivity_melted=melt(Productivity,id.vars = "Connectance")
      
      #Plot of Productivity for each network, community size and strength
      #############################################################################################
      # ggplot(data=Productivity_melted,aes(x=Connectance,y=value, group=variable))+
      # geom_point(aes(shape=variable),size=2.5)+
      # #scale_shape_manual(values = c(1,2,5))+
      # geom_smooth(method = "loess",span=0.85,se=T,aes(linetype=variable,fill=variable),size=1.5, color="black")+
      # scale_linetype_manual(values=c("solid", "dotted","dashed"),
      #                       breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
      #                      labels=c("Nonswitching","Elimination", "Optimization"))+
      # scale_shape_manual(values = c(1,2,5),breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
      #                      labels=c("Nonswitching","Elimination", "Optimization"))+
      # scale_fill_grey(start = 0.3, end=0.9,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
      #                   labels=c("Nonswitching","Elimination", "Optimization")) +
      # theme_few()+
      # labs(y="Total Biomass",x="Connectance",aspect.ratio = 1)+
      # theme(legend.position = c(0.8,0.85), legend.title=element_blank(),aspect.ratio = 1,
      #       legend.key.width = unit(1.5,"cm"))+
      # theme(axis.title=element_text(face="plain"),
      #       axis.text =element_text(face = "plain",colour = "black"))
      # 
      # #Save Total biomass for each network,community size and competition strength
      # #Ntwork: Network ;S : Community size; SD: Competition strength
      # ggsave(paste0("Total Biomass_",toString(c(Ntwork,S,SD)),".pdf"),width = 7.5,height = 7.5)
      # ggsave(paste0("Total Biomass_",toString(c(Ntwork,S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
      ##################################################################################################### 
      ##NESTEDNESS  
      
      Nestedness= cbind(BIOM_NEST_STAB_CONNECTANCE[1],BIOM_NEST_STAB_CONNECTANCE[5:7])
      Nestedness_melted=melt(Nestedness,id.vars = "Connectance")
      
      #Plot of Nestedness for each network, community size and strength
      ######################################################################################################
      # ggplot(data=Nestedness_melted,aes(x=Connectance,y=value, group=variable))+
      #   geom_point(aes(shape=variable),size=2.5)+
      #   geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=1.5,color="black")+
      #   scale_linetype_manual(values=c("solid", "dotted","dashed"),
      #                         breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
      #                         labels=c("Nonswitching","Elimination", "Optimization"))+
      #   scale_shape_manual(values = c(1,2,5),breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
      #                         labels=c("Nonswitching","Elimination", "Optimization"))+
      #   scale_fill_grey(start = 0.3, end=0.9,breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
      #                   labels=c("Nonswitching","Elimination", "Optimization")) +
      #   theme_few()+
      #   labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
      #   theme(legend.position = c(0.15,0.85),legend.title = element_blank(),aspect.ratio = 1,
      #         legend.key.width = unit(1.5,"cm"))+
      #   theme(axis.title=element_text(face="plain"),
      #         axis.text =element_text(face = "plain",colour = "black"))
      # #Save Nestedness for each network,community size and competition strength
      # #Ntwork: Network ;S : Community size; SD: Competition strength
      # ggsave(paste0("Nestedness_",toString(c(Ntwork,S,SD)),".pdf"),width = 7.5,height = 7.5)
      # ggsave(paste0("Nestedness_",toString(c(Ntwork,S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
      #
      
      #################################################################################################
      ## LOCAL STABILITY
      
      Local_Stab= cbind(BIOM_NEST_STAB_CONNECTANCE[1],BIOM_NEST_STAB_CONNECTANCE[8:10])
      Local_Stab_melted=melt(Local_Stab,id.vars = "Connectance")
      
      #Plot of Productivity for each network, community size and strength
      
      ###################################################################################################
      # ggplot(data=Local_Stab_melted,aes(x=Connectance,y=value, group=variable))+
      #   geom_hline(aes(yintercept=0),size=0.5,colour='black',linetype="dotted")+
      #   geom_point(aes(shape=variable),size=2.5)+
      #   geom_smooth(method = "loess",se=T, aes(linetype=variable,fill=variable),size=1.5, color="black")+
      #   scale_linetype_manual(values=c("solid", "dotted","dashed"),
      #                         breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
      #                         labels=c("Nonswitching","Elimination", "Optimization"))+
      #    scale_shape_manual(values = c(1,2,5),breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
      #                         labels=c("Nonswitching","Elimination", "Optimization"))+
      #   scale_fill_grey(start = 0.3, end=0.9,breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
      #                   labels=c("Nonswitching","Elimination", "Optimization")) +
      #   theme_few()+
      #   labs(y=expression("Re"~(lambda["Max"])),x="Connectance",aspect.ratio = 1)+
      #   theme(legend.position = c(0.85,0.15),legend.title = element_blank(), aspect.ratio = 1,
      #         legend.key.width = unit(1.5,"cm"))+
      #   theme(axis.title=element_text(face="plain"),
      #         axis.text =element_text(face = "plain",colour = "black"))
      # #Save stability for each network,community size and competition strength
      # #Ntwork: Network ;S : Community size; SD: Competition strength
      # ggsave(paste0("LocalStabilityNetwork_",toString(c(Ntwork,S,SD)),".pdf"),width = 7.5,height = 7.5)
      # ggsave(paste0("LocalStabilityNetwork_",toString(c(Ntwork,S,SD)),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
      #  
      ############################################################################################
      
      # Species competitiveness computation. For each connectance values stored in a dataframe.
      #Competitiveness computed by taking the sum of the row values of comunity matrix for species competition from species and
      #   column sum for species competion effect to other species.
      ############################################################################################
      #Combinining Nonswitching, elimination and optimization average population abundances at for each connectance value           
      
      Compet_Abunda= do.call(rbind, Map(data.frame, Non_Switching_Competitiveness = colMeans(connect_compet_nonswitch_competitiveness),
                                        Non_Switching_Endurance = colMeans(connect_compet_nonswitch_endurance),
                                        Elimination_Competitiveness = colMeans(connect_compet_elim_competitiveness),
                                        Elimination_Endurance = colMeans(connect_compet_elim_endurance),
                                        Optimization_Competitiveness = colMeans(connect_compet_opt_competitiveness),
                                        Optimization_Endurance = colMeans(connect_compet_opt_endurance),
                                        NonSw_abundance=colMeans(con_abund_NonSw),
                                        Nswitch_elim_abundance=colMeans(con_abund_elimination),Nswitch_opt_abundance=colMeans(con_abund_optimization) ))
      
      
      require("stats")
      
      Competitiveness = dplyr::select(Compet_Abunda ,Non_Switching_Competitiveness,Elimination_Competitiveness,Optimization_Competitiveness,
                                      NonSw_abundance,Nswitch_elim_abundance,Nswitch_opt_abundance)
      
      minX = min(Competitiveness$Non_Switching_Competitiveness, Competitiveness$Elimination_Competitiveness,Competitiveness$Optimization_Competitiveness)
      minY = min(Competitiveness$NonSw_abundance, Competitiveness$Nswitch_elim_abundance,Competitiveness$Nswitch_opt_abundance)
      maxX = max(Competitiveness$Non_Switching_Competitiveness, Competitiveness$Elimination_Competitiveness,Competitiveness$Optimization_Competitiveness)
      maxY = max(Competitiveness$NonSw_abundance, Competitiveness$Nswitch_elim_abundance,Competitiveness$Nswitch_opt_abundance)
      
      
      Compe1 <- loess(NonSw_abundance ~ Non_Switching_Competitiveness,data=Competitiveness)
      Compe2<- loess(Nswitch_elim_abundance~Elimination_Competitiveness, data = Competitiveness)
      Compe3<-loess(Nswitch_opt_abundance~Optimization_Competitiveness,data=Competitiveness)
      
      ##############################################################################################
      # tiff(paste0("Competitivenessplot_",toString(c(Ntwork,S,SD)),".tiff"), width = 7.5, height = 7.5, units = 'in', res = 300)
      # #pdf(paste0("Competitivenessplot_",toString(c(Ntwork,S,SD)),".pdf"),width=7.5,height=7.5)
      # par(mar=c(5, 4, 4, 2) + 0.1,pty='s')
      # plot(NonSw_abundance ~Non_Switching_Competitiveness , data=Competitiveness,pch=1,cex=1.5, 
      #      xlim= c(minX,maxX), ylim=c(minY,maxY),xlab = "",ylab = "")
      # j <- order(Competitiveness$Non_Switching_Competitiveness)
      # lines(Competitiveness$Non_Switching_Competitiveness[j],Compe1$fitted[j],lwd=3, lty=1)
      # points(Competitiveness$Elimination_Competitiveness,Competitiveness$Nswitch_elim_abundance,cex=1.5,pch=2)
      # j <- order(Competitiveness$Elimination_Competitiveness)
      # lines(Competitiveness$Elimination_Competitiveness[j],Compe2$fitted[j],lwd=3,lty=3)
      # points(Competitiveness$Optimization_Competitiveness,Competitiveness$Nswitch_opt_abundance,cex=1.5,pch=5)
      # j <- order(Competitiveness$Optimization_Competitiveness)
      # lines(Competitiveness$Optimization_Competitiveness[j],Compe3$fitted[j],lwd=3,lty=5)
      # legend("topleft", legend = c("NonSwitching", "Elimination","Optimization"), pch = c(1, 2,5), lty = c(1, 3,5))
      # title(ylab="Biomass",xlab = "Competitiveness", line=2.2, cex.axis=2.5,cex.lab=1.5)
      # par(mar=c(5, 4, 4, 2) + 0.1)
      # dev.off()
      
      ################################################################################################
      Endurance = dplyr::select(Compet_Abunda ,Non_Switching_Endurance,Elimination_Endurance,Optimization_Endurance,
                                NonSw_abundance,Nswitch_elim_abundance,Nswitch_opt_abundance)
      
      minX = min(Endurance$Non_Switching_Endurance, Endurance$Elimination_Endurance,Endurance$Optimization_Endurance)
      minY = min(Endurance$NonSw_abundance, Endurance$Nswitch_elim_abundance,Endurance$Nswitch_opt_abundance)
      maxX = max(Endurance$Non_Switching_Endurance, Endurance$Elimination_Endurance,Endurance$Optimization_Endurance)
      maxY = max(Endurance$NonSw_abundance, Endurance$Nswitch_elim_abundance,Endurance$Nswitch_opt_abundance)
      
      
      endure1 <- loess(NonSw_abundance ~ Non_Switching_Endurance,data=Endurance)
      endure2<- loess(Nswitch_elim_abundance~Elimination_Endurance, data = Endurance)
      endure3<-loess(Nswitch_opt_abundance~Optimization_Endurance,data=Endurance)
      
      ################################################################################################
      # tiff(paste0("Enduranceplot_",toString(c(Ntwork,S,SD)),".tiff"), width = 7.5, height = 7.5, units = 'in', res = 300)
      # #pdf(paste0("Enduranceplot_",toString(c(Ntwork,S,SD)),".pdf"),width=7.5,height=7.5)
      # par(mar=c(5, 4, 4, 2) + 0.1,pty='s')
      # plot(NonSw_abundance ~Non_Switching_Endurance , data=Endurance,pch=1,cex=1.5, 
      #      xlim= c(minX,maxX), ylim=c(minY,maxY),xlab = "",ylab = "")
      # j <- order(Endurance$Non_Switching_Endurance)
      # lines(Endurance$Non_Switching_Endurance[j],endure1$fitted[j],lwd=3, lty=1)
      # points(Endurance$Elimination_Endurance,Endurance$Nswitch_elim_abundance,cex=1.5,pch=2)
      # j <- order(Endurance$Elimination_Endurance)
      # lines(Endurance$Elimination_Endurance[j],endure2$fitted[j],lwd=3,lty=3)
      # points(Endurance$Optimization_Endurance,Endurance$Nswitch_opt_abundance,cex=1.5,pch=5)
      # j <- order(Endurance$Optimization_Endurance)
      # lines(Endurance$Optimization_Endurance[j],endure3$fitted[j],lwd=3,lty=5)
      # legend("topleft", legend = c("NonSwitching", "Elimination","Optimization"), pch = c(1, 2,5), lty = c(1, 3,5))
      # title(ylab="Biomass",xlab = "Endurance", line=2.2, cex.axis=2.5,cex.lab=1.5)
      # par(mar=c(5, 4, 4, 2) + 0.1)
      # dev.off()
      #               
      ###################################################################################################
      
      
      #Dataframes saved uniquely for each competition strength values
      
      Competitiveness$SD=SD
      Compet_sigma_datalist1[[SD]]=Competitiveness
      
      Endurance$SD=SD
      Compet_sigma_datalist2[[SD]]=Endurance
      
      ############################################################################################              
      ##SPECIES DIVERSITY MEASURE USING RANK ABUNDANCE CURVE
      ############################################################################################   
      
      conn_per_abund_NOnSW=data.frame(t(con_abund_NonSw))
      conn_per_abund_SW=data.frame(t(con_abund_elimination))
      conn_per_abund_Opt=data.frame(t(con_abund_optimization))
      
      Species_rankConnect=list()
      
      # Rank abundance curves. And ranking the switch based on non switching ranking
      
      for (C in 1:length(connectance_sequence)) {
        
        rakingNonswitch=order(-conn_per_abund_NOnSW[,C])
        sortedNonswitch=conn_per_abund_NOnSW[,C][rakingNonswitch]
        rankingSwitch=order(-conn_per_abund_SW[,C])
        ranking_Opt_Switch=order(-conn_per_abund_Opt[,C])
        sortedElimSwitch=conn_per_abund_SW[,C][rankingSwitch]
        sortedOptSwitch=conn_per_abund_Opt[,C][ranking_Opt_Switch]
        
        Species_rank= data.frame(cbind(speciesrank=rank(-sortedNonswitch),sortedNonswitch,sortedElimSwitch,sortedOptSwitch))
        
        Species_rank$C=C
        Species_rankConnect[[C]]=Species_rank
        Species_rank_melted=melt(Species_rank,id.vars = c("speciesrank","C"))
        
      }
      
      Species_rankndata=bind_rows(Species_rankConnect)
      
      # Species rank data saved uniquely for each species competition strength
      Species_rankndata$SD=SD
      
      Species_rankdatalist[[SD]]=Species_rankndata
      
    } # End of standard deviation(SD) loop
    
    ######################################################################
    # SAVING DATA UNIQUELY FOR EACH SPECIES COMMUNITY SIZE
    # Combining the dataframes produced by each  community size(S) uniquely into one dataframe
    # Abundsigma_data     : species abundances  
    # Compet_sigma_data1  : competitveness along the column
    # Compet_sigma_data2  : competitveness along the row
    ######################################################################
    
    Species_rankdata=bind_rows(Species_rankdatalist)
    Species_rankdata$S=S
    Overall_Speciesranking[[S]]=Species_rankdata
    
    Abundsigma_data=bind_rows(Abundsigma_datalist)
    Abundsigma_data$S=S
    CombinedAbund_datalist[[S]]=Abundsigma_data
    SpeciesAbunda_data=bind_rows(CombinedAbund_datalist)
    
    Compet_sigma_data1=bind_rows(Compet_sigma_datalist1)
    Compet_sigma_data1$S=S
    CombinedCompet_datalist1[[S]]=Compet_sigma_data1
    SpeciesCompet_data1=bind_rows(CombinedCompet_datalist1)
    
    Compet_sigma_data2=bind_rows(Compet_sigma_datalist2)
    Compet_sigma_data2$S=S
    CombinedCompet_datalist2[[S]]=Compet_sigma_data2
    SpeciesCompet_data2=bind_rows(CombinedCompet_datalist2)
    
    
  } # End of Community size (S) loop.
  
  #########################################################################################
  # Combining the dataframes into one at end of community size(S) loop
  # SpeciesAbunda_data          : species abundances  
  # SpeciesCompet_data1         : competitveness along the column
  # SpeciesCompet_data2         : competitveness along the row
  #Overall_CombinedSpecies_rank : Species ranking 
  #########################################################################################
  
  SpeciesAbunda_data=bind_rows(CombinedAbund_datalist)
  
  SpeciesAbunda_data$Ntwork=Ntwork
  Network_Abundsigma_datalist[[Ntwork]]=SpeciesAbunda_data
  Network_SpeciesAbunda_data=bind_rows(Network_Abundsigma_datalist)
  
  SpeciesCompet_data1=bind_rows(CombinedCompet_datalist1)
  SpeciesCompet_data1$Ntwork = Ntwork
  Network_Compet_sigma_datalist1[[Ntwork]] = SpeciesCompet_data1
  Network_CombinedCompet_data = bind_rows(Network_Compet_sigma_datalist1)
  
  SpeciesCompet_data2=bind_rows(CombinedCompet_datalist2)
  SpeciesCompet_data2$Ntwork = Ntwork 
  Network_Compet_sigma_datalist2[[Ntwork]] = SpeciesCompet_data2
  Network_SpeciesCompet_Data2 = bind_rows(Network_Compet_sigma_datalist2)
  
  Overall_COmbinedSpecies_rank=bind_rows(Overall_Speciesranking)
  Overall_COmbinedSpecies_rank$Ntwork = Ntwork
  Network_Species_rankdatalist[[Ntwork]] = Overall_COmbinedSpecies_rank
  Network_Species_rankdata = bind_rows(Network_Species_rankdatalist)
  
  
  ################################################################################################
  
  SpeciesCommunity= SpeciesAbunda_data
  
  
  #################################################################################################
  #Combined Community plots for different values of sigma and S (number of species)
  # Community plots comparison under different values of  standard deviation and species population
  # S along the rows and SD( for community standard deviation) along the columns on the partitions
  ################################################################################################
  
  # Total Biomass versus connectance plots for each network
  
  CommunityProductivity= dplyr::select(SpeciesCommunity,Connectance,Abundance_Nonswitch,Abund_elim_switch,Abund_opt_switch,SD,S)
  
  CommunityProductivity_melted=melt(CommunityProductivity,id.vars = c("Connectance","SD","S"))
  ggplot(data=CommunityProductivity_melted,aes(x=Connectance,y=value, group=variable))+
    geom_point(aes(shape=variable),size=2.5)+
    facet_grid(S~SD, scales = "free",labeller = label_bquote(S == .(S), sigma ==.(SD)))+  ## partitioning S-rows and SD-columns
    ## line of best fit using LOESS
    geom_smooth(method = "loess",se=T,aes(linetype=variable, fill=variable),size=1.5,color="black")+
    scale_linetype_manual(values=c("solid", "dotted","dashed"),
                          breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                          labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_shape_manual(values = c(1,2,5),breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                       labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_fill_grey(start=0.3, end=0.9,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                    labels=c("Nonswitching","Elimination", "Optimization"))+
    theme_few()+
    
    labs(y="Total Biomass",x="Connectance",aspect.ratio = 1)+
    theme(legend.position =c(0.85,0.25),legend.title = element_blank(),aspect.ratio = 1,
          legend.key.width = unit(1.25,"cm"))+
    theme(axis.title=element_text(face="bold"),axis.text =element_text(face = "plain",colour = "black"))
  #Total Biomass plot saved uniquely for each network size
  ggsave(paste0("TotalCommunityBiomassNetwork_",toString(Ntwork),".pdf"), width = 7.5,height = 7.5)
  ggsave(paste0("TotalCommunityBiomassNetwork_",toString(Ntwork),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300) 
  
  ### Community stability versus connectance for each network
  
  CommunityLocal_Stab= cbind(SpeciesCommunity[1],SpeciesCommunity[8:10],SpeciesCommunity[11:12])
  CommunityLocal_Stab_melted=melt(CommunityLocal_Stab,id.vars = c("Connectance","SD","S"))
  ggplot(data=CommunityLocal_Stab_melted,aes(x=Connectance,y=value, group=variable))+
    geom_hline(aes(yintercept=0),size=0.5,colour='red',linetype="dashed")+  ## dashed line at y=0.
    geom_point(aes(shape=variable),size=2.5)+
    facet_grid(S~SD, scales = "free",labeller = label_bquote(S == .(S), sigma ==.(SD)))+
    ## line of best fit using LOESS 
    geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=1.5, color="black")+
    scale_linetype_manual(values = c("solid","dotted","dashed"),
                          breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                          labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_shape_manual(values = c(1,2,5),breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                       labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_fill_grey(start=0.3, end=0.9,breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                    labels=c("Nonswitching","Elimination", "Optimization"))+
    theme_few()+
    labs(y=expression("Re"~(lambda["Max"])),x="Connectance",aspect.ratio = 1)+
    theme(legend.position = c(0.14,0.3),legend.title = element_blank(),aspect.ratio = 1,
          legend.key.width = unit(1.5,"cm"))+
    theme(axis.title=element_text(face="plain"),axis.text =element_text(face = "plain",colour = "black"))
  ggsave(paste0("LocalStabilityNetwork_",toString(Ntwork),".pdf"),width = 7.5,height = 7.5)
  ggsave(paste0("LocalStabilityNetwork_",toString(Ntwork),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
  
  
  ### Nestedness versus connectance for each network
  
  CommunityNestedness= cbind(SpeciesCommunity[1],SpeciesCommunity[5:7],SpeciesCommunity[11:12])
  CommunityNestedness_melted=melt(CommunityNestedness,id.vars = c("Connectance","SD","S"))
  ggplot(data=CommunityNestedness_melted,aes(x=Connectance,y=value, group=variable))+
    geom_point(aes(shape=variable),size=2.5)+
    facet_grid(S~SD, scales = "free",labeller = label_bquote(S == .(S), sigma ==.(SD)))+
    geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable), size=1.5, color="black")+
    scale_linetype_manual(values = c("solid","dotted","dashed"),
                          breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                          labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_shape_manual(values = c(1,2,5),breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                       labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_fill_grey(start=0.3, end=0.9,breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                    labels=c("Nonswitching","Elimination", "Optimization"))+
    theme_few()+
    labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
    theme(legend.position = c(0.22, 0.72),legend.title = element_blank(),aspect.ratio = 1,
          legend.key.width = unit(1.5,"cm"))+
    theme(axis.title=element_text(face="plain"),
          axis.text =element_text(face = "plain",colour = "black"))
  ggsave(paste0("NestednessNetwork_",toString(Ntwork),".pdf"),width = 7.5,height = 7.5)
  ggsave(paste0("NestednessNetwork_",toString(Ntwork),".tiff"), width = 7.5, height = 7.5, units = 'in', dpi = 300)
  
  
  
  ## Species competitiveness plot per network
  SpeciesCompetitivenessData= SpeciesCompet_data1
  
  minX = min(SpeciesCompetitivenessData$Non_Switching_Competitiveness, SpeciesCompetitivenessData$Elimination_Competitiveness,SpeciesCompetitivenessData$Optimization_Competitiveness)
  minY = min(SpeciesCompetitivenessData$NonSw_abundance, SpeciesCompetitivenessData$Nswitch_elim_abundance,SpeciesCompetitivenessData$Nswitch_opt_abundance)
  maxX = max(SpeciesCompetitivenessData$Non_Switching_Competitiveness, SpeciesCompetitivenessData$Elimination_Competitiveness,SpeciesCompetitivenessData$Optimization_Competitiveness)
  maxY = max(SpeciesCompetitivenessData$NonSw_abundance, SpeciesCompetitivenessData$Nswitch_elim_abundance,SpeciesCompetitivenessData$Nswitch_opt_abundance)
  
  
  NetwCompe1 <- loess(NonSw_abundance ~ Non_Switching_Competitiveness,data=SpeciesCompetitivenessData)
  NetwCompe2<- loess(Nswitch_elim_abundance~Elimination_Competitiveness, data = SpeciesCompetitivenessData)
  NetwCompe3<-loess(Nswitch_opt_abundance~Optimization_Competitiveness,data=SpeciesCompetitivenessData)
  
  tiff(paste0("NetworkCompetitivenessplot_",toString(Ntwork),".tiff"), width = 7.5, height = 7.5, units = 'in', res = 300)
  #pdf(paste0("NetworkCompetitivenessplot_",toString(Ntwork),".pdf"),width=7.5,height=7.5)
  par(mar=c(5, 4, 4, 2) + 0.1,pty='s')
  plot(NonSw_abundance ~Non_Switching_Competitiveness , data=SpeciesCompetitivenessData,pch=1,cex=1.5, 
       xlim= c(minX,maxX), ylim=c(minY,maxY),xlab = "",ylab = "")
  j <- order(SpeciesCompetitivenessData$Non_Switching_Competitiveness)
  lines(SpeciesCompetitivenessData$Non_Switching_Competitiveness[j],NetwCompe1$fitted[j],lwd=3, lty=1)
  points(SpeciesCompetitivenessData$Elimination_Competitiveness,SpeciesCompetitivenessData$Nswitch_elim_abundance,cex=1.5,pch=2)
  j <- order(SpeciesCompetitivenessData$Elimination_Competitiveness)
  lines(SpeciesCompetitivenessData$Elimination_Competitiveness[j],NetwCompe2$fitted[j],lwd=3,lty=3)
  points(SpeciesCompetitivenessData$Optimization_Competitiveness,SpeciesCompetitivenessData$Nswitch_opt_abundance,cex=1.5,pch=5)
  j <- order(SpeciesCompetitivenessData$Optimization_Competitiveness)
  lines(SpeciesCompetitivenessData$Optimization_Competitiveness[j],NetwCompe3$fitted[j],lwd=3,lty=5)
  legend("topleft", legend = c("NonSwitching", "Elimination","Optimization"), pch = c(1, 2,5), lty = c(1, 3,5))
  title(ylab="Biomass",xlab = "Competitiveness", line=2.2, cex.axis=2.5,cex.lab=1.5)
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  
  
  # Endurance plot per network
  SpeciesEnduranceData= SpeciesCompet_data2
  
  
  minX = min(SpeciesEnduranceData$Non_Switching_Endurance, SpeciesEnduranceData$Elimination_Endurance,SpeciesEnduranceData$Optimization_Endurance)
  minY = min(SpeciesEnduranceData$NonSw_abundance, SpeciesEnduranceData$Nswitch_elim_abundance,SpeciesEnduranceData$Nswitch_opt_abundance)
  maxX = max(SpeciesEnduranceData$Non_Switching_Endurance, SpeciesEnduranceData$Elimination_Endurance,SpeciesEnduranceData$Optimization_Endurance)
  maxY = max(SpeciesEnduranceData$NonSw_abundance, SpeciesEnduranceData$Nswitch_elim_abundance,SpeciesEnduranceData$Nswitch_opt_abundance)
  
  
  Netwendure1 <- loess(NonSw_abundance ~ Non_Switching_Endurance,data=SpeciesEnduranceData)
  Netwendure2<- loess(Nswitch_elim_abundance~Elimination_Endurance, data = SpeciesEnduranceData)
  Netwendure3<-loess(Nswitch_opt_abundance~Optimization_Endurance,data=SpeciesEnduranceData)
  
  tiff(paste0("NetworkEnduranceplot_",toString(Ntwork),".tiff"), width = 7.5, height = 7.5, units = 'in', res = 300)
  #pdf(paste0("NetworkEnduranceplot_",toString(Ntwork),".pdf"),width=7.5,height=7.5)
  par(mar=c(5, 4, 4, 2) + 0.1,pty='s')
  plot(NonSw_abundance ~Non_Switching_Endurance , data=SpeciesEnduranceData,pch=1,cex=1.5, 
       xlim= c(minX,maxX), ylim=c(minY,maxY),xlab = "",ylab = "")
  j <- order(SpeciesEnduranceData$Non_Switching_Endurance)
  lines(SpeciesEnduranceData$Non_Switching_Endurance[j],Netwendure1$fitted[j],lwd=3, lty=1)
  points(SpeciesEnduranceData$Elimination_Endurance,SpeciesEnduranceData$Nswitch_elim_abundance,cex=1.5,pch=2)
  j <- order(SpeciesEnduranceData$Elimination_Endurance)
  lines(SpeciesEnduranceData$Elimination_Endurance[j],Netwendure2$fitted[j],lwd=3,lty=3)
  points(SpeciesEnduranceData$Optimization_Endurance,SpeciesEnduranceData$Nswitch_opt_abundance,cex=1.5,pch=5)
  j <- order(SpeciesEnduranceData$Optimization_Endurance)
  lines(SpeciesEnduranceData$Optimization_Endurance[j],Netwendure3$fitted[j],lwd=3,lty=5)
  legend("topleft", legend = c("NonSwitching", "Elimination","Optimization"), pch = c(1, 2,5), lty = c(1, 3,5))
  title(ylab="Biomass",xlab = "Endurance", line=2.2, cex.axis=2.5,cex.lab=1.5)
  par(mar=c(5, 4, 4, 2) + 0.1)
  dev.off()
  
  
  
} # End of Networks loop.

#########################################################################################
#SAVING NETWORK DATA 

Network_SpeciesAbunda_data=bind_rows(Network_Abundsigma_datalist)
Network_SpeciesCompet_data1=bind_rows(Network_Compet_sigma_datalist1)
Network_SpeciesCompet_Data2 = bind_rows(Network_Compet_sigma_datalist2)
Network_Species_rankdata = bind_rows(Network_Species_rankdatalist)

###############################################################

##Comparison switching to May's Criterion (Sigma* Sqrt((S-1)C))

##############################################################

NetworksCommunity=Network_SpeciesAbunda_data

#############################################################

NetworksCommunity$SD[NetworksCommunity$SD==1] =std_dvtn[1] 
NetworksCommunity$SD[NetworksCommunity$SD==2] =std_dvtn[2]
NetworksCommunity$SD[NetworksCommunity$SD==3] =std_dvtn[3]

############################################################ 
#COMPUTING  MAY'S CRITERIA (Sigma* Sqrt((S-1)C))
############################################################ 
NetworksCommunity$MaysRule <- with(NetworksCommunity, SD* sqrt(S*Connectance))


#Total Biomass Vs May's Criterion
for (Netwks in 1:CommNetworks){
  MayProductivity=dplyr::select(NetworksCommunity,Ntwork,MaysRule,Abundance_Nonswitch,Abund_elim_switch,Abund_opt_switch)
  MayProductivity_melted=melt(MayProductivity,id.vars = c("MaysRule","Ntwork"))
  ggplot(subset(MayProductivity_melted,Ntwork%in%Netwks),aes(x=MaysRule,y=value, group=variable))+
    geom_point(aes(shape=variable),size=5)+
    ylim(3,20) +
    # stat_smooth_func(geom = "text",method = "lm",hjust=-3,parse=T)+
    geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=2, color="black")+
    scale_linetype_manual(values=c("solid", "dotted","dashed")
                          ,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                          labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_shape_manual(values = c(1,2,5),breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                       labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_fill_grey(start=0.3, end=0.9,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                    labels=c("Nonswitching","Elimination", "Optimization"))+
    theme_few()+
    labs(y="Total Biomass",x=expression("Log"(~sigma*sqrt("(S-1)C"))),aspect.ratio = 1)+
    theme(legend.position = c(0.15,0.15),legend.title = element_blank(), aspect.ratio = 1,
          legend.text = element_text(size=14),legend.key.width = unit(1.5, "cm"))+
    theme(axis.title=element_text(face="plain",size = 20.5),
          axis.text =element_text(face = "plain",colour = "black",size=14))
  ggsave(paste0("May'sRuleProductivity_",toString(Netwks),".pdf"),width = 7.5,height = 7.5)
  ggsave(paste0("May'sRuleProductivity_",toString(Netwks),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
  
  
  # # Stability vs May's Criterion
  
  MayLocal_Stab= dplyr::select(NetworksCommunity,Stability_NonSW,Stability_Elim_sw,Stability_Opt_sw,Ntwork,MaysRule)
  MayLocal_Stab_melted=melt(MayLocal_Stab,id.vars = c("MaysRule","Ntwork"))
  ggplot(subset(MayLocal_Stab_melted,Ntwork%in%Netwks),aes(x=MaysRule,y=value, group=variable))+
    geom_hline(aes(yintercept=0),size=0.9,colour='red',linetype="dashed")+
    geom_point(aes(shape=variable),size=5)+
    #scale_shape_manual(values=c(1,2,5))+
    ylim(-0.5, 0.1)+
    geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=2,color="black")+
    scale_linetype_manual(values=c("solid", "dotted","dashed"),
                          breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                          labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_shape_manual(values= c(1,2,5),breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                       labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_fill_grey(start=0.3, end=0.9,breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                    labels=c("Nonswitching","Elimination", "Optimization"))+
    theme_few()+
    labs(y=expression("Re"~(lambda["Max"])),x=expression("Log"(~sigma*sqrt("(S-1)C"))),aspect.ratio = 1)+
    theme(legend.position = c(0.85,0.2),legend.title = element_blank(),aspect.ratio = 1,
          legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
    theme(axis.title=element_text(face="plain",size = 20.5),
          axis.text =element_text(face = "plain",colour = "black",size = 14))
  ggsave(paste0("May'sRulestability_",toString(Netwks),".pdf"),width = 7.5,height = 7.5)
  ggsave(paste0("May'sRulestability_",toString(Netwks),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
  
  # Nestedness vs May's criterion
  
  MayNestedness= dplyr::select(NetworksCommunity,Nested_Nonswitch,Nest_elim_switch,Nest_opt_switch,Ntwork,MaysRule)
  MayNestedness_melted=melt(MayNestedness,id.vars = c("MaysRule","Ntwork"))
  ggplot(subset(MayNestedness_melted,Ntwork%in%Netwks),aes(x=MaysRule,y=value, group=variable))+
    geom_point(aes(shape=variable),size=5)+
    #scale_shape_manual(values = c(1,2,5))+
    geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=2,color="black")+
    scale_linetype_manual(values=c("solid", "dotted","dashed"),
                          breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                          labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_shape_manual(values = c(1,2,5),breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                       labels=c("Nonswitching","Elimination", "Optimization"))+
    scale_fill_grey(start=0.3,end=0.9,breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                    labels=c("Nonswitching","Elimination", "Optimization"))+
    theme_few()+
    labs(y="Nestedness",x=expression("Log"~(sigma*sqrt("(S-1)C"))),aspect.ratio = 1)+
    theme(legend.position = c(0.85,0.15),legend.title = element_blank(), aspect.ratio = 1,
          legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
    theme(axis.title=element_text(face="plain", size = 20.5),
          axis.text =element_text(face = "plain",colour = "black",size = 14))
  ggsave(paste0("May'sRule_vs_Nestedness_",toString(Netwks),".pdf"),width = 7.5,height = 7.5)
  ggsave(paste0("May'sRule_vs_Nestedness_",toString(Netwks),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)
  
}
##################################################################################

##COMBINED TOTAL BIOMASS, NESTEDNESS & STABILITY PLOTS FOR ALL NETWORKS

##############################################################################

Combined_NetworkData= Network_SpeciesAbunda_data

##Network Biomass

Network_Biomass = dplyr::select(Combined_NetworkData,Connectance,Abundance_Nonswitch,Abund_elim_switch,Abund_opt_switch,SD,S,Ntwork)
Network_Biomass_melted=melt(Network_Biomass,id.vars = c("Connectance","SD","S","Ntwork"))
ggplot(subset(Network_Biomass_melted,SD%in%2 ),aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable),size=5)+
  #scale_shape_manual(values = c(1,2,5))+
  geom_smooth(method = "loess",span=0.85,se=T,aes(linetype=variable,fill=variable),size=2, color="black")+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual( values = c(1,2,5),
                      breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                      labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.9,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Total Biomass",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.8,0.85), legend.title=element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))
ggsave(paste0("NetworkTotal Biomass.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("NetworkTotal Biomass.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)

#Nestedness Plot

Network_Nestedness=dplyr::select(Combined_NetworkData,Connectance,Nested_Nonswitch,Nest_elim_switch,Nest_opt_switch,SD,S,Ntwork)
Network_Nestedness_melted=melt(Network_Nestedness,id.vars = c("Connectance","SD","S","Ntwork"))
ggplot(subset(Network_Nestedness_melted, SD%in%1),aes(x=Connectance,y=value, group=variable))+
  geom_point(aes(shape=variable),size=5)+
  #scale_shape_manual(values = c(1,2,5))+
  geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=2,color="black")+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual(values = c(1,2,5),
                     breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                     labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.9,breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.2,0.85),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))

# Saving the plots to Folder (pdf & Tiff files )

ggsave(paste0("NetworkNestedness.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("NetworkNestedness.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)

# Community Resilience .

Network_LocalStability = dplyr::select(Combined_NetworkData,Connectance,Stability_NonSW,Stability_Elim_sw,Stability_Opt_sw,SD,S,Ntwork)
Network_LocalStability_melted=melt(Network_LocalStability,id.vars = c("Connectance","SD","S","Ntwork"))
ggplot(subset(Network_LocalStability_melted,SD%in%1),aes(x=Connectance,y=value, group=variable))+
  geom_hline(aes(yintercept=0),size=1,colour='black',linetype="dotted")+
  geom_point(aes(shape=variable),size=5)+
  #scale_shape_manual(values = c(1,2,5))+
  geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=2, color="black")+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual( values = c(1,2,5),
                      breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                      labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.9,breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y=expression("Re"~(lambda["Max"])),x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.8,0.15),legend.title = element_blank(), aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size = 20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))

# Saving the plots to Folder (pdf & Tiff files )

ggsave(paste0("NetworkLocalStability.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("NetworkLocalStability.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)

##################################################################################################

##SPECIES RANKING USING RANK ABUNDANCE CURVES

##################################################################################################


##CHOOSING LEAST, MEDIUM, AND HIGH CONNECTANCE

Overall_COmbinedSpecies_rank$C[Overall_COmbinedSpecies_rank$C==which(connectance_sequence==min(connectance_sequence),arr.ind = T)]=connectance_sequence[which(connectance_sequence==min(connectance_sequence),arr.ind = T)] = connectance_sequence[1]
Overall_COmbinedSpecies_rank$C[Overall_COmbinedSpecies_rank$C==which(connectance_sequence==median(connectance_sequence),arr.ind = T)]=connectance_sequence[which(connectance_sequence==median(connectance_sequence),arr.ind = T)] = connectance_sequence[9]
Overall_COmbinedSpecies_rank$C[Overall_COmbinedSpecies_rank$C==which(connectance_sequence==max(connectance_sequence),arr.ind = T)]=connectance_sequence[which(connectance_sequence==max(connectance_sequence),arr.ind = T)]=connectance_sequence[17]


Overall_COmbinedSpecies_rank$SD[Overall_COmbinedSpecies_rank$SD==1]=std_dvtn[1]
Overall_COmbinedSpecies_rank$SD[Overall_COmbinedSpecies_rank$SD==2]=std_dvtn[2]
Overall_COmbinedSpecies_rank$SD[Overall_COmbinedSpecies_rank$SD==3]=std_dvtn[3]

Overall_COmbinedSpecies_rankmelt=melt(Overall_COmbinedSpecies_rank,id.vars = c("speciesrank","C","SD","S","Ntwork"))

# Species ranking in 20 species community

ggplot(subset(Overall_COmbinedSpecies_rankmelt,C%in%c(0.1,0.5,0.9) & S%in%20),aes(x=speciesrank,y=value, group=variable))+
  geom_point(aes(shape=variable),size=2.5)+
  facet_grid(SD~S+C, scales = "free",labeller = label_bquote(sigma == .(SD), cols = list(S ==.(S),C==.(C))))+
  geom_line(aes(linetype=variable),size=0.75)+
  scale_linetype_manual(values=c("solid","dotted","dashed"),
                        breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                        labels=c("Nonswitching","Elimination","Optimization"))+
  scale_shape_manual(values=c(1,2,5),
                     breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                     labels=c("Nonswitching","Elimination","Optimization"))+
  theme_few()+
  labs(y="Abundance",x="Species rank")+
  theme(legend.position = "bottom",legend.title=element_blank(),aspect.ratio = 1,
        legend.key.width = unit(1.5, "cm"))+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))
ggsave("Species ranking25.pdf",width = 10,height = 10)
ggsave(paste0("Species ranking25.tiff"), width = 8, height = 8, units = 'in', dpi = 300)


# Species ranking in 50 species community

ggplot(subset(Overall_COmbinedSpecies_rankmelt,C%in%c(0.1,0.5,0.9) & S%in%50),aes(x=speciesrank,y=value, group=variable))+
  geom_point(aes(shape=variable),size=2.5)+
  facet_grid(SD~S+C, scales = "free", labeller = label_bquote(sigma == .(SD), cols=list(S ==.(S),C==.(C))))+
  geom_line(aes(linetype=variable),size=0.75)+
  scale_shape_manual(values=c(1,2,5),breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                     labels=c("Nonswitching","Elimination","Optimization"))+
  scale_linetype_manual(values=c("solid","dotted","dashed"),
                        breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                        labels=c("Nonswitching","Elimination","Optimization"))+
  theme_few()+
  labs(y="Abundance",x="Species rank")+
  theme(legend.position = "bottom",legend.title = element_blank(),aspect.ratio = 1,
        legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain"),
        axis.text =element_text(face = "plain",colour = "black"))
ggsave("Species ranking50.pdf",width = 10,height = 10)
ggsave(paste0("Species ranking50.tiff"), width = 8, height = 8, units = 'in', dpi = 300)
#######################################################################



##SPECIES RANK CURVE

Overall_COmbinedSpecies_rankmelt=melt(Overall_COmbinedSpecies_rank,id.vars = c("speciesrank","C","SD","S","Ntwork"))

###############################################           
##Rank abundance plot         
###################################################################################################       

ggplot(subset(Overall_COmbinedSpecies_rankmelt, C%in% 5.0 & SD %in% 0.1 & S%in%20),aes(x=speciesrank,y=value, group=variable))+
  geom_point(aes(shape=variable),size=5)+
  #scale_shape_manual(values = c(1,2,5))+
  #facet_grid(SD~S+C, scales = "free", labeller = label_bquote(sigma == .(SD), cols=list(S ==.(S),C==.(C))))+
  geom_line(aes(linetype=variable),size=2)+
  scale_shape_manual(values=c(1,2,5),breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                     labels=c("Nonswitching","Elimination","Optimization"))+
  scale_linetype_manual(values=c("solid","dotted","dashed"),
                        breaks=c("sortedNonswitch","sortedElimSwitch","sortedOptSwitch"),
                        labels=c("Nonswitching","Elimination","Optimization"))+
  theme_few()+
  labs(y="Abundance",x="Species rank")+
  theme(legend.position =c(0.2,0.1),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))
#ggsave("Species ranking50.pdf",width = 10,height = 10)
#(paste0("Species ranking50.tiff"), width = 8, height = 8, units = 'in', dpi = 300)
ggsave(paste0("ranking_",toString(Ntwork),".pdf"),width = 7.5,height = 7.5)
ggsave(paste0("ranking_",toString(Ntwork),".tiff"),width = 7.5,height = 7.5,units = 'in', dpi=300)



###COMPARISON OF MAY'S & NEW STABILITY CRITERION

SpeciesCommunity_Comparison = NetworksCommunity
SpeciesCommunity_Comparison= SpeciesCommunity_Comparison[SpeciesCommunity_Comparison$Connectance=="0.15",]
SpeciesCommunity_Comparison$SigmaMDEr<- with(SpeciesCommunity_Comparison, SD*sqrt(pi*(0.5- (Connectance/3))))

SpeciesCommunity_Comparison$FullyConnected=with(SpeciesCommunity_Comparison,  Connectance*SD*sqrt(S-1))
SpeciesCommunity_Comparison$MaysRule <- with(SpeciesCommunity_Comparison, SD* sqrt((S-1)*Connectance))
SpeciesCommunity_Comparison$Switching <-with(SpeciesCommunity_Comparison,  (1/Connectance)*SigmaMDEr*sqrt((S-1)*Connectance))


MayStability_OurRule= dplyr::select(SpeciesCommunity_Comparison,Connectance,SD,S,Ntwork,MaysRule,FullyConnected,Switching,Stability_NonSW,Stability_Elim_sw,Stability_Opt_sw)

## averaging the all networks' values per may's rule group.
MayStability_OurRulePerNetwork = MayStability_OurRule %>% group_by(MaysRule) %>% summarise_all(funs(mean(.,na.rm=T)))

MayStability_OurRule_Melt= melt(MayStability_OurRulePerNetwork, id.vars =c("Connectance","SD","S","Ntwork","Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"))
ggplot(MayStability_OurRule_Melt ,aes(x=value,y=Stability_NonSW, group=variable))+
  geom_hline(aes(yintercept=0),size=1.5,colour='red',linetype="dashed")+
  geom_vline(aes(xintercept =1), size= 1.5, linetype="dotted",color=2)+
  geom_vline(aes(xintercept=sqrt(0.1)),size=1.5,linetype="dotted",color=3)+
  geom_vline(aes(xintercept=(1/0.15)*sqrt(pi*(0.5-1/45))),size=1.5,linetype="dotted",color=4)+
  geom_point(aes(shape=variable),size=3.5)+
  scale_shape_manual(values = c(1,2,5))+
  geom_smooth(method = "loess",se=F,aes(linetype=variable,color=variable),size=2)+
  theme_few()+
  labs(y=expression("Re"~(lambda["Max"])),x=expression(~sigma*sqrt("(S-1)")),aspect.ratio = 1)+
  theme(legend.position = c(0.55,0.15),legend.title = element_blank(), aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20),
        axis.text =element_text(face = "plain",colour = "black",size=15))
ggsave(paste("NewCriteria.pdf"), width = 7.5,height = 7.5)
ggsave(paste0("NewCriteria.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


## Maximum real eigenvalues vs the MAY stability criterion
MayLocal_Stab= dplyr::select(NetworksCommunity,Stability_NonSW,Stability_Elim_sw,Stability_Opt_sw,Ntwork,MaysRule)
Summarised_MayLocal_Stab <- MayLocal_Stab %>% group_by(MaysRule)%>% summarise_all(funs(mean))
MayLocal_Stab_melted=melt(Summarised_MayLocal_Stab,id.vars = c("MaysRule","Ntwork"))
ggplot(MayLocal_Stab_melted,aes(x=MaysRule,y=value, group=variable))+
  geom_hline(aes(yintercept=0),size=0.9,colour='red',linetype="dashed")+
  geom_point(aes(shape=variable),size=5)+
  #scale_shape_manual(values=c(1,2,5))+
  ylim(-0.5, 0.1)+
  geom_smooth(method = "loess",se=T,aes(linetype=variable),size=2,color="black")+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual(values= c(1,2,5),breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                     labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y=expression("Re"~(lambda["Max"])),x=expression("Log"(~sigma*sqrt("(S-1)C"))),aspect.ratio = 1)+
  theme(legend.position = c(0.85,0.2),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size = 20.5),
        axis.text =element_text(face = "plain",colour = "black",size = 14))
ggsave(paste0("May'sRulestability_",toString(Netwks),".pdf"),width = 7.5,height = 7.5)
ggsave(paste0("May'sRulestability_",toString(Netwks),".tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


##COMBINED NETWORK AVERAGES FOR STABILITY, TOTAL BIOMASS & NESTEDNESS
Network_Dataset = Combined_NetworkData[1:10]

pd= position_dodge(width = 0.2)
## Averaged STABILITY for all the networks for each connectance value

ALL_Network_Stability =dplyr::select(Network_Dataset,Connectance,Stability_NonSW,Stability_Elim_sw,Stability_Opt_sw)
ALL_Network_Stability_Melt = melt(ALL_Network_Stability,id.vars ="Connectance")
Summarised_Networks_Stability= ALL_Network_Stability_Melt %>% group_by(Connectance, variable) %>%summarise(mean=mean(value), sd=sd(value))
ggplot(Summarised_Networks_Stability, aes(x=Connectance, y=mean, group=variable))+
  geom_hline(aes(yintercept=0),size=0.9,colour='red',linetype="dashed")+
  geom_point(position=pd,aes(shape=variable),size=5)+
  ylim(-0.5, 0.1)+
  geom_smooth(method = "loess",se=F,aes(linetype=variable),size=2,color="black",position=pd)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.05,position=pd)+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual(values= c(1,2,5),breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                     labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y=expression("Re"~(lambda["Max"])),x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.8,0.2),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size = 20.5),
        axis.text =element_text(face = "plain",colour = "black",size = 14))

ggsave(paste0("AveragedNetworkStability.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("AveragedNetworkStability.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)

##plot 2 without errorbars

ggplot(Summarised_Networks_Stability, aes(x=Connectance, y=mean, group=variable))+
  geom_hline(aes(yintercept=0),size=0.9,colour='red',linetype="dashed")+
  geom_point(aes(shape=variable),size=5)+
  ylim(-0.5, 0.1)+
  geom_smooth(data=ALL_Network_Stability_Melt,method = "loess",se=T,aes(linetype=variable),size=2,color="black")+
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.05,position=pd)+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual(values= c(1,2,5),breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                     labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3,end=0.8,breaks=c("Stability_NonSW","Stability_Elim_sw","Stability_Opt_sw"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y=expression("Re"~(lambda["Max"])),x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.8,0.2),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size = 20.5),
        axis.text =element_text(face = "plain",colour = "black",size = 14))

ggsave(paste0("AveragedNetworkStabilityGrey.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("AveragedNetworkStabilityGrey.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


## Averaged TOTAL BIOMASS plots for all the community networks
Selected_Network_Biomass =dplyr::select(Network_Dataset,Connectance,Abundance_Nonswitch,Abund_elim_switch,Abund_opt_switch)
Selected_Network_Biomass_melt = melt(Selected_Network_Biomass, id.vars = "Connectance")
Summarised_Networks_Biomass= Selected_Network_Biomass_melt %>% group_by(Connectance, variable) %>%summarise(mean=mean(value), sd=sd(value))

##Plot 1 with error bars
ggplot(Summarised_Networks_Biomass,aes(x=Connectance,y=mean, group=variable))+
  geom_point(position=pd,aes(shape=variable),size=5)+
  geom_smooth(method = "loess",span=0.85, se=F, aes(linetype=variable,fill=variable),size=2, color="black",position = pd)+
  geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.05,position=pd)+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual( values = c(1,2,5),
                      breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                      labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.8,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Total Biomass",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.8,0.85), legend.title=element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))

ggsave(paste0("AveragedNetworkBiomass.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("AveragedNetworkBiomass.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)

##Plot 2 with grey smoothing 
ggplot(Summarised_Networks_Biomass,aes(x=Connectance,y=mean, group=variable))+
  geom_point(aes(shape=variable),size=5)+
  geom_smooth(method = "loess",span=0.85, se=T, aes(linetype=variable,fill=variable),size=2, color="black")+
  #geom_errorbar(aes(ymin=mean-sd,ymax=mean+sd),width=0.05,position=pd)+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual( values = c(1,2,5),
                      breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                      labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.8,breaks=c("Abundance_Nonswitch","Abund_elim_switch","Abund_opt_switch"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Total Biomass",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.8,0.85), legend.title=element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))

ggsave(paste0("AveragedNetworkBiomassGRey.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("AveragedNetworkBiomassGrey.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)



#Averaged NESTEDNESS plot for all community networks

##Plot 1 with error bar
Selected_Network_Nestedness = dplyr::select(Network_Dataset,Connectance,Nested_Nonswitch,Nest_elim_switch,Nest_opt_switch)
Selected_Network_Nestedness_melt = melt(Selected_Network_Nestedness, id.vars = "Connectance")
Summarised_Networks_Nestedness = Selected_Network_Nestedness_melt %>% group_by(Connectance,variable) %>% summarise(mean=mean(value),sd=sd(value))
ggplot(Summarised_Networks_Nestedness,aes(x=Connectance,y=mean, group=variable))+
  geom_point(position=pd,aes(shape=variable),size=5)+
  geom_smooth(method = "loess",se=F,aes(linetype=variable,fill=variable),size=2,color="black",position=pd)+
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.05, position = pd)+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual(values = c(1,2,5),
                     breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                     labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.9,breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.2,0.85),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))

ggsave(paste0("AveragedNetworkNestedness.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("AveragedNetworkNestedness.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)


##plot 2 without errorbars
ggplot(Summarised_Networks_Nestedness,aes(x=Connectance,y=mean, group=variable))+
  geom_point(aes(shape=variable),size=5)+
  geom_smooth(method = "loess",se=T,aes(linetype=variable,fill=variable),size=2,color="black")+
  #geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd),width=0.05, position = pd)+
  scale_linetype_manual(values=c("solid", "dotted","dashed"),
                        breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                        labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_shape_manual(values = c(1,2,5),
                     breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                     labels=c("Nonswitching","Elimination", "Optimization"))+
  scale_fill_grey(start=0.3, end=0.7,breaks=c("Nested_Nonswitch","Nest_elim_switch","Nest_opt_switch"),
                  labels=c("Nonswitching","Elimination", "Optimization"))+
  theme_few()+
  labs(y="Nestedness",x="Connectance",aspect.ratio = 1)+
  theme(legend.position = c(0.2,0.85),legend.title = element_blank(),aspect.ratio = 1,
        legend.text = element_text(size=14),legend.key.width = unit(1.5,"cm"))+
  theme(axis.title=element_text(face="plain",size=20.5),
        axis.text =element_text(face = "plain",colour = "black",size=14))

ggsave(paste0("AveragedNetworkNestednessGrey.pdf"),width = 7.5,height = 7.5)
ggsave(paste0("AveragedNetworkNestednessGrey.tiff"), width = 7.5,height = 7.5, units = 'in', dpi = 300)



###########################################################################################
#SAVING ALL DATA SETS TO A FOLDER
###########################################################################################

save(SpeciesCommunity,SpeciesAbunda_data,BIOM_NEST_STAB_CONNECTANCE,Network_SpeciesAbunda_data,NetworksCommunity,
     Network_SpeciesCompet_data1,Network_SpeciesCompet_Data2,Network_Species_rankdata, file="Community.RData")


# ##################################
####################
##
##
##/** END OF SINGLE COMMUNITY SIMULATION 

#############################################################################################################
##########################################################################################################
#load("Community.RData")
