knitr::opts_chunk$set(echo = TRUE)
library(mrgsolve)
library(knitr)
library(PKPDmisc)
library(PKNCA)
library(knitr)
library(parallel)
library(kableExtra)
library(future.apply)
library(tidyverse)
library(officedown)
library(gridExtra)
library(tidyr)
library(dplyr)
library(readr)



#MFE Work
models <- "/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models"

setwd(models)

CLhep<- rnorm(10000,10950,5600)
CLhep<-CLhep[CLhep>=200 &  CLhep<=12000]

CLrenal<- rnorm(10000,121,20) 
CLrenal<-CLrenal[CLrenal<140]
Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.8 & Kpf<1.25]

out_0<-data.frame(matrix(1:21, ncol=21))

names(out_0)<-c("ID","time","amt","x","IRa","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva","CortTis","CortTis2")

c <-detectCores(all.tests=FALSE, logical =TRUE) #count cores on the CPU
cores <- c-1                                    #use all cores minus 1


options(future.fork.endable = TRUE)  #### from package future.apply, using forked approach
plan(multiprocess, workers=cores)    #### assign num

#####################################################
#####################################################
#####################################################
#####################################################
#####################################################
#####################################################

repeat{
  for(i in 0:500){ 
    it<-i+1
    
    names(it)<-"x"
    idataCLh<- expand.idata(CLhep = (sample(CLhep,33,replace=FALSE)))
    idataCLr<- expand.idata(Clrenal=(sample(CLrenal,33, replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,33,replace=FALSE)))
    itx<- merge(it,dem)
    idata1<- merge(itx,idataCLh)
    idata2<- merge(idata1,idataCLr)
    idata3<- merge(idata2,idataKpf)
    
    
    
    sp <- split(idata3,idata3$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod_pt<-loadso(mread("phenytoin_091121b"))
      as.data.frame(
        mod_pt %>%
          data_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva, CortTis,CortTis2) %>%
          carry_out(amt,x,IRa)%>%
          mrgsim(maxsteps=100000, tgrid=t_w))})%>%
      bind_rows()
    
    
    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    
    if (it==500) {
      return(out_0); break}
    
  }}


out_0<-out_0[-c(1),]


#ROI1 Cortex
#ROI2 Basal Ganglia
#ROI3 Thalamus
#ROI4 ROB
 

out_tiss<- out_0%>% group_by(time)%>%
  summarise(
    CxEc10=quantile(Cec1, probs=0.10),
    CxEc50=quantile(Cec1, probs=0.50),
    CxEc90=quantile(Cec1, probs=0.90),
    CxIc10=quantile(Icf1, probs=0.10),
    CxIc50=quantile(Icf1, probs=0.50),
    CxIc90=quantile(Icf1, probs=0.90),
    
    BGEc10=quantile(Cec2, probs=0.10),
    BGEc50=quantile(Cec2, probs=0.50),
    BGEc90=quantile(Cec2, probs=0.90),
    BGIc10=quantile(Icf2, probs=0.10),
    BGIc50=quantile(Icf2, probs=0.50),
    BGIc90=quantile(Icf2, probs=0.90),
    
    THEc10=quantile(Cec3, probs=0.10),
    THEc50=quantile(Cec3, probs=0.50),
    THEc90=quantile(Cec3, probs=0.90),
    THIc10=quantile(Icf3, probs=0.10),
    THIc50=quantile(Icf3, probs=0.50),
    THIc90=quantile(Icf3, probs=0.90),
    
    
    ROBEc10=quantile(Cec4, probs=0.10),
    ROBEc50=quantile(Cec4, probs=0.50),
    ROBEc90=quantile(Cec4, probs=0.90),
    ROBIc10=quantile(Icf4, probs=0.10),
    ROBIc50=quantile(Icf4, probs=0.50),
    ROBIc90=quantile(Icf4, probs=0.90),
    
  ) 


tiss_NCA<- out_0 %>% group_by(ID,x)%>%
  summarise(
    #AUC0-t 
    Cec1AUC=auc_partial(time,Cec1),
    Cec2AUC=auc_partial(time,Cec2),
    Cec3AUC=auc_partial(time,Cec3),
    Cec4AUC=auc_partial(time,Cec4),
    Icf1AUC=auc_partial(time,Icf1),
    Icf2AUC=auc_partial(time,Icf2),
    Icf3AUC=auc_partial(time,Icf3),
    Icf4AUC=auc_partial(time,Icf4),
    
    #Cavg 0-24
    Cec1CAV=Cec1AUC/24,
    Cec2CAV=Cec2AUC/24,
    Cec3CAV=Cec3AUC/24,
    Cec4CAV=Cec4AUC/24,
    Icf1CAV=Icf1AUC/24,
    Icf2CAV=Icf2AUC/24,
    Icf3CAV=Icf3AUC/24,
    Icf4CAV=Icf4AUC/24,
    
    #Cmax
    Cec1Cmax=max(Cec1),
    Cec2Cmax=max(Cec2),
    Cec3Cmax=max(Cec3),
    Cec4Cmax=max(Cec4),
    Icf1Cmax=max(Icf1),
    Icf2Cmax=max(Icf2),
    Icf3Cmax=max(Icf3),
    Icf4Cmax=max(Icf4),
    
    #tmax
    Cec1Tmax=time[which(Cec1==Cec1Cmax)],
    
  )

NCA_sum<- tiss_NCA %>% group_by(x)%>%
  summarise(
    
    #CAV
    Cec150=quantile(Cec1CAV,probs=0.50),
    Cec110=quantile(Cec1CAV,probs=0.10),
    Cec190=quantile(Cec1CAV,probs=0.90),
    Icf150=quantile(Icf1CAV,probs=0.50),
    Icf110=quantile(Icf1CAV,probs=0.10),
    Icf190=quantile(Icf1CAV,probs=0.90),
    
    Cec250=quantile(Cec2CAV,probs=0.50),
    Cec210=quantile(Cec2CAV,probs=0.10),
    Cec290=quantile(Cec2CAV,probs=0.90),
    Icf250=quantile(Icf2CAV,probs=0.50),
    Icf210=quantile(Icf2CAV,probs=0.10),
    Icf290=quantile(Icf2CAV,probs=0.90),
    
    Cec350=quantile(Cec3CAV,probs=0.50),
    Cec310=quantile(Cec3CAV,probs=0.10),
    Cec390=quantile(Cec3CAV,probs=0.90),
    Icf350=quantile(Icf3CAV,probs=0.50),
    Icf310=quantile(Icf3CAV,probs=0.10),
    Icf390=quantile(Icf3CAV,probs=0.90),
    
    Cec450=quantile(Cec4CAV,probs=0.50),
    Cec410=quantile(Cec4CAV,probs=0.10),
    Cec490=quantile(Cec4CAV,probs=0.90),
    Icf450=quantile(Icf4CAV,probs=0.50),
    Icf410=quantile(Icf4CAV,probs=0.10),
    Icf490=quantile(Icf4CAV,probs=0.90),
    
    
    
    
    
    
    
    
  )


NCA_ratio<- NCA_sum %>% 
  summarize(
        Cec190a=mean(Cec190),
        Cec110a=mean(Cec110),
        
        Cec290a=mean(Cec290),
        Cec210a=mean(Cec210),
        
        Cec390a=mean(Cec390),
        Cec310a=mean(Cec310),
        
        Cec490a=mean(Cec490),
        Cec410a=mean(Cec410),
        
        Icf110a=mean(Icf110),
        Icf190a=mean(Icf190),
        
        Icf210a=mean(Icf210),
        Icf290a=mean(Icf290),
        
        Icf310a=mean(Icf310),
        Icf390a=mean(Icf390),
        
        Icf410a=mean(Icf410),
        Icf490a=mean(Icf490)
        
        
  )%>%
  select("Cec190a","Cec110a",
         "Cec290a","Cec210a",
         "Cec390a","Cec310a",
         "Cec490a","Cec410a",
         
         "Icf110a","Icf190a",
         "Icf210a","Icf290a",
         "Icf310a","Icf390a",
         "Icf410a","Icf490a"
         )


NCA_ratio2<- NCA_ratio %>%
  
  #ROI1 Cortex
  #ROI2 Basal Ganglia
  #ROI3 Thalamus
  #ROI4 ROB
#CXBG,CXTH,CXROB
#BGTH,BGROB
#THROB
  
  mutate(
    rCx_BgEc10=Cec110a/Cec210a,
    rCx_BgEc90=Cec190a/Cec290a,
    
    rCx_THEc10=Cec110a/Cec310a,
    rCx_THEc90=Cec190a/Cec390a,
    
    rCx_ROBEc10=Cec110a/Cec410a,
    rCx_ROBEc90=Cec190a/Cec490a,
    
    rBG_THEc10=Cec210a/Cec310a,
    rBG_THEc90=Cec290a/Cec390a,
    
    rBG_ROBEc10=Cec210a/Cec410a,
    rBG_ROBEc90=Cec290a/Cec490a,
    
    rTH_ROBEc10=Cec310a/Cec410a,
    rTH_ROBEc90=Cec390a/Cec490a,
    
    
    
    rCx_BgIc10=Icf110a/Icf210a,
    rCx_BgIc90=Icf190a/Icf290a,
    
    rCx_THIc10=Icf110a/Icf310a,
    rCx_THIc90=Icf190a/Icf390a,
    
    rCx_ROBIc10=Icf110a/Icf410a,
    rCx_ROBIc90=Icf190a/Icf490a,
    
    rBG_THIc10=Icf210a/Icf310a,
    rBG_THIc90=Icf290a/Icf390a,
    
    rBG_ROBIc10=Icf210a/Icf410a,
    rBG_ROBIc90=Icf290a/Icf490a,
    
    rTH_ROBIc10=Icf310a/Icf410a,
    rTH_ROBIc90=Icf390a/Icf490a,
    
#EC/ICratios  
#ROI1 Cortex
#ROI2 Basal Ganglia
#ROI3 Thalamus
#ROI4 ROB
    rCx10=Icf110a/Cec110a,
    rCx90=Icf190a/Cec190a,
    
    rBG10=Icf210a/Cec210a,
    rBG90=Icf290a/Cec290a,

    rTH10=Icf310a/Cec310a,
    rTH90=Icf390a/Cec390a,
    
    rROB10=Icf410a/Cec410a,
    rROB90=Icf490a/Cec490a

  )%>%
  select(starts_with("r"))

NCA_tables<-NCA_ratio2 %>% 
    pivot_longer(cols=starts_with("r"),
                 names_to="RegionConcPtile",
                 values_to="Value")


NCA_tables<-NCA_tables %>%
    mutate(
        Regions=(ifelse(startsWith(RegionConcPtile,"rCx"),"Cortex",
                        ifelse(startsWith(RegionConcPtile,"rBG"),"Basal Ganglia","Thalamus"))
        )
    )
write.csv(NCA_tables,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/outputs/Tables for Pub/phenytoin_ratios.csv",
          row.names=FALSE)

CX<-NCA_sum %>% 
  mutate(
    group="Cortex",
    EC10=Cec110,
    EC50=Cec150,
    EC90=Cec190,
    IC10=Icf110,
    IC50=Icf150,
    IC90=Icf190
  )%>%
  select("x","group","EC10","EC90","EC50","IC10","IC50","IC90")

BG<-NCA_sum %>% 
  mutate(
    group="Basal Ganglia",
    EC10=Cec210,
    EC50=Cec250,
    EC90=Cec290,
    IC10=Icf210,
    IC50=Icf250,
    IC90=Icf290
  )%>%
  select("x","group","EC10","EC90","EC50","IC10","IC50","IC90")

TH<- NCA_sum %>% 
  mutate(
    group="Thalamus",
    EC10=Cec310,
    EC50=Cec350,
    EC90=Cec390,
    IC10=Icf310,
    IC50=Icf350,
    IC90=Icf390
  )%>%
  select("x","group","EC10","EC90","EC50","IC10","IC50","IC90")


ROB<- NCA_sum %>% 
  mutate(
    group="ROB",
    EC10=Cec410,
    EC50=Cec450,
    EC90=Cec490,
    IC10=Icf410,
    IC50=Icf450,
    IC90=Icf490
  )%>%
  select("x","group","EC10","EC90","EC50","IC10","IC50","IC90")


All<- rbind(CX,BG,TH,ROB)


boxEC<- ggplot(All, aes(x=group, y=EC50))+
  geom_boxplot()

boxEC

boxIC<- ggplot(All, aes(x=group,y=IC50))+
  geom_boxplot()
boxIC



boxp<- ggplot(All, aes(x=group,y=))


CxEc<-ggplot() +
  geom_ribbon(data=out_tiss, aes(time,ymin=CxEc10, ymax=CxEc90), color="blue",alpha=0.5) +
  geom_line(data=out_tiss,aes(time,CxEc50), color="black") 

CxIc<-ggplot() +
  geom_ribbon(data=out_tiss, aes(time,ymin=CxIc10, ymax=CxIc90), color="blue",alpha=0.5) +
  geom_line(data=out_tiss,aes(time,CxIc50), color="black") 




BGEc<-ggplot() +
  geom_ribbon(data=out_tiss, aes(time,ymin=BGEc10, ymax=BGEc90), color="blue",alpha=0.5) +
  geom_line(data=out_tiss,aes(time,BGEc50), color="black") 

BGIc<-ggplot() +
  geom_ribbon(data=out_tiss, aes(time,ymin=BGIc10, ymax=BGIc90), color="blue",alpha=0.5) +
  geom_line(data=out_tiss,aes(time,BGIc50), color="black") 



ALLec<-ggplot() +
  geom_ribbon(data=out_tiss, aes(time,ymin=CxEc10, ymax=CxEc90), color="blue",alpha=0.5) +
  #geom_line(data=out_tiss,aes(time,CxEc50), color="blue4")+
  geom_ribbon(data=out_tiss, aes(time,ymin=BGEc10, ymax=BGEc90), color="yellow",alpha=0.25) 
  #geom_line(data=out_tiss,aes(time,BGEc50), color="yellow4") 

ALLec