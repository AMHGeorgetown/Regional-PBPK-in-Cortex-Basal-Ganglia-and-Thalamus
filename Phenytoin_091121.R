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
 

path<- "/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/outputs/"
wd  <- "/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models"

setwd(wd)
options(knitr.table.format = "pipe")



#import data from Wilder 1977

dem <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Phenytoin/Data Files/wilder_data.csv") #IDs, IV rates, Dose, BW and groupings

demworal <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Phenytoin/Data Files//wilder_data_woral.csv") #Includes oral dosing as well for repeated doses

pt_pl <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Phenytoin/Data Files/wilder_plasma_Naremoved.csv") #plasma values
pt_cs <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Phenytoin/Data Files/wilder_CSFgroup_csf.csv")
pt_ti <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Phenytoin/Data Files/wilder_TissueGroup_Tissue.csv")
pt_pl <-transform(pt_pl,
                  conc=as.numeric(conc),
                  timechar =as.character(time)
)

t_w<-c(0.08333333, 0.116666667,0.1333333,0.133333334,0.166666667,0.183333334,0.2,0.233333333,0.25,0.266666667,0.316666666,
        0.333333333,0.333333337,0.35,0.366666667,0.416666666,0.416666667,0.499999997,0.5,0.516666667,
        0.583333333,0.583333334,0.616666667,0.65,0.666666667,0.683333334,0.75,0.833333333,0.833333337,
        0.9,0.916666666,0.999999997,1.083333333,1.116666667,1.166666667,1.249999997,1.25,1.316666663,
        1.333333333,1.333333334,1.5,1.999999997,2.166666667,2.25,2.333333333,2.5,6.166666667,6.25,
        6.333333333,6.5,12.16666667,12.25,12.33333333,12.5,24.16666667,24.25,24.33333333,24.5)

plas_pats<- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14)
tisu_pats<- c(25,27)
csf_pats <- c(28,29,30,31,32,33)

 
mod_pt<-mread("phenytoin_091121b")

set.seed(10950)
CLhep<- rnorm(10000,10950,5600)
CLhep<-CLhep[CLhep>=200 &  CLhep<=12000]

set.seed(121)
CLrenal<- rnorm(10000,121,20) 
CLrenal<-CLrenal[CLrenal<140]

set.seed(1)
Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.9 & Kpf<1.15]

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
      mod<-loadso(mod_pt)
      as.data.frame(
        mod %>%
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

out_0<-out_0%>%
              mutate(
                CSFPlr=Cp/Csa
              )


#MFE work
#All Plasma
pdose<- out_0%>% 
          filter(amt !=0)%>%  mutate(dose=amt)%>%
                            select("ID","dose","IRa","x")



out_1<- merge(out_0,pdose, by=c("ID","IRa","x"))        


psum<- out_1 %>% group_by(x,ID,dose,IRa,time)%>%
                summarise(
                  conc_mean = mean(Cp)
                )

psum2<- psum %>% group_by(x,dose,IRa,time)%>%
  summarise(
  avgconc = mean(conc_mean)
)


dem2<- dem %>% select(-("time"))%>% 
                                  filter(ID %in% plas_pats)
 

osum<- merge(dem2,pt_pl, by="ID")%>% 
  mutate( dose=amt)%>%
  select("ID","IRa","time","TAI","conc","dose")


mfe<-merge(psum2,osum, by=c("dose","IRa","time"))

mfe21<- mfe %>% 
  mutate(
    FE=ifelse(conc>0 & avgconc>0,
              (ifelse(conc>avgconc,conc/avgconc,avgconc/conc)),
              1),
    LFE=log10(FE))


avgMFE2<-mfe21%>%
summarise(
    MFE_avg=mean(LFE)
  )

Exp <- function(x) {
  a <- 10^x
}




WilderP2<- avgMFE2%>%  summarise(
    LFE=(signif(Exp(MFE_avg), digits=4))
  )
WilderP2

  


p_pt_2<- out_1%>%group_by(time)%>% 
  summarise(
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95)
  )



o_pt1_2<-osum%>% 
   group_by(time)%>%
  summarize(
    median_p=median(conc,na.rm=TRUE),
    mean_p  =mean(conc, na.rm=TRUE),
    Cp05=quantile(conc, probs=0.05),
    Cp50=quantile(conc, probs=0.50),
    Cp95=quantile(conc, probs=0.95)
  )



p12<- ggplot() + 
  geom_ribbon(data=p_pt_2, aes(time,ymin=Cp05, ymax=Cp95), color="black", alpha=0.5)+
  geom_line(data=p_pt_2, aes(time,Cp50), color="black")+
  geom_point(data=pt_pl, aes(TAI,conc))

p12<- p12+labs( x="Time (hrs)", y="Concentration (ug/mL)",
                 title="Phenytoin IV Plasma", subtitle="All patients")+theme_bw() 

p12L<-p12+labs( x="Time (hrs)", y="LOG Concentration (ug/mL)",
        title="Phenytoin IV Plasma", subtitle="All patients")+theme_bw()+
  scale_y_continuous(trans="log10")


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_AllPlasma.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_AllPlasma_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12L)
dev.off()


write.csv(WilderP2,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_MFE_allplasma.csv",row.names=TRUE)

p12+xlim(0,2)+ylim(0,58)
p12+xlim(0,5)
p12+coord_cartesian(xlim=c(0,2), ylim=c(0,60), expand=TRUE)
p12+xlim(0,2)



#CSF
#plasma for CSF subjects

pdose<- out_0%>% 
 
    filter(amt !=0)%>%  
      mutate(dose=amt)%>%
        select("ID","dose","IRa","x")
 



out_1<- merge(out_0,pdose, by=c("ID","IRa","x"))        


psum<- out_1 %>% group_by(x,ID,dose,IRa,time)%>%
  summarise(
    conc_mean = mean(Cp)
  )

psum2<- psum %>% group_by(x,dose,IRa,time)%>%filter(ID %in% csf_pats)%>%
  summarise(
    avgconc = mean(conc_mean)
  ) 
  


dem2<- dem %>% 
  select(-("time"))%>% 
    filter(ID %in% csf_pats)


osum<- merge(dem2,pt_pl, by="ID")%>%
  mutate( dose=amt)%>%
    filter(ID %in% csf_pats)%>%
        select("ID","IRa","time","conc","dose","TAI")


mfe<-merge(psum2,osum, by=c("dose","IRa","time"))

mfe21<- mfe %>% 
  mutate(
    FE=ifelse(conc>0 & avgconc>0,
              (ifelse(conc>avgconc,conc/avgconc,avgconc/conc)),
              1),
    LFE=log10(FE))


avgMFE2<-mfe21%>%
  summarise(
    MFE_avg=mean(LFE)
  )

Exp <- function(x) {
  a <- 10^x
}




WilderP2c<- avgMFE2%>%  summarise(
  LFE=(signif(Exp(MFE_avg), digits=4))
)
WilderP2c


write.csv(WilderP2,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/MFE_plasma4CSF.csv",row.names=TRUE)

#####


p_pt_2CSF<- out_1%>%group_by(time)%>%  filter(ID %in% csf_pats)%>%
  summarise(
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95),
    min = min(Cp),
    max = max(Cp)
  ) 



 


p12<- ggplot() + 
  geom_ribbon(data=p_pt_2CSF, aes(time,ymin=min, ymax=max), color="black", alpha=0.5)+
  geom_line(data=p_pt_2CSF, aes(time,Cp50), color="black")+
  geom_point(data=osum, aes(time,conc)) 

p12<-p12 +xlim(0,2)+labs( x="Time (hrs)", y="Concentration (ug/mL)",
                          title="Phenytoin IV Plasma", subtitle="CSF patients")+theme_bw() 
p12L<-p12+labs( x="Time (hrs)", y="LOG Concentration (ug/mL)",
                title="Phenytoin IV Plasma", subtitle="All patients")+theme_bw()+
                scale_y_continuous(trans="log10") 




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_Plasma4CSF.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_Plasma4CSF_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12L)
dev.off()


write.csv(WilderP2c,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_MFE_Plasma4CSF.csv",row.names=TRUE)



outmfec<- filter(out_0,time >0 &ID %in% csf_pats )

outmfe2c<- merge(pt_cs,outmfec, by=c("ID","time"))%>% select(x,time,conc,Csa,ID)

outmfe3c<- outmfe2c %>% 
  group_by(x) %>% 
  mutate(
    FE=ifelse(conc>Csa,conc/Csa,Csa/conc),
    LFE=log10(FE))%>%
  summarise(
    MFE=mean(LFE)
  )

avgMFEc<- mean(outmfe3c$MFE)
Exp <- function(x) {
  a <- 10^x
}



WilderC<- signif(Exp(avgMFEc), digits=4)
WilderC



p_pt_2c<- out_1%>%group_by(time)%>% filter(ID %in% csf_pats)%>%
  summarise(
    Cs05=quantile(Csa, probs=0.05),
    Cs50=quantile(Csa, probs=0.50),
    Cs95=quantile(Csa, probs=0.95)
  )


osum2<- merge(dem2,pt_cs, by="ID")%>% 
  mutate( dose=amt)%>%
  select("ID","IRa","time","conc","dose","tai")%>%
  mutate(
    timeg= ifelse(tai<0,0,tai)
  )

 


p12c<- ggplot() + 
  geom_ribbon(data=p_pt_2c, aes(time,ymin=Cs05, ymax=Cs95), color="black", alpha=0.5)+
  geom_line(data=p_pt_2c, aes(time,Cs50), color="black")+
  geom_point(data=pt_cs, aes(time,conc))


p12cL<-p12c+xlim(0,1.5)+labs( x="Time (hrs)", y="LOG Concentration (ug/mL)",
            title="Phenytoin IV CSF", subtitle="CSF patients")+theme_bw()+
              scale_y_continuous(trans="log10") 
p12cL


p12c<-p12c+xlim(0,1.5)+ylim(0,2.5)+
      labs( x="Time (hrs)", y="Concentration (ug/mL)",
          title="Phenytoin IV CSF", subtitle="CSF patients")+theme_bw() 
p12c




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_CSF.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12c)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_CSF_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12cL)
dev.off()


write.csv(WilderC,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_MFE_CSF.csv",row.names=TRUE)









#Plasma:CSF Ratios

csr<- out_1 %>% filter(time>0 & ID %in% csf_pats)%>% select("ID","time","Cp","Csa","CSFPlr") %>%
  mutate(
    tai=ifelse(ID==28,time-(7/60),
               ifelse(ID==29,time-(10/60),
                      ifelse(ID==30,time-(15/60),
                             ifelse(ID==31,time-(20/60),
                                    ifelse(ID==32, time-(15/60),
                                        ifelse(ID==33,time-(20/60),time))))))
   
  )

Csr_5<- csr %>% filter(tai>=0.08 &tai<=0.085 ) %>%
  summarise (
    CP_rmean=mean(CSFPlr),
    CP_rmedian=median(CSFPlr),
    CP_rmin =min(CSFPlr),
    CP_rmax= max(CSFPlr),
    CP_sd = sqrt(var(CSFPlr))
  )


Csr_1015<- csr %>% filter(tai>=.15 & tai<=0.26 )  %>%
  summarise (
    CP_rmean=mean(CSFPlr),
    CP_rmedian=median(CSFPlr),
    CP_rmin =min(CSFPlr),
    CP_rmax= max(CSFPlr),
    CP_sd = sqrt(var(CSFPlr))
  )



Csr_2530<- csr %>% filter(tai ==0.5) %>%
  summarise (
    CP_rmean=mean(CSFPlr),
    CP_rmedian=median(CSFPlr),
    CP_rmin =min(CSFPlr),
    CP_rmax= max(CSFPlr),
    CP_sd = sqrt(var(CSFPlr))
  )


Csr_50<- csr %>% filter(tai > 0.8 & tai<1.1 )  %>%
  summarise (
    CP_rmean=mean(CSFPlr),
    CP_rmedian=median(CSFPlr),
    CP_rmin =min(CSFPlr),
    CP_rmax= max(CSFPlr),
    CP_sd = sqrt(var(CSFPlr))
  )


write.csv(Csr_60,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/CSF_PlasmaR_60.csv",row.names=TRUE)



#Brain tissue 

#Plasma of brain tissue patients

pdose<- out_0%>% 
  
  filter(amt !=0)%>%  
  mutate(dose=amt)%>%
  select("ID","dose","IRa","x")




out_1<- merge(out_0,pdose, by=c("ID","IRa","x"))        


psum<- out_1 %>% group_by(x,ID,dose,IRa,time)%>%
  summarise(
    conc_mean = mean(Cp)
  )

psum2<- psum %>% group_by(x,dose,IRa,time)%>%
  summarise(
    avgconc = mean(conc_mean)
  )


dem2<- dem %>% 
  select(-("time"))%>% 
  filter(time>0 & ID %in% tisu_pats)


osum<- merge(dem2,pt_pl, by="ID")%>%
  mutate( dose=amt)%>%
  filter(ID %in% tisu_pats)%>%
  select("ID","IRa","time","conc","dose","TAI")


mfe<-merge(psum2,osum, by=c("dose","time"))

mfe21<- mfe %>% 
  mutate(
    FE=ifelse(conc>0 & avgconc>0,
              (ifelse(conc>avgconc,conc/avgconc,avgconc/conc)),
              1),
    LFE=log10(FE))


avgMFE2<-mfe21%>%
  summarise(
    MFE_avg=mean(LFE)
  )

Exp <- function(x) {
  a <- 10^x
}




WilderP2t<- avgMFE2%>%  summarise(
  LFE=(signif(Exp(MFE_avg), digits=4))
)
WilderP2t


write.csv(WilderP2t,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_MFE_Plasma4Tissue.csv",row.names=TRUE)




pdose<- out_0%>% 
  
  filter(amt !=0)%>%  
  mutate(dose=amt)%>%
  select("ID","dose","IRa","x")




out_1<- merge(out_0,pdose, by=c("ID","IRa","x"))        


psum<- out_1 %>% group_by(x,ID,dose,IRa,time)%>%
  summarise(
    conc_mean = mean(CortTis2)
  )

psum2<- psum %>% group_by(x,dose,IRa,time)%>%filter(ID %in% csf_pats)%>%
  summarise(
    avgconc = mean(conc_mean)
  )


dem2<- dem %>% 
  select(-("time"))%>% 
  filter(ID %in% tisu_pats)


osum<- merge(dem2,pt_ti, by="ID")%>%
  mutate( dose=amt)%>%
  filter(ID %in% tisu_pats)%>%
  select("ID","IRa","time","conc","dose")


mfe<-merge(psum2,osum, by=c("dose","time"))

mfe21<- mfe %>% 
  mutate(
    FE=ifelse(conc>0 & avgconc>0,
              (ifelse(conc>avgconc,conc/avgconc,avgconc/conc)),
              1),
    LFE=log10(FE))


avgMFE2<-mfe21%>%
  summarise(
    MFE_avg=mean(LFE)
  )

Exp <- function(x) {
  a <- 10^x
}




WilderT<- avgMFE2%>%  summarise(
  LFE=(signif(Exp(MFE_avg), digits=4))
)
WilderT

write.csv(WilderT,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_MFE_Tissue.csv",row.names=TRUE)



p_pt_2b<- out_1%>%group_by(time)%>%  filter(ID %in% tisu_pats) %>% 
  summarise(
    Ctis05=quantile((CortTis2), probs=0.05),
    Ctis50=quantile((CortTis2), probs=0.50),
    Ctis95=quantile((CortTis2), probs=0.95)
  )

obs_graph<-pt_ti%>% 
  filter(ID %in% tisu_pats)


p12b<- ggplot() + 
  geom_ribbon(data=p_pt_2b, aes(time,ymin=Ctis05, ymax=Ctis95), color="black", alpha=0.5)+
  geom_line(data=p_pt_2b, aes(time,Ctis50), color="black")+
  geom_point(data=obs_graph, aes(time, conc) )

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  
p12t<-p12b+xlim(0,2.5)+
   labs( x="Time (hrs)", y="Concentration (ug/mL)",
        title="Phenytoin IV Tissue", subtitle="Tissue patients")+theme_bw() 



p12tl<-p12b+scale_y_continuous(trans="log10") +xlim(0,2.5)+
  labs( x="Time (hrs)", y="Concentration (ug/mL)",title="Phenytoin IV Tissue", subtitle="Tissue patients")+
  theme_bw() 




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_tissue.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12t)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_tissue_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12tl)
dev.off()



p_pt_2t<- out_1%>%group_by(time)%>%   filter(ID %in% tisu_pats)%>%
  summarise(
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95),
    min = min(Cp),
    max = max(Cp)
  )




obs_graph<-pt_pl%>% 
  filter(ID %in% tisu_pats)


p13b<- ggplot() + 
  geom_ribbon(data=p_pt_2t, aes(time,ymin=Cp05, ymax=Cp95), color="black", alpha=0.5)+
  geom_line(data=p_pt_2t, aes(time,Cp50), color="black")+
  geom_point(data=obs_graph, aes(time, conc) )


p13tp<-p13b+xlim(0,2.5)+
  labs( x="Time (hrs)", y="Concentration (ug/mL)",
        title="Phenytoin IV Plasma", subtitle="Tissue patients")+theme_bw() 



p13tl<-p13b+scale_y_continuous(trans="log10") +xlim(0,2.5)+
  labs( x="Time (hrs)", y="Concentration (ug/mL)",title="Phenytoin IV Tissue", subtitle="Tissue patients")+
  theme_bw() 


p13tl
p13tp

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_Plasma4tissue.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p13tp)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_tPlasma4tissue_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p13tl)
dev.off()




p12tp +xlim(0,2)+
  labs( x="Time (hrs)", y="LOG Concentration (ug/mL)",
                      title="Phenytoin IV CSF", subtitle="CSF patients")+theme_bw()++scale_y_continuous(trans="log10")





p12c<- ggplot() + 
  geom_ribbon(data=p_pt_2c, aes(time,ymin=Cs05, ymax=Cs95), color="black", alpha=0.5)+
  geom_line(data=p_pt_2c, aes(time,Cs50), color="black")+
  geom_point(data=pt_cs, aes(time,conc))


p12cL<-p12c+xlim(0,1.5)+labs( x="Time (hrs)", y="LOG Concentration (ug/mL)",
                              title="Phenytoin IV CSF", subtitle="CSF patients")+theme_bw()+
  scale_y_continuous(trans="log10") 
p12cL


p12c<-p12c+xlim(0,1.5)+ylim(0,2.5)+
  labs( x="Time (hrs)", y="Concentration (ug/mL)",
        title="Phenytoin IV CSF", subtitle="CSF patients")+theme_bw() 
p12c




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_CSF.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12c)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_CSF_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p12cL)
dev.off()


write.csv(WilderC,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Phenytoin/Wilder_MFE_CSF.csv",row.names=TRUE)




###Increased BBB perm


CLhep<- rnorm(10000,10950,5600)
CLhep<-CLhep[CLhep>=200 &  CLhep<=12000]

CLrenal<- rnorm(10000,121,20) 
CLrenal<-CLrenal[CLrenal<140]
Kpf<- rnorm(10000,1,0.3)
Kpf<-Kpf[Kpf>0.4 & Kpf<1.6]

Fpbbb<-runif(10000,0.00006,0.00036)
 

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
  for(i in 0:100){ 
    it<-i+1
    
    names(it)<-"x"
    idataCLh<- expand.idata(CLhep = (sample(CLhep,33,replace=FALSE)))
    idataCLr<- expand.idata(Clrenal=(sample(CLrenal,33, replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,33,replace=FALSE)))
    idataFpbbb<- expand.idata(Fpbbb=(sample(Fpbbb,33,replace=FALSE)))

    
    itx<- merge(it,dem)
    idata1<- merge(itx,idataCLh)
    idata2<- merge(idata1,idataCLr)
    idata3<- merge(idata2,idataKpf)
    idata4<- merge(idata3,idataFpbbb)
    
    
    sp <- split(idata4,idata4$ID)#### population data set split into individual IDs to allow for parallel processing
    
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
    
    if (it==100) {
      return(out_0); break}
    
  }}


out_0<-out_0[-c(1),]

pdose<- out_0%>% 
  filter(amt !=0)%>%  mutate(dose=amt)%>%
  select("ID","dose","IRa","x")



out_1<- merge(out_0,pdose, by=c("ID","IRa","x"))        


p_pt_2b<- out_1%>%group_by(time)%>% 
  summarise(
    Cp05=quantile((CortTis2), probs=0.05),
    Cp50=quantile((CortTis2), probs=0.50),
    Cp95=quantile((CortTis2), probs=0.95)
  ) 


obs_graph<-pt_ti%>% 
  filter(ID %in% tisu_pats)


p12b<- ggplot() + 
  geom_ribbon(data=p_pt_2b, aes(time,ymin=Cp05, ymax=Cp95), color="black", alpha=0.5)+
  geom_line(data=p_pt_2b, aes(time,Cp50), color="black")+
  geom_point(data=obs_graph, aes(time, conc), color="red")

p12b+xlim(0,2)+scale_y_continuous(trans="log10")