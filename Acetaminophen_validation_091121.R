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



#MFE Work
models <- "/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models"

setwd(models)






  
  
  
  
  #dosing parameters for all Acetaminophen studies
  names_d <- as.list(c("oral", "IV","IRa"))
  ac_d_ib <- as.list(c(0,1,20)) #Data 1 Bannwarth - 1000 mg IV  3 minute infusion
  ac_d_o <- as.list(c(1,0,0))  #Data 2 - > 1000 mg/Oral /6 hour - Langford CSF + Plasma 1 
                               #Data 3 - > 1500 mg/Oral /6 hour - Langford CSF + Plasma 2 
  
  ac_d_i <- as.list(c(0,1,4)) #Data 4 - > 1000 mg/IV   /6 hour - Singla CSF + Plasma 15 minute infusion
  ac_d_i2 <- as.list(c(0,1,6)) #Data 5 - > 1000 mg/IV   /6 hour - Langford CSF + Plasma 3 10 minute infusion
  names(ac_d_o) <- names_d
  names(ac_d_i) <- names_d
  names(ac_d_i2) <- names_d
  names(ac_d_ib) <- names_d
  


######################################################################
#                                                                    #
#                      Langford                                      #
#                      CSF/Plasma                                    #
#                                                                    #
######################################################################

################
#   IV         #
################

###IV Data Import###

LangIV  <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Langford_1g_IV_plasma.csv")
LangIVC <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Langford_1g_IV_CSF.csv")


ltime<-c(0,0.25,0.5,0.75,1,1.5,2,3,6)

mod_ac <- mread("acetaminophen_091121")


set.seed(12349)
bwlan <-rnorm(10000,77000,17000) #BW Langford 
bwlan<- bwlan[bwlan>45000 & bwlan<120000]

Renal <-rnorm(10000, 780, 180 )  #renal clearance from Gaohua 2016 with ~30% CV
Renal <-Renal[Renal > 500 & Renal < 1060]
 
Hep   <-rnorm(10000,22800, 6840) #hep clearance from Gaouhia 2016 with ~30% CV 
Hep   <-Hep[Hep>500 & Hep <45200  ]

Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.25 & Kpf<1.75]

out_0<-data.frame(matrix(1:17, ncol=17))

names(out_0)<-c("ID","time","x","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva")

c <-detectCores(all.tests=FALSE, logical =TRUE) #count cores on the CPU
c<- c-1                                    #use all cores minus 1


options(future.fork.endable = TRUE)  #### from package future.apply, using forked approach
plan(multiprocess, workers=c)    #### assign num
repeat{
  for(i in 0:1000 ){ 
    it<-i+1
    
    names(it)<-"x"
    
    idataCLr<- expand.idata(CLrenal=(sample(Renal,7, replace = FALSE)))
    idataCLh<- expand.idata(CLhep = (sample(Hep,7,replace=FALSE)))
    idataBW <- expand.idata(BW = (sample(bwlan, 7, replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,7,replace=FALSE)))
    
    met<- merge(idataCLh,idataCLr)
    met2<- merge(met,idataBW)
    
    idata1<- merge(it,met2)

    idata2<- merge(idata1, ac_d_i2)
    idata3<- merge(idata2,idataKpf)
    
    
    sp <- split(idata3,idata3$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod_ac<-loadso(mread("acetaminophen_091121b"))
      as.data.frame(
        mod_ac %>%
          ev(amt=1e6)%>%
          idata_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva) %>%
          carry_out(x)%>%
          mrgsim(maxsteps=100000, tgrid=ltime))})%>%
      bind_rows()

    
    

    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    if (it==1000  ) {return(out_0);
      break}
    
  }
  
}
out_0<-out_0[-c(1),]

#write.csv(out_0,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/IV1000mg.csv",row.names=TRUE)


qplot<-out_0%>%
  group_by(time)%>% 
  summarise(
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95)
  )


p<-ggplot() +
  geom_ribbon(data=qplot, aes(time,ymin=Cp05, ymax=Cp95), color="blue",alpha=0.5) +
  geom_line(data=qplot,aes(time,Cp50), color="black")
p

 
#MB#


sum_mb<- out_0%>% group_by(x,ID)%>%
  summarise(
    Max = max(Mball),
    
    auc = auc_partial(time,Mball),
  )


expected<-as.list(c(1e6,6e6))
enms<-c("expmax","expauc")
names(expected)<- enms

sum_mb2<-as.data.frame(
  merge(sum_mb, expected)
)


pcterr_max<- with(sum_mb2, 100*((sum_mb2$Max-sum_mb2$expmax)/sum_mb2$expmax)) %>% na.omit()

Max<- mean(pcterr_max)

pcterr_auc<-with(sum_mb2, 100*((sum_mb2$auc-sum_mb2$expauc)/sum_mb2$expauc))%>% na.omit()
Auc<- mean(pcterr_auc)

graph<-as.data.frame(data1)

p1<- ggplot()+
  geom_line(data=graph, aes(time,Mball))+labs(title="Whole system mass balance", subtitle = "Acetaminophen- Langfor - IV 1000 mg dose", x="Time (hour)", y="Mass (ug)" )+
  annotate("text",x=5,y=2.5e5,label=Max)+annotate("text",x=5,y=3.1e5,label="Maximum Value Error Percentage")+
  annotate("text",x=5,y=5e5,label=Auc)+annotate("text",x=5,y=5.6e5,label="AUC Error Percentage")


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/MB1000mgIV.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p1)
dev.off()


#####INPUT LANGFORD IV NCA PARAMETERS######

###Plasma
#Literature values
cmaxL<-46.1
cmaxLmin<-21.7
cmaxLmax<-99.7


CmaxObs<-cbind(cmaxL,cmaxLmax,cmaxLmin)
cat<-"Obs"
CmaxObs<- as.data.frame(merge(CmaxObs,cat))
names(CmaxObs)<- c("Median Cmax", "Upper Cmax", "Lower Cmax","type")


auc1L<-1688
auc1Lmin<-880
auc1Lmax<-2992


Auc1Obs<-cbind(auc1L,auc1Lmax,auc1Lmin)
Auc1Obs<- as.data.frame(merge(Auc1Obs,cat))
names(Auc1Obs)<- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")


auc2L<-973
auc2Lmin<- 437
auc2Lmax<- 1647


Auc2Obs<-cbind(auc2L,auc2Lmax,auc2Lmin)
Auc2Obs<- as.data.frame(merge(Auc2Obs,cat))
names(Auc2Obs)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")
  



aucL<- 3924
aucLmin <- 2937
aucLmax <- 7323



AucObs<-cbind(aucL,aucLmax,auc2Lmin)
AucObs<- as.data.frame(merge(AucObs,cat))
names(AucObs)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")



tmaxL<-15
tmaxmin<-15
tmaxmax<-15


#Predicted values - Plasma

NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Cp),
    auc = auc_partial(time,Cp),
    auc_inf= auc_inf(time,Cp),
    tmax  = (time[which(Cp==Cmax)]
    )
  ) 

out1h<- filter(out_0, time<=1)
out2h<- filter(out_0, time>=1 & time<=2)

NCA_p1h1 <- out1h %>% group_by(x,ID)%>%
        summarise(
          auc1h=auc_partial(time,Cp)
        )
NCA_p1h2 <- out2h %>% group_by(x,ID)%>%
  summarise(
    auc2h=auc_partial(time,Cp)
  )

NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmed=median(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmed =median(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  tmax=   median(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

NCA_1h<-NCA_p1h1%>%group_by(x)%>%
                          summarise(
                            auc1hmed=median(auc1h),
                            auc1hmin=min(auc1h),
                            auc1hmax=max(auc1h)
                          )

NCA_2h<-NCA_p1h2%>%group_by(x)%>%
  summarise(
    auc2hmed=median(auc2h),
    auc2hmin=min(auc2h),
    auc2hmax=max(auc2h)
  )


cmaxP<-mean(NCA_0$Cmax)  #ng/mL = mg/L
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmed) * 60  #convert ng*hr/mL to mg*min/hr
minauc<-mean(NCA_0$minauc) * 60
maxauc<-mean(NCA_0$maxauc) * 60

tmaxP<- mean(NCA_0$tmax)*60  #convert to min
mintmax<- mean(NCA_0$tmax) * 60 
maxtmax<- mean(NCA_0$tmax) * 60

auc1hp<- mean(NCA_1h$auc1hmed) * 60 
auc1min<- mean(NCA_1h$auc1hmin) * 60 
auc1max <- mean(NCA_1h$auc1hmax) * 60



auc2hp<- mean(NCA_2h$auc2hmed) * 60 
auc2min<- mean(NCA_2h$auc2hmin) * 60 
auc2max <- mean(NCA_2h$auc2hmax) * 60



outmfe3<- outmfe2 %>% 
  group_by(x) %>% 
  mutate(
    FE=ifelse(conc>Cp,conc/Cp,Cp/conc),
    LFE=log10(FE))%>%
  summarise(
    MFE=mean(LFE)
  )

cmax_FE<- signif(ifelse(cmaxL>cmaxP,cmaxL/cmaxP,cmaxP/cmaxL), digits=4)
names(cmax_FE)<- "Cmax"


auc_FE<- signif(ifelse(aucP>aucL,aucP/aucL,aucL/aucP), digits=4) 
names(auc_FE)<- "AUC"



auc1_FE<- signif(ifelse(auc1L>auc1hp,auc1L/auc1hp,auc1hp/auc1L), digits=4)
names(auc1_FE)<-"AUC1hour"

auc2_FE<- signif(ifelse(auc2L>auc2hp,auc2L/auc2hp,auc2hp/auc2L), digits=4)
names(auc2_FE)<- "AUC2hour"

tmax_FE<- signif(ifelse(tmaxL>tmaxP,tmaxL/tmaxP,tmaxP/tmaxL),digits=4) 
names(tmax_FE)<-"Tmax"


#datasets for graphs of NCA 
pred<-"Pred"

CmaxObs #
CmaxPred<- cbind(cmaxP,cmaxmax,cmaxmin)
CmaxPred<-merge(CmaxPred,pred)
names(CmaxPred)<- c("Median Cmax", "Upper Cmax", "Lower Cmax", "type")

CmaxG<- rbind(CmaxObs,CmaxPred)
CmaxG


t_cmax<- ggplot(CmaxG, aes(x=type, y=cmaxP, color=type))+
        geom_pointrange(aes(ymin=cmaxmin, ymax=cmaxmax))+
        labs(title="Cmax", x="Source", y="Cmax (ng/mL)")

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgIVtmax.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_cmax)
dev.off()




AucObs
AucPred<- cbind(aucP,maxauc,minauc)
AucPred<- merge(AucPred,pred)
names(AucPred)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


AucG <- rbind(AucObs,AucPred)  
AucG


t_auc<- ggplot(AucG, aes(x=type, y=aucP, color=type))+
  geom_pointrange(aes(ymin=minauc, ymax=maxauc))+
  labs(title="AUC0-6", x="Source", y="AUC (ng*min/mL)")

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgIVtauc.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc)
dev.off()
 


Auc1Obs
Auc1Pred<- cbind(auc1hp,auc1max,auc1min)
Auc1Pred<- merge(Auc1Pred,pred)
names(Auc1Pred) <- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")

Auc1hG<- rbind(Auc1Obs,Auc1Pred)
Auc1hG

t_auc1h <- ggplot(Auc1hG, aes(x=type, y=auc1hp, color=type))+
  geom_pointrange(aes(ymin=auc1min, ymax=auc1max))+
  labs(title="AUC0-1", x="Source", y="AUC (ng*min/mL)")

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgIVtauc1h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc1h)
dev.off()


Auc2Obs
Auc2Pred<- cbind(auc2hp, auc2max, auc2min)
Auc2Pred<- merge(Auc2Pred,pred)
names(Auc2Pred)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")

Auc2hG<- rbind(Auc2Obs,Auc2Pred)
Auc2hG


t_auc2h <- ggplot(Auc2hG, aes(x=type, y=auc2hp, color=type))+
  geom_pointrange(aes(ymin=auc2min, ymax=auc2max))+
  labs(title="AUC2-6", x="Source", y="AUC (ng*min/mL)")

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgIVtauc2h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc2h)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgIVNCA.png",width= 9 * ppi, height=4*ppi, res=ppi)

grid.arrange(t_cmax,t_auc,t_auc1h,t_auc2h)
dev.off()





all<-cbind(cmax_FE,auc_FE,auc1_FE,auc2_FE,tmax_FE)

 
LFE<-log10(all)

mean_LFE<-mean(LFE)







Exp <- function(x) {
  a <- 10^x
}

MFE_L<-Exp(mean_LFE)
names(MFE_L)<-("Mean Log-Fold Error")

 
MFE_L
all


write.csv(MFE_L,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/MFELIV1000mg.csv",row.names=TRUE)
write.csv(all,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/ALLFEIV1000mg.csv",row.names=TRUE)



#CSF analysis


#Literature values
cmaxLc<-12.6
cmaxLcmin<-4.2
cmaxLcmax<-15.5


CmaxObsc<-cbind(cmaxLc,cmaxLcmax,cmaxLcmin)
 
CmaxObsc<- as.data.frame(merge(CmaxObsc,cat))
names(CmaxObsc)<- c("Median Cmax", "Upper Cmax", "Lower Cmax","type")


auc1Lc<-110
auc1Lcmin<-45
auc1Lcmax<-275


Auc1Obsc<-cbind(auc1Lc,auc1Lcmax,auc1Lcmin)
Auc1Obsc<- as.data.frame(merge(Auc1Obsc,cat))
names(Auc1Obsc)<- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")


auc2Lc<-547
auc2Lcmin<- 183
auc2Lcmax<- 806


Auc2Obsc<-cbind(auc2Lc,auc2Lcmax,auc2Lcmin)
Auc2Obsc<- as.data.frame(merge(Auc2Obsc,cat))
names(Auc2Obsc)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")




aucLc<- 3488
aucLcmin <- 1083
aucLcmax <- 3724



AucObsc<-cbind(aucLc,aucLcmax,aucLcmin)
AucObsc<- as.data.frame(merge(AucObsc,cat))
names(AucObsc)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")



tmaxLc<-180
tmaxLcmin<-90
tmaxLcmax<-360
TmaxObsc<- cbind(tmaxLc,tmaxLcmin,tmaxLcmax)
TmaxObsc<- as.data.frame(merge(TmaxObsc,cat))

#Predicted values - CSF

NCA_c1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Csa),
    auc = auc_partial(time,Csa),
    auc_inf= auc_inf(time,Csa),
    tmax  = (time[which(Csa==Cmax)]
    )
  ) 

out1h<- filter(out_0, time<=1)
out2h<- filter(out_0, time>=1 & time<=2)

NCA_c1h1 <- out1h %>% group_by(x,ID)%>%
  summarise(
    auc1h=auc_partial(time,Csa)
  )
NCA_c1h2 <- out2h %>% group_by(x,ID)%>%
  summarise(
    auc2h=auc_partial(time,Csa)
  )

NCA_0c<- as.data.frame(NCA_c1%>%
                        group_by(x)%>%
                        summarise(Cmaxmed=median(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmed =median(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  tmax=   median(tmax))
)

NCA_1hc<-NCA_c1h1%>%group_by(x)%>%
  summarise(
    auc1hmed=median(auc1h),
    auc1hmin=min(auc1h),
    auc1hmax=max(auc1h)
  )

NCA_2hc<-NCA_c1h2%>%group_by(x)%>%
  summarise(
    auc2hmed=median(auc2h),
    auc2hmin=min(auc2h),
    auc2hmax=max(auc2h)
  )

cmaxPc<-mean(NCA_0c$Cmax)  #ng/mL = mg/L
cmaxminc<-mean(NCA_0c$cmaxmin)
cmaxmaxc<-mean(NCA_0c$cmaxmax)


aucPc<- mean(NCA_0c$aucmed) * 60  #convert ng*hr/mL to mg*min/hr
minaucc<-mean(NCA_0c$minauc) * 60
maxaucc<-mean(NCA_0c$maxauc) * 60

tmaxPc<- mean(NCA_0c$tmax)*60  #convert to min
mintmaxc<- mean(NCA_0c$tmax) * 60 
maxtmaxc<- mean(NCA_0c$tmax) * 60

auc1hpc<- mean(NCA_1hc$auc1hmed) * 60 
auc1minc<- mean(NCA_1hc$auc1hmin) * 60 
auc1maxc <- mean(NCA_1hc$auc1hmax) * 60



auc2hpc<- mean(NCA_2hc$auc2hmed) * 60 
auc2minc<- mean(NCA_2hc$auc2hmin) * 60 
auc2maxc <- mean(NCA_2hc$auc2hmax) * 60

cmax_FEc<-signif(ifelse(cmaxLc>cmaxPc,cmaxLc/cmaxPc,cmaxPc/cmaxLc),digits=4)
names(cmax_FEc)<- "Cmax"


auc_FEc<-signif(ifelse(aucLc>aucPc,aucLc/aucPc,aucPc/aucLc),digits=4)
names(auc_FEc)<- "AUC"


auc1_FEc<- signif(ifelse(auc1Lc>auc1hpc,auc1Lc/auc1hpc,auc1hp/auc1Lc),digits=4)
names(auc1_FEc)<-"AUC1hour"

auc2_FEc<- signif(ifelse(auc2Lc>auc2hpc,auc2Lc/auc2hpc,auc2hpc/auc2Lc),digits=4)
names(auc2_FEc)<- "AUC2hour"

tmax_FEc<- signif(ifelse(tmaxPc>tmaxLc,tmaxPc/tmaxLc,tmaxLc/tmaxPc),digits=4)
names(tmax_FEc)<-"Tmax"


#datasets for graphs of NCA 
pred<-"Pred"

CmaxObsc #
CmaxPredc<- cbind(cmaxPc,cmaxmaxc,cmaxminc)
CmaxPredc<-merge(CmaxPredc,pred)
names(CmaxPredc)<- c("Median Cmax", "Upper Cmax", "Lower Cmax", "type")

CmaxGc<- rbind(CmaxObsc,CmaxPredc)
CmaxGc


t_cmaxc<- ggplot(CmaxGc, aes(x=type, y=cmaxPc, color=type))+
  geom_pointrange(aes(ymin=cmaxminc, ymax=cmaxmaxc))+
  labs(title="Cmax", x="Source", y="Cmax (ng/mL)")+
  annotate("text",x=1.5,y=10,label="Mean Fold-Error")+annotate("text", x=1.5, y=8, label=cmax_FEc)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgIVcmax.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_cmaxc)
dev.off()




AucObsc
AucPredc<- cbind(aucPc,maxaucc,minaucc)
AucPredc<- merge(AucPredc,pred)
names(AucPredc)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


AucGc <- rbind(AucObsc,AucPredc)  
AucGc


t_aucc<- ggplot(AucGc, aes(x=type, y=aucPc, color=type))+
  geom_pointrange(aes(ymin=minaucc, ymax=maxaucc))+
  labs(title="AUC0-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=4300,label="Mean Fold-Error")+annotate("text", x=1.5, y=4000, label=auc_FEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgIVtauc.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_aucc)
dev.off()



Auc1Obsc
Auc1Predc<- cbind(auc1hpc,auc1maxc,auc1minc)
Auc1Predc<- merge(Auc1Predc,pred)
names(Auc1Predc) <- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")

Auc1hGc<- rbind(Auc1Obsc,Auc1Predc)
Auc1hGc

t_auc1hc <- ggplot(Auc1hGc, aes(x=type, y=auc1hpc, color=type))+
  geom_pointrange(aes(ymin=auc1minc, ymax=auc1maxc))+
  labs(title="AUC0-1", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=300,label="Mean Fold-Error")+annotate("text", x=1.5, y=275, label=auc1_FEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgIVtauc1h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc1hc)
dev.off()


Auc2Obsc
Auc2Predc<- cbind(auc2hpc, auc2maxc, auc2minc)
Auc2Predc<- merge(Auc2Predc,pred)
names(Auc2Predc)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")

Auc2hGc<- rbind(Auc2Obsc,Auc2Predc)
Auc2hGc


t_auc2hc <- ggplot(Auc2hGc, aes(x=type, y=auc2hpc, color=type))+
  geom_pointrange(aes(ymin=auc2minc, ymax=auc2maxc))+
  labs(title="AUC2-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=400,label="Mean Fold-Error")+annotate("text", x=1.5, y=375, label=auc2_FEc)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgIVtauc2h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc2hc)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgIVNCA.png",width= 9 * ppi, height=4*ppi, res=ppi)

grid.arrange(t_cmaxc,t_aucc,t_auc1hc,t_auc2hc)
dev.off()



allc<-cbind(cmax_FEc,auc_FEc,auc1_FEc,auc2_FEc,tmax_FEc)


LFEc<-log10(allc)

mean_LFEc<-mean(LFEc)







 
MFE_Lc<-Exp(mean_LFEc)
allc
MFE_Lc

write.csv(MFE_Lc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSFMFELIV1000mg.csv",row.names=TRUE)
write.csv(allc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSFALLFEIV1000mg.csv",row.names=TRUE)






################
# Oral 1000    #
################

LangO1k  <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Langford_1g_PO_Plasma.csv")
LangO1kc <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Langford_1g_PO_CSF.csv")




mod_ac <- mread("acetaminophen_091121b")


set.seed(12349)
bwlan <-rnorm(10000,69500,14000) #BW Langford 
bwlan<-bwlan[bwlan>60000 & bwlan<120000]



Renal <-rnorm(10000, 1400, 180 )  #renal clearance from Gaohua 2016 with ~30% CV
Renal <-Renal[Renal > 500]

Hep   <-rnorm(10000,44000, 6840) #hep clearance from Gaouhia 2016 with ~30% CV 
Hep   <-Hep[Hep>500  ]
Abs   <- rnorm(10000,2,0.5) #Calculated from abs half life Divoll 1982
Abs <-Abs[Abs>0.5 & Abs<3.5]
Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.25 & Kpf<1.75]

c <-detectCores(all.tests=FALSE, logical =TRUE) #count cores on the CPU
#use all cores minus 1


options(future.fork.endable = TRUE)  #### from package future.apply, using forked approach
plan(multiprocess, workers=c)    #### assign num


out_0<-data.frame(matrix(1:17, ncol=17))

names(out_0)<-c("ID","time","x","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva")


repeat{
  for(i in 0:1000){ 
    it<-i+1
    
    names(it)<-"x"
    
    idataCLr<- expand.idata(CLrenal=(sample(Renal,7, replace = FALSE)))
    idataCLh<- expand.idata(CLhep = (sample(Hep,7,replace=FALSE)))
    idataBW <- expand.idata(BW = (sample(bwlan, 7, replace=FALSE)))
    idataAbs <- expand.idata(Abs=(sample(Abs,7, replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,7,replace=FALSE)))
    
    met<- merge(idataCLh,idataCLr)
    met2<- merge(met,idataBW)
    
    idata1<- merge(it,met2)
    idata2<- merge(idata1,idataAbs)
    idata3<- merge(idata2,idataKpf)
    idata4<- merge(idata3, ac_d_o)
    
    
    
    sp <- split(idata4,idata4$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod_ac<-loadso(mread("acetaminophen_091121b"))
      as.data.frame(
        mod_ac %>%
          ev(amt=1e6)%>%
          idata_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva) %>%
          carry_out(x)%>%
          mrgsim(maxsteps=100000, tgrid=ltime))})%>%
      bind_rows()
    
     
    

    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    if (it==1000) {
      return(out_0); break}
    
  }
  
}
out_0<-out_0[-c(1),]


###PLOTS 

#MB#


sum_mb<- out_0%>% group_by(x,ID)%>%
  summarise(
    Max = max(Mball),
    
    auc = auc_partial(time,Mball),
  )


expected<-as.list(c(1e6,6e6))
enms<-c("expmax","expauc")
names(expected)<- enms

sum_mb2<-as.data.frame(
  merge(sum_mb, expected)
)


pcterr_max<- with(sum_mb2, 100*((sum_mb2$Max-sum_mb2$expmax)/sum_mb2$expmax)) %>% na.omit()

Max<- mean(pcterr_max)

pcterr_auc<-with(sum_mb2, 100*((sum_mb2$auc-sum_mb2$expauc)/sum_mb2$expauc))%>% na.omit()
Auc<- mean(pcterr_auc)

graph<-as.data.frame(data1)

p1<- ggplot()+
  geom_line(data=graph, aes(time,Mball))+labs(title="Whole system mass balance", subtitle = "Acetaminophen- Langfor - Oral 1000 mg dose", x="Time (hour)", y="Mass (ug)" )+
  annotate("text",x=5,y=2.5e5,label=Max)+annotate("text",x=5,y=3.1e5,label="Maximum Value Error Percentage")+
  annotate("text",x=5,y=5e5,label=Auc)+annotate("text",x=5,y=5.6e5,label="AUC Error Percentage")


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/MB1000mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p1)
dev.off()



#####INPUT LANGFORD Oral NCA PARAMETERS######

#write.csv(out_0,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/Oral1000mg.csv",row.names=TRUE)


qplot<-out_0%>%
    group_by(time)%>% 
      summarise(
        Cp05=quantile(Cp, probs=0.05),
        Cp50=quantile(Cp, probs=0.50),
        Cp95=quantile(Cp, probs=0.95)
      )


p<-ggplot() +
  geom_ribbon(data=qplot, aes(time,ymin=Cp05, ymax=Cp95), color="blue",alpha=0.5) +
  geom_line(data=qplot,aes(time,Cp50), color="black")
p
#NCA 

#####INPUT LANGFORD Oral 1000 mg NCA PARAMETERS######

###Plasma
#Literature values
cmaxL<-18
cmaxLmin<-2.8
cmaxLmax<-30.8


CmaxObs<-cbind(cmaxL,cmaxLmax,cmaxLmin)
cat<-"Obs"
CmaxObs<- as.data.frame(merge(CmaxObs,cat))
names(CmaxObs)<- c("Median Cmax", "Upper Cmax", "Lower Cmax","type")


auc1L<-87
auc1Lmin<-0
auc1Lmax<-907


Auc1Obs<-cbind(auc1L,auc1Lmax,auc1Lmin)
Auc1Obs<- as.data.frame(merge(Auc1Obs,cat))
names(Auc1Obs)<- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")


auc2L<-283
auc2Lmin<- 53
auc2Lmax<- 1775


Auc2Obs<-cbind(auc2L,auc2Lmax,auc2Lmin)
Auc2Obs<- as.data.frame(merge(Auc2Obs,cat))
names(Auc2Obs)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")




aucL<- 2659
aucLmin <- 527
aucLmax <- 5616



AucObs<-cbind(aucL,aucLmax,auc2Lmin)
AucObs<- as.data.frame(merge(AucObs,cat))
names(AucObs)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")



tmaxL<-120
tmaxLmin<- 30
tmaxLmax<- 360

#Predicted values - Plasma

NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Cp),
    auc = auc_partial(time,Cp),
    auc_inf= auc_inf(time,Cp),
    tmax  = (time[which(Cp==Cmax)]
    )
  ) 

out1h<- filter(out_0, time<=1)
out2h<- filter(out_0, time>=1 & time<=2)

NCA_p1h1 <- out1h %>% group_by(x,ID)%>%
  summarise(
    auc1h=auc_partial(time,Cp)
  )
NCA_p1h2 <- out2h %>% group_by(x,ID)%>%
  summarise(
    auc2h=auc_partial(time,Cp)
  )

NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmed=median(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmed =median(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  tmax=   median(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

NCA_1h<-NCA_p1h1%>%group_by(x)%>%
  summarise(
    auc1hmed=median(auc1h),
    auc1hmin=min(auc1h),
    auc1hmax=max(auc1h)
  )

NCA_2h<-NCA_p1h2%>%group_by(x)%>%
  summarise(
    auc2hmed=median(auc2h),
    auc250  =quantile(auc2h,probs=0.5),
    auc2hmin=min(auc2h),
    auc2hmax=max(auc2h)
  )


cmaxP<-mean(NCA_0$Cmaxmed)  #ng/mL = mg/L
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmed) * 60  #convert ng*hr/mL to mg*min/L
minauc<-mean(NCA_0$minauc) * 60
maxauc<-mean(NCA_0$maxauc) * 60

tmaxP<- mean(NCA_0$tmax)*60  #convert to min
mintmax<- mean(NCA_0$mintmax) * 60 
maxtmax<- mean(NCA_0$maxtmax) * 60

auc1hp<- mean(NCA_1h$auc1hmed) * 60 
auc1min<- mean(NCA_1h$auc1hmin) * 60 
auc1max <- mean(NCA_1h$auc1hmax) * 60



auc2hp<- mean(NCA_2h$auc250) * 60 
auc2min<- mean(NCA_2h$auc2hmin) * 60 
auc2max <- mean(NCA_2h$auc2hmax) * 60

cmax_FE<- signif(ifelse(cmaxL>cmaxP,cmaxL/cmaxP,cmaxP/cmaxL), digits=4)
names(cmax_FE)<- "Cmax"


auc_FE<- signif(ifelse(aucP>aucL,aucP/aucL,aucL/aucP), digits=4) 
names(auc_FE)<- "AUC"


auc1_FE<- signif(ifelse(auc1hp>auc1L,auc1hp/auc1L, auc1L/auc1hp), digits=4)
names(auc1_FE)<-"AUC1hour"

auc2_FE<- signif(ifelse(auc2hp>auc2L,auc2hp/auc2L,auc2L/auc2hp), digits=4)
names(auc2_FE)<- "AUC2hour"

tmax_FE<- signif(ifelse(tmaxL>=tmaxP,tmaxL/tmaxP,tmaxP/tmaxL),digits=4)
names(tmax_FE)<-"Tmax"


#datasets for graphs of NCA 
pred<-"Pred"

CmaxObs #
CmaxPred<- cbind(cmaxP,cmaxmax,cmaxmin)
CmaxPred<-merge(CmaxPred,pred)
names(CmaxPred)<- c("Median Cmax", "Upper Cmax", "Lower Cmax", "type")

CmaxG<- rbind(CmaxObs,CmaxPred)
CmaxG


t_cmax<- ggplot(CmaxG, aes(x=type, y=cmaxP, color=type))+
  geom_pointrange(aes(ymin=cmaxmin, ymax=cmaxmax))+
  labs(title="Cmax", x="Source", y="Cmax (ng/mL)")+
  annotate("text",x=1.5,y=25,label="Mean Fold-Error")+annotate("text", x=1.5, y=20, label=cmax_FE)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/Cmax1000mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_cmax)
dev.off()


AucObs
AucPred<- cbind(aucP,maxauc,minauc)
AucPred<- merge(AucPred,pred)
names(AucPred)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


AucG <- rbind(AucObs,AucPred)  
AucG


t_auc<- ggplot(AucG, aes(x=type, y=aucP, color=type))+
  geom_pointrange(aes(ymin=minauc, ymax=maxauc))+
  labs(title="AUC0-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=3900,label="Mean Fold-Error")+annotate("text", x=1.5, y=3700, label=auc_FE)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/AUC1000mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc)
dev.off()





Auc1Obs
Auc1Pred<- cbind(auc1hp,auc1max,auc1min)
Auc1Pred<- merge(Auc1Pred,pred)
names(Auc1Pred) <- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")

Auc1hG<- rbind(Auc1Obs,Auc1Pred)
Auc1hG

t_auc1h <- ggplot(Auc1hG, aes(x=type, y=auc1hp, color=type))+
  geom_pointrange(aes(ymin=auc1min, ymax=auc1max))+
  labs(title="AUC0-1", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=1500,label="Mean Fold-Error")+annotate("text", x=1.5, y=1400, label=auc1_FE)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/AUC011000mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc1h)
dev.off()


Auc2Obs
Auc2Pred<- cbind(auc2hp, auc2max, auc2min)
Auc2Pred<- merge(Auc2Pred,pred)
names(Auc2Pred)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")

Auc2hG<- rbind(Auc2Obs,Auc2Pred)
Auc2hG


t_auc2h <- ggplot(Auc2hG, aes(x=type, y=auc2hp, color=type))+
  geom_pointrange(aes(ymin=auc2min, ymax=auc2max))+
  labs(title="AUC2-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=1000,label="Mean Fold-Error")+annotate("text", x=1.5, y=800, label=auc2_FE)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/AUC011000mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc2h)
dev.off()




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/NCA1000mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
grid.arrange(t_cmax,t_auc,t_auc1h,t_auc2h)

dev.off()


alloral<-cbind(cmax_FE,auc_FE,auc1_FE,auc2_FE,tmax_FE)


LFE<-log10(all)

mean_LFE<-mean(LFE)







Exp <- function(x) {
  a <- 10^x
}

MFE_L1000oral<-Exp(mean_LFE)
alloral
MFE_L1000oral


write.csv(MFE_L1000oral,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/ORALMFELIV1000mg.csv",row.names=TRUE)
write.csv(alloral1000,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/ALLFEORAL1000mg.csv",row.names=TRUE)







#CSF analysis


###CSF
#Literature values
cmaxLc<-8.8
cmaxLcmin<-4.1
cmaxLcmax<-20.2


CmaxObsc<-cbind(cmaxLc,cmaxLcmax,cmaxLcmin)

CmaxObsc<- as.data.frame(merge(CmaxObsc,cat))
names(CmaxObsc)<- c("Median Cmax", "Upper Cmax", "Lower Cmax","type")


auc1Lc<-11
auc1Lcmin<-0
auc1Lcmax<-52


Auc1Obsc<-cbind(auc1Lc,auc1Lcmax,auc1Lcmin)
Auc1Obsc<- as.data.frame(merge(Auc1Obsc,cat))
names(Auc1Obsc)<- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")


auc2Lc<-175
auc2Lcmin<- 0
auc2Lcmax<- 321


Auc2Obsc<-cbind(auc2Lc,auc2Lcmax,auc2Lcmin)
Auc2Obsc<- as.data.frame(merge(Auc2Obsc,cat))
names(Auc2Obsc)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")




aucLc<- 3488
aucLcmin <- 45
aucLcmax <- 275



AucObsc<-cbind(aucLc,aucLcmax,auc2Lcmin)
AucObsc<- as.data.frame(merge(AucObsc,cat))
names(AucObsc)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")



tmaxLc<-180
tmaxLcmin<-90
tmaxLcmax<-180
TmaxObsc<- cbind(tmaxLc,tmaxLcmax,tmaxLcmin)
TmaxObsc<- as.data.frame(merge(TmaxObsc,cat))
names(TmaxObsc)<- c("Median Tmax", "Upper Tmax", "Lower Tmax", "type")
#Predicted values - CSF

NCA_c1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Csa),
    auc = auc_partial(time,Csa),
    auc_inf= auc_inf(time,Csa),
    tmax  = (time[which(Csa==Cmax)]
    )
  ) 

out1h<- filter(out_0, time<=1)
out2h<- filter(out_0, time>=1 & time<=2)

NCA_c1h1 <- out1h %>% group_by(x,ID)%>%
  summarise(
    auc1h=auc_partial(time,Csa)
  )
NCA_c1h2 <- out2h %>% group_by(x,ID)%>%
  summarise(
    auc2h=auc_partial(time,Csa)
  )

NCA_0c<- as.data.frame(NCA_c1%>%
                         group_by(x)%>%
                         summarise(Cmaxmed=median(Cmax),
                                   cmaxmin=min(Cmax),
                                   cmaxmax=max(Cmax),
                                   aucmed =median(auc),
                                   minauc=min(auc),
                                   maxauc=max(auc),
                                   tmaxmed=   median(tmax),
                                   tmaxmin=min(tmax),
                                   tmaxmax=max(tmax)
))

NCA_1hc<-NCA_c1h1%>%group_by(x)%>%
  summarise(
    auc1hmed=median(auc1h),
    auc1hmin=min(auc1h),
    auc1hmax=max(auc1h)
  )

NCA_2hc<-NCA_c1h2%>%group_by(x)%>%
  summarise(
    auc2hmed=median(auc2h),
    auc2hmin=min(auc2h),
    auc2hmax=max(auc2h)
  )


cmaxPc<-mean(NCA_0c$Cmax)  #ng/mL = mg/L
aucPc<- mean(NCA_0c$auc) * 60  #convert ng*hr/mL to mg*min/hr


tmaxPc<- mean(NCA_0c$tmaxmed)*60  #convert to min
tmaxminPc<-mean(NCA_0c$tmaxmin)* 60
tmaxmaxPc<- mean(NCA_0c$tmaxmax) * 60


auc1hc<- mean(NCA_1hc$auc1hmed) * 60 
auc2hc<- mean(NCA_2hc$auc2hmed) * 60 

cmax_FEc<- signif(ifelse(cmaxLc>cmaxPc,cmaxLc/cmaxPc,cmaxPc/cmaxLc), digits=4)
names(cmax_FEc)<- "Cmax"


auc_FEc<- signif(ifelse(aucLc>aucPc,aucLc/aucPc,aucPc/aucLc), digits=4)
names(auc_FEc)<- "AUC"


auc1_FEc<- signif(ifelse(auc1hc>auc1Lc,auc1hc/auc1Lc,auc1Lc/auc1hc), digits=4)
names(auc1_FEc)<-"AUC1hour"

auc2_FEc<- signif(ifelse(auc2Lc>auc2hc,auc2Lc/auc2hc,auc2hc/auc2Lc), digits=4)
names(auc2_FEc)<- "AUC2hour"

tmax_FEc<- signif(ifelse(tmaxPc>tmaxLc,tmaxPc/tmaxLc, tmaxLc/tmaxPc),digits=4)
names(tmax_FEc)<-"Tmax"


#datasets for graphs of NCA 
pred<-"Pred"


TmaxObsc #
TmaxPredc<- cbind(tmaxPc,tmaxmaxPc,tmaxminPc)
TmaxPredc<-merge(TmaxPredc,pred)
names(TmaxPredc)<- c("Median Tmax", "Upper Tmax", "Lower Tmax", "type")

TmaxGc<- rbind(TmaxObsc,TmaxPredc)
TmaxGc


t_tmaxc<- ggplot(TmaxGc, aes(x=type, y=tmaxPc, color=type))+
  geom_pointrange(aes(ymin=tmaxminPc, ymax=tmaxmaxPc))+
  labs(title="Tmax", x="Source", y="Tmax (min)")+
  annotate("text",x=1.5,y=325,label="Mean Fold-Error")+annotate("text", x=1.5, y=300, label=tmax_FEc)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgOralTmax.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_tmaxc)
dev.off()



CmaxObsc #
CmaxPredc<- cbind(cmaxPc,cmaxmaxc,cmaxminc)
CmaxPredc<-merge(CmaxPredc,pred)
names(CmaxPredc)<- c("Median Cmax", "Upper Cmax", "Lower Cmax", "type")

CmaxGc<- rbind(CmaxObsc,CmaxPredc)
CmaxGc


t_cmaxc<- ggplot(CmaxGc, aes(x=type, y=cmaxPc, color=type))+
  geom_pointrange(aes(ymin=cmaxminc, ymax=cmaxmaxc))+
  labs(title="Cmax", x="Source", y="Cmax (ng/mL)")+
  annotate("text",x=1.5,y=8,label="Mean Fold-Error")+annotate("text", x=1.5, y=6, label=cmax_FEc)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgOralcmax.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_cmaxc)
dev.off()




AucObsc
AucPredc<- cbind(aucPc,maxaucc,minaucc)
AucPredc<- merge(AucPredc,pred)
names(AucPredc)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


AucGc <- rbind(AucObsc,AucPredc)  
AucGc


t_aucc<- ggplot(AucGc, aes(x=type, y=aucPc, color=type))+
  geom_pointrange(aes(ymin=minaucc, ymax=maxaucc))+
  labs(title="AUC0-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=4300,label="Mean Fold-Error")+annotate("text", x=1.5, y=4000, label=auc_FEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1000mgOralauc.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_aucc)
dev.off()



Auc1Obsc
Auc1Predc<- cbind(auc1hpc,auc1maxc,auc1minc)
Auc1Predc<- merge(Auc1Predc,pred)
names(Auc1Predc) <- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")

Auc1hGc<- rbind(Auc1Obsc,Auc1Predc)
Auc1hGc

t_auc1hc <- ggplot(Auc1hGc, aes(x=type, y=auc1hpc, color=type))+
  geom_pointrange(aes(ymin=auc1minc, ymax=auc1maxc))+
  labs(title="AUC0-1", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=300,label="Mean Fold-Error")+annotate("text", x=1.5, y=275, label=auc1_FEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgOraltauc1h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc1hc)
dev.off()


Auc2Obsc
Auc2Predc<- cbind(auc2hpc, auc2maxc, auc2minc)
Auc2Predc<- merge(Auc2Predc,pred)
names(Auc2Predc)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")

Auc2hGc<- rbind(Auc2Obsc,Auc2Predc)
Auc2hGc


t_auc2hc <- ggplot(Auc2hGc, aes(x=type, y=auc2hpc, color=type))+
  geom_pointrange(aes(ymin=auc2minc, ymax=auc2maxc))+
  labs(title="AUC2-6", x="Source", y="AUC (ng*min/mL)") 

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgORaltauc2h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc2hc)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1000mgOralNCA.png",width= 9 * ppi, height=4*ppi, res=ppi)

grid.arrange(t_cmaxc,t_aucc,t_auc1hc,t_auc2hc,t_tmaxc, ncol=3,nrow=2)
dev.off()



allc<-cbind(cmax_FEc,auc_FEc,auc1_FEc,auc2_FEc,tmax_FEc)


LFEc<-log10(allc)

mean_LFEc<-mean(LFEc)








MFE_Lc<-Exp(mean_LFEc)
allc
MFE_Lc

write.csv(MFE_Lc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSFMFEORAL1000mg.csv",row.names=TRUE)
write.csv(allc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSFALLFEORAL1000mg.csv",row.names=TRUE)









################
# Oral 1500    #
################

LangO15k  <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Langford_1.5g_PO_Plasma.csv")
LangO15kc <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Langford_1.5g_PO_CSF.csv")

mod_ac <- mread("acetaminophen_091121b")


set.seed(12349)
bwlan <-rnorm(10000,98000,23000) #BW Langford 
 
bwlan<- bwlan[bwlan>70000 & bwlan<120000]



Renal <-rnorm(10000, 1400, 180 )  #renal clearance from Gaohua 2016 with ~30% CV
Renal <-Renal[Renal > 500]

Hep   <-rnorm(10000,44000, 6840) #hep clearance from Gaouhia 2016 with ~30% CV 
Hep   <-Hep[Hep>500  ]
Abs   <- rnorm(10000,2,0.5) #Calculated from abs half life Divoll 1982
Abs <-Abs[Abs>0.5 & Abs<3.5]

Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.25 & Kpf<1.75]


out_0<-data.frame(matrix(1:17, ncol=17))

names(out_0)<-c("ID","time","x","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva")


repeat{
  for(i in 0:1000){ 
    it<-i+1
    
    names(it)<-"x"
    
    idataCLr<- expand.idata(CLrenal=(sample(Renal,7, replace = FALSE)))
    idataCLh<- expand.idata(CLhep = (sample(Hep,7,replace=FALSE)))
    idataBW <- expand.idata(BW = (sample(bwlan, 7, replace=FALSE)))
    idataAbs <- expand.idata(Abs=(sample(Abs,7, replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,7,replace=FALSE)))
    met<- merge(idataCLh,idataCLr)
    met2<- merge(met,idataBW)
    
    idata1<- merge(it,met2)
    idata2<- merge(idata1,idataAbs)
    idata3<- merge(idata2, ac_d_o)
    idata4<- merge(idata3,idataKpf)
    
    
    sp <- split(idata4,idata4$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod_ac<-loadso(mread("acetaminophen_091121b"))
      as.data.frame(
        mod_ac %>%
          ev(amt=1.5e6)%>%
          idata_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva) %>%
          carry_out(x)%>%
          mrgsim(maxsteps=100000, tgrid=ltime))})%>%
      bind_rows()
    
    
    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    if (it==1000) {
      return(out_0); break}
    
  }
  
}
out_0<-out_0[-c(1),]




#####INPUT LANGFORD Oral NCA PARAMETERS 1500 mg######

###Plasma
cmaxL<-20.8
cmaxLmin<-11.6
cmaxLmax<-36

CmaxObs<-cbind(cmaxL,cmaxLmax,cmaxLmin)
cat<-"Obs"
CmaxObs<- as.data.frame(merge(CmaxObs,cat))
names(CmaxObs)<- c("Median Cmax", "Upper Cmax", "Lower Cmax","type")

aucL<- 3878
aucLmin <- 2823
aucLmax<- 6582


AucObs<-cbind(aucL,aucLmax,auc2Lmin)
AucObs<- as.data.frame(merge(AucObs,cat))
names(AucObs)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


  tmaxL <- 90
minTmaxL<- 45
maxTmaxL<- 240

TmaxObs<- cbind(tmaxL,minTmaxL,maxTmaxL)
TmaxObs<- as.data.frame(merge(TmaxObs,cat))

auc1L<- 253
auc1Lmin <- 3
auc1Lmax <- 928

Auc1Obs<-cbind(auc1L,auc1Lmax,auc1Lmin)
Auc1Obs<- as.data.frame(merge(Auc1Obs,cat))
names(Auc1Obs)<- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")



auc2L <- 954
auc2Lmin <- 353
auc2Lmax <- 1731


Auc2Obs<-cbind(auc2L,auc2Lmax,auc2Lmin)
Auc2Obs<- as.data.frame(merge(Auc2Obs,cat))
names(Auc2Obs)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")


NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Cp),
    auc = auc_partial(time,Cp),
    auc_inf= auc_inf(time,Cp),
    tmax  = (time[which(Cp==Cmax)]
    )
  ) 

out1h<- filter(out_0, time<=1)
out2h<- filter(out_0, time>=1 & time<=2)

NCA_p1h1 <- out1h %>% group_by(x,ID)%>%
  summarise(
    auc1h=auc_partial(time,Cp)
  )
NCA_p1h2 <- out2h %>% group_by(x,ID)%>%
  summarise(
    auc2h=auc_partial(time,Cp)
  )

NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmed=median(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmed =median(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  tmax=   median(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

NCA_1h<-NCA_p1h1%>%group_by(x)%>%
  summarise(
    auc1hmed=median(auc1h),
    auc1hmin=min(auc1h),
    auc1hmax=max(auc1h)
  )

NCA_2h<-NCA_p1h2%>%group_by(x)%>%
  summarise(
    auc2hmed=median(auc2h),
    auc2hmin=min(auc2h),
    auc2hmax=max(auc2h)
  )


cmaxP<-mean(NCA_0$Cmaxmed)  #ng/mL = mg/L
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmed) * 60  #convert ng*hr/mL to mg*min/hr
minauc<-mean(NCA_0$minauc) * 60
maxauc<-mean(NCA_0$maxauc) * 60

tmaxP<- mean(NCA_0$tmax)*60  #convert to min
mintmax<- mean(NCA_0$tmax) * 60 
maxtmax<- mean(NCA_0$tmax) * 60

auc1hp<- mean(NCA_1h$auc1hmed) * 60 
auc1min<- mean(NCA_1h$auc1hmin) * 60 
auc1max <- mean(NCA_1h$auc1hmax) * 60



auc2hp<- mean(NCA_2h$auc2hmed) * 60 
auc2min<- mean(NCA_2h$auc2hmin) * 60 
auc2max <- mean(NCA_2h$auc2hmax) * 60

cmax_FE<- signif(ifelse(cmaxL>cmaxP,cmaxL/cmaxP,cmaxP/cmaxL), digits=4)
names(cmax_FE)<- "Cmax"


auc_FE<- signif(ifelse(aucP>aucL,aucP/aucL,aucL/aucP), digits=4) 
names(auc_FE)<- "AUC"


auc1_FE<- signif(ifelse(auc1L>auc1hp,auc1L/auc1hp,auc1hp/auc1L), digits=4)
names(auc1_FE)<-"AUC1hour"

auc2_FE<- signif(ifelse(auc2hp>auc2L,auc2hp/auc2L,auc2L/auc2hp), digits=4)
names(auc2_FE)<- "AUC2hour"

tmax_FE<- signif(ifelse(tmaxL>tmaxP,tmaxL/tmaxP,tmaxP/tmaxL), digits=4)
names(tmax_FE)<-"Tmax"


#datasets for graphs of NCA 
pred<-"Pred"

CmaxObs #
CmaxPred<- cbind(cmaxP,cmaxmax,cmaxmin)
CmaxPred<-merge(CmaxPred,pred)
names(CmaxPred)<- c("Median Cmax", "Upper Cmax", "Lower Cmax", "type")

CmaxG<- rbind(CmaxObs,CmaxPred)
CmaxG


t_cmax<- ggplot(CmaxG, aes(x=type, y=cmaxP, color=type))+
  geom_pointrange(aes(ymin=cmaxmin, ymax=cmaxmax))+
  labs(title="Cmax", x="Source", y="Cmax (ng/mL)")+
  annotate("text",x=1.5,y=25,label="Mean Fold-Error")+annotate("text", x=1.5, y=20, label=cmax_FE)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/Cmax1500mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_cmax)
dev.off()


AucObs
AucPred<- cbind(aucP,maxauc,minauc)
AucPred<- merge(AucPred,pred)
names(AucPred)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


AucG <- rbind(AucObs,AucPred)  
AucG


t_auc<- ggplot(AucG, aes(x=type, y=aucP, color=type))+
  geom_pointrange(aes(ymin=minauc, ymax=maxauc))+
  labs(title="AUC0-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=3900,label="Mean Fold-Error")+annotate("text", x=1.5, y=3700, label=auc_FE)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/AUC1500mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc)
dev.off()





Auc1Obs
Auc1Pred<- cbind(auc1hp,auc1max,auc1min)
Auc1Pred<- merge(Auc1Pred,pred)
names(Auc1Pred) <- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")

Auc1hG<- rbind(Auc1Obs,Auc1Pred)
Auc1hG

t_auc1h <- ggplot(Auc1hG, aes(x=type, y=auc1hp, color=type))+
  geom_pointrange(aes(ymin=auc1min, ymax=auc1max))+
  labs(title="AUC0-1", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=1500,label="Mean Fold-Error")+annotate("text", x=1.5, y=1400, label=auc1_FE)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/AUC011500mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc1h)
dev.off()


Auc2Obs
Auc2Pred<- cbind(auc2hp, auc2max, auc2min)
Auc2Pred<- merge(Auc2Pred,pred)
names(Auc2Pred)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")

Auc2hG<- rbind(Auc2Obs,Auc2Pred)
Auc2hG


t_auc2h <- ggplot(Auc2hG, aes(x=type, y=auc2hp, color=type))+
  geom_pointrange(aes(ymin=auc2min, ymax=auc2max))+
  labs(title="AUC2-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=1000,label="Mean Fold-Error")+annotate("text", x=1.5, y=800, label=auc2_FE)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/AUC011500mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc2h)
dev.off()




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/NCA1500mgORAL.png",width= 9 * ppi, height=4*ppi, res=ppi)
grid.arrange(t_cmax,t_auc,t_auc1h,t_auc2h)

dev.off()


all<-cbind(cmax_FE,auc_FE,auc1_FE,auc2_FE,tmax_FE)


LFE<-log10(all)

mean_LFE<-mean(LFE)







Exp <- function(x) {
  a <- 10^x
}

MFE_L<-Exp(mean_LFE)
all
MFE_L


write.csv(MFE_L,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/ORALMFELIV1500mg.csv",row.names=TRUE)
write.csv(all,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/ALLFEORAL1500mg.csv",row.names=TRUE)



#CSF analysis
###CSF
#Literature values
cmaxLc<-14.8
cmaxLcmin<-5.2
cmaxLcmax<-18.1


CmaxObsc<-cbind(cmaxLc,cmaxLcmax,cmaxLcmin)

CmaxObsc<- as.data.frame(merge(CmaxObsc,cat))
names(CmaxObsc)<- c("Median Cmax", "Upper Cmax", "Lower Cmax","type")


auc1Lc<-7
auc1Lcmin<-0
auc1Lcmax<-371


Auc1Obsc<-cbind(auc1Lc,auc1Lcmax,auc1Lcmin)
Auc1Obsc<- as.data.frame(merge(Auc1Obsc,cat))
names(Auc1Obsc)<- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")


auc2Lc<-203
auc2Lcmin<- 33
auc2Lcmax<- 872


Auc2Obsc<-cbind(auc2Lc,auc2Lcmax,auc2Lcmin)
Auc2Obsc<- as.data.frame(merge(Auc2Obsc,cat))
names(Auc2Obsc)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")




aucLc<- 2715
aucLcmin <- 962
aucLcmax <- 4761



AucObsc<-cbind(aucLc,aucLcmax,auc2Lcmin)
AucObsc<- as.data.frame(merge(AucObsc,cat))
names(AucObsc)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")



tmaxLc<-240
tmaxLcmin<-180
tmaxLcmax<-360
TmaxObsc<- cbind(tmaxLc,tmaxLcmin,tmaxLcmax)
TmaxObsc<- as.data.frame(merge(TmaxObsc,cat))



#Predicted values - CSF

NCA_c1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Csa),
    auc = auc_partial(time,Csa),
    auc_inf= auc_inf(time,Csa),
    tmax  = (time[which(Csa==Cmax)]
    )
  ) 

out1h<- filter(out_0, time<=1)
out2h<- filter(out_0, time>=1 & time<=2)

NCA_c1h1 <- out1h %>% group_by(x,ID)%>%
  summarise(
    auc1h=auc_partial(time,Csa)
  )
NCA_c1h2 <- out2h %>% group_by(x,ID)%>%
  summarise(
    auc2h=auc_partial(time,Csa)
  )

NCA_0c<- as.data.frame(NCA_c1%>%
                         group_by(x)%>%
                         summarise(Cmaxmed=median(Cmax),
                                   cmaxmin=min(Cmax),
                                   cmaxmax=max(Cmax),
                                   aucmed =median(auc),
                                   minauc=min(auc),
                                   maxauc=max(auc),
                                   tmaxmed=   median(tmax),
                                   tmaxmin=  min(tmax),
                                   tmaxmax=  max(tmax))
)

NCA_1hc<-NCA_c1h1%>%group_by(x)%>%
  summarise(
    auc1hmed=median(auc1h),
    auc1hmin=min(auc1h),
    auc1hmax=max(auc1h)
  )

NCA_2hc<-NCA_c1h2%>%group_by(x)%>%
  summarise(
    auc2hmed=median(auc2h),
    auc2hmin=min(auc2h),
    auc2hmax=max(auc2h)
  )


cmaxPc<-mean(NCA_0c$Cmaxmed)  #ng/mL = mg/L
aucPc<- mean(NCA_0c$aucmed) * 60  #convert ng*hr/mL to mg*min/hr
tmaxPc<- mean(NCA_0c$tmaxmed)*60  #convert to min

auc1hc<- mean(NCA_1hc$auc1hmed) * 60 
auc2hc<- mean(NCA_2hc$auc2hmed) * 60 

cmax_FEc<-signif(ifelse(cmaxLc>cmaxPc,cmaxLc/cmaxPc,cmaxPc/cmaxLc),digits=4)
names(cmax_FEc)<- "Cmax"


auc_FEc<- signif(ifelse(aucLc>aucPc,aucLc/aucPc,aucPc/aucLc),digits=4)
names(auc_FEc)<- "AUC"


auc1_FEc<- signif(ifelse(auc1hc>auc1Lc,auc1hc/auc1Lc,auc1Lc/auc1hc),digits=4)
names(auc1_FEc)<-"AUC1hour"

auc2_FEc<- signif(ifelse(auc2Lc>auc2hc,auc2Lc/auc2hc,auc2hc/auc2Lc), digits=4) 
names(auc2_FEc)<- "AUC2hour"

tmax_FEc<- signif(ifelse(tmaxPc>tmaxLc,tmaxPc/tmaxLc,tmaxLc/tmaxPc), digits=4)
names(tmax_FEc)<-"Tmax"



#datasets for graphs of NCA 
pred<-"Pred"

CmaxObsc #
CmaxPredc<- cbind(cmaxPc,cmaxmaxc,cmaxminc)
CmaxPredc<-merge(CmaxPredc,pred)
names(CmaxPredc)<- c("Median Cmax", "Upper Cmax", "Lower Cmax", "type")

CmaxGc<- rbind(CmaxObsc,CmaxPredc)
CmaxGc


t_cmaxc<- ggplot(CmaxGc, aes(x=type, y=cmaxPc, color=type))+
  geom_pointrange(aes(ymin=cmaxminc, ymax=cmaxmaxc))+
  labs(title="Cmax", x="Source", y="Cmax (ng/mL)")+
  annotate("text",x=1.5,y=10,label="Mean Fold-Error")+annotate("text", x=1.5, y=8, label=cmax_FEc)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1500mgIVcmax.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_cmaxc)
dev.off()




AucObsc
AucPredc<- cbind(aucPc,maxaucc,minaucc)
AucPredc<- merge(AucPredc,pred)
names(AucPredc)<- c("Median AUC0-6", "Upper AUC0-6", "Lower AUC0-6","type")


AucGc <- rbind(AucObsc,AucPredc)  
AucGc


t_aucc<- ggplot(AucGc, aes(x=type, y=aucPc, color=type))+
  geom_pointrange(aes(ymin=minaucc, ymax=maxaucc))+
  labs(title="AUC0-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=4300,label="Mean Fold-Error")+annotate("text", x=1.5, y=4000, label=auc_FEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/1500mgoralauc.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_aucc)
dev.off()



Auc1Obsc
Auc1Predc<- cbind(auc1hpc,auc1maxc,auc1minc)
Auc1Predc<- merge(Auc1Predc,pred)
names(Auc1Predc) <- c("Median AUC1hour", "Upper AUC1hour", "Lower AUC1hour","type")

Auc1hGc<- rbind(Auc1Obsc,Auc1Predc)
Auc1hGc

t_auc1hc <- ggplot(Auc1hGc, aes(x=type, y=auc1hpc, color=type))+
  geom_pointrange(aes(ymin=auc1minc, ymax=auc1maxc))+
  labs(title="AUC0-1", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=300,label="Mean Fold-Error")+annotate("text", x=1.5, y=275, label=auc1_FEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1500mgoraltauc1h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc1hc)
dev.off()


Auc2Obsc
Auc2Predc<- cbind(auc2hpc, auc2maxc, auc2minc)
Auc2Predc<- merge(Auc2Predc,pred)
names(Auc2Predc)<- c("Median AUC2hour", "Upper AUC2hour", "Lower AUC2hour","type")

Auc2hGc<- rbind(Auc2Obsc,Auc2Predc)
Auc2hGc


t_auc2hc <- ggplot(Auc2hGc, aes(x=type, y=auc2hpc, color=type))+
  geom_pointrange(aes(ymin=auc2minc, ymax=auc2maxc))+
  labs(title="AUC2-6", x="Source", y="AUC (ng*min/mL)")+
  annotate("text",x=1.5,y=400,label="Mean Fold-Error")+annotate("text", x=1.5, y=375, label=auc2_FEc)


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1500mgoraltauc2h.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(t_auc2hc)
dev.off()


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSF1500mgoralNCA.png",width= 9 * ppi, height=4*ppi, res=ppi)

grid.arrange(t_cmaxc,t_aucc,t_auc1hc,t_auc2hc)
dev.off()



allc<-cbind(cmax_FEc,auc_FEc,auc1_FEc,auc2_FEc,tmax_FEc)


LFEc<-log10(allc)

mean_LFEc<-mean(LFEc)

Exp <- function(x) {
  a <- 10^x
}

MFE_Lc<-Exp(mean_LFEc)
allc

MFE_Lc








write.csv(MFE_Lc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSFMFELoral1500mg.csv",row.names=TRUE)
write.csv(allc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Langford/CSFALLFEoral1500mg.csv",row.names=TRUE)








######################################################################
#                                                                    #
#                      Singla                                        #
#                      CSF/Plasma                                    #
#                                                                    #
######################################################################


singIVp<- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Singla IV plasma.csv")




singIVc<- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/SinglaIVCSF.csv")
singOrc<- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/SinglaOralCSF.csv")


singIVp2<- singIVp%>%
                      group_by(time)%>%
                                        summarise( 
                                          concavg=mean(conc),
                                          concmed=median(conc),
                                          concsd =sqrt(var(conc)))


###########IV Singla 1000 mg over 15 minutes
 

bwsing<- rnorm(10000,80000,21000)
bwsing<- bwsing[bwsing>=50000 & bwsing<=110000]



Renal <-rnorm(10000, 1400, 180 )  #renal clearance from Gaohua 2016 with ~30% CV
Renal <-Renal[Renal > 500]

Hep   <-rnorm(10000,44000, 6840) #hep clearance from Gaouhia 2016 with ~30% CV 
Hep   <-Hep[Hep>500  ]


Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.25 & Kpf<1.75]
 

out_0<-data.frame(matrix(1:17, ncol=17))

names(out_0)<-c("ID","time","x","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva")

stime<-c(0,.25,.5,.75,1,1.5,2,3,4,6)


c <-detectCores(all.tests=FALSE, logical =TRUE) #count cores on the CPU
cores <- c-1                                    #use all cores minus 1


options(future.fork.endable = TRUE)  #### from package future.apply, using forked approach
plan(multiprocess, workers=cores)    #### assign num

repeat{
  for(i in 0:1000){ 
    it<-i+1
    
    names(it)<-"x"
    idataBW <- expand.idata(BW=(sample(bwsing,6,replace=FALSE)))
    idataCLr<- expand.idata(CLrenal=(sample(Renal,6, replace = FALSE)))
    idataCLh<- expand.idata(CLhep = (sample(Hep,6,replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,7,replace=FALSE)))
    
    met<- merge(idataCLh,idataCLr)
    met2<- merge(met,idataBW)
    
    idata1<- merge(it,met2)
    idata2<- merge(idata1, ac_d_i)
    idata3<- merge(idata2,idataKpf)
    
    
    sp <- split(idata3,idata3$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod_ac<-loadso(mread("acetaminophen_091121b"))
      as.data.frame(
        mod_ac %>%
          ev(amt=1e6)%>%
          idata_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva) %>%
          carry_out(x)%>%
          mrgsim(maxsteps=100000, tgrid=stime))})%>%
      bind_rows()
    

    
    
    
    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    if (it==1000 ) {return(out_0);
      break}
    
  }
  
}
out_0<-out_0[-c(1),]



#MFE by point
mfetimes<- c(.25,.5,.75,1,1.5,2,3,4,6)
#Plasma

outmfe<- filter(out_0,time %in% mfetimes )

outmfe2<- merge(singIVp2,outmfe, by="time")

outmfe3<- outmfe2 %>% 
  group_by(x) %>% 
  mutate(
    FE=ifelse(concavg>Cp,concavg/Cp,Cp/concavg),
    LFE=log10(FE))%>%
  summarise(
    MFE=mean(LFE)
  )

avgMFE<- mean(outmfe3$MFE)
Exp <- function(x) {
  a <- 10^x
}



SingVLFEp<- signif(Exp(avgMFE), digits=4)
SingVLFEp

out_p<-out_0%>% group_by(time)%>%
  summarise(
    mean=mean(Cp),
    med=median(Cp),
     
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95)
  ) 

 


#CSF



outmfec<- filter(out_0,time %in% mfetimes )


outmfe2c<- merge(singIVc,outmfec, by="time")

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



SingIVLFEc<- signif(Exp(avgMFEc), digits=4)
SingIVLFEc


###PLOTS 

out_p<- out_0%>% group_by(time)%>%
  summarise(
    conc=mean(Cp),
    neduab=median(Cp),
    sd  = sqrt(var(Cp)),
    up  = max(Cp),
    lo  = min(Cp),
    iqr = IQR(Cp),
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95)
  ) %>%
  mutate(up2 =conc+iqr,
         lo2 =conc-iqr)



p1<-ggplot() +
  geom_line(data=out_p, aes(time,Cp05), color="red")+
  geom_line(data=out_p, aes(time,Cp50), color="red")+
  geom_line(data=out_p, aes(time,Cp95), color="red")+
  geom_line(data=out_p,aes(time,conc), color="black") +
  geom_point(data=singIVp2,aes(time,concavg, color="Observed"))

p1

p<-p+labs(colour="Median", x="Time (hrs)", y="Concentration (mg/L)",
          title="Singla IV Plasma Comparison", subtitle="1000 mg")+
  annotate("text",x=3,y=50,label="Mean LOG Fold-Error")+annotate("text", x=3, y=40, label=SingVLFEp)



p<-ggplot() +
  geom_ribbon(data=out_p, aes(time,ymin=lo, ymax=up), color="blue",alpha=0.5) +
  geom_line(data=out_p,aes(time,conc), color="black") +
  geom_point(data=singIVp2,aes(time,concavg, color="Observed"))

p<-p+labs(colour="Median", x="Time (hrs)", y="Concentration (mg/L)",
          title="Singla IV Plasma Comparison", subtitle="1000 mg")+
  annotate("text",x=3,y=50,label="Mean LOG Fold-Error")+annotate("text", x=3, y=40, label=SingVLFEp)
ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaIV1000mg.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p1)
dev.off()



out_c<- out_0%>% group_by(time)%>%
  summarise(
    conc=median(Csa),
    mean=mean(Csa),
    sd  = sqrt(var(Csa)),
    up  = max(Csa),
    lo  = min(Csa),
    Cs05=quantile(Csa, probs=0.05),
    Cs50=quantile(Csa,probs=0.50),
    Cs95=quantile(Csa,probs=0.95)
  )

c<-ggplot() +
  geom_ribbon(data=out_c, aes(time,ymin=Cs05, ymax=Cs95), color="blue",alpha=0.5) +
  geom_line(data=out_c,aes(time,Cs50)) + 
  geom_point(data=LangIVC,aes(time,conc, color="Obs"))

c<-c+labs(colour="Observed Median", x="Time (hrs)", y="Concentration (mg/L)",
          title="Singla IV CSF Comparison", subtitle="1000 mg")+
  annotate("text",x=0.5,y=9,label="Mean Fold-Error")+annotate("text", x=0.5, y=6, label=SingIVLFEc)

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/CSFIV1000mg.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(c)
dev.off()

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaCSFIV1000mg.png",width= 9 * ppi, height=4*ppi, res=ppi)
grid.arrange(p,c)
dev.off()



###Plasma - IV
cmaxS<-21.6
cmaxSsd<- 17.9

aucS06<- 42.5
aucS06sd<-16.5

aucSinf<- 50.0
aucSinfsd<- 18.7

tmaxS<- 0.25




NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Cp),
    auc = auc_partial(time,Cp),
    auc_inf= auc_inf(time,Cp),
    tmax  = (time[which(Cp==Cmax)]
    )
  ) 



NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmean=mean(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmean =mean(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  aucImean=mean(auc_inf),
                                  aucImax =max(auc_inf),
                                  aucImin =min(auc_inf),
                                  tmax=   mean(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

cmaxP<-mean(NCA_0$Cmaxmean)  #ng/mL = mg/L
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmean)  
minauc<-mean(NCA_0$minauc)  
maxauc<-mean(NCA_0$maxauc)  

tmaxP<- mean(NCA_0$tmax) 
mintmax<- mean(NCA_0$tmax)  
maxtmax<- mean(NCA_0$tmax)  

aucIP<-mean(NCA_0$aucImean)
minaucI<-mean(NCA_0$aucImin)
maxaucI<-mean(NCA_0$aucImax)


cmax_FE<- signif(ifelse(cmaxS>cmaxP,cmaxS/cmaxP,cmaxP/cmaxS), digits=4)
names(cmax_FE)<- "Cmax"


auc_FE<- signif(ifelse(aucP>aucS,aucP/aucS,aucS/aucP), digits=4) 
names(auc_FE)<- "AUC"

aucI_FE<- signif(ifelse(aucIP>aucSinf, aucIP/aucSinf,aucSinf/aucIP), digits=4)
names(aucI_FE)<-"AUCinf"

tmax_FE<- signif(ifelse(tmaxS>tmaxP,tmaxS/tmaxP,tmaxP/tmaxS), digits=4)
names(tmax_FE)<-"Tmax"




all<-cbind(cmax_FEc,auc_FEc,auc1_FEc,auc2_FEc,tmax_FEc)


LFE<-log10(all)

mean_LFE<-mean(LFE)

write.csv(MFE_L,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaLFELIV1000mg.csv",row.names=TRUE)
write.csv(all,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaLFEIV1000mg.csv",row.names=TRUE)

###CSF - IV
cmaxSc<-5.94
cmaxScsd<- 18.4

aucS06c<- 24.9
aucS06csd<-17.4


tmaxSc<- 2
tmaxScmin<-1
tmaxScmax<-4




NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Csa),
    auc = auc_partial(time,Csa),
    tmax  = (time[which(Csa==Cmax)]
    )
  ) 



NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmean=mean(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmean =mean(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  tmax=   mean(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

cmaxP<-mean(NCA_0$Cmaxmean)  #ng/mL = mg/L
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmean)  
minauc<-mean(NCA_0$minauc)  
maxauc<-mean(NCA_0$maxauc)  

tmaxP<- mean(NCA_0$tmax) 
mintmax<- mean(NCA_0$tmax)  
maxtmax<- mean(NCA_0$tmax)  

 


cmax_FEc<- signif(ifelse(cmaxSc>cmaxP,cmaxSc/cmaxP,cmaxP/cmaxSc), digits=4)
names(cmax_FEc)<- "Cmax"


auc_FEc<- signif(ifelse(aucP>auc06Sc,aucP/auc06Sc,auc06Sc/aucP), digits=4) 
names(auc_FEc)<- "AUC"

 

tmax_FEc<- signif(ifelse(tmaxSc>tmaxP,tmaxSc/tmaxP,tmaxP/tmaxSc), digits=4)
names(tmax_FEc)<-"Tmax"




allc<-cbind(cmax_FEc,auc_FEc,tmax_FEc)


LFEc<-log10(allc)

mean_LFEc<-mean(LFEc)


######################################################################
#                                                                    #
#                      Singla Oral                                   #
#                      CSF/Plasma                                    #
#                                                                    #
######################################################################
singOrp<- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/Singla Oral plasma.csv")
singOrc<- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Acetaminophen/Data Files/SinglaOralCSF.csv")



singOrp2<- singOrp%>%
  group_by(time)%>%
  summarise( 
    concavg=mean(conc),
    concmed=median(conc),
    concsd =sqrt(var(conc)))



bwsing<- rnorm(10000,80000,21000)
bwsing<- bwsing[bwsing>=50000 & bwsing<=110000]



Renal <-rnorm(10000, 1400, 180 )  #renal clearance from Gaohua 2016 with ~30% CV
Renal <-Renal[Renal > 500]

Hep   <-rnorm(10000,44000, 6840) #hep clearance from Gaouhia 2016 with ~30% CV 
Hep   <-Hep[Hep>500  ]


Abs   <- rnorm(10000,2,0.5) #Calculated from abs half life Divoll 1982
Abs <-Abs[Abs>0.5 & Abs<3.5]

Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.5 & Kpf<1.25]


out_0<-data.frame(matrix(1:17, ncol=17))

names(out_0)<-c("ID","time","x","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva")

stime<-c(0,.25,.5,.75,1,1.5,2,3,4,6)


c <-detectCores(all.tests=FALSE, logical =TRUE) #count cores on the CPU
cores <- c-1                                    #use all cores minus 1


options(future.fork.endable = TRUE)  #### from package future.apply, using forked approach
plan(multiprocess, workers=cores)    #### assign num

repeat{
  for(i in 0:1000){ 
    it<-i+1
    
    
    names(it)<-"x"
    idataCLr<- expand.idata(CLrenal=(sample(Renal,7, replace = FALSE)))
    idataCLh<- expand.idata(CLhep = (sample(Hep,7,replace=FALSE)))
    idataBW <- expand.idata(BW = (sample(bwlan, 7, replace=FALSE)))
    idataAbs <- expand.idata(Abs=(sample(Abs,7, replace=FALSE)))
    idataKpf<- expand.idata(Kpf=(sample(Kpf,7,replace=FALSE)))
    
    met<- merge(idataCLh,idataCLr)
    met2<- merge(met,idataBW)
    
    idata1<- merge(it,met2)
    idata2<- merge(idata1,idataAbs)
    idata3<- merge(idata2, ac_d_o)
    idata4<- merge(idata3,idataKpf)
    
    
    sp <- split(idata4,idata4$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod_ac<-loadso(mread("acetaminophen_091121b"))
      as.data.frame(
        mod_ac %>%
          ev(amt=1e6)%>%
          idata_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva) %>%
          carry_out(x)%>%
          mrgsim(maxsteps=100000, tgrid=stime))})%>%
      bind_rows()
    
    
    
    
    
    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    if (it==1000 ) {return(out_0);
      break}
    
  }
  
}
out_0<-out_0[-c(1),]



#MFE by point
mfetimes<- c(.25,.5,.75,1,1.5,2,3,4,6)
#Plasma

outmfe<- filter(out_0,time %in% mfetimes )

outmfe2<- merge(singOrp2,outmfe, by="time")

outmfe3<- outmfe2 %>% 
  group_by(x) %>% 
  mutate(
    FE=ifelse(concavg>Cp,concavg/Cp,Cp/concavg),
    LFE=log10(FE))%>%
  summarise(
    MFE=mean(LFE)
  )

avgMFE<- mean(outmfe3$MFE)
Exp <- function(x) {
  a <- 10^x
}



SingORLFEp<- signif(Exp(avgMFE), digits=4)
SingORLFEp

 



#CSF



outmfec<- filter(out_0,time %in% mfetimes )


outmfe2c<- merge(singOrc,outmfec, by="time")

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



SingORLFEc<- signif(Exp(avgMFEc), digits=4)
SingORLFEc


###PLOTS 

out_p<- out_0%>% group_by(time)%>%
  summarise(
    conc=mean(Cp),
    neduab=median(Cp),
    sd  = sqrt(var(Cp)),
    up  = max(Cp),
    lo  = min(Cp),
    iqr = IQR(Cp),
    Cp05=quantile(Cp, probs=0.05),
    Cp50=quantile(Cp, probs=0.50),
    Cp95=quantile(Cp, probs=0.95)
  ) %>%
  mutate(up2 =conc+iqr,
         lo2 =conc-iqr)



p1<-ggplot() +
  geom_line(data=out_p, aes(time,Cp05), color="red")+
  geom_line(data=out_p, aes(time,Cp50), color="red")+
  geom_line(data=out_p, aes(time,Cp95), color="red")+
  geom_line(data=out_p,aes(time,conc), color="black") +
  geom_point(data=singIVp2,aes(time,concavg, color="Observed"))

p1

p<-p+labs(colour="Median", x="Time (hrs)", y="Concentration (mg/L)",
          title="Singla Oral Comparison", subtitle="1000 mg")+
  annotate("text",x=3,y=50,label="Mean LOG Fold-Error")+annotate("text", x=3, y=40, label=SingORLFEp)



p<-ggplot() +
  geom_ribbon(data=out_p, aes(time,ymin=lo, ymax=up), color="blue",alpha=0.5) +
  geom_line(data=out_p,aes(time,conc), color="black") +
  geom_point(data=singIVp2,aes(time,concavg, color="Observed"))

p<-p+labs(colour="Median", x="Time (hrs)", y="Concentration (ug/mL)",
          title="Singla Oral Plasma Comparison", subtitle="1000 mg")+
  annotate("text",x=3,y=50,label="Mean LOG Fold-Error")+annotate("text", x=3, y=40, label=SingORLFEp)
p
ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaOral1000mg.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(p1)
dev.off()



out_c<- out_0%>% group_by(time)%>%
  summarise(
    conc=median(Csa),
    mean=mean(Csa),
    sd  = sqrt(var(Csa)),
    up  = max(Csa),
    lo  = min(Csa),
    Cs05=quantile(Csa,probs=0.05),
    Cs50=quantile(Csa,probs=0.50),
    Cs95=quantile(Csa,probs=0.95)
  ) 


c<-ggplot() +
  geom_ribbon(data=out_c, aes(time,ymin=Cs05, ymax=Cs95), color="blue",alpha=0.5) +
  geom_line(data=out_c,aes(time,Cs50)) + 
  geom_point(data=singOrc,aes(time,conc, color="Obs"))

c<-c+labs(colour="Observed Median", x="Time (hrs)", y="Concentration (ug/mL)",
          title="Singla Oral CSF Comparison", subtitle="1000 mg")#+
 # annotate("text",x=0.5,y=9,label="Mean Fold-Error")+annotate("text", x=0.5, y=6, label=SingORLFEc)
c+scale_y_continuous(trans="log10")

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/CSFOral1000mg.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(c)
dev.off()

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaCSFOral1000mg.png",width= 9 * ppi, height=4*ppi, res=ppi)
grid.arrange(p,c)
dev.off()



###Plasma - Oral
cmaxS<-12.3
cmaxSsd<- 45.2

aucS06<- 29.4
aucS06sd<-52.3

aucSinf<- 44.4
aucSinfsd<- 35.4

tmaxS<- 1
tmaxSmin<- 0.5
tmaxSmax<- 2.0




NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Cp),
    auc = auc_partial(time,Cp),
    auc_inf= auc_inf(time,Cp),
    tmax  = (time[which(Cp==Cmax)]
    )
  ) 



NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmean=mean(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmean =mean(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  aucImean=mean(auc_inf),
                                  aucImax =max(auc_inf),
                                  aucImin =min(auc_inf),
                                  tmax=   mean(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

cmaxP<-mean(NCA_0$Cmaxmean)  
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmean)  
minauc<-mean(NCA_0$minauc)  
maxauc<-mean(NCA_0$maxauc)  

tmaxP<- mean(NCA_0$tmax) 
mintmax<- mean(NCA_0$tmax)  
maxtmax<- mean(NCA_0$tmax)  

aucIP<-mean(NCA_0$aucImean)
minaucI<-mean(NCA_0$aucImin)
maxaucI<-mean(NCA_0$aucImax)


cmax_FE<- signif(ifelse(cmaxS>cmaxP,cmaxS/cmaxP,cmaxP/cmaxS), digits=4)
names(cmax_FE)<- "Cmax"


auc_FE<- signif(ifelse(aucP>aucS,aucP/aucS,aucS/aucP), digits=4) 
names(auc_FE)<- "AUC"

aucI_FE<- signif(ifelse(aucIP>aucSinf, aucIP/aucSinf,aucSinf/aucIP), digits=4)
names(aucI_FE)<-"AUCinf"

tmax_FE<- signif(ifelse(tmaxS>tmaxP,tmaxS/tmaxP,tmaxP/tmaxS), digits=4)
names(tmax_FE)<-"Tmax"




all<-cbind(cmax_FEc,auc_FEc,auc1_FEc,auc2_FEc,tmax_FEc)


LFE<-log10(all)

mean_LFE<-mean(LFE)

write.csv(MFE_L,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaLFELOral1000mg.csv",row.names=TRUE)
write.csv(all,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/PlasmaLFEOral1000mg.csv",row.names=TRUE)

###CSF - Oral
cmaxSc<-3.72
cmaxScsd<- 39.1

aucS06c<- 14.2
aucS06csd<-52.1


tmaxSc<- 4
tmaxScmin<-0.75
tmaxScmax<-6




NCA_p1<- out_0%>% group_by(x,ID)%>%
  summarise(
    Cmax = max(Csa),
    auc = auc_partial(time,Csa),
    tmax  = (time[which(Csa==Cmax)]
    )
  ) 



NCA_0<- as.data.frame(NCA_p1%>%
                        group_by(x)%>%
                        summarise(Cmaxmean=mean(Cmax),
                                  cmaxmin=min(Cmax),
                                  cmaxmax=max(Cmax),
                                  aucmean =mean(auc),
                                  minauc=min(auc),
                                  maxauc=max(auc),
                                  tmax=   mean(tmax),
                                  mintmax=min(tmax),
                                  maxtmax=max(tmax))
)

cmaxP<-mean(NCA_0$Cmaxmean)  #ng/mL = mg/L
cmaxmin<-mean(NCA_0$cmaxmin)
cmaxmax<-mean(NCA_0$cmaxmax)


aucP<- mean(NCA_0$aucmean)  
minauc<-mean(NCA_0$minauc)  
maxauc<-mean(NCA_0$maxauc)  

tmaxP<- mean(NCA_0$tmax) 
mintmax<- mean(NCA_0$tmax)  
maxtmax<- mean(NCA_0$tmax)  




cmax_FEc<- signif(ifelse(cmaxSc>cmaxP,cmaxSc/cmaxP,cmaxP/cmaxSc), digits=4)
names(cmax_FEc)<- "Cmax"


auc_FEc<- signif(ifelse(aucP>auc06Sc,aucP/auc06Sc,auc06Sc/aucP), digits=4) 
names(auc_FEc)<- "AUC"



tmax_FEc<- signif(ifelse(tmaxSc>tmaxP,tmaxSc/tmaxP,tmaxP/tmaxSc), digits=4)
names(tmax_FEc)<-"Tmax"




allc<-cbind(cmax_FEc,auc_FEc,tmax_FEc)


LFEc<-log10(allc)

mean_LFEc<-mean(LFEc)

write.csv(mean_LFEc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/CSFLFELOral1000mg.csv",row.names=TRUE)
write.csv(allc,"/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Acetaminophen/Singla/CSFLFEOral1000mg.csv",row.names=TRUE)


