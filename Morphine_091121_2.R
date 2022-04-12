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



mod_mo<-mread("Morphine_091121b")


##############################################################
#       Ederoth and Bouw Combined Data                       #
#                                                            #
#                                                            #
##############################################################

EB_PL <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Morphine/ALL_EB_Plasma.csv")
EB_EC_B <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Data Extractions/Morphine/ALL_BETTER_ECF.csv")



EB_TIME<-c(0,	
5,	2.5,
10,	7.5,
15,	12.5,
20,	17.5,
25,	22.5,
30,	27.5,
40,	35,
50,	45,
60,	55,
70,	65,
80,	75,
90,	85,
100,	95,
110,	105,
120,	115,
130,	125,
140,	135,
150,	145,
160,	155,
170,	165,
180,	175)
EB_TIME2<- EB_TIME/60

#Ederoth/Bouw infusion parameters - 10 mg (1e7 ng) in 10 min ()
names_d <- as.list(c("time","oral", "IV","IRa","amt","evid","cmt"))
b <- as.list(c(0,0,1,6,1e7,1,1))
names(b)<-names_d



set.seed(72000)
BW<-rnorm(10000, 72000,20000)
BW<- BW[BW>45000 & BW<120000]

set.seed(300)
CLr <- rnorm(10000,300,100)
CLr<- CLr[CLr>100  ]

set.seed(1)
Kpf<- rnorm(10000,1,0.5)
Kpf<-Kpf[Kpf>0.25 & Kpf<1.75]


out_0<-data.frame(matrix(1:17, ncol=17))

names(out_0)<-c("ID","time","x","Cp","Mball", "Csa", "Cec1", "Cec2", "Cec3", "Cec4", "Icf1", 
                "Icf2", "Icf3", "Icf4","Ceca","Bra","Cmva")

cores <- 6                                   #Use 6 cores

options(future.fork.endable = TRUE)  #### from package future.apply, using forked approach
plan(multiprocess, workers=cores)    #### assign num



repeat{
  for(i in 0:1000){ 
    it<-i+1
    
    names(it)<-"x"
    
    idataCL<- expand.idata(CLrange=(sample(CLr,3,replace=FALSE)))
    idataBW<- expand.idata(BW=(sample(BW,3,replace=FALSE)))
    idataKpf<-expand.idata(Kpf=(sample(Kpf,3,replace=FALSE)))
    idata1<-merge(idataCL,idataBW, by="ID")
    
    data<- merge(idata1,b)
    data2<-merge(data,idataKpf)
    
    data3<-merge(it,data2)
    
 
    
    
    sp <- split(data3,data3$ID)#### population data set split into individual IDs to allow for parallel processing
    
    out<- future_lapply(sp, function(chunk) {  
      mod<-loadso(mod_mo)
      as.data.frame(
        mod_mo %>%
          data_set(chunk) %>%
          Req(Cp, Mball, Csa, Ceca,Cec1, Cec2, Cec3, Cec4, Icf1, Icf2, Icf3, Icf4, Bra, Cmva) %>%
          carry_out(x)%>%
          mrgsim(maxsteps=100000, tgrid=EB_TIME2))})%>%
      bind_rows()
    
    
    
    out<-as.data.frame(out)
    out_0<-rbind(out_0,out)
    if (it==1000) {
      return(out_0); break}
    
  }
  
}



out_0<-out_0[-c(1),]

sum<-out_0%>% group_by(x,time)%>%
  summarise(
 ConcP=mean(Cp),
 ConcE=mean(Cec1),
 Cp05 =quantile(Cp, probs=0.05),
 Cp50 =quantile(Cp, probs=0.50),
 Cp95 =quantile(Cp, probs=0.95),
 Cec05 =quantile(Cec1, probs=0.05),
 Cec50 =quantile(Cec1, probs=0.50),
 Cec95 =quantile(Cec1, probs=0.95)
  )%>% mutate(
    timehours=time
  )
sum2<-sum%>%
            group_by(timehours)%>%
  summarise(
    ConcPm=mean(ConcP),
    ConcEm=mean(ConcE),
    Cp05x=mean(Cp05),
    Cp50x=mean(Cp50),
    Cp95x=mean(Cp95),
    ConcP05=quantile(ConcP,probs = 0.05),
    ConcP50=quantile(ConcP,probs = 0.50),
    ConcP95=quantile(ConcP,probs = 0.95),
    ConEP05=quantile(ConcE,probs = 0.05),
    ConcEP50=quantile(ConcE,probs = 0.50),
    ConcEP95=quantile(ConcE,probs = 0.95),
    Cec05x=mean(Cec05),
    Cec50x=mean(Cec50),
    Cec95x=mean(Cec95),
  )

 


MFE_PL<- merge(sum2,EB_PL, by="timehours")
 

outmfe2<- MFE_PL %>% 
 
  mutate(
    FE=ifelse(conc>Cp50x,conc/Cp50x,Cp50x/conc),
    LFE=log10(FE))%>%
  summarise(
    MFE=mean(LFE)
  )

avgMFE<- mean(outmfe2$MFE)
Exp <- function(x) {  
  a <- 10^x
}



EB_FEp<- signif(Exp(avgMFE), digits=4)
EB_FEp


p<-ggplot() +
  geom_ribbon(data=sum2, aes(timehours,ymin=ConcP05, ymax=ConcP95), color="black",alpha=0.5) +
  geom_line(data=sum2,aes(timehours,ConcP50), color="black") +
  geom_point(data=EB_PL,aes(timehours,conc))

plin<-p+labs( x="Time (hrs)", y="Concentration (ng/mL)",
      title="Morphine IV 10 mg - Plasma", subtitle="10 minute infusion")+theme_bw()

plin


ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Morphine/EB/PEB_Plasma_Linear.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(plin)
dev.off()

pl<-p+labs( x="Time (hrs)", y="LOG Concentration (ng/mL)",
            title="Morphine IV 10 mg - Plasma", subtitle="10 minute infusion")+theme_bw()+
            scale_y_continuous(trans="log10")

pl

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Morphine/EB/PEB_Plasma_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(pl)
dev.off()

MFE_EC<- merge(sum2,EB_EC_B, by="timehours")


outmfe2c<- MFE_EC %>% 
  
  mutate(
    FE=ifelse(conc>Cec50x,conc/Cec50x,Cec50x/conc),
    LFE=log10(FE))%>%
  summarise(
    MFE=mean(LFE)
  )

avgMFEc<- mean(outmfe2c$MFE)
Exp <- function(x) {  
  a <- 10^x
}



EB_FEc<- signif(Exp(avgMFEc), digits=4)
EB_FEc


c<-ggplot() +
  geom_ribbon(data=sum2, aes(timehours,ymin=Cec05x, ymax=Cec95x), color="black",alpha=0.5) +
  geom_line(data=sum2,aes(timehours,Cec50x), color="black") +
  geom_point(data=EB_EC_B,aes(timehours,conc))

clin<-c+labs( x="Time (hrs)", y="Concentration (ng/mL)",
       title="Morphine IV 10 mg - ECF", subtitle="10 minute infusion")+theme_bw()




ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Morphine/EB/PEB_ECF_Linear.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(clin)
dev.off()

cl<-c+labs( x="Time (hrs)", y="LOG Concentration (ng/mL)",
            title="Morphine IV 10 mg - ECF", subtitle="10 minute infusion")+theme_bw()+
             scale_y_continuous(trans="log10")

cl

ppi<-300
png("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/Simulation Exports/Morphine/EB/PEB_ECF_LOG.png",width= 9 * ppi, height=4*ppi, res=ppi)
print(cl)
dev.off()






