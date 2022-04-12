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
###This code was taken from Utsey 2020 - it generates the partition coefficient based on Rogers/Rowland 2015-2016

calcKp_RR <- function(logP, pKa=0, fup, BP=1, type=1, dat){
  
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs", "Adipose", "Plasma"))  #df for all tissues except for adipose, RBCs, and plasma
  dat_ad <- dat %>% filter(tissue == "Adipose")  #df for adipose
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for aplasma
  
  
  pH_IW <- 7       #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  Ka <- 10^(-pKa)
  HCT <- 0.45 #hematocrit
  
  
  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT)
  
  X <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_IW-pKa),
              #3-monoprotic base
              10^(pKa-pH_IW),
              #4-diprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_IW)+10^(pKa[1]+pKa[2]-2*pH_IW), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_IW)+10^(pH_IW-pKa[1]),  
              #7-triprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2])+10^(3*pH_IW-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_IW)+10^(pKa[3]+pKa[2]-2*pH_IW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_IW)+10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IW-pKa[1])+10^(pKa[3]-pH_IW)+10^(pKa[2]+pKa[3]-2*pH_IW))       
  
  Y <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_P-pKa),
              #3-monoprotic base
              10^(pKa-pH_P), 
              #4-diprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_P)+10^(pKa[1]+pKa[2]-2*pH_P), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_P)+10^(pH_P-pKa[1]),  
              #7-triprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2])+10^(3*pH_P-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_P)+10^(pKa[3]+pka[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_P)+10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_P-pKa[1])+10^(pKa[3]-pH_P)+10^(pKa[2]+pKa[3]-2*pH_P))       
  
  Z <- switch(type,
              #1-neutral
              1,   
              #2-monoprotic acid
              1,
              #3-monoprotic base
              10^(pKa-pH_RBC), 
              #4-diprotic acid
              1,
              #5-diprotic base
              10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_RBC)+10^(pH_RBC-pKa[1]),  
              #7-triprotic acid
              1,  
              #8-triprotic base
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pka[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC)) 
  
  
  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
  
  
  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)
  
  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals
  
  
  # Multiply by fup to get Kp rather than Kpu
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_AP*dat_all$f_a_pl*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_AP*dat_ad$f_a_pl*X)/(1 + Y))*fup  #lipid
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$AR*X)/(1 + Y))*fup #lipid
  }else{  #neutrals
    Kp_all <- (dat_all$f_ew + ((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$f_n_l + (0.3*P + 0.7)*dat_all$f_n_pl))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))*fup  #non lipid
    Kp_ad <- (dat_ad$f_ew + ((1 + X)/(1 + Y))*dat_ad$f_iw + ((P_OW*dat_ad$f_n_l + (0.3*P_OW + 0.7)*dat_ad$f_n_pl))/(1 + Y) + (Ka_PR*dat_ad$LR*X)/(1 + Y))*fup  #lipid
  }
  
  
  nms_all <- dat_all$tissue %>% substr(1,2) %>% tolower()
  nms_all <- paste("Kp", nms_all, sep="")
  nms <- c("Kpad",nms_all)
  Kp <- as.list(c(Kp_ad,Kp_all))
  names(Kp) <- nms
  
  
  return(Kp)
}

#The below data was also taken from the Utsey 2020 publication. It is a unified database of human values used in the calculation of Kps

dat1 <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models/unifiedtissue.csv")

#Run the macro with the following inputs, in order
#Logp pka fup bp type #dat


#Acetaminophen
Kp_ac <- calcKp_RR(0.46,pKa=9.46,0.855,BP=1.09,type=1,dat1) #Saleh 2021

#Morphine

Kp_mo<- calcKp_RR(0.87,pKa=10.26,0.65,BP=1.02, type=3, dat1) #Ball 2012 

#Nilotinib


Kp_nl <-   calcKp_RR(4.5,pKa=c(5.4,3.9),0.016,BP=0.68,type=5,dat1) #Drugbank/Heimbach 


#Phenytoin

Kp_pt <-  calcKp_RR(2.47,pKa=9.47,0.1,0.61,type=1,dat1) #Saleh 2021/Polasek 2009




#Assembly of Kps into table 

#Acetaminophen

ac_kp <- as.data.frame(Kp_ac)
drug2<- as.data.frame("Acetaminophen")
ac_kp <- cbind(drug2,ac_kp)
colnames(ac_kp) <- c('Drug','Adipose','Bone','Brain','Heart','Kidney','Gut','Liver','Lung','Muscle','Skin','Spleen')


#Morphine



mo_kp <- as.data.frame(Kp_mo) 
druga <- as.data.frame ("Morphine")
mo_kp <- cbind(druga,mo_kp)
colnames(mo_kp) <- c('Drug','Adipose','Bone','Brain','Heart','Kidney','Gut','Liver','Lung','Muscle','Skin','Spleen')


#Nilotinib



nl_kp <- as.data.frame(Kp_nl)
drug3<- as.data.frame("Nilotinib")
nl_kp <- cbind(drug3,nl_kp)
colnames(nl_kp) <- c('Drug','Adipose','Bone','Brain','Heart','Kidney','Gut','Liver','Lung','Muscle','Skin','Spleen')





#Phenytoin



pt_kp <- as.data.frame(Kp_pt) 
drug1 <- as.data.frame ("Phenytoin")
pt_kp <- cbind(drug1,pt_kp)
colnames(pt_kp) <- c('Drug','Adipose','Bone','Brain','Heart','Kidney','Gut','Liver','Lung','Muscle','Skin','Spleen')





names_dat<-c("Drug","Log P","pKa1","pKa2","Fup","B:P","Charge")
 

ac_dat<-c("Acetaminophen",0.46,9.46,NA,0.855,1.09,"Neutral")



mo_dat<-c("Morphine",0.87,10.26,NA,0.65,1.02,"monoprotic base")




nl_dat<-c("Nilotinib",4.51,5.4,3.9,0.016,0.68,"diprotic base")



#Phenytoin

 
pt_dat<-c("Phenytoin",2.47,9.47,NA,0.1,0.61,"Neutral")


all_dat<-as.data.frame(rbind(ac_dat,mo_dat,nl_dat,pt_dat))

colnames(all_dat)<-names_dat
rownames(all_dat)<-c()



##TABLE FOR PUBLICATION

kp_wb <- rbind(ac_kp, mo_kp, nl_kp, pt_kp)
knitr::kable(all_dat,"pipe")
knitr::kable(kp_wb, "pipe",digits=3)




#CALCULATE Plasma/EC partition using RR derived equation

calcKp_EC <- function(logP, pKa=0, fup, BP=1, type=1, dat){
  
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs","Plasma","Adipose","Bone","Brain","Heart","Kidney","Gut","Liver",
                                           "Lung","Muscle","Skin","Spleen","RBCs","Plasma"))  #df for brain regions
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for aplasma
  
  
  pH_EW <- dat_all$pHew     #pH of intracellular tissue water
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  Ka <- 10^(-pKa)
  HCT <- 0.45 #hematocrit
  
  
  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT)
  
  X <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_EW-pKa),
              #3-monoprotic base
              10^(pKa-pH_EW),
              #4-diprotic acid
              10^(pH_EW-pKa[1])+10^(2*pH_EW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_EW)+10^(pKa[1]+pKa[2]-2*pH_EW), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_EW)+10^(pH_EW-pKa[1]),  
              #7-triprotic acid
              10^(pH_EW-pKa[1])+10^(2*pH_EW-pKa[1]-pKa[2])+10^(3*pH_EW-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_EW)+10^(pKa[3]+pKa[2]-2*pH_EW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_EW),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_EW)+10^(pH_EW-pKa[1])+10^(2*pH_EW-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_EW-pKa[1])+10^(pKa[3]-pH_EW)+10^(pKa[2]+pKa[3]-2*pH_EW))       
  
  Y <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_P-pKa),
              #3-monoprotic base
              10^(pKa-pH_P), 
              #4-diprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_P)+10^(pKa[1]+pKa[2]-2*pH_P), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_P)+10^(pH_P-pKa[1]),  
              #7-triprotic acid
              10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2])+10^(3*pH_P-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_P)+10^(pKa[3]+pka[2]-2*pH_P)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_P),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_P)+10^(pH_P-pKa[1])+10^(2*pH_P-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_P-pKa[1])+10^(pKa[3]-pH_P)+10^(pKa[2]+pKa[3]-2*pH_P))       
  
  Z <- switch(type,
              #1-neutral
              1,   
              #2-monoprotic acid
              1,
              #3-monoprotic base
              10^(pKa-pH_RBC), 
              #4-diprotic acid
              1,
              #5-diprotic base
              10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_RBC)+10^(pH_RBC-pKa[1]),  
              #7-triprotic acid
              1,  
              #8-triprotic base
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pka[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC)) 
  
  
  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
  
  
  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)
  
  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals
  
  
  
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (((1 + X)/(1 + Y))*dat_all$f_ew + ((P*dat_all$fnl_o + (0.3*P + 0.7)*dat_all$fnpl_o))/(1 + Y) + (Ka_AP*dat_all$AP_o*X)/(1 + Y))*fup  #non lipid
    
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (((1 + X)/(1 + Y))*dat_all$f_ew + ((P*dat_all$fnl_o + (0.3*P + 0.7)*dat_all$fnpl_o))/(1 + Y) + (Ka_PR*dat_all$AR*X)/(1 + Y))*fup  #non lipid
    
  }else{  #neutrals
    Kp_all <- (((1 + X)/(1 + Y))*dat_all$f_ew + ((P*dat_all$fnl_o + (0.3*P + 0.7)*dat_all$fnpl_o))/(1 + Y) + (Ka_PR*dat_all$LR*X)/(1 + Y))*fup  #non lipid
    
  }
  
  
  nms_all <- dat_all$tissue %>% substr(1,3) %>% tolower()
  nms_all <- paste("Kp_ec", nms_all, sep="")
  
  Kp_ec <- as.list(c(Kp_all))
  names(Kp_ec) <- nms_all
  
  
  return(Kp_ec)
}

#IC

calcKp_IC <- function(logP, pKa=0, fup,BP=1, type=1, dat){
  
  
  dat_all <- dat %>% filter(!tissue %in% c("RBCs","Plasma","Adipose","Bone","Brain","Heart","Kidney","Gut","Liver",
                                           "Lung","Muscle","Skin","Spleen","RBCs","Plasma"))  #df for brain regions
  dat_rbc <- dat %>% filter(tissue == "RBCs") #df for RBCs
  dat_plas <- dat %>% filter(tissue == "Plasma") #df for aplasma
  
  
  pH_EW <- dat_all$pHew     #pH of intracellular tissue water
  pH_IW <- dat_all$pHiw
  pH_P <- 7.4      #pH of plasma
  pH_RBC <- 7.22    #pH of blood cells
  P <- 10^(logP)   # octonal:water partition coeff
  logP_OW <- 1.115*logP - 1.35 #oil:water partition coeff
  P_OW <- 10^(logP_OW) 
  Ka <- 10^(-pKa)
  HCT <- 0.45 #hematocrit
  
  #[EW->IW, P->EW]
  #Calculate Kp values
  Kpu_bc <- (HCT - 1 + BP)/(HCT)
  
  X <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_IW-pKa),
              #3-monoprotic base
              10^(pKa-pH_IW),
              #4-diprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_IW)+10^(pKa[1]+pKa[2]-2*pH_IW), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_IW)+10^(pH_IW-pKa[1]),  
              #7-triprotic acid
              10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2])+10^(3*pH_IW-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_IW)+10^(pKa[3]+pKa[2]-2*pH_IW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_IW),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_IW)+10^(pH_IW-pKa[1])+10^(2*pH_IW-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_IW-pKa[1])+10^(pKa[3]-pH_IW)+10^(pKa[2]+pKa[3]-2*pH_IW))       
  
  Y <- switch(type,
              #1-neutral
              0,   
              #2-monoprotic acid
              10^(pH_EW-pKa),
              #3-monoprotic base
              10^(pKa-pH_EW), 
              #4-diprotic acid
              10^(pH_EW-pKa[1])+10^(2*pH_EW-pKa[1]-pKa[2]),
              #5-diprotic base
              10^(pKa[2]-pH_EW)+10^(pKa[1]+pKa[2]-2*pH_EW), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_EW)+10^(pH_EW-pKa[1]),  
              #7-triprotic acid
              10^(pH_EW-pKa[1])+10^(2*pH_EW-pKa[1]-pKa[2])+10^(3*pH_EW-pKa[1]-pKa[2]-pKa[3]),  
              #8-triprotic base
              10^(pKa[3]-pH_EW)+10^(pKa[3]+pka[2]-2*pH_EW)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_EW),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_EW)+10^(pH_EW-pKa[1])+10^(2*pH_EW-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_EW-pKa[1])+10^(pKa[3]-pH_EW)+10^(pKa[2]+pKa[3]-2*pH_EW))       
  
  Z <- switch(type,
              #1-neutral
              1,   
              #2-monoprotic acid
              1,
              #3-monoprotic base
              10^(pKa-pH_RBC), 
              #4-diprotic acid
              1,
              #5-diprotic base
              10^(pKa[2]-pH_RBC)+10^(pKa[1]+pKa[2]-2*pH_RBC), 
              #6-monoprotic acid monoprotic base (acid comes first)
              10^(pKa[2]-pH_RBC)+10^(pH_RBC-pKa[1]),  
              #7-triprotic acid
              1,  
              #8-triprotic base
              10^(pKa[3]-pH_RBC)+10^(pKa[3]+pka[2]-2*pH_RBC)+10^(pKa[1]+pKa[2]+pKa[3]-3*pH_RBC),  
              #9-diprotic acid monoprotic base (first two are acid)
              10^(pKa[3]-pH_RBC)+10^(pH_RBC-pKa[1])+10^(2*pH_RBC-pKa[1]-pKa[2]), 
              #10-diprotic base monoprotic acid (first one is acid)
              10^(pH_RBC-pKa[1])+10^(pKa[3]-pH_RBC)+10^(pKa[2]+pKa[3]-2*pH_RBC)) 
  
  
  Ka_PR <- (1/fup - 1 - (P*dat_plas$f_n_l + (0.3*P + 0.7)*dat_plas$f_n_pl)/(1+Y))
  Ka_AP <- (Kpu_bc - (1 + Z)/(1 + Y)*dat_rbc$f_iw - (P*dat_rbc$f_n_l + (0.3*P + 0.7)*dat_rbc$f_n_pl)/(1 + Y)) * (1 + Y)/dat_rbc$f_a_pl/Z
  
  
  # Assign the moderate to strong bases type_calc=1 and everything else type_calc=2 
  type_calc <- ifelse((type==3 & pKa[1]>7) | (type==5 & pKa[1] >7) | (type==6 & pKa[2] > 7) | (type==8 & pKa[1] > 7) | (type==9 & pKa[3]>7) | (type==10 & pKa[2]>7), 1,2)
  
  # Re-assign the neutrals type_calc=3
  if(type==1){type_calc=3}  #neutrals
  
  
  if(type_calc==1){  #moderate to strong bases
    Kp_all <- (((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$fnl_i + (0.3*P + 0.7)*dat_all$fnpl_i))/(1 + Y) /(Ka_AP*dat_all$AP_i*X)/(1 + Y))*fup
    
  }else if(type_calc==2){   #acidic and zwitterions
    Kp_all <- (((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$fnl_i + (0.3*P + 0.7)*dat_all$fnpl_i))/(1 + Y))*fup   
    
  }else{  #neutrals
    Kp_all <- (((1 + X)/(1 + Y))*dat_all$f_iw + ((P*dat_all$fnl_i + (0.3*P + 0.7)*dat_all$fnpl_i))/(1 + Y))*fup   #non lipid
    
  }
  
  
  nms_all <- dat_all$tissue %>% substr(1,3) %>% tolower()
  nms_all <- paste("Kp_ic", nms_all, sep="")
  
  Kp_ic <- as.list(c(Kp_all))
  names(Kp_ic) <- nms_all
  
  
  return(Kp_ic)
}

 
dat5 <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models/ROIvalues_final_628.csv")
dat6 <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models/ROIvalues_final_628_pt1.csv")
dat7 <- read.csv("/Volumes/GoogleDrive/My Drive/PBPK/MRG solve PBPK/models/ROIvalues_final_628_pt2.csv")
#Acetaminophen
#Logp pka fup bp type #dat

Kp_ac_ec <- calcKp_EC(0.46, #logp
                      pKa=9.38, #pka
                      0.855, #fup
                      1.09, #bp
                      type=1, #acid/base type - See RR 
                      dat5 #tissue data
)



Kp_ac_ic <- calcKp_IC(0.46, #logp
                      pKa=9.38, #pka
                     0.855,
                      1.09, #bp
                      type=1, #acid/base type - See RR 
                      dat5 #tissue data
)


ac_kpic <- as.data.frame(Kp_ac_ic) 
drug1 <- as.data.frame ("Acetaminophen")
ac_kpic <- cbind(drug1,ac_kpic)
colnames(ac_kpic) <- c('Drug','Thalamus Kp IC', 'Basal Ganglia KP IC','Cortex KP IC','ROB KP IC')

ac_kpec <- as.data.frame(Kp_ac_ec) 
drug1 <- as.data.frame ("Acetaminophen")
ac_kpec <- cbind(drug1,ac_kpec)
colnames(ac_kpec) <- c('Drug','Thalamus Kp EC', 'Basal Ganglia KP EC','Cortex KP EC','ROB KP EC')



#Morphine


Kp_mo_ec <- calcKp_EC(0.87,pKa=10.96,0.65,1.02,type=3,dat5)
Kp_mo_ic <- calcKp_IC(0.87,pKa=10.96,0.65,1.02,type=3,dat5)



mo_kpic <- as.data.frame(Kp_mo_ic) 
drug1 <- as.data.frame ("Morphine")
mo_kpic <- cbind(drug1,mo_kpic)
colnames(mo_kpic) <- c('Drug','Thalamus Kp IC', 'Basal Ganglia KP IC','Cortex KP IC','ROB KP IC')

mo_kpec <- as.data.frame(Kp_mo_ec) 
drug1 <- as.data.frame ("Morphine")
mo_kpec <- cbind(drug1,mo_kpec)
colnames(mo_kpec) <- c('Drug','Thalamus Kp EC', 'Basal Ganglia KP EC','Cortex KP EC','ROB KP EC')

#Phenytoin


Kp_pt_ec <-  calcKp_EC(2.47,pKa=9.47,0.1,0.61,type=2,dat5)
Kp_pt_ic <-  calcKp_IC(2.47,pKa=9.47,0.1,0.61,type=2,dat7)



pt_kpic <- as.data.frame(Kp_pt_ic) 
drug1 <- as.data.frame ("Phenytoin")
pt_kpic <- cbind(drug1,pt_kpic)
colnames(pt_kpic) <- c('Drug','Thalamus Kp IC', 'Basal Ganglia KP IC','Cortex KP IC','ROB KP IC')

pt_kpec <- as.data.frame(Kp_pt_ec) 
drug1 <- as.data.frame ("Phenytoin")
pt_kpec <- cbind(drug1,pt_kpec)
colnames(pt_kpec) <- c('Drug','Thalamus Kp EC', 'Basal Ganglia KP EC','Cortex KP EC','ROB KP EC')




#nilotinib

Kp_nl_ec <-   calcKp_EC(4.51,pKa=c(5.4,3.9),0.016,0.68,type=1,dat5)
Kp_nl_ic <-   calcKp_IC(4.51,pKa=c(5.4,3.9),0.016,0.68,type=1,dat5)
 
nl_kpic <- as.data.frame(Kp_nl_ic) 
drug1 <- as.data.frame ("Nilotinib")
nl_kpic <- cbind(drug1,nl_kpic)
colnames(nl_kpic) <- c('Drug','Thalamus Kp IC', 'Basal Ganglia KP IC','Cortex KP IC','ROB KP IC')

nl_kpec <- as.data.frame(Kp_nl_ec) 
drug1 <- as.data.frame ("Nilotinib")
nl_kpec <- cbind(drug1,nl_kpec)
colnames(nl_kpec) <- c('Drug','Thalamus Kp EC', 'Basal Ganglia KP EC','Cortex KP EC','ROB KP EC')





kp_ec <- rbind(ac_kpec, mo_kpec, nl_kpec, pt_kpec)



kp_ic <- rbind(ac_kpic,mo_kpic, nl_kpic, pt_kpic)

knitr::kable(kp_ec,"pipe", digits=3)
knitr::kable(kp_ic, "pipe",digits=3)


