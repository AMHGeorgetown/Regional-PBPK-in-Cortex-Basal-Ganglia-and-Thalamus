
[ PROB ] 

# CNS PBPK MODEL

This BASE model implements a WBPBPK with CNS model. 
CNS model comprises 4 ROIs in the Brain and 3CSF compartments
Perfusion compartments are applied to each brain ROI and 1 CSF compartment
All Brain compartments have identical parameters split through the 4 regions 
including flow, SA, ECF/ICF volume, and partition coefficient


Drug specific parameters in all sections are set to 1 and included in idata set for simulations


[ SET ] end = 12, delta = 0.1

[ CMT ] 
D     //; dose
Aad   //; adipose
Agu   //; gut
Aki   //; kidney
Ali   //; liver
Amu   //; muscle
Asp   //; spleen
Ave   //; venous blood
Aar   //; arterial blood
Alp   //; rest of body
Aou   //; Out
//{CNS COMPARTMENTS}
//{Brain}
//; ROI1 - Cortex
Aecf1 //; ECF 
Aicf1  //; ICF 
Acb1   //; cerebral blood Brain
Alyso1 //; Lyso 

//; ROI2 - Basal Ganglia
Aecf2  //; ECF 
Aicf2  //; ICF 
Acb3   //; cerebral blood Brain
Alyso2 //; Lyso 


//; ROI3 - Thalamus
Aecf3 //; ECF 
Aicf3  //; ICF 
Acb4   //; cerebral blood Brain
Alyso3 //; Lyso 


//; ROI4 - ROB
Aecf4 //; ECF 
Aicf4  //; ICF 
Acb5   //; cerebral blood Brain
Alyso4 //; Lyso 

//{CSF}
Avent //; Ventricles
Acis   //; Basal Cisterns
Asas    //; Subarachnoid space
Acb2   //;  cerebral blood CSF (CP)

[ PARAM ] 

BW = 70000	//; BW (g)

// {Fractional tissue volumes}

FVad = 0.213     //; adipose
FVbr = 0.02      //; brain
FVgu = 0.0171    //; gut
FVki = 0.0044    //; kidney
FVli = 0.021     //; liver 
FVmu = 0.4       //; muscle
FVsp = 0.0026    //; spleen 
FVve = 0.0514    //; venous
FVar = 0.0257    //; arterial 



//{CNS}


//{Brain}
Bden   =  1.035    //; density of brain tissue
Fvcb   = 0.10    //; CB/ROI
SAf    = 120       //; cm2/g tissue of ROI
Fvcx   = 0.70      //; fraction of total brain volume Cortex
FVbg   = 0.035     //; fraction of total brain volume BG
FVth   = 0.013     //; fraction of total brain volume Thalamus
FVicf  = 0.686     //; fraction on ROI volume for ICF
FVcyt  = 0.677     //; fraction on ROI volume for Cytosol (ICF less Lyso)
FVecf  = 0.215     //; fraction on ROI volume for ECF
FVlyso = 0.008     //; fraction volume on ROI ICF 


//{CSF}
FVcb2 =   0.000028 //; cerebral bood - CSF
FVcis  =  0.000071  //; basal cisterns
FVvent =  0.0005 //; ventricles
FVsas  =  0.00157   //; subarachnoid spaces




// {Fractional tissue blood flows}



FQad = 0.05      //; adipose 
FQgu = 0.146462  //; gut
FQki = 0.19      //; kidney 
FQh  = 0.215385  //; hepatic (venous side) 
FQrt = 1         //; lung
FQmu = 0.17      //; muscle 
FQsp = 0.017231  //; spleen
 
//{CNS}

tort = 1.5 //; tortousity of brain mass

//{Brain as a function of CO}
FQbr1 = 0.1039       //; Cortex
FQbr3 = 0.0059       //; BG
FQbr4 = 0.0022       //; TH
FQbr5 = 0.0087       //; ROB

//{CSF as a function of CO}
FQbr2 = 0.00138      //; CB  - CSF


//{ASymmetry factors}
//{BBB} If net influx, EF is set to 1 and IF is adjusted; If net efflux,
//I     IF is set to 1 and EF is adjusted

IF1  = 1; //; Asymmetry factor from BBB to MV Influx
EF1  = 1; //; Asymmetry factor from MV to BBB Efflux
//{BCSFB}
IF2  = 1 ; // Asymmetry factor from BCSFB to MV
EF2  = 1; // Asymmetry factor from MV to BCSFB


//{CNS specific parameter}
//{Surface Areas have to be in cm}

Qcsf   =  21 //; CSF flow 
Fsacp  = 8333 //; Surface area of CP as cm2/g 
Masscp = 0.000028 //; mass of CP as function of BW

Ftm =    0.998 //; fraction available for TM 
Fpbbb  = 0.00006   //; fraction available for PC bbb
Fpbcsfb= 0.00016   // fraction available for PC bcsfb

WIroi1 = 0.000005  //; Width Paracellular RO1 1
WIroi2 = 0.000005  //; Width Paracellular RO1 2
WIroi3 = 0.000005  //; Width Paracellular RO1 3
WIroi4 = 0.000005  //; Width Paracellular RO1 4




WIbcsfb = 0.000005//; Width BCSFB/BBB -> distance to travel, not gaps b/w cells


Fsabcm1 = 0.2375 // fraction of ICF volume m2/mL -> convert to cm2 by 10000; 
Fsabcm2 = 0.2375 //;
Fsabcm3 = 0.2375 //;
Fsabcm4 = 0.2375//;  

Fsalyso1 = 0.1875 //; Fraction of ICF volume
Fsalyso2 = 0.1875 //; 
Fsalyso3 = 0.1875 //; 
Fsalyso4 = 0.1875 //; 


// {COMPOUND SPECIFIC PARAMETERS}


// {Drug Specific Parameters}

Fucsf = 0.999   //; fraction unbound in CSF
fup   = 0.855  //; fraction unbound in plasma
BP    = 1.09  //; blood to plasma ratio
mw    = 151.2   //; Molecular weight
logp  = 0.46   //; n-octanol lipophilicity

acid   = 0//; Drug is monoprotic acide (1= yes, 0= no)  
base   = 0 //; Drug is monoprotic base  (1=yes, 0 = no)
dpbase = 0 //; Drug is diprotic base  ( 1= yes, 0 = no)
neut   = 1

pka1 = 9.46 //;  
pka2 = 1   //; 

 
pHp    = 7.4 //; pH plasma
pHecf  = 7.4//;   pH ECF
pHicf  = 7.0//;   pH ICF
pHcsf  = 7.3//;   pH CSF
pHlyso = 5.0//;   pH Lyso


// {Clearances}



CLhep= 0  // Hepatic clearance
CLrenal   =  0  //; Renal clearance 
CLbrain = 0  

// {Absorption}
oral = 0  // 1 for yes, 0 for no
IV   = 0 //; 1 for yes , 0 for no
IRa  = 0//; infusion rate Dose/Time (hr-1) a
Abs = 0  //; Ka (hr-1)
Kpf=0
 

[ MAIN  ]

// {Tissue to plasma partition coefficients}

double Kpad = 0.260 * Kpf; //; adipose

double Kpgu = 0.784 * Kpf; //; gut
double Kpki = 0.757 * Kpf; //; kidney
double Kpli = 0.734 * Kpf;  //; liver
double Kpmu = 0.684 * Kpf; //; muscle
double Kpsp = 0.720 * Kpf;  //; spleen
double Kplp = 1 * Kpf; //; low perfusion


double Kpbrec1 = 0.157 * Kpf;//; Cortex cb to ecf 
double Kpbric1 = 0.515 * Kpf;//; Ecf to ICF


double Kpbrec3 = 0.189* Kpf;//; BG cb to ecf 
double Kpbric3 = 0.589 * Kpf;//; Ecf to ICF



double Kpbrec4 = 0.192* Kpf;//; TH cb to ecf 
double Kpbric4 = 0.544 * Kpf;  //; Ecf to ICF

double Kpbrec5 = 0.192* Kpf;//; ROB cb to ecf 
double Kpbric5 = 0.591 * Kpf; //; Ecf to ICF


// {Total tissue volumes - mL}



double CO =  15 * double(pow((BW/1000),0.74));

double FVlp = 1 -  (FVad + FVgu + FVki + FVli + FVmu + 
				   FVsp + FVve + FVar + FVbr  +
				   FVcb2  + FVvent + FVcis  + FVsas) ;



double Vad = BW*FVad;  // adipose 

double Vbr = BW*FVbr;  // brain 
double Vgu = BW*FVgu;  // gut 

double Vki = BW*FVki;  // kidney
double Vli = BW*FVli;  // liver 
 
double Vmu = BW*FVmu;  // muscle

double Vsp = BW*FVsp;  // spleen
 
double Vve = BW*FVve;  // venous blood
double Var = BW*FVar;  // arterial blood

double Vlp = BW*FVlp;  // rest of body
double Vou = 1      ;  //out







//{CNS}


//{Brain}


double TBM = BW * FVbr; // total brain MASS
double TBV = TBM * (1/Bden);// total brain VOLUME





//{ROI 1 - Cortex}

double Vcx    = Fvcx * TBV       ; // Volume CX
double SAroi1 = Vcx * Bden * SAf ; //Surface area CX MV
double SAroi1t = SAroi1 * Ftm    ; //; adjusted for SA available for TC
double SAroi1p = SAroi1 *Fpbbb   ; //; Adjusted for SA available for PC
double Vcb1   = Fvcb   * Vcx     ; //; Cerebral blood supply - Brain ROI1
double Vlyso1 = FVlyso * Vcx     ; //; Lyso ROI 1
double Vicf1  = FVcyt *Vcx       ; //: ICF ROI 1/Cytosol only
double Vecf1  = FVecf  * Vcx     ; //; Volume ECF relative to Brain
double Ticf1  = FVicf   * Vcx    ; // Total volume of ICF
double SAbcm1 = Fsabcm1 * Ticf1 * 10000  ; //SA brain cell membrane
double SAlyso1 = Fsalyso1 * Ticf1 * 10000; //SA lysosomal membrane


//{ROI 2 - Basal Ganglia}
double Vbg    = FVbg * TBV  ; //Volume BG
double SAroi2 = Vbg  * Bden * SAf; //Surface area BG MV
double SAroi2t = SAroi2 * Ftm; 
double SAroi2p = SAroi2 * Fpbbb; 
double Vlyso2 = FVlyso  * Vbg; 
double Vcb3   = Fvcb   * Vbg;  //; Cerebral blood supply - Brain ROI2 
double Vicf2  = FVcyt * Vbg; 
double Vecf2  = FVecf  * Vbg;
double Ticf2  = FVicf   * Vbg   ; // Total volume of ICF
double SAbcm2 = Fsabcm2 * Ticf2 * 10000 ; //SA brain cell membrane
double SAlyso2 = Fsalyso2 * Ticf2 * 10000; //SA lysosomal membrane

//{ROI 3 - Thalamus}
double Vth    = FVth * TBV; // Volume Thalamus
double SAroi3 = Vth * Bden * SAf; 
double SAroi3t = SAroi3* Ftm; 
double SAroi3p = SAroi3 * Fpbbb; 
double Vlyso3 = FVlyso * Vth ; 
double Vcb4   = Fvcb   *  Vth;  //; Cerebral blood supply - Brain ROI2 
double Vicf3  = FVcyt*Vth;  
double Vecf3  = FVecf  * Vth;
double Ticf3  = FVicf * Vth ; 
double SAbcm3 = Fsabcm3 * Ticf3 * 10000 ; //SA brain cell membrane
double SAlyso3 = Fsalyso3 * Ticf3 * 10000; //SA lysosomal membrane





//{ROI 4 - ROB}
double Vrob = TBV - (Vcx + Vth + Vbg); 
double SAroi4 = Vrob * Bden * SAf; 
double SAroi4t = SAroi4 * Ftm  ; 
double SAroi4p = SAroi4 * Fpbbb; 
double Vlyso4  = FVlyso * Vth ; 
double Vcb5   = Fvcb * Vrob ; 
double Vicf4  = FVcyt * Vrob; 
double Vecf4 = FVecf * Vrob ; 
double Ticf4 = FVicf * Vrob;  
double SAbcm4 = Fsabcm4 * Ticf4 * 10000 ; //SA brain cell membrane
double SAlyso4 = Fsalyso4 * Ticf4 * 10000; //SA lysosomal membrane




//{CSF}
double Vcb2  = FVcb2  * BW ; //; Cerebral blood supply - CSF
double Vvent = FVvent * BW ; //; Volume Ventricles to Whole body
double Vcis  = FVcis  * BW ; //; basal cisterns volume to whole body
double Vsas  = FVsas  * BW ; //; SAS volume
double Mcp   = Masscp * BW ; //; mass cp 
double SAbcsfb = Mcp * Fsacp ; // SAbcsbf
double SAbcsfbt = SAbcsfb * Ftm;
double SAbcsfbp = SAbcsfb * Fpbcsfb;
// {Total tissue blood flows - mL/hr}
double FQlp = 1 - (FQad + FQki + FQh + 
				   FQmu + 
			   FQbr1 + FQbr2 + FQbr3 + FQbr4 + FQbr5) ; 


double QC  = CO*1000 ;  // cardiac output (mL/hr)
double Qad = QC*FQad       ;  // adipose 


double Qgu = QC*FQgu       ;  // gut

double Qki = QC*FQki       ;  // kidney 
double Qh  = QC*FQh        ;  // hepatic (venous side)
double Qha = Qh - Qgu - Qsp;  // hepatic artery 
double Qrt = QC*FQrt       ;  // Return flow 
double Qmu = QC*FQmu       ;  // muscle 

double Qsp = QC*FQsp       ;  // spleen 
 
double Qlp = QC * FQlp     ;  // rest of body

//{CNS}

//{Brain}
double Qbr1 = QC * FQbr1    ;  // Brain ROI 1
double Qbr3 = QC * FQbr3    ;  // Brain ROI 2
double Qbr4 = QC * FQbr4    ;  // Brain ROI 3
double Qbr5 = QC * FQbr5    ;  // Brain ROI 4





//{CSF}
double Qbr2 = QC * FQbr2      ;  //  CSF


//{CNS parameters}

double PHF11 = acid * ((1+(pow(10,(pHp - pka1)))) / (1+(pow(10, (pHecf-pka1))))) + 
               base * ((1+(pow(10,(pHp - pka1)))) / (1+(pow(10, (pka1-pHecf))))) +
             dpbase * ((1 +(pow(10, pka2 - pHp) + pow(10,(pka1 + pka2 - 2 * pHp))))/
                        1 +(pow(10, pka2 - pHecf) + pow(10,(pka1 + pka2 - 2 * pHecf))))+
                        neut*1; //; ionization factor from ECF to CB
double PHF12 = PHF11 ; //; ionization factor from ECF to CB
double PHF13 = PHF11 ; //; ionization factor from ECF to CB
double PHF14 = PHF11 ; //; ionization factor from ECF to CB





double PHF2 = acid * ((1+(pow(10,(pHp - pka1)))) / (1+(pow(10, (pHcsf-pka1))))) + 
               base * ((1+(pow(10,(pHp - pka1)))) / (1+(pow(10, (pka1-pHcsf))))) +
             dpbase * ((1 +(pow(10, pka2 - pHp) + pow(10,(pka1 + pka2 - 2 * pHp))))/
                        1 +(pow(10, pka2 - pHcsf) + pow(10,(pka1 + pka2 - 2 * pHcsf))))+
                        neut*1 ; //; ionization factor from CSF to CB


double PHF31 = acid * ((1+(pow(10,(pHecf - pka1)))) / (1+(pow(10, (pHicf-pka1))))) + 
               base * ((1+(pow(10,(pHecf - pka1)))) / (1+(pow(10, (pka1-pHicf))))) +
             dpbase * ((1 +(pow(10, pka2 - pHecf) + pow(10,(pka1 + pka2 - 2 * pHecf))))/
                        1 +(pow(10, pka2 - pHicf) + pow(10,(pka1 + pka2 - 2 * pHicf)))) + 
                        neut*1; //; Ionization factor from ICF to ECF
double PHF32 = PHF31 ; 
double PHF33 = PHF31 ; 
double PHF34 = PHF31 ; 




double PHF41 = acid * ((1+(pow(10,(pHicf - pka1)))) / (1+(pow(10, (pHlyso-pka1))))) + 
               base * ((1+(pow(10,(pHicf - pka1)))) / (1+(pow(10, (pka1-pHlyso))))) +
             dpbase * ((1 +(pow(10, pka2 - pHicf) + pow(10,(pka1 + pka2 - 2 * pHicf))))/
                        1 +(pow(10, pka2 - pHlyso) + pow(10,(pka1 + pka2 - 2 * pHlyso))))+
                        neut*1 ; //; Ionization factor from Lyso to ICF
double PHF42 = PHF41 ; 
double PHF43 = PHF41 ; 
double PHF44 = PHF41 ; 



//{CNS permeability and diffusion}

double MWl = log10(mw) ;
double DaQ = (double (pow(10, double (-4.113 - 0.4609 * MWl)))) ;
double PoT = (double (pow(10, double(0.939 * logp - 6.210))))   ;

double Qecf = DaQ / (pow(2,tort)) ; // calculated ECF diffusion

double QTroi1 = 0.5 * PoT * SAroi1t  * 60; 
double QTroi2 = 0.5 * PoT * SAroi2t  * 60; 
double QTroi3 = 0.5 * PoT * SAroi3t  * 60; 
double QTroi4 = 0.5 * PoT * SAroi4t  * 60; 

double QProi1 = (DaQ/WIroi1) * SAroi1p * 60; 
double QProi2 = (DaQ/WIroi2) * SAroi2p * 60; 
double QProi3 = (DaQ/WIroi3) * SAroi3p * 60; 
double QProi4 = (DaQ/WIroi4) * SAroi4p * 60; 

double QTbcsfb =  (0.5 * PoT * SAbcsfbt) * 60 ; 
double QPbcsfb = (((DaQ/WIbcsfb) * SAbcsfbp)) * 60; 

double cl11 = (QTroi1 + QProi1)*IF1; //From CB to ECF 
double cl21 = cl11  * EF1 ; //From ECF to CB


double cl12 = (QTroi2+ QProi2)* IF1; //From CB to ECF  
double cl22 = cl12  * EF1 ; //From ECF to CB


double cl13 = (QTroi3 + QProi3) * IF1; //From CB to ECF
double cl23 = cl13  * EF1 ; //From ECF to CB



double cl14 = (QTroi4 + QProi4) * IF1; //From CB to ECF  
double cl24 = cl14  * EF1 ; //From ECF to CB




double cl3 = (QTbcsfb + QPbcsfb) * IF2 ;  //; from choroid plexus blood to CSF
double cl4 = cl3 * PHF2 * EF2  ; //; from CSF to choroid plexus blood


double Qbcmin1 = SAbcm1 * PoT * (60);  //BCM permeability
double Qbcmou1 = Qbcmin1  ;

double Qbcmin2 = SAbcm2 * PoT * (60);  //BCM permeability
double Qbcmou2 = Qbcmin2  ;

double Qbcmin3 = SAbcm3 * PoT * (60);  //BCM permeability
double Qbcmou3 = Qbcmin3  ;

double Qbcmin4 = SAbcm4 * PoT * (60);  //BCM permeability
double Qbcmou4 = Qbcmin4  ;





double Qlysoi1 = SAlyso1 * PoT * (60);
double Qlysoo1 = Qlysoi1 * PHF41; 

double Qlysoi2 = SAlyso2 * PoT * (60);
double Qlysoo2 = Qlysoi2 * PHF42;

double Qlysoi3 = SAlyso3 * PoT * (60);
double Qlysoo3 = Qlysoi3 * PHF43;

double Qlysoi4 = SAlyso4 * PoT * (60);
double Qlysoo4 = Qlysoi4 * PHF44; 
 

 


[ ODE ]

double Cadipose  = Aad/Vad;  // adipose 


double Cgut      = Agu/Vgu;  // gut

double Ckidney   = Aki/Vki;  // kidney 
double Cliver    = Ali/Vli;  // liver 
 
double Cmu   = Amu/Vmu;  // muscle

double Cspleen   = Asp/Vsp;  // spleen 
 
double Cvenous   = Ave/Vve;  // venous blood
double Carterial = Aar/Var;  // arterial blood
double Clp     = Alp/Vlp;  // rest of body


//{CNS}

//{Brain}
double Ccb1 = Acb1 / Vcb1       ;
double Cecf1 = Aecf1 / Vecf1    ; 
double Cicf1  = Aicf1 / Vicf1   ;
double Clyso1 = Alyso1 / Vlyso1 ;


double Ccb3 = Acb3 / Vcb3       ;
double Cecf2 = Aecf2 / Vecf2    ; 
double Cicf2  = Aicf2 / Vicf2   ;
double Clyso2 = Alyso2 / Vlyso2 ;


double Ccb4   = Acb4 / Vcb4     ;
double Cecf3  = Aecf3 / Vecf3   ; 
double Cicf3  = Aicf3 / Vicf3   ;
double Clyso3 = Alyso3 / Vlyso3 ;

double Ccb5   = Acb5 / Vcb5     ;
double Cecf4  = Aecf4 / Vecf4   ; 
double Cicf4  = Aicf4 / Vicf4   ;
double Clyso4 = Alyso4 / Vlyso4 ;




 
//{CSF}
double Ccb2 = Acb2 / Vcb2       ; 
double Cvent = Avent / Vvent    ; 
double Ccis  = Acis / Vcis      ; 
double Csas = Asas / Vsas       ; 


// {Calculation of free concentrations - mg/L}

double Cliverfree  = Cliver*fup;  // liver 
double Ckidneyfree = Ckidney*fup; // kidney 
 
// {Clearance calculations}



double Venous = 
  Qad*(Cadipose/Kpad*BP)  + 
  Qki*(Ckidney/Kpki*BP) + 
  Qh*(Cliver/Kpli*BP)    + Qmu*(Cmu/Kpmu*BP) +
  Qlp*(Clp/Kplp*BP) +
  Qbr1*(Ccb1 * BP) + Qbr2 * (Ccb2*BP ) + Qbr3 * (Ccb3*BP) +
  Qbr4*(Ccb4*BP) +   Qbr5*(Ccb5*BP) + (IV * Absorption) + Qcsf*Csas

   ;

double Absorption = (oral*D*Abs) + (IV * D*IRa);

dxdt_Aad = Qad*(Carterial - Cadipose/Kpad*BP);    // adipose


dxdt_Agu = (oral* Absorption) + 
           Qgu*(Carterial - Cgut/Kpgu*BP);        // gut

dxdt_Aki = Qki*(Carterial - Ckidney/Kpki*BP) - 
           CLrenal*Ckidneyfree;                   // kidney
dxdt_Ali = Qha*Carterial + 
           Qgu*(Cgut/Kpgu*BP) + 
           Qsp*(Cspleen/Kpsp*BP) - 
           Qh*(Cliver/Kpli*BP) - 
           Cliverfree*CLhep;                      // liver 

dxdt_Amu = Qmu*(Carterial - Cmu/Kpmu*BP);     // muscle

dxdt_Asp = Qsp*(Carterial - Cspleen/Kpsp*BP);     // spleen
 
dxdt_Ave = Venous - Qrt*Cvenous ;                  // venous blood
dxdt_Aar = Qrt*(Cvenous) - Qrt*Carterial;   // arterial blood
dxdt_Alp = Qlp*(Carterial - Clp/Kplp*BP);       // rest of body
dxdt_D   = - Absorption;                          // oral dosing
dxdt_Aou =  Cliverfree*CLhep + CLrenal*Ckidneyfree + 
            (Cicf1/Kpbric1)*CLbrain + (Cicf2/Kpbric3)*CLbrain + 
            (Cicf3/Kpbric4)*CLbrain + (Cicf4/Kpbric5)*CLbrain; 

//{CNS model}


//{Brain}
dxdt_Acb1 = Qbr1 * (Carterial - Ccb1*BP) - (cl11 *Ccb1*BP) + (cl21 *(Cecf1/Kpbrec1)) ; // Cerebral blood
dxdt_Aecf1 = cl11*Ccb1*BP - cl21*(Cecf1/Kpbrec1) - Qecf * (Cecf1) - Qbcmin1*(Cecf1/Kpbrec1) + Qbcmou1*(Cicf1/Kpbric1);   //ECF of the brain
dxdt_Aicf1 = Qbcmin1 * (Cecf1/Kpbrec1) -(Qbcmou1 * (Cicf1/Kpbric1)) - (Cicf1/Kpbric1)*CLbrain - Qlysoi1 * (Cicf1/Kpbric1) + Qlysoo1 * Clyso1 ; //ICF of brain
dxdt_Alyso1 = (Qlysoi1 * Cicf1/Kpbric1)  - Qlysoo1 * Clyso1;


dxdt_Acb3 = Qbr3 * (Carterial - (Ccb3*BP)) - (cl12 * Ccb3*BP) + (cl22 *(Cecf2/Kpbrec3)) ; // Cerebral blood
dxdt_Aecf2 = cl12*Ccb3*BP - cl22*(Cecf2/Kpbrec3) - Qecf * (Cecf2) - Qbcmin2*(Cecf2/Kpbrec3) + Qbcmou2*(Cicf2/Kpbric3);   //ECF of the brain
dxdt_Aicf2 = Qbcmin2 * (Cecf2/Kpbrec3) - Qbcmou2 * (Cicf2/Kpbric3) - Qlysoi2 * (Cicf2/Kpbric3) -(Cicf2/Kpbric3)*CLbrain + Qlysoo2 * Clyso2 ; //ICF of brain
dxdt_Alyso2 = (Qlysoi2 * Cicf2/Kpbric3)  - Qlysoo2 * Clyso2;


dxdt_Acb4 = Qbr4 * (Carterial - (Ccb4*BP)) - (cl13 * Ccb4*BP) + (cl23 * (Cecf3/Kpbrec4)) ; // Cerebral blood
dxdt_Aecf3 = cl13*Ccb4*BP - cl23*(Cecf3/Kpbrec4) - Qecf * (Cecf3) - Qbcmin3*(Cecf3/Kpbrec4) + Qbcmou3*(Cicf3/Kpbric4);   //ECF of the brain
dxdt_Aicf3 = Qbcmin3 * (Cecf3/Kpbrec4) - Qbcmou3 * (Cicf3/Kpbric4) - Qlysoi3 * (Cicf3/Kpbric4) -(Cicf3/Kpbric4)*CLbrain + Qlysoo3 * Clyso3 ; //ICF of brain
dxdt_Alyso3 = (Qlysoi3 * Cicf3/Kpbric4)  - Qlysoo3 * Clyso3;



dxdt_Acb5 = Qbr5 * (Carterial - (Ccb5*BP)) - (cl14 * Ccb5*BP ) + (cl24 * Cecf4/Kpbrec5) ; // Cerebral blood
dxdt_Aecf4 = cl14*Ccb5*BP - cl24*(Cecf4/Kpbrec5) - Qecf * (Cecf4) - Qbcmin4*(Cecf4/Kpbrec5) + Qbcmou4*(Cicf4/Kpbric5);   //ECF of the brain
dxdt_Aicf4 = Qbcmin4 * (Cecf4/Kpbrec5) - Qbcmou4 * (Cicf4/Kpbric5) - Qlysoi4 * (Cicf4/Kpbric5) -(Cicf4/Kpbric5)*CLbrain + Qlysoo4 * Clyso4 ; //ICF of brain
dxdt_Alyso4 = (Qlysoi4 * Cicf4/Kpbric5)  - Qlysoo4 * Clyso4 ;






//{CSF}
dxdt_Acb2 = Qbr2 * (Carterial - (Ccb2*BP)) - (cl3 * Ccb2*BP) + (cl4 * Cvent)  ; // CB to ventricles
dxdt_Avent = cl3 * Ccb2*BP - (cl4* Cvent) + Qecf * (Cecf1 + Cecf2 + Cecf3 + Cecf4)-  Qcsf * Cvent; 
dxdt_Acis = Qcsf * Cvent - Qcsf*Ccis; 
dxdt_Asas = Qcsf * Ccis - Qcsf*Csas; 



[ TABLE ] 
Cvenous = Ave/Vve;
capture Cp = (Cvenous/BP) ; // venous plasma in ng/mL
capture Cpm = (Cp/mw) * 1000; // concentration in nmol/l
 
capture Mball = (Aad + Agu + Aki + Ali + 
                Amu + Asp + Ave + Aar +
                Alp + Aou + D + 
                Acb1 + Aecf1 + Aicf1 + Alyso1 + 
                Acb3 + Aecf2 + Aicf2 + Alyso2 +
                Acb4 + Aecf3 + Aicf3 + Alyso3 +
                Acb5 + Aecf4 + Aicf4 + Alyso4 +
                Acb2 + Avent + Acis + Asas); 

capture amt = D; 

capture Csa = (Asas/Vsas) * Fucsf ; // SAS in mcg/mL 
capture Csam = (Csa/mw) *1000; // 

capture Cec1 = Aecf1 / Vecf1 ;
capture Cec2 = Aecf2 / Vecf2 ;
capture Cec3 = Aecf3 / Vecf3 ;   
capture Cec4 = Aecf4 / Vecf4 ;

capture Icf1  = Aicf1 / Vicf1 ; 
capture Icf2  = Aicf2 / Vicf2 ; 
capture Icf3  = Aicf3 / Vicf3 ; 
capture Icf4  = Aicf4 / Vicf4 ; 

capture Cmva = (Acb1+Acb3+Acb4+Acb5)/(Vcb1+Vcb3+Vcb4+Vcb5);




capture Ceca = (Aecf1 + Aecf2 + Aecf3 + Aecf4 ) / (Vecf1 + Vecf2 + Vecf3 + Vecf4) ;
capture Cicfa = (Aicf1 + Aicf2 + Aicf3 + Aicf4) / (Vicf1 + Vicf2 + Vicf3 + Vicf4) ; 
capture Clysa = (Alyso1+Alyso2+Alyso3+Alyso4) / (Vlyso1 + Vlyso2 + Vlyso3 + Vlyso4);
capture Bra = ((Aecf1 + Aecf2 + Aecf3 + Aecf4 ) + (Aicf1 + Aicf2 + Aicf3 + Aicf4)+
			(Alyso1+Alyso2+Alyso3+Alyso4)) /
			((Vecf1 + Vecf2 + Vecf3 + Vecf4)+ ((Vicf1 + Vicf2 + Vicf3 + Vicf4)+
			(Vlyso1 + Vlyso2 + Vlyso3 + Vlyso4)));
 
[ CAPTURE ] Cvenous Cp Cpm Mball  amt Csa Csam Ceca Cec1 Cec2 Cec3 Cec4 Icf1 Icf2 Icf2 
			Icf4 Cicfa cl11 cl12 Qecf QProi1 QTroi1 SAroi1 cl21
			Ccb1 Ccb2 Ccb3 Ccb4 Ccb5 Bra Cmva
  