
#This is the code necessary to include 15N mass balance constraints.
#Input arguments are the matrix of approximate equalities, a vector of flows, a vector of 
#the weighting factors to be applied to each flow, a vector holding the
#initial d15N values for each group, and each of the known 15N values.
ResetRN15 <- function(Aa,flows,wts,d15N,del15NMeso,del15NMacro,del15NSalp){
  
  flows <- flows/wts
  
  #Assumed fractionation constants from the supplement of Stukel 2018 CRD
  Eps_NO3up <- -5      #Uptake of nitrate, note this was -6 in CRD
  Eps_NH4up <- -10     #Uptake of ammonia
  Eps_Nit <- -14       #Nitrification (NH4->NO3)
  Eps_ExcL <- -5       #Excretion of Mesozoo and Fish
  Eps_EgL <- -1.5      #Egestion of Mesozoo and Fish
  Eps_ExcS <- -1       #Excretion of Protistan Zoos and salps
  Eps_EgS <- -1        #Egestion of Protistan Zoos and salps
  Eps_ExcB <- -5       #Excretion of Bacteria
  Eps_UpB <- -1        #Bacterial uptake of DOM
  Eps_Remin <- 0       #Solubilization of Det to DOM

  #Setting the appropriate entries in d15N0 that correspond to each starting 15N
  col0 <- 36  #The entry before the first 15N entry
  NO3_col <- 1
  NH4_col <- 2
  Pico_col <- 3
  Dtm_col <- 4
  Flag_col <- 5
  HNF_col <- 6
  Mic_col <- 7
  Meso_col <- 8
  Macro_col <- 9
  VMMacro_col <- 10
  Gel_col <- 11
  Salp_col <- 12
  Amph_col <- 13
  PiscivFish_col <- 14
  GelativFish_col <- 15
  Myct_col <- 16
  bac_col <- 17
  sdet_col <- 18
  mdet_col <- 19
  ldet_col <- 20
  SalpDet_col <- 21
  Dom_col <- 22
  SalpMort_col <- 23
  Export_col <- 24
  upNO3_col <- 25

  #Withholding columns that we have 15N data for and setting the rest to a starting point
  #of 0
  RNO3 <- d15N[NO3_col]
  RNH4 <- d15N[NH4_col]
  RPico <- d15N[Pico_col]
  RDtm <- d15N[Dtm_col]
  RFlag <- d15N[Flag_col]
  RHNF <- d15N[HNF_col]
  RMic <- d15N[Mic_col]
  #RMeso <- d15N[Meso_col]
  #RMacro <- d15N[Macro_col]
  #RVMMacro <- d15N[VMMacro_col]
  RGel <- d15N[Gel_col]
  #RSalp <- d15N[Salp_col]
  RAmph <- d15N[Amph_col]
  RPiscivFish <- d15N[PiscivFish_col]
  RGelativFish <- d15N[GelativFish_col]
  RMyct <- d15N[Myct_col]
  Rbac <- d15N[bac_col]
  Rsdet <- d15N[sdet_col]
  Rmdet <- d15N[mdet_col]
  Rldet <- d15N[ldet_col]
  RSalpDet <- d15N[SalpDet_col]
  RDom <- d15N[Dom_col]
  RSalpMort <- d15N[SalpMort_col]
  RExport <- d15N[Export_col]
  RupNO3 <- d15N[upNO3_col]
  

  RMeso <- del15NMeso
  RMacro <- del15NMacro
  RVMMacro <- del15NMacro
  # Note the following code is necessary because no salp 15N is available for nonsalp cycles 3 and 5.
  #Ergo if we have the data, use it. If not, set it to the default.
  if (identical(del15NSalp,NaN)){
    RSalp <- d15N[Salp_col]
  }else{
    RSalp <- del15NSalp
  }
  
  
  # RN2 <- 0.0036765
  # RupNO3 <- del15NupNO3/1000*RN2+RN2
  # RNO3s <- del15NNO3s/1000*RN2+RN2
  # RNH4s <- del15NNH4s/1000*RN2+RN2
  # RCyas <- del15NCyas/1000*RN2+RN2
  # RTris <- del15NTris/1000*RN2+RN2
  # RDtms <- del15NDtms/1000*RN2+RN2
  # RFlags <- del15NFlags/1000*RN2+RN2
  # RHNFs <- del15NHNFs/1000*RN2+RN2
  # RMICs <- del15NMICs/1000*RN2+RN2
  # RHerbNVMs <- del15NHerbNVMs/1000*RN2+RN2
  # RApps <- del15NApps/1000*RN2+RN2
  # RClads <- del15NClads/1000*RN2+RN2
  # RnvmCals <- del15NnvmCals/1000*RN2+RN2
  # RChaetos <- del15NChaetos/1000*RN2+RN2
  # RPoecils <- del15NPoecils/1000*RN2+RN2
  # RGelPreds <- del15NGelPreds/1000*RN2+RN2
  # RPlankfishs <- del15NPlankfishs/1000*RN2+RN2
  # RPreflexs <- del15NPreflexs/1000*RN2+RN2
  # RPostflexs <- del15NPostflexs/1000*RN2+RN2
  # RBacs <- del15NBacs/1000*RN2+RN2
  # RDets <- del15NDets/1000*RN2+RN2
  # RLdets <- del15NLdets/1000*RN2+RN2
  # RDons <- del15NDons/1000*RN2+RN2
  # RHerbVM <- del15HerbVM/1000*RN2+RN2
  # RvmCal <- del15NvmCal/1000*RN2+RN2


  #This is listing the columns for each flow in the A matrices. Copy from excel
  Upwelling	=1
  NH4toNO3	=2
  NH4toPico	=3
  NH4TOdtm	=4
  NH4TOflag	=5
  NO3toPico	=6
  NO3TOdtm	=7
  NO3TOflag	=8
  bacTOflag	=9
  bacTOhnf	=10
  bacTOmic	=11
  bacTOnh4	=12
  PicoTOhnf	=13
  PicoTOmic	=14
  PicoTOflag	=15
  PicoTOSalp	=16
  PicoTOdom	=17
  PicoTOsdet	=18
  dtmTOmic	=19
  dtmTOMeso	=20
  dtmTOMacro	=21
  dtmTOVMMacro	=22
  dtmTOSalp	=23
  dtmTOdom	=24
  dtmTOmdet	=25
  flagTOmic	=26
  flagTOMeso	=27
  flagTOMacro	=28
  flagTOVMMacro	=29
  flagTOSalp	=30
  flagTOdom	=31
  flagTOsdet	=32
  hnfTOmic	=33
  hnfTOMeso	=34
  hnfTOMacro	=35
  hnfTOVMMacro	=36
  hnfTOSalp	=37
  hnfTOnh4	=38
  hnfTOdom	=39
  hnfTOsdet	=40
  micTOMeso	=41
  micTOMacro	=42
  micTOVMMacro	=43
  micTOSalp	=44
  micTOnh4	=45
  micTOdom	=46
  micTOsdet	=47
  MesoTOMeso	=48
  MesoTOMacro	=49
  MesoTOVMMacro	=50
  MesoTOgel	=51
  MesoTOMyct	=52
  MesoTOHTL	=53
  MesoTOnh4	=54
  MesoTOdom	=55
  MesoTOmdet	=56
  MacroTOgel	=57
  MacroTOMyct	=58
  MacroTOHTL	=59
  MacroTOnh4	=60
  MacroTOdom	=61
  MacroTOldet	=62
  VMMacroTOPiscivFish	=63
  VMMacroTOgel	=64
  VMMacroTOMyct	=65
  VMMacroTOHTL	=66
  VMMacroTOnh4	=67
  VMMacroTOdom	=68
  VMMacroTOldet	=69
  VMMacroTOdnh4	=70
  VMMacroTOddom	=71
  gelTOAmph	=72
  gelTOMyct	=73
  gelTOGelativFish	=74
  gelTOHTL	=75
  gelTOnh4	=76
  gelTOdom	=77
  gelTOldet	=78
  MyctTOPiscivFish	=79
  MyctTOnh4	=80
  MyctTOdom	=81
  MyctTOddom =82
  MyctTOdnh4	=83
  MyctTOpoop	=84
  MyctTObiomass	=85
  SalpTOAmph	=86
  SalpTOMyct	=87
  SalpTOnh4	=88
  SalpTOdom	=89
  SalpTOSalpDet	=90
  SalpTOSalpMort	=91
  SalpTOHTL	=92
  SalpTOdnh4	=93
  SalpTOddom	=94
  SalpMortTOGelativFish	=95
  SalpMortTOSink	=96
  AmphTOPiscivFish	=97
  AmphTOMyct	=98
  AmphTOHTL	=99
  AmphTOnh4	=100
  AmphTOdom	=101
  AmphTOldet	=102
  PiscivFishTOdnh4	=103
  PiscivFishTOddom	=104
  PiscivFishTOpoop	=105
  PiscivFishTObiomass	=106
  GelativFishTOdnh4	=107
  GelativFishTOddom	=108
  GelativFishTOpoop	=109
  GelativFishTObiomass	=110
  domTObac	=111
  sdetTOhnf	=112
  sdetTOmic	=113
  sdetTOdom	=114
  mdetTOmic	=115
  mdetTOMeso	=116
  mdetTOdom	=117
  mdetTOSink	=118
  ldetTOdom	=119
  ldetTOSink	=120
  SalpDetTOdom	=121
  SalpDetTOSink	=122
  
  
  
  
  #Adding up the total amount consumed by each group
  HNFcons <- flows[PicoTOhnf] + flows[bacTOhnf] + flows[sdetTOhnf]
  MICcons <- flows[PicoTOmic] + flows[dtmTOmic] +flows[flagTOmic] + flows[hnfTOmic] + flows[bacTOmic] + flows[sdetTOmic]
  Mesocons <- flows[dtmTOMeso] + flows[flagTOMeso] + flows[hnfTOMeso] + flows[micTOMeso] + flows[MesoTOMeso] + flows[mdetTOMeso]
  Macrocons <- flows[dtmTOMacro] + flows[hnfTOMacro] + flows[micTOMacro] + flows[MesoTOMacro]
  VMMacrocons <- flows[dtmTOVMMacro] + flows[hnfTOVMMacro] + flows[micTOVMMacro] + flows[MesoTOVMMacro] 
  Gelcons <- flows[MesoTOgel]  + flows[MacroTOgel] + flows[VMMacroTOgel]
  Salpcons <- flows[PicoTOSalp] + flows[dtmTOSalp] + flows[flagTOSalp] + flows[hnfTOSalp] + flows[micTOSalp]
  Amphcons <- flows[gelTOAmph] + flows[SalpTOAmph] 
  PiscivFishcons <- flows[VMMacroTOPiscivFish] + flows[MyctTOPiscivFish] + flows[AmphTOPiscivFish]
  GelativFishcons <- flows[gelTOGelativFish] + flows[SalpMortTOGelativFish] 
  Myctcons <- flows[MesoTOMyct]  + flows[MacroTOMyct] + flows[VMMacroTOMyct] + flows[SalpTOMyct] + flows[AmphTOMyct]
  
  #Calculating the isotopic composition of the food of each group
  RHNFfood <- (RPico*flows[PicoTOhnf] + Rbac*flows[bacTOhnf] + Rsdet*flows[sdetTOhnf])/HNFcons
  RMICfood <- (RPico*flows[PicoTOmic] + RDtm*flows[dtmTOmic] + RFlag*flows[flagTOmic] + RHNF*flows[hnfTOmic] + Rbac*flows[bacTOmic] + Rsdet*flows[sdetTOmic])/MICcons
  RMesofood <- (RDtm*flows[dtmTOMeso] + RFlag*flows[flagTOMeso] + RHNF*flows[hnfTOMeso] + RMic*flows[micTOMeso] + RMeso*flows[MesoTOMeso] + Rmdet*flows[mdetTOMeso])/Mesocons
  RMacrofood <- (RDtm*flows[dtmTOMacro] + RHNF*flows[hnfTOMacro] + RMic*flows[micTOMacro] + RMeso*flows[MesoTOMacro])/Macrocons
  RVMMacrofood <- (RDtm*flows[dtmTOVMMacro] + RHNF*flows[hnfTOVMMacro] + RMic*flows[micTOVMMacro] + RMeso*flows[MesoTOVMMacro] )/VMMacrocons
  RGelfood <- (RMeso*flows[MesoTOgel] + RMacro*flows[MacroTOgel] + RVMMacro*flows[VMMacroTOgel])/Gelcons
  RSalpfood <- (RPico*flows[PicoTOSalp] + RDtm*flows[dtmTOSalp] + RFlag*flows[flagTOSalp] + RHNF*flows[hnfTOSalp] + RMic*flows[micTOSalp])/Salpcons
  RAmphfood <- (RGel*flows[gelTOAmph] + RSalp*flows[SalpTOAmph])/Amphcons
  RPiscivFishfood <- (RVMMacro*flows[VMMacroTOPiscivFish] + RMyct*flows[MyctTOPiscivFish] + RAmph*flows[AmphTOPiscivFish])/PiscivFishcons
  RGelativFishfood <- (RGel*flows[gelTOGelativFish] + RSalp*flows[SalpMortTOGelativFish])/GelativFishcons
  RMyctfood <- (RMeso*flows[MesoTOMyct] + RMacro*flows[MacroTOMyct] + RVMMacro*flows[VMMacroTOMyct] + RSalp*flows[SalpTOMyct] + RAmph*flows[AmphTOMyct])/Myctcons
  
 #Computing the sum of each flow for normalization purposes
  totflowNO3 <- sum(flows[c(Upwelling,NH4toNO3)])
  totflowNH4 <- sum(flows[c(bacTOnh4,hnfTOnh4,micTOnh4,MesoTOnh4,MacroTOnh4,VMMacroTOnh4,gelTOnh4,MyctTOnh4,SalpTOnh4,AmphTOnh4)])
  totflowPico <- sum(flows[c(NH4toPico,NO3toPico)])
  totflowDtm <- sum(flows[c(NH4TOdtm,NO3TOdtm)])
  totflowFlag <- sum(flows[c(NH4TOflag,NO3TOflag,bacTOflag,PicoTOflag)])
  totflowHNF <- sum(flows[c(bacTOhnf,PicoTOhnf,sdetTOhnf)])
  totflowMic <- sum(flows[c(bacTOmic,PicoTOmic,dtmTOmic,flagTOmic,hnfTOmic,sdetTOmic)])
  totflowMeso <- sum(flows[c(dtmTOMeso,flagTOMeso,hnfTOMeso,MesoTOMeso,mdetTOMeso)])
  totflowMacro <- sum(flows[c(dtmTOMacro,flagTOMacro,hnfTOMacro,micTOMacro,MesoTOMacro)])
  totflowVMMacro <- sum(flows[c(dtmTOVMMacro,flagTOVMMacro,hnfTOVMMacro,micTOVMMacro,MesoTOVMMacro)])
  totflowGel <- sum(flows[c(MesoTOgel,MacroTOgel,VMMacroTOgel)])
  totflowSalp <- sum(flows[c(PicoTOSalp,dtmTOSalp,flagTOSalp,hnfTOSalp,micTOSalp)])
  totflowAmph <- sum(flows[c(gelTOAmph,SalpTOAmph)])
  totflowPiscivFish <- sum(flows[c(VMMacroTOPiscivFish,MyctTOPiscivFish,AmphTOPiscivFish)])
  totflowGelativFish <- sum(flows[c(gelTOGelativFish,SalpMortTOGelativFish)])
  totflowMyct <- sum(flows[c(MesoTOMyct,MacroTOMyct,SalpTOMyct,VMMacroTOMyct,AmphTOMyct)])
  totflowBac <- sum(flows[c(domTObac)])
  totflowSdet <- sum(flows[c(PicoTOsdet,flagTOsdet,hnfTOsdet,micTOsdet)])
  totflowMdet <- sum(flows[c(dtmTOmdet,MesoTOmdet)])
  totflowLdet <- sum(flows[c(MacroTOldet,VMMacroTOldet,gelTOldet,AmphTOldet)])
  totflowSalpDet <- sum(flows[c(SalpTOSalpDet)])
  totflowDOM <- sum(flows[c(PicoTOdom,dtmTOdom,flagTOdom,hnfTOdom,micTOdom,MesoTOdom,MacroTOdom,VMMacroTOdom,gelTOdom,MyctTOdom,SalpTOdom,AmphTOdom,sdetTOdom,mdetTOdom,ldetTOdom,SalpDetTOdom)])
  totflowSalpMort <-sum(flows[c(SalpTOSalpMort)])
  totflowExport <- sum(flows[c(MesoTOHTL,MacroTOHTL,VMMacroTOHTL,VMMacroTOdnh4,VMMacroTOddom,gelTOHTL,MyctTOdnh4,MyctTOddom,MyctTOpoop,MyctTObiomass,SalpTOHTL,SalpTOdnh4,SalpTOddom,SalpMortTOSink,AmphTOHTL,PiscivFishTOdnh4,PiscivFishTOddom,PiscivFishTOpoop,PiscivFishTObiomass,GelativFishTOdnh4,GelativFishTOddom,GelativFishTOpoop,GelativFishTObiomass,mdetTOSink,ldetTOSink,SalpDetTOSink)])

#This codifies every entry in the d15N portion of the excel matrix from top left
#to bottom right, column by column
  Aa[col0+NO3_col,Upwelling] <- RupNO3
  Aa[col0+NO3_col,NH4toNO3] <- RNH4+Eps_Nit
  Aa[col0+NO3_col,NO3toPico] <- -(RNO3+Eps_NO3up)
  Aa[col0+NO3_col,NO3TOdtm] <- -(RNO3+Eps_NO3up)
  Aa[col0+NO3_col,NO3TOflag] <- -(RNO3+Eps_NO3up)
  Aa[col0+NH4_col,NH4toNO3] <- -(RNH4+Eps_Nit)
  Aa[col0+NH4_col,NH4toPico] <- -(RNH4+Eps_NH4up)
  Aa[col0+NH4_col,NH4TOdtm] <- -(RNH4+Eps_NH4up)
  Aa[col0+NH4_col,NH4TOflag] <- -(RNH4+Eps_NH4up)
  Aa[col0+NH4_col,bacTOnh4] <- Rbac+Eps_ExcB
  Aa[col0+NH4_col,hnfTOnh4] <- RHNF+Eps_ExcS
  Aa[col0+NH4_col,micTOnh4] <- RMic+Eps_ExcS
  Aa[col0+NH4_col,MesoTOnh4] <- RMeso+Eps_ExcL
  Aa[col0+NH4_col,MacroTOnh4] <- RMacro+Eps_ExcL
  Aa[col0+NH4_col,VMMacroTOnh4] <- RVMMacro+Eps_ExcL
  Aa[col0+NH4_col,gelTOnh4] <- RGel+Eps_ExcL
  Aa[col0+NH4_col,MyctTOnh4] <- RMyct+Eps_ExcL
  Aa[col0+NH4_col,SalpTOnh4] <- RSalp+Eps_ExcS
  Aa[col0+NH4_col,AmphTOnh4] <- RAmph+Eps_ExcL
  Aa[col0+Pico_col,NH4toPico] <- RNH4+Eps_NH4up
  Aa[col0+Pico_col,NO3toPico] <- RNO3+Eps_NO3up
  Aa[col0+Pico_col,PicoTOhnf] <- -RPico
  Aa[col0+Pico_col,PicoTOmic] <- -RPico
  Aa[col0+Pico_col,PicoTOflag] <- -RPico
  Aa[col0+Pico_col,PicoTOSalp] <- -RPico
  Aa[col0+Pico_col,PicoTOdom] <- -RPico
  Aa[col0+Pico_col,PicoTOsdet] <- -RPico
  Aa[col0+Dtm_col,NH4TOdtm] <- RNH4+Eps_NH4up
  Aa[col0+Dtm_col,NO3TOdtm] <- RNO3+Eps_NO3up
  Aa[col0+Dtm_col,dtmTOmic] <- -RDtm
  Aa[col0+Dtm_col,dtmTOMeso] <- -RDtm
  Aa[col0+Dtm_col,dtmTOMacro] <- -RDtm
  Aa[col0+Dtm_col,dtmTOVMMacro] <- -RDtm
  Aa[col0+Dtm_col,dtmTOSalp] <- -RDtm
  Aa[col0+Dtm_col,dtmTOdom] <- -RDtm
  Aa[col0+Dtm_col,dtmTOmdet] <- -RDtm
  Aa[col0+Flag_col,NH4TOflag] <- RNH4+Eps_NH4up
  Aa[col0+Flag_col,NO3TOflag] <- RNO3+Eps_NO3up
  Aa[col0+Flag_col,bacTOflag] <- Rbac
  Aa[col0+Flag_col,PicoTOflag] <- RPico
  Aa[col0+Flag_col,flagTOmic] <- -RFlag
  Aa[col0+Flag_col,flagTOMeso] <- -RFlag
  Aa[col0+Flag_col,flagTOMacro] <- -RFlag
  Aa[col0+Flag_col,flagTOVMMacro] <- -RFlag
  Aa[col0+Flag_col,flagTOSalp] <- -RFlag
  Aa[col0+Flag_col,flagTOdom] <- -RFlag
  Aa[col0+Flag_col,flagTOsdet] <- -RFlag
  Aa[col0+HNF_col,bacTOhnf] <- Rbac
  Aa[col0+HNF_col,PicoTOhnf] <- RPico
  Aa[col0+HNF_col,hnfTOmic] <- -RHNF
  Aa[col0+HNF_col,hnfTOMeso] <- -RHNF
  Aa[col0+HNF_col,hnfTOMacro] <- -RHNF
  Aa[col0+HNF_col,hnfTOVMMacro] <- -RHNF
  Aa[col0+HNF_col,hnfTOSalp] <- -RHNF
  Aa[col0+HNF_col,hnfTOnh4] <- -(RHNF+Eps_ExcS)
  Aa[col0+HNF_col,hnfTOdom] <- -(RHNF+Eps_ExcS)
  Aa[col0+HNF_col,hnfTOsdet] <- -(RHNFfood+Eps_EgS)
  Aa[col0+HNF_col,sdetTOhnf] <- Rsdet
  Aa[col0+Mic_col,bacTOmic] <-  Rbac
  Aa[col0+Mic_col,PicoTOmic] <- RPico
  Aa[col0+Mic_col,dtmTOmic] <- RDtm
  Aa[col0+Mic_col,flagTOmic] <- RFlag
  Aa[col0+Mic_col,hnfTOmic] <- RHNF
  Aa[col0+Mic_col,micTOMeso] <- -RMic
  Aa[col0+Mic_col,micTOMacro] <- -RMic
  Aa[col0+Mic_col,micTOVMMacro] <- -RMic
  Aa[col0+Mic_col,micTOSalp] <- -RMic
  Aa[col0+Mic_col,micTOnh4] <- -(RMic+Eps_ExcS)
  Aa[col0+Mic_col,micTOdom] <- -(RMic+Eps_ExcS)
  Aa[col0+Mic_col,micTOsdet] <- -(RMICfood+Eps_EgS)
  Aa[col0+Mic_col,sdetTOmic] <-  Rsdet
  Aa[col0+Mic_col,mdetTOmic] <-  Rmdet
  Aa[col0+Meso_col,dtmTOMeso] <- RDtm
  Aa[col0+Meso_col,flagTOMeso] <- RFlag
  Aa[col0+Meso_col,hnfTOMeso] <- RHNF
  Aa[col0+Meso_col,micTOMeso] <- RMic
  Aa[col0+Meso_col,MesoTOMacro] <- -RMeso
  Aa[col0+Meso_col,MesoTOVMMacro] <- -RMeso
  Aa[col0+Meso_col,MesoTOgel] <- -RMeso
  Aa[col0+Meso_col,MesoTOMyct] <- -RMeso
  Aa[col0+Meso_col,MesoTOHTL] <- -RMeso
  Aa[col0+Meso_col,MesoTOnh4] <- -(RMeso+Eps_ExcL)
  Aa[col0+Meso_col,MesoTOdom] <- -(RMeso+Eps_ExcL)
  Aa[col0+Meso_col,MesoTOmdet] <- -(RMesofood+Eps_EgL)
  Aa[col0+Meso_col,mdetTOMeso] <- Rmdet
  Aa[col0+Macro_col,dtmTOMacro] <- RDtm
  Aa[col0+Macro_col,flagTOMacro] <- RFlag
  Aa[col0+Macro_col,hnfTOMacro] <- RHNF
  Aa[col0+Macro_col,micTOMacro] <- RMic
  Aa[col0+Macro_col,MesoTOMacro] <- RMeso
  Aa[col0+Macro_col,MacroTOgel] <- -RMacro
  Aa[col0+Macro_col,MacroTOMyct] <- -RMacro
  Aa[col0+Macro_col,MacroTOHTL] <- -RMacro
  Aa[col0+Macro_col,MacroTOnh4] <- -(RMacro+Eps_ExcL)
  Aa[col0+Macro_col,MacroTOdom] <- -(RMacro+Eps_ExcL)
  Aa[col0+Macro_col,MacroTOldet] <- -(RMacrofood+Eps_EgL)
  Aa[col0+VMMacro_col,dtmTOVMMacro] <- RDtm
  Aa[col0+VMMacro_col,flagTOVMMacro] <- RFlag
  Aa[col0+VMMacro_col,hnfTOVMMacro] <- RHNF
  Aa[col0+VMMacro_col,micTOVMMacro] <- RMic
  Aa[col0+VMMacro_col,MesoTOVMMacro] <- RMeso
  Aa[col0+VMMacro_col,VMMacroTOPiscivFish] <- -RVMMacro
  Aa[col0+VMMacro_col,VMMacroTOgel] <- -RVMMacro
  Aa[col0+VMMacro_col,VMMacroTOMyct] <- -RVMMacro
  Aa[col0+VMMacro_col,VMMacroTOHTL] <- -RVMMacro
  Aa[col0+VMMacro_col,VMMacroTOnh4] <- -(RVMMacro+Eps_ExcL)
  Aa[col0+VMMacro_col,VMMacroTOdom] <- -(RVMMacro+Eps_ExcL)
  Aa[col0+VMMacro_col,VMMacroTOldet] <- -(RVMMacrofood+Eps_EgL)
  Aa[col0+VMMacro_col,VMMacroTOdnh4] <- -(RVMMacro+Eps_ExcL)
  Aa[col0+VMMacro_col,VMMacroTOddom] <- -(RVMMacro+Eps_ExcL)
  Aa[col0+Gel_col,MesoTOgel] <- RMeso
  Aa[col0+Gel_col,MacroTOgel] <- RMacro
  Aa[col0+Gel_col,VMMacroTOgel] <- RVMMacro
  Aa[col0+Gel_col,gelTOAmph] <- -RGel
  Aa[col0+Gel_col,gelTOMyct] <- -RGel
  Aa[col0+Gel_col,gelTOGelativFish] <- -RGel
  Aa[col0+Gel_col,gelTOHTL] <- -RGel
  Aa[col0+Gel_col,gelTOnh4] <- -(RGel+Eps_ExcL)
  Aa[col0+Gel_col,gelTOdom] <- -(RGel+Eps_ExcL)
  Aa[col0+Gel_col,gelTOldet] <- -(RGelfood+Eps_EgL)
  Aa[col0+Salp_col,PicoTOSalp] <- RPico
  Aa[col0+Salp_col,dtmTOSalp] <- RDtm
  Aa[col0+Salp_col,flagTOSalp] <- RFlag
  Aa[col0+Salp_col,hnfTOSalp] <- RHNF
  Aa[col0+Salp_col,micTOSalp] <- RMic
  Aa[col0+Salp_col,SalpTOAmph] <- -RSalp
  Aa[col0+Salp_col,SalpTOMyct] <- -RSalp
  Aa[col0+Salp_col,SalpTOnh4] <- -(RSalp+Eps_ExcS)
  Aa[col0+Salp_col,SalpTOdom] <- -(RSalp+Eps_ExcS)
  Aa[col0+Salp_col,SalpTOSalpDet] <- -(RSalpfood+Eps_EgS)
  Aa[col0+Salp_col,SalpTOSalpMort] <- -RSalp
  Aa[col0+Salp_col,SalpTOHTL] <- -RSalp
  Aa[col0+Salp_col,SalpTOdnh4] <- -(RSalp+Eps_ExcS)
  Aa[col0+Salp_col,SalpTOddom] <- -(RSalp+Eps_ExcS)
  Aa[col0+Amph_col,gelTOAmph] <- RGel
  Aa[col0+Amph_col,SalpTOAmph] <- RSalp
  Aa[col0+Amph_col,AmphTOPiscivFish] <- -RAmph
  Aa[col0+Amph_col,AmphTOMyct] <- -RAmph
  Aa[col0+Amph_col,AmphTOHTL] <- -RAmph
  Aa[col0+Amph_col,AmphTOnh4] <- -(RAmph+Eps_ExcL)
  Aa[col0+Amph_col,AmphTOdom] <- -(RAmph+Eps_ExcL)
  Aa[col0+Amph_col,AmphTOldet] <- -(RAmphfood+Eps_EgL)
  Aa[col0+PiscivFish_col,VMMacroTOPiscivFish] <- RVMMacro
  Aa[col0+PiscivFish_col,MyctTOPiscivFish] <- RMyct
  Aa[col0+PiscivFish_col,AmphTOPiscivFish] <- RAmph
  Aa[col0+PiscivFish_col,PiscivFishTOdnh4] <- -(RPiscivFish+Eps_ExcL)
  Aa[col0+PiscivFish_col,PiscivFishTOddom] <- -(RPiscivFish+Eps_ExcL)
  Aa[col0+PiscivFish_col,PiscivFishTOpoop] <- -(RPiscivFishfood+Eps_EgL)
  Aa[col0+PiscivFish_col,PiscivFishTObiomass] <- -RPiscivFish
  Aa[col0+GelativFish_col,gelTOGelativFish] <- RGel
  Aa[col0+GelativFish_col,SalpMortTOGelativFish] <- RSalp
  Aa[col0+GelativFish_col,GelativFishTOdnh4] <- -(RGelativFish+Eps_ExcL)
  Aa[col0+GelativFish_col,GelativFishTOddom] <- -(RGelativFish+Eps_ExcL)
  Aa[col0+GelativFish_col,GelativFishTOpoop] <- -(RGelativFishfood+Eps_EgL)
  Aa[col0+GelativFish_col,GelativFishTObiomass] <- -RGelativFish
  Aa[col0+Myct_col,MesoTOMyct] <- RMeso
  Aa[col0+Myct_col,MacroTOMyct] <- RMacro
  Aa[col0+Myct_col,VMMacroTOMyct] <- RVMMacro
  Aa[col0+Myct_col,gelTOMyct] <- RGel
  Aa[col0+Myct_col,MyctTOPiscivFish] <- -RMyct
  Aa[col0+Myct_col,MyctTOnh4] <- -(RMyct+Eps_ExcL)
  Aa[col0+Myct_col,MyctTOdom] <- -(RMyct+Eps_ExcL)
  Aa[col0+Myct_col,MyctTOdnh4] <- -(RMyct+Eps_ExcL)
  Aa[col0+Myct_col,MyctTOddom] <- -(RMyct+Eps_ExcL)
  Aa[col0+Myct_col,MyctTOpoop] <- -(RMyctfood+Eps_EgL)
  Aa[col0+Myct_col,MyctTObiomass] <- -RMyct
  Aa[col0+Myct_col,SalpTOMyct] <- RSalp
  Aa[col0+Myct_col,AmphTOMyct] <- RAmph
  Aa[col0+bac_col,bacTOflag] <- -Rbac
  Aa[col0+bac_col,bacTOhnf] <- -Rbac
  Aa[col0+bac_col,bacTOmic] <- -Rbac
  Aa[col0+bac_col,bacTOnh4] <- -(Rbac+Eps_ExcB)
  Aa[col0+bac_col,domTObac] <- RDom+Eps_UpB
  Aa[col0+sdet_col,PicoTOsdet] <- RPico
  Aa[col0+sdet_col,flagTOsdet] <- RFlag
  Aa[col0+sdet_col,hnfTOsdet] <- (RHNFfood+Eps_EgS)
  Aa[col0+sdet_col,micTOsdet] <- (RMICfood+Eps_EgS)
  Aa[col0+sdet_col,sdetTOhnf] <- -Rsdet
  Aa[col0+sdet_col,sdetTOmic] <- -Rsdet
  Aa[col0+sdet_col,sdetTOdom] <- -(Rsdet+Eps_Remin)
  Aa[col0+mdet_col,dtmTOmdet] <- RDtm
  Aa[col0+mdet_col,MesoTOmdet] <- (RMesofood+Eps_EgL)
  Aa[col0+mdet_col,mdetTOmic] <- -Rmdet
  Aa[col0+mdet_col,mdetTOMeso] <- -Rmdet
  Aa[col0+mdet_col,mdetTOdom] <- -(Rmdet+Eps_Remin)
  Aa[col0+mdet_col,mdetTOSink] <- -Rmdet
  Aa[col0+ldet_col,MacroTOldet] <- (RMacrofood+Eps_EgL)
  Aa[col0+ldet_col,VMMacroTOldet] <- (RVMMacrofood+Eps_EgL)
  Aa[col0+ldet_col,gelTOldet] <- (RGelfood+Eps_EgL)
  Aa[col0+ldet_col,AmphTOldet] <- (RAmphfood+Eps_EgL)
  Aa[col0+ldet_col,ldetTOdom] <- -(Rldet+Eps_Remin)
  Aa[col0+ldet_col,ldetTOSink] <- -Rldet
  Aa[col0+SalpDet_col,SalpTOSalpDet] <- (RSalpfood+Eps_EgL)
  Aa[col0+SalpDet_col,SalpDetTOdom] <- -(RSalpDet+Eps_Remin)
  Aa[col0+SalpDet_col,SalpDetTOSink] <- -RSalpDet
  Aa[col0+Dom_col,PicoTOdom] <- RPico
  Aa[col0+Dom_col,dtmTOdom] <- RDtm
  Aa[col0+Dom_col,flagTOdom] <- RFlag
  Aa[col0+Dom_col,hnfTOdom] <- (RHNF+Eps_ExcS)
  Aa[col0+Dom_col,micTOdom] <- (RMic+Eps_ExcS)
  Aa[col0+Dom_col,MesoTOdom] <- (RMeso+Eps_ExcL)
  Aa[col0+Dom_col,MacroTOdom] <- (RMacro+Eps_ExcL)
  Aa[col0+Dom_col,VMMacroTOdom] <- (RVMMacro+Eps_ExcL)
  Aa[col0+Dom_col,gelTOdom] <- (RGel+Eps_ExcL)
  Aa[col0+Dom_col,MyctTOdom] <- (RMyct+Eps_ExcL)
  Aa[col0+Dom_col,SalpTOdom] <- (RSalp+Eps_ExcS)
  Aa[col0+Dom_col,AmphTOdom] <- (RAmph+Eps_ExcL)
  Aa[col0+Dom_col,domTObac] <- -(RDom+Eps_UpB)
  Aa[col0+Dom_col,sdetTOdom] <- (Rsdet+Eps_Remin)
  Aa[col0+Dom_col,mdetTOdom] <- (Rmdet+Eps_Remin)
  Aa[col0+Dom_col,ldetTOdom] <- (Rldet+Eps_Remin)
  Aa[col0+Dom_col,SalpDetTOdom] <- (RSalpDet+Eps_Remin)
  Aa[col0+SalpMort_col,SalpTOSalpMort] <- RSalp
  Aa[col0+SalpMort_col,SalpMortTOGelativFish] <- -RSalp
  Aa[col0+SalpMort_col,SalpMortTOSink] <- -RSalp
  Aa[col0+Export_col,MesoTOHTL] <- RMeso
  Aa[col0+Export_col,MacroTOHTL] <- RMacro
  Aa[col0+Export_col,VMMacroTOHTL] <- RVMMacro
  Aa[col0+Export_col,VMMacroTOdnh4] <- (RVMMacro+Eps_ExcL)
  Aa[col0+Export_col,VMMacroTOddom] <- (RVMMacro+Eps_ExcL)
  Aa[col0+Export_col,gelTOHTL] <- RGel
  Aa[col0+Export_col,MyctTOdnh4] <- (RMyct+Eps_ExcL)
  Aa[col0+Export_col,MyctTOddom] <- (RMyct+Eps_ExcL)
  Aa[col0+Export_col,MyctTOpoop] <- (RMyctfood+Eps_EgL)
  Aa[col0+Export_col,MyctTObiomass] <- RMyct
  Aa[col0+Export_col,SalpTOHTL] <- RSalp
  Aa[col0+Export_col,SalpTOdnh4] <- (RSalp+Eps_ExcS)
  Aa[col0+Export_col,SalpTOddom] <- (RSalp+Eps_ExcS)
  Aa[col0+Export_col,SalpMortTOSink] <- RSalp
  Aa[col0+Export_col,AmphTOHTL] <- RAmph
  Aa[col0+Export_col,PiscivFishTOdnh4] <- (RPiscivFish+Eps_ExcL)
  Aa[col0+Export_col,PiscivFishTOddom] <- (RPiscivFish+Eps_ExcL)
  Aa[col0+Export_col,PiscivFishTOpoop] <- (RPiscivFishfood+Eps_EgL)
  Aa[col0+Export_col,PiscivFishTObiomass] <- RPiscivFish
  Aa[col0+Export_col,GelativFishTOdnh4] <- (RGelativFish+Eps_ExcL)
  Aa[col0+Export_col,GelativFishTOddom] <- (RGelativFish+Eps_ExcL)
  Aa[col0+Export_col,GelativFishTOpoop] <- (RGelativFishfood+Eps_EgL)
  Aa[col0+Export_col,GelativFishTObiomass] <- RGelativFish
  Aa[col0+Export_col,mdetTOSink] <- Rmdet
  Aa[col0+Export_col,ldetTOSink] <- Rldet
  Aa[col0+Export_col,SalpDetTOSink] <- RSalpDet
  
  
  
      
  Aa[col0+NO3_col,] <- Aa[col0+NO3_col,]/totflowNO3
  Aa[col0+NH4_col,] <- Aa[col0+NH4_col,]/totflowNH4
  Aa[col0+Pico_col,] <- Aa[col0+Pico_col,]/totflowPico
  Aa[col0+Dtm_col,] <- Aa[col0+Dtm_col,]/totflowDtm
  Aa[col0+Flag_col,] <- Aa[col0+Flag_col,]/totflowFlag
  Aa[col0+HNF_col,] <- Aa[col0+HNF_col,]/totflowHNF
  Aa[col0+Mic_col,] <- Aa[col0+Mic_col,]/totflowMic
  Aa[col0+Meso_col,] <- Aa[col0+Meso_col,]/totflowMeso
  Aa[col0+Macro_col,] <- Aa[col0+Macro_col,]/totflowMacro
  Aa[col0+VMMacro_col,] <- Aa[col0+VMMacro_col,]/totflowVMMacro
  Aa[col0+Gel_col,] <- Aa[col0+Gel_col,]/totflowGel
  Aa[col0+Salp_col,] <- Aa[col0+Salp_col,]/totflowSalp
  Aa[col0+SalpMort_col,] <- Aa[col0+SalpMort_col,]/totflowSalpMort
  Aa[col0+Amph_col,] <- Aa[col0+Amph_col,]/totflowAmph
  Aa[col0+PiscivFish_col,] <- Aa[col0+PiscivFish_col,]/totflowPiscivFish
  Aa[col0+GelativFish_col,] <- Aa[col0+GelativFish_col,]/totflowGelativFish
  Aa[col0+Myct_col,] <- Aa[col0+Myct_col,]/totflowMyct
  Aa[col0+bac_col,] <- Aa[col0+bac_col,]/totflowBac
  Aa[col0+sdet_col,] <- Aa[col0+sdet_col,]/totflowSdet
  Aa[col0+mdet_col,] <- Aa[col0+mdet_col,]/totflowMdet
  Aa[col0+ldet_col,] <- Aa[col0+ldet_col,]/totflowLdet
  Aa[col0+SalpDet_col,] <- Aa[col0+SalpDet_col,]/totflowSalpDet
  Aa[col0+Dom_col,] <- Aa[col0+Dom_col,]/totflowDOM
  Aa[col0+Export_col,] <- Aa[col0+Export_col,]/totflowExport
  
  
  
  for (i in 1:length(wts))  {
    Aa[,i] <- Aa[,i]/wts[i]
  }
  #print(Aa[13,10])
  return(Aa)
}



### Notify Function ###
pb = txtProgressBar(0, 1, style=3)
notify = function(i, n) {
  setTxtProgressBar(pb, i/n)
}
