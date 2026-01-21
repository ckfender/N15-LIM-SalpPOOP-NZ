function [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,Macro_col,...
    VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
    mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumnsC1(Aa,flows,d15N,del15NMeso,del15NMacro,del15NSalp)
  %Assumed fractionation constants from the supplement of Stukel 2018 CRD
  Eps_NO3up = -5;      %Uptake of nitrate, note this was -6 in CRD
  Eps_NH4up = -10;     %Uptake of ammonia
  Eps_Nit = -14;       %Nitrification (NH4->NO3)
  Eps_ExcL = -5;       %Excretion of Mesozoo and Fish
  Eps_EgL = -1.5;      %Egestion of Mesozoo and Fish
  Eps_ExcS = -1;       %Excretion of Protistan Zoos
  Eps_EgS = -1;        %Egestion of Protistan Zoos
  Eps_ExcB = -5;       %Excretion of Bacteria
  Eps_UpB = -1;        %Bacterial uptake of DOM
  Eps_Remin = 0;       %Solubilization of Det to DOM

  %Setting the appropriate entries in d15N0 that correspond to each starting 15N
col0 =  30;  %The entry before the first 15N entry;
NO3_col =  1;
NH4_col =  2;
Pico_col =  3;
Dtm_col =  4;
Flag_col =  5;
HNF_col =  6;
Mic_col =  7;
Meso_col =  8;
Macro_col =  9;
VMMacro_col =  10;
Gel_col =  11;
Salp_col =  12;
Amph_col =  13;
PiscivFish_col =  14;
GelativFish_col =  15;
Myct_col =  16;
bac_col =  17;
sdet_col =  18;
mdet_col =  19;
ldet_col =  20;
SalpDet_col =  21;
Dom_col =  22;
SalpMort_col =  23;
Export_col =  24;
upNO3_col =  25;
