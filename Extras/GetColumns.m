function [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,...
    VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
    mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumns(Aa,flows,d15N,del15NMeso,del15NMacro,del15NSalp)
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
  
col0 =  39;  %The entry before the first 15N entry;
NO3_col =  1;
NH4_col =  2;
Pico_col =  3;
Dtm_col =  4;
Flag_col =  5;
HNF_col =  6;
Mic_col =  7;
Meso_col =  8;
VMMeso_col =  9;
Macro_col =  10;
VMMacro_col =  11;
Gel_col =  12;
Salp_col =  13;
Amph_col =  14;
PiscivFish_col =  15;
GelativFish_col =  16;
Myct_col =  17;
bac_col =  18;
sdet_col =  19;
mdet_col =  20;
ldet_col =  21;
SalpDet_col =  22;
Dom_col =  23;
SalpMort_col =  24;
Export_col =  25;
upNO3_col =  26;

