clear all
close all

%%Initialization of excel mass balance/constraint files and input files into
%%readable matrices to be run in R

%%Set input file to be read in, in this case which cycle
Cycle=5;

%%first, get the name of spreadsheet with data
sheet=['N15InverseModelNZ.xlsx'];
workbook='NZ';

%%now start reading in the data, starting with sizes
datsize=xlsread(sheet,workbook,'B2:E2');
neeq=datsize(1);   %Number of exact equalities (mass balances)
naeq=datsize(2);  %Number of approximate equalities (well known Inputs with associated uncertainty + 15N data)
ngt0=datsize(3);  %Number of inequalities (> < relationships like GGE, AE, or max ingestion rate) 
nvar=datsize(4);  %Number of flows (i.e. rows in matrix)

ngt=ngt0+nvar;		% total no. ineqs, including all flows >0
neq=neeq+naeq;      % total number of equalities

rskip=5;    %rows to skip before data
cskip=8;    %cols to skip before data
ccl=cskip+1;	% col of upper left number to be read
crl=rskip+1;    % row of upper left number to be read
crr=rskip+nvar+1; %row of bottom right number to be read
ccr=cskip+neq+ngt0+2;  %col of bottom right number to be read
crr2=rskip+nvar+2;
cread=[char(ExcelCol(ccl)),num2str(crl),':',char(ExcelCol(ccr)),num2str(crr)]; %Excel code for data to be read in
M0=xlsread(sheet,workbook,cread); %The actual input matrix
if Cycle==2
    cread=[char(ExcelCol(ccl-4)),num2str(crl),':',char(ExcelCol(ccl-4)),num2str(crr-1)];
elseif Cycle==3
    cread=[char(ExcelCol(ccl-3)),num2str(crl),':',char(ExcelCol(ccl-3)),num2str(crr-1)];
elseif Cycle==4
    cread=[char(ExcelCol(ccl-2)),num2str(crl),':',char(ExcelCol(ccl-2)),num2str(crr-1)];
elseif Cycle==5
    cread=[char(ExcelCol(ccl-1)),num2str(crl),':',char(ExcelCol(ccl-1)),num2str(crr-1)];
end
wts=xlsread(sheet,workbook,cread); %The weights column of the input matrix

% % now, sort out
% % do some minimal processing
M=[M0(:,1:(neq+ngt0)),eye(nvar+1,nvar)]; %Expands M0 by adding columns to hold uncertainties. Adds a diagonal line of 1s to the M0 matrix of size=rows
A=M(1:nvar,1:neq)'; %The exact and approximate equalities, transposed
G=M(1:nvar,neq+1:neq+ngt)'; %Just the inequalities, transposed
b=M(nvar+1,1:neq)'; %The b values to which each equality or inequality is equal (0 for mass balance for example)
Ae=A(1:neeq,:); %The exact equalities, transposed
Aa=A(neeq+1:neeq+naeq,:); %The approximate equalities, transposed
be=b(1:neeq); %The values to which the exact equalities are equal (all 0)
ba=b(neeq+1:neeq+naeq); %The values to which the approximate equalities (rate measurements) are equal
sdba=NaN(size(ba)); %Empty vector to hold the sd of the ba

Inputs=xlsread('InputsNZ.xlsx','InputsNZ','D3:M41');    %loading ba and sdba
col=Cycle*2-1;  %column of Inputs (in matlab) corresponding to Cycle run

CN=106/16; %Redfield C:N ratio

%Approximate equalities from Inputs file
ba(1)=Inputs(1,Cycle*2-1);  %NPP
ba(2)=Inputs(2,Cycle*2-1);  %Sedtrap Flux (70m)
ba(3)=Inputs(3,Cycle*2-1);  %Pico NPP
ba(4)=Inputs(4,Cycle*2-1);  %Flag NPP
ba(5)=Inputs(5,Cycle*2-1);  %Diatom NPP
ba(6)=Inputs(6,Cycle*2-1);  %Pico grazing mort
ba(7)=Inputs(7,Cycle*2-1);  %Flag grazing mort
ba(8)=Inputs(8,Cycle*2-1);  %Dtm grazing mort
ba(9)=Inputs(9,Cycle*2-1);  %Protistan grazing
ba(10)=Inputs(10,Cycle*2-1);  %Meso grazing
ba(11)=Inputs(11,Cycle*2-1);  %VMMeso grazing
ba(12)=Inputs(12,Cycle*2-1);  %Macro grazing
ba(13)=Inputs(13,Cycle*2-1);  %VMMacro grazing
ba(14)=Inputs(14,Cycle*2-1);  %Salp Auto grazing
ba(15)=Inputs(15,Cycle*2-1);  %Salp Hetero grazing
ba(16)=Inputs(16,Cycle*2-1);   %Salp Pico grazing
ba(17)=Inputs(17,Cycle*2-1);   %Salp dtm grazing
ba(18)=Inputs(18,Cycle*2-1);   %Other salp grazing
ba(19)=Inputs(19,Cycle*2-1);   %Salp fecal pellet production
ba(20)=Inputs(20,Cycle*2-1);   %Bacterial production
ba(21)=Inputs(21,Cycle*2-1);   %Pinkerton Myct-Meso Consumption
ba(22)=Inputs(22,Cycle*2-1);   %Pinkerton Myct-VMMeso Contumption
ba(23)=Inputs(23,Cycle*2-1);   %Pinkerton Myct-VMMacro Consumption
ba(24)=Inputs(24,Cycle*2-1);   %Pinkerton Myct-Amph Consumption
ba(25)=Inputs(25,Cycle*2-1);   %Pinkerton Myct-Salp Consumption
ba(26)=Inputs(26,Cycle*2-1);   %Pinkerton Myct-Gel Consumption
ba(27)=Inputs(27,Cycle*2-1);   %Pinkerton Pisc-Myct Consumption
ba(28)=Inputs(28,Cycle*2-1);   %Pinkerton Pisc-Macro Consumption
ba(29)=Inputs(29,Cycle*2-1);   %Pinkerton Pisc-VMMacro Consumption
ba(30)=Inputs(30,Cycle*2-1);   %Pinkerton Pisc-Amph Consumption
ba(31)=Inputs(31,Cycle*2-1);   %Pinkerton Gelativ-Salp Consumption
ba(32)=Inputs(32,Cycle*2-1);   %Pinkerton Gelativ-Gel Consumption
ba(33)=Inputs(33,Cycle*2-1);   %Pinkerton HTL-Meso Consumption
ba(34)=Inputs(34,Cycle*2-1);   %Pinkerton HTL-VMMeso Consumption
ba(35)=Inputs(35,Cycle*2-1);   %Pinkerton HTL-Macro Consumption
ba(36)=Inputs(36,Cycle*2-1);   %Pinkerton HTL-VMMacro Consumption
ba(37)=Inputs(37,Cycle*2-1);   %Pinkerton HTL-Amph Consumption
ba(38)=Inputs(38,Cycle*2-1);   %Pinkerton HTL-Salp Consumption
ba(39)=Inputs(39,Cycle*2-1);   %Pinkerton HTL-Gel Consumption


sdba(1)=Inputs(1,Cycle*2);  %NPP
sdba(2)=Inputs(2,Cycle*2);  %Sedtrap Flux (70m)
sdba(3)=Inputs(3,Cycle*2);  %Pico NPP
sdba(4)=Inputs(4,Cycle*2);  %Flag NPP
sdba(5)=Inputs(5,Cycle*2);  %Diatom NPP
sdba(6)=Inputs(6,Cycle*2);  %Pico grazing mort
sdba(7)=Inputs(7,Cycle*2);  %Flag grazing mort
sdba(8)=Inputs(8,Cycle*2);  %Dtm grazing mort
sdba(9)=Inputs(9,Cycle*2);  %Protistan grazing
sdba(10)=Inputs(10,Cycle*2);  %Meso grazing
sdba(11)=Inputs(11,Cycle*2);  %VMMeso grazing
sdba(12)=Inputs(12,Cycle*2);  %Macro grazing
sdba(13)=Inputs(13,Cycle*2);  %VMMacro grazing
sdba(14)=Inputs(14,Cycle*2);  %Salp Auto grazing
sdba(15)=Inputs(15,Cycle*2);  %Salp Hetero grazing
sdba(16)=Inputs(16,Cycle*2);   %Salp Pico grazing
sdba(17)=Inputs(17,Cycle*2);   %Salp dtm grazing
sdba(18)=Inputs(18,Cycle*2);   %Other salp grazing
sdba(19)=Inputs(19,Cycle*2);   %Salp fecal pellet production
sdba(20)=Inputs(20,Cycle*2);   %Bacterial production
sdba(21)=Inputs(21,Cycle*2);   %Pinkerton Myct-Meso Consumption
sdba(22)=Inputs(22,Cycle*2);   %Pinkerton Myct-VMMeso Contumption
sdba(23)=Inputs(23,Cycle*2);   %Pinkerton Myct-Macro Consumption
sdba(24)=Inputs(24,Cycle*2);   %Pinkerton Myct-VMMacro Consumption
sdba(25)=Inputs(25,Cycle*2);   %Pinkerton Myct-Amph Consumption
sdba(26)=Inputs(26,Cycle*2);   %Pinkerton Myct-Salp Consumption
sdba(27)=Inputs(27,Cycle*2);   %Pinkerton Myct-Gel Consumption
sdba(28)=Inputs(28,Cycle*2);   %Pinkerton Pisc-Myct Consumption
sdba(29)=Inputs(29,Cycle*2);   %Pinkerton Pisc-VMMacro Consumption
sdba(30)=Inputs(30,Cycle*2);   %Pinkerton Pisc-Amph Consumption
sdba(31)=Inputs(31,Cycle*2);   %Pinkerton Gelativ-Salp Consumption
sdba(32)=Inputs(32,Cycle*2);   %Pinkerton Gelativ-Gel Consumption
sdba(33)=Inputs(33,Cycle*2);   %Pinkerton HTL-Meso Consumption
sdba(34)=Inputs(34,Cycle*2);   %Pinkerton HTL-VMMeso Consumption
sdba(35)=Inputs(35,Cycle*2);   %Pinkerton HTL-Macro Consumption
sdba(36)=Inputs(36,Cycle*2);   %Pinkerton HTL-VMMacro Consumption
sdba(37)=Inputs(37,Cycle*2);   %Pinkerton HTL-Amph Consumption
sdba(38)=Inputs(38,Cycle*2);   %Pinkerton HTL-Salp Consumption
sdba(39)=Inputs(39,Cycle*2);   %Pinkerton HTL-Gel Consumption
sdba(40:naeq)=NaN; %The d15N approximate equalities have unknown uncertainty
 
b=[be;ba];   %Recombining b vectors now that they have been set up

InEqualInputs = xlsread('InputsNZ.xlsx','InputsNZ','D44:M76');

h=M(nvar+1,neq+1:neq+ngt)'; %h is the equivalent of b for the equalities, the value they are greater or less than
%%Note that positive values denote minimums, negative maximums. The
%%following indexing should correspond to column the inequality falls in the
%%excel sheet (so h(1) is Salp FP Flux and h(2) is NH4 Excretion->20% ingestion)
Temp1=InEqualInputs(1,col);
Temp2=InEqualInputs(2,col);
Temp3=InEqualInputs(3,col);
BotDepth=InEqualInputs(4,col);
SalpFPFlux=InEqualInputs(5,col);
Picobiomass=InEqualInputs(6,col);
Diatombiomass=InEqualInputs(7,col);
Flagbiomass=InEqualInputs(8,col);
HNFbiomass=InEqualInputs(9,col);
MICbiomass=InEqualInputs(10,col);
MesoDaybiomass=InEqualInputs(11,col);
MesoNightbiomass=InEqualInputs(12,col);
MacroDaybiomass=InEqualInputs(13,col);
MacroNightbiomass=InEqualInputs(14,col);
Amphbiomass=InEqualInputs(15,col);
Gelbiomass=InEqualInputs(16,col);
SalpDaybiomass=InEqualInputs(17,col);
SalpNightbiomass=InEqualInputs(18,col);
Myctbiomass=InEqualInputs(19,col);
PiscivFishbiomass=InEqualInputs(20,col);
GelativFishbiomass=InEqualInputs(21,col);
HNFsize=InEqualInputs(22,col);
Flagsize=InEqualInputs(23,col);
Microsize=InEqualInputs(24,col);
Mesosize=InEqualInputs(25,col);
Macrosize=InEqualInputs(26,col);
Amphsize=InEqualInputs(27,col);
GelPredsize=InEqualInputs(28,col);
Salpsize=InEqualInputs(29,col);
Myctsize=InEqualInputs(30,col);
PiscivFishsize=InEqualInputs(31,col);
GelativFishsize=InEqualInputs(32,col);
SalpResp=InEqualInputs(33,col);

%Minimum salp fecal pellet flux as that identified via microscopy
h(1)=SalpFPFlux;
%The following are minimum respiration (ug N ind-1 h-1) from Ikeda (2014) - they
%parameterize ammonium excretion from weight (in units of mg ind-1), depth,
%and temperature. Final units are (mmol N m-2 d-1)
%Copepods, Euphausiids, and Amphipods were found not to be significantly
%different. Salps are not currently included due to plenty of other
%constraints but they easily could be.
Depth1=70; %Euphotic zone depth, depth at which daytime excretion is done
SurfTemp=Temp1+273.15; %All Ikeda temps are in Kelvin
CN=6.625;
%Non-DVM biomass experiences surface temps the entire day, DVM portion of
%biomass experiences surface temp half the day and a different temp the
%other half (depends on migration depth). Night BM-DayBM=DVMBiomass
a0=17.945;
a1=0.771;
a2=-5.528;
a3=-0.082;
CW=Mesosize;
Depth2=200; %VM depth
%For nonmigrating mesozoos
DeepTemp=Temp2+273.15;
DVMBiomass=MesoNightbiomass-MesoDaybiomass; %For DVMers, only nighttime respiration contributes to NH4 flow into system, daytime respiration is a sink
NonDVMBiomass=MesoDaybiomass;   %All hours are spent in the euphotic and contribute NH4
NonDVMExc=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  %NH4 exc in ug ind-1 h-1
NonDVMfactor=24/14*NonDVMExc/CW;  %Multiply by 24 hours for h-1 to d-1, divide by 14 g/mol for ugN to umolN, divide by carbon weight, so this has units of umol N / day / mg C
h(14)=NonDVMfactor*(NonDVMBiomass*CN*12)/1000; %Need to multiply by CN*12 for per mg C to mmolN because biomass is in units of mmol N (so umol N / day / mmolN).  Need to divide by 1000 to go from umol N / day to mmol N / day
h(15) = -6*h(14); %max resp is 3x the above
DVMExcNight=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  
DVMfactorNight=12/14*DVMExcNight/CW;  
DVMfinalNight=DVMfactorNight*(DVMBiomass*CN*12)/1000; 
DVMExcDay=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/DeepTemp))*exp(a3*log(Depth2)); 
DVMfactorDay=12/14*DVMExcDay/CW;  
DVMfinalDay=DVMfactorDay*(DVMBiomass*CN*12)/1000; 
h(16)=DVMfinalDay+DVMfinalNight;
h(17) = -6*h(16);
%For Macro and VMMacro
CW=Macrosize;
Depth2=200;
DeepTemp=Temp2+273.15;
DVMBiomass=MacroNightbiomass-MacroDaybiomass; 
NonDVMBiomass=MacroDaybiomass;   
NonDVMExc=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  
NonDVMfactor=24/14*NonDVMExc/CW;  
h(18)=NonDVMfactor*(NonDVMBiomass*CN*12)/1000; 
h(19) = -6*h(18);
DVMExcNight=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  
DVMfactorNight=12/14*DVMExcNight/CW;  
DVMfinalNight=DVMfactorNight*(DVMBiomass*CN*12)/1000; 
DVMExcDay=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/DeepTemp))*exp(a3*log(Depth2)); 
DVMfactorDay=12/14*DVMExcDay/CW;  
DVMfinalDay=DVMfactorDay*(DVMBiomass*CN*12)/1000; 
h(20)=DVMfinalDay+DVMfinalNight;
h(21) = -6*h(20);
%For Amphipods
CW=Amphsize;
NonDVMBiomass=Amphbiomass;   
NonDVMExc=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  
NonDVMfactor=24/14*NonDVMExc/CW;  
h(22)=NonDVMfactor*(NonDVMBiomass*CN*12)/1000; 
h(23) = -6*h(22);
%GelPred, cnidarians and ctenophores have additional coefficients
a10=0.655;
a11=0.878;
CW=GelPredsize;
NonDVMBiomass=Gelbiomass;   
NonDVMExc=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1))*exp(a10)*exp(a11);  
NonDVMfactor=24/14*NonDVMExc/CW;  
NonDVMfinal=NonDVMfactor*(NonDVMBiomass*CN*12)/1000; 
h(24)=NonDVMfactor*(NonDVMBiomass*CN*12)/1000; 
h(25) = -6*h(24);

%Minimum respiration for fish is derived from Ikeda (2016)'s oxygen
%respiration equation using wet weight and then converting to ammonia
%excretion using his median O:N relationship
%Model 1 for all fish
a0=19.491;
a1=0.885;
a2=-5.770;
a3=-0.261;
DeepTemp=Temp2+273.15;
WW=Myctsize;
DVMRespNight=exp(a0)*exp(a1*log(WW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  %O2 resp in uL ind-1 h-1
DVMExcNight=DVMRespNight/22.4/24.2;   %uL to umols, umols O2 to umols N using Ikeda 2016 mean O:N of 24.2
DVMfactorNight=12*DVMExcNight/(WW*.1);  %Multiply by 12 hours, divide by wet weight for umol N / day/ mg C, assume WW is 10% of CW (Pinkerton , 2011)
DVMfinalNight=DVMfactorNight*(Myctbiomass*CN*12)/1000; %Need to multiply by CN*12 for per mg C to mmolN because biomass is in units of mmol N (so umol N / day / mmolN).  Need to divide by 1000 to go from umol N / day to mmol N / day
DVMRespDay=exp(a0)*exp(a1*log(WW))*exp(a2*(1000/DeepTemp))*exp(a3*log(Depth2)); 
DVMExcDay=DVMRespDay/22.4/24.2;
DVMfactorDay=12*DVMExcDay/(WW*.1);  
DVMfinalDay=DVMfactorDay*(Myctbiomass*CN*12)/1000; 
h(26)=DVMfinalDay+DVMfinalNight;

%Piscivorous fish
DeepTemp=Temp3+273.15;
WW=PiscivFishsize;
Resp=exp(a0)*exp(a1*log(WW))*exp(a2*(1000/DeepTemp))*exp(a3*log(BotDepth));  
Exc=Resp/22.4/24.2;  
factor=24*Exc/(WW*.1);  
h(27)=factor*(PiscivFishbiomass*CN*12)/1000;

%Salp minimum respiration is calculated using Ikeda 2014
a0=17.945;
a1=0.771;
a2=-5.528;
a3=-0.082;
a13=0.871;
CW=Salpsize;
Depth2=200;
DeepTemp=Temp2+273.15;
DVMBiomass=SalpNightbiomass-SalpDaybiomass; 
NonDVMBiomass=SalpDaybiomass;   
NonDVMExc=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1))*exp(a13);    
NonDVMfactor=24/14*NonDVMExc/CW;  
NonDVMfinal=NonDVMfactor*(NonDVMBiomass*CN*12)/1000; 
DVMExcNight=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/SurfTemp))*exp(a3*log(Depth1));  
DVMfactorNight=12/14*DVMExcNight/CW;  
DVMfinalNight=DVMfactorNight*(DVMBiomass*CN*12)/1000; 
DVMExcDay=exp(a0)*exp(a1*log(CW))*exp(a2*(1000/DeepTemp))*exp(a3*log(Depth2)); 
DVMfactorDay=12/14*DVMExcDay/CW;  
DVMfinalDay=DVMfactorDay*(DVMBiomass*CN*12)/1000; 
h(28)=NonDVMfinal+DVMfinalDay+DVMfinalNight;
h(29) = -6*h(28);

%Constraining respiration of vertical migrators to be higher in the surface
%than at depth, effectively inputting the typed cells in the Excel matrix
VMMesoexc_col=30;
VMMacroexc_col=31;
Salpexc_col=32;
Myctexc_col=33;
VMMesoTOnh4 = 69;
VMMesoTOdnh4 = 72;
VMMacroTOnh4 = 84;
VMMacroTOdnh4 = 87;
SalpTOnh4 = 105;
SalpTOdnh4 = 110;
MyctTOnh4 = 97;
MyctTOdnh4 = 99;
a2=-5.528;
G(VMMesoexc_col,VMMesoTOnh4)=exp(a2*(1000/(Temp2+273.15)));
G(VMMesoexc_col,VMMesoTOdnh4)=-exp(a2*(1000/(Temp1+273.15)));
G(VMMacroexc_col,VMMacroTOnh4)=exp(a2*(1000/(Temp2+273.15)));
G(VMMacroexc_col,VMMacroTOdnh4)=-exp(a2*(1000/(Temp1+273.15)));
G(Salpexc_col,SalpTOnh4)=exp(a2*(1000/(Temp2+273.15)));
G(Salpexc_col,SalpTOdnh4)=-exp(a2*(1000/(Temp1+273.15)));
a2=-5.770;
G(Myctexc_col,MyctTOnh4)=exp(a2*(1000/(Temp2+273.15)));
G(Myctexc_col,MyctTOdnh4)=-exp(a2*(1000/(Temp1+273.15)));

VMMesoexcdom_col=67;
VMMacroexcdom_col=68;
Salpexcdom_col=69;
Myctexcdom_col=70;
VMMesoTOdom = 70;
VMMesoTOddom = 73;
VMMacroTOdom = 85;
VMMacroTOddom = 88;
SalpTOdom = 106;
SalpTOddom = 111;
MyctTOdom = 98;
MyctTOddom = 100;
a2=-5.528;
G(VMMesoexcdom_col,VMMesoTOdom)=exp(a2*(1000/(Temp2+273.15)));
G(VMMesoexcdom_col,VMMesoTOddom)=-exp(a2*(1000/(Temp1+273.15)));
G(VMMacroexcdom_col,VMMacroTOdom)=exp(a2*(1000/(Temp2+273.15)));
G(VMMacroexcdom_col,VMMacroTOddom)=-exp(a2*(1000/(Temp1+273.15)));
G(Salpexcdom_col,SalpTOdom)=exp(a2*(1000/(Temp2+273.15)));
G(Salpexcdom_col,SalpTOddom)=-exp(a2*(1000/(Temp1+273.15)));
a2=-5.770;
G(Myctexcdom_col,MyctTOdom)=exp(a2*(1000/(Temp2+273.15)));
G(Myctexcdom_col,MyctTOddom)=-exp(a2*(1000/(Temp1+273.15)));


AllPhytoNPP = Inputs(1,col);
PicoNPP = Inputs(3,col);
FlagNPP = Inputs(4,col);
DtmNPP = Inputs(5,col);

%The following are max ingestion constraints from Hansen. Note that the
%temps are all for the surface, I believe under the assumption that the
%majority of feeding happens within the euphotic zone. If you knew how much
%time was spent at depth you could split these into two equations with that
%time instead of multiplying by 24 hours.
a6 = 0.1+0.24;  %THis is the value for protozooplankton from Hansen et al (1997, p. 11) plus one s.e.
b6 = -0.20;   %This is the value for protozooplankton from Hansen et al. (1997, p. 11) )

%Mixotrophic Flagellates
PredVol=Flagsize;  %um^3,   Average volume
%Multiplying by 24 to go from h-1 to d-1. The rest comes from the equation
%log(Imax)=log(I0)+log(Q10)*((t-t0)/10) where t is the temperature, and t0
%is the reference temperature for the Q10 being used. Using the 
%overall average Q10 of 2.8 at 20 C.
%Since Imax is specific ingestion it described the prey mass ingested per 
%predator mass per hour, thus it just needs to be multiplied by the observed 
%predator biomass for realized prey ingestion.
h(123)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*Flagbiomass;
%Heterotrophic Nanoflagellates
PredVol=HNFsize;
h(124)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*HNFbiomass;
%Microzoos
PredVol=Microsize;
h(125)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*MICbiomass;

a6 = -1.78+0.47;  %This is the value for metazooplankton from Hansen et al (1997, p. 11) plus one s.e.
b6 = 0.08;   %This is the value for metazooplankton from Hansen et al. (1997, p. 11))
%Meso
DWC=0.4; %Dry weight to carbon conversion
WWDW = 0.2;  %Wet weight to Dryweight conversion
PredVol=Mesosize/DWC/WWDW*10^9;   %10^9 goes from mm^3 to um^3. 1 gram of water has a volume of 1 mL or 1 cm^3 so 1 mg = 1 mm^3
h(126)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*MesoDaybiomass;
%VMMeso
DWC=0.4; %Dry weight to carbon conversion
WWDW = 0.2;  %Wet weight to Dryweight conversion
PredVol=Mesosize/DWC/WWDW*10^9;   
h(127)=-12*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*(MesoNightbiomass-MesoDaybiomass);
%Macro
DWC=0.4; %Dry weight to carbon conversion
WWDW = 0.2;  %Wet weight to Dryweight conversion
PredVol=Macrosize/DWC/WWDW*10^9;   %10^9 goes from mm^3 to um^3
h(128)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*MacroDaybiomass;
%VMMacro
DWC=0.4; %Dry weight to carbon conversion
WWDW = 0.2;  %Wet weight to Dryweight conversion
PredVol=Mesosize/DWC/WWDW*10^9;   
h(129)=-12*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*(MacroNightbiomass-MacroDaybiomass);
%Amphipod
DWC=0.4; %Dry weight to carbon conversion
WWDW = 0.2;  %Wet weight to Dryweight conversion
PredVol=Amphsize/DWC/WWDW*10^9;   %10^9 goes from mm^3 to um^3
h(130)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*Amphbiomass;
%Gelatinous Predators
DWC=0.4; %Dry weight to carbon conversion
WWDW = 0.2;  %Wet weight to Dryweight conversion
PredVol=GelPredsize/DWC/WWDW*10^9;   %10^9 goes from mm^3 to um^3
h(131)=-24*10^(a6+b6*log10(PredVol)+log10(2.8)*(Temp1-20)/10)*Gelbiomass;


%Capping max autotroph production to double the biomass
h(132)=-2*Picobiomass;  %Cyano max production at double biomass
h(133)=-2*Diatombiomass;  %Diatom max production at double biomass
h(134)=-2*Flagbiomass;  %Flag max production at double biomass


%Removing remaining NaN entries with missing equalities
inds=find(isnan(h));
G(inds,:)=[];
h(inds)=[];
inds=find(isnan(ba));
Aa(inds,:)=[];
ba(inds)=[];
sdba(inds)=[];



A0 = A;
G0 = G;
Ae0=Ae;
Aa0=Aa;
for i=1:length(wts)
    A(:,i)=A(:,i)/wts(i);
    Aa(:,i)=Aa(:,i)/wts(i);
    Ae(:,i)=Ae(:,i)/wts(i);
    G(:,i)=G(:,i)/wts(i);
end


d15NInputs=xlsread('InputsNZ.xlsx','InputsNZ','D79:M81');



%Just clearing variables that we don't need to see
clear Araw
clear Awt
clear Graw
clear M
clear M0
clear captions
clear cc2
clear cc3
clear ccl
clear ccr
clear chsheet
clear comment
clear constraintout
clear constraints
clear cread
clear cread2
clear cread3
clear cread4
clear crl
clear crr
clear crr2
clear crr3
clear crr4
clear cskip
clear datsize
clear equalities
clear equalitiesout
clear linsum
clear neq
clear naeq
clear neeq
clear ngt
clear ngt0
clear norm2
clear npar
clear nvar
clear outgate
clear rc1
clear cnorm
clear rr
clear rrout
clear rrtrans
clear rskip
clear sheet
clear sols
clear sqrterror
clear temp
clear time1
clear weight
clear SD


InputCol = col;

%Saving matrices to .mat file
if length(ba)>0
    save(['N15NZInverseCycle',num2str(Cycle),'.mat'],'A','Ae','Aa','G','b','be','ba','h','Inputs','sdba','InputCol','d15NInputs','wts','A0','Ae0','Aa0','G0')
else
    save([sheetp,'.',Model,'.mat'],'A','G','b','h','wts')
end