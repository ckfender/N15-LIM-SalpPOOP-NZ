clear all
close all

%%Script to create a figure summarizing organismal size, trophic levels, 
%%and magnitude ofpredation
%Loading trophic level file produced by PlotTL as well as Approximate
%Equality Inputs file
TL2 = readtable('TL.xlsx','Sheet','Mean');
Biom=xlsread('InputsNZ.xlsx','InputsNZ','D54:M57');
nvmMesoBiom = Biom(1,[1,3,5,7,9])';
vmMesoBiom = Biom(2,[1,3,5,7,9])'-Biom(1,[1,3,5,7,9])';
nvmMacroBiom = Biom(3,[1,3,5,7,9])';
vmMacroBiom = Biom(4,[1,3,5,7,9])'-Biom(3,[1,3,5,7,9])';
TL2(:,["NO3" "NH4" "sdet" "mdet" "ldet" "SalpDet" "Dom" "SalpMort"])=[];
%Representing TL for Meso and Macro as mean of the migrant and resident
%populations scaled to their biomasses
TL2{2:5,"nvmMeso"} = (TL2{2:5,"nvmMeso"}.*nvmMesoBiom(2:5) +...
    TL2{2:5,"VMMeso"}.*vmMesoBiom(2:5))./(nvmMesoBiom(2:5)+vmMesoBiom(2:5));
TL2{1:5,"nvmMacro"} = (TL2{1:5,"nvmMacro"}.*nvmMacroBiom +...
    TL2{1:5,"vmMacro"}.*vmMacroBiom)./(nvmMacroBiom+vmMacroBiom);
TL2(:,["VMMeso","vmMacro"])=[];
TL2 = renamevars(TL2,["nvmMeso","nvmMacro"],["Meso","Macro"]);

SalpTL = table2array(mean(TL2([1 2 4],:),1));
NonSalpTL = table2array(mean(TL2([3 5],:),1));
%Vector of average sizes for each entry in um
Size=[1.5,20,10,10,40,2500,5000,10000,40000,25000,1300000,350000,50000,0.5];
Size=log10(Size);
Labs=TL2.Properties.VariableNames;

load("prodmean.mat"); %Loading production matrix produced by PlotIndirect.m
IngS = mean(prodmean(:,:,[1 2 4]),3);
IngS(:,8)=IngS(:,8)+IngS(:,9);
IngS(:,10)=IngS(:,10)+IngS(:,11);
IngS(8,:)=IngS(8,:)+IngS(8,:);
IngS(10,:)=IngS(10,:)+IngS(11,:);
IngS(:,[9 11])=[];
IngS([9 11],:)=[];
IngS(:,[1 2 17:(end-1)])=[]; %removing non-organisms 
IngS([1 2 17:(end-1)],:)=[]; %removing non-organisms 
IngN = mean(prodmean(:,:,[3 5]),3);
IngN(:,8)=IngN(:,8)+IngN(:,9);
IngN(:,10)=IngN(:,10)+IngN(:,11);
IngN(8,:)=IngN(8,:)+IngN(8,:);
IngN(10,:)=IngN(10,:)+IngN(11,:);
IngN(:,[9 11])=[];
IngN([9 11],:)=[];
IngN(:,[1 2 17:(end-1)])=[]; %removing non-organisms 
IngN([1 2 17:(end-1)],:)=[]; %removing non-organisms 
ProdS = log10(sum(IngS(1:14,:),2));
ProdN = log10(sum(IngN(1:14,:),2));
%%Transformation to make arrows more distinct
% shift=min([IngS(find(IngS>0)); IngN(find(IngN>0))])*.95;
% IngS=log10(IngS./shift);
% IngN=log10(IngN./shift);
% IngS = log10(IngS./shift)./log10(50); %log base 50
% IngN = log10(IngN./shift)./log10(50);
%%Setting salp values to be relative to nonsalp 'control'
IngSnorm=log10(IngS./IngN)./log10(50);

%Generating a vector of color values for Prod sampled from a continous palette
% Normalize values to the range [0, 1], where 0 is the lowest log10 value
% the colorbar will represent and 1 is the largest
Min = floor(min([ProdS; ProdN]));
Max = ceil(max([ProdS; ProdN]));
ProdSNorm = (ProdS - Min) / (Max - Min);
ProdNNorm = (ProdN - Min) / (Max - Min);
% Define a colormap and number of colors
num_colors = 1000; % Adjust as needed for smoothness
cmap = parula(num_colors); % Use your desired colormap
% Map normalized values to the colormap
color_indices1 = round(ProdSNorm * (num_colors - 1)) + 1;
color_indices2 = round(ProdNNorm * (num_colors - 1)) + 1;
colorsS = cmap(color_indices1, :);
colorsN = cmap(color_indices2, :);
% (Optional) Plot to visualize colors
figure;
scatter(1:length(ProdS), ProdS, 100, colorsS, 'filled');
xlabel('Index');
ylabel('Value');
title('Values with Corresponding Colors');
cb=colorbar;
colormap(cmap);
% Customize tick locations and labels
num_ticks = 6; % Number of ticks you want
tick_values = linspace(Min, Max, num_ticks); % Tick values
tick_positions = (tick_values - Min) / (Max - Min); % Normalized positions
% Update colorbar properties
cb.Ticks = tick_positions; % Set tick locations
cb.TickLabels = {'10^-^4', '10^-^3', '10^-^2', '10^-^1', '10^0', '10^1'}; % Set tick labels
cb.Label.String = 'Value'; % Add a label to the colorbar

Pico_col=1;
Dtm_col=2;
Flag_col=3;
HNF_col=4;
Mic_col=5;
Meso_col=6;
Macro_col=7;
Gel_col=8;
Salp_col=9;
Amph_col=10;
Pisciv_col=11;
Gelativ_col=12;
Myct_col=13;
Bac_col=14;
HTL_col=15;

%Saving x-y locations of each organism
PicoS=[Size(1) SalpTL(1)];
DtmS=[Size(2) SalpTL(2)];
FlagS=[Size(3) SalpTL(3)];
HNFS=[Size(4) SalpTL(4)];
MicS=[Size(5) SalpTL(5)];
MesoS=[Size(6) SalpTL(6)];
MacroS=[Size(7) SalpTL(7)];
GelS=[Size(8) SalpTL(8)];
SalpS=[Size(9) SalpTL(9)];
AmphS=[Size(10) SalpTL(10)];
PiscivS=[Size(11) SalpTL(11)];
GelativS=[Size(12) SalpTL(12)];
MyctS=[Size(13) SalpTL(13)];
BacS=[Size(14) SalpTL(14)];

PicoN=[Size(1) NonSalpTL(1)];
DtmN=[Size(2) NonSalpTL(2)];
FlagN=[Size(3) NonSalpTL(3)];
HNFN=[Size(4) NonSalpTL(4)];
MicN=[Size(5) NonSalpTL(5)];
MesoN=[Size(6) NonSalpTL(6)];
MacroN=[Size(7) NonSalpTL(7)];
GelN=[Size(8) NonSalpTL(8)];
SalpN=[Size(9) NonSalpTL(9)];
AmphN=[Size(10) NonSalpTL(10)];
PiscivN=[Size(11) NonSalpTL(11)];
GelativN=[Size(12) NonSalpTL(12)];
MyctN=[Size(13) NonSalpTL(13)];
BacN=[Size(14) NonSalpTL(14)];

f1=figure('Position',[100 100 1000 500]);
figure(f1)
subplot(1,2,1)
xlim([-1 8])
ylim([0.5 6])
XTick = 10.^(0:7);
XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));
set(gca,'XTick',log10(XTick))
set(gca,'XTickLabel',XTickLabels)
set(gca,'Ygrid','on')
set(gca,'Xgrid','on')
set(gca,'GridLineStyle','--')
set(gca,'GridAlpha',0.1)
xlabel('Size (um)')
ylabel('Trophic Level')
hold on
%Pico arrows
dp=FlagS-PicoS;
quiver(PicoS(1),PicoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Pico_col,Flag_col))
dp=HNFS-PicoS;
quiver(PicoS(1),PicoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Pico_col,HNF_col))
dp=MicS-PicoS;
quiver(PicoS(1),PicoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Pico_col,Mic_col))
dp=SalpS-PicoS;
quiver(PicoS(1),PicoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Pico_col,Salp_col))
%Dtm arrows
dp=MicS-DtmS;
quiver(DtmS(1),DtmS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Dtm_col,Mic_col))
dp=MesoS-DtmS;
quiver(DtmS(1),DtmS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Dtm_col,Meso_col))
dp=MacroS-DtmS;
quiver(DtmS(1),DtmS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Dtm_col,Macro_col))
dp=SalpS-DtmS;
quiver(DtmS(1),DtmS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Dtm_col,Salp_col))
%Flag arrows
dp=MicS-FlagS;
quiver(FlagS(1),FlagS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Flag_col,Mic_col))
dp=MesoS-FlagS;
quiver(FlagS(1),FlagS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Flag_col,Meso_col))
dp=MacroS-FlagS;
quiver(FlagS(1),FlagS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Flag_col,Macro_col))
dp=SalpS-FlagS;
quiver(FlagS(1),FlagS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Flag_col,Salp_col))
%HNF arrows
dp=MicS-HNFS;
quiver(HNFS(1),HNFS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(HNF_col,Mic_col))
dp=MesoS-HNFS;
quiver(HNFS(1),HNFS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(HNF_col,Meso_col))
dp=MacroS-HNFS;
quiver(HNFS(1),HNFS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(HNF_col,Macro_col))
dp=SalpS-HNFS;
quiver(HNFS(1),HNFS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(HNF_col,Salp_col))
%Mic arrows
dp=MesoS-MicS;
quiver(MicS(1),MicS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Mic_col,Meso_col))
dp=MacroS-MicS;
quiver(MicS(1),MicS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Mic_col,Macro_col))
dp=SalpS-MicS;
quiver(MicS(1),MicS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Mic_col,Salp_col))
%Meso arrows
dp=MacroS-MesoS;
quiver(MesoS(1),MesoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Meso_col,Macro_col))
dp=GelS-MesoS;
quiver(MesoS(1),MesoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Meso_col,Gel_col))
dp=MyctS-MesoS;
quiver(MesoS(1),MesoS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Meso_col,Myct_col))
%Macro arrows
dp=GelS-MacroS;
quiver(MacroS(1),MacroS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Macro_col,Gel_col))
dp=PiscivS-MacroS;
quiver(MacroS(1),MacroS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Macro_col,Pisciv_col))
dp=MyctS-MacroS;
quiver(MacroS(1),MacroS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Macro_col,Myct_col))
%Gel arrows
dp=AmphS-GelS;
quiver(GelS(1),GelS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Gel_col,Amph_col))
dp=GelativS-GelS;
quiver(GelS(1),GelS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Gel_col,Gelativ_col))
dp=MyctS-GelS;
quiver(GelS(1),GelS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Gel_col,Myct_col))
%Salp arrows
dp=AmphS-SalpS;
quiver(SalpS(1),SalpS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Salp_col,Amph_col))
dp=GelativS-SalpS;
quiver(SalpS(1),SalpS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Salp_col,Gelativ_col))
dp=MyctS-SalpS;
quiver(SalpS(1),SalpS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Salp_col,Myct_col))
%Amph arrows
dp=PiscivS-AmphS;
quiver(AmphS(1),AmphS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Amph_col,Pisciv_col))
dp=MyctS-AmphS;
quiver(AmphS(1),AmphS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Amph_col,Myct_col))
%Myct arrows
dp=PiscivS-MyctS;
quiver(MyctS(1),MyctS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Myct_col,Pisciv_col))
%Bac arrows
dp=FlagS-BacS;
quiver(BacS(1),BacS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Bac_col,Flag_col))
dp=HNFS-BacS;
quiver(BacS(1),BacS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Bac_col,HNF_col))
dp=MicS-BacS;
quiver(BacS(1),BacS(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngS(Bac_col,Mic_col))
%Labels
%text(Size,SalpTL,Labs,'Color','red')
%Circles around taxa with color as productivity
scatter(Size, SalpTL,400,colorsS,'LineWidth',2);


subplot(1,2,2)%Pico arrows
xlim([-1 8])
ylim([0.5 6])
XTick = 10.^(0:7);
XTickLabels = cellstr(num2str(round(log10(XTick(:))), '10^%d'));
set(gca,'XTick',log10(XTick))
set(gca,'XTickLabel',XTickLabels)
set(gca,'Ygrid','on')
set(gca,'Xgrid','on')
set(gca,'GridLineStyle','--')
set(gca,'GridAlpha',0.1)
xlabel('Size (um)')
ylabel('Trophic Level')
hold on
%Pico arrows
dp=FlagN-PicoN;
quiver(PicoN(1),PicoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Pico_col,Flag_col))
dp=HNFN-PicoN;
quiver(PicoN(1),PicoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Pico_col,HNF_col))
dp=MicN-PicoN;
quiver(PicoN(1),PicoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Pico_col,Mic_col))
dp=SalpN-PicoN;
quiver(PicoN(1),PicoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Pico_col,Salp_col))
%Dtm arrows
dp=MicN-DtmN;
quiver(DtmN(1),DtmN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Dtm_col,Mic_col))
dp=MesoN-DtmN;
quiver(DtmN(1),DtmN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Dtm_col,Meso_col))
dp=MacroN-DtmN;
quiver(DtmN(1),DtmN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Dtm_col,Macro_col))
dp=SalpN-DtmN;
quiver(DtmN(1),DtmN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Dtm_col,Salp_col))
%Flag arrows
dp=MicN-FlagN;
quiver(FlagN(1),FlagN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Flag_col,Mic_col))
dp=MesoN-FlagN;
quiver(FlagN(1),FlagN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Flag_col,Meso_col))
dp=MacroN-FlagN;
quiver(FlagN(1),FlagN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Flag_col,Macro_col))
dp=SalpN-FlagN;
quiver(FlagN(1),FlagN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Flag_col,Salp_col))
%HNF arrows
dp=MicN-HNFN;
quiver(HNFN(1),HNFN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(HNF_col,Mic_col))
dp=MesoN-HNFN;
quiver(HNFN(1),HNFN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(HNF_col,Meso_col))
dp=MacroN-HNFN;
quiver(HNFN(1),HNFN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(HNF_col,Macro_col))
dp=SalpN-HNFN;
quiver(HNFN(1),HNFN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(HNF_col,Salp_col))
%Mic arrows
dp=MesoN-MicN;
quiver(MicN(1),MicN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Mic_col,Meso_col))
dp=MacroN-MicN;
quiver(MicN(1),MicN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Mic_col,Macro_col))
dp=SalpN-MicN;
quiver(MicN(1),MicN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Mic_col,Salp_col))
%Meso arrows
dp=MacroN-MesoN;
quiver(MesoN(1),MesoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Meso_col,Macro_col))
dp=GelN-MesoN;
quiver(MesoN(1),MesoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Meso_col,Gel_col))
dp=MyctN-MesoN;
quiver(MesoN(1),MesoN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Meso_col,Myct_col))
%Macro arrows
dp=GelN-MacroN;
quiver(MacroN(1),MacroN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Macro_col,Gel_col))
dp=PiscivN-MacroN;
quiver(MacroN(1),MacroN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Macro_col,Pisciv_col))
dp=MyctN-MacroN;
quiver(MacroN(1),MacroN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Macro_col,Myct_col))
%Gel arrows
dp=AmphN-GelN;
quiver(GelN(1),GelN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Gel_col,Amph_col))
dp=GelativN-GelN;
quiver(GelN(1),GelN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Gel_col,Gelativ_col))
dp=MyctN-GelN;
quiver(GelN(1),GelN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Gel_col,Myct_col))
%Salp arrows
dp=AmphN-SalpN;
quiver(SalpN(1),SalpN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Salp_col,Amph_col))
dp=GelativN-SalpN;
quiver(SalpN(1),SalpN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Salp_col,Gelativ_col))
dp=MyctN-SalpN;
quiver(SalpN(1),SalpN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Salp_col,Myct_col))
%Amph arrows
dp=PiscivN-AmphN;
quiver(AmphN(1),AmphN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Amph_col,Pisciv_col))
dp=MyctN-AmphN;
quiver(AmphN(1),AmphN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Amph_col,Myct_col))
%Myct arrows
dp=PiscivN-MyctN;
quiver(MyctN(1),MyctN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Myct_col,Pisciv_col))
%Bac arrows
dp=FlagN-BacN;
quiver(BacN(1),BacN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Bac_col,Flag_col))
dp=HNFN-BacN;
quiver(BacN(1),BacN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Bac_col,HNF_col))
dp=MicN-BacN;
quiver(BacN(1),BacN(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',IngN(Bac_col,Mic_col))
%Legend arrow
Start=[6 1];
End=[7.5 1];
dp=End-Start;
quiver(Start(1),Start(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',1)
% quiver(Start(1),Start(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',log10(1/shift))
% quiver(Start(1),Start(2),dp(1),dp(2),0,'Color',[120/255 120/255 120/255],'LineWidth',log10(1/shift)/log10(100))
%Labels
%text(Size,NonSalpTL,Labs,'Color','red')
%Circles around taxa with color as productivity
%Circles around taxa with color as productivity
scatter(Size, NonSalpTL,400,colorsN,'LineWidth',2);

f2=figure('Position',[100 100 1000 500]);
figure(f2)
cb=colorbar;
colormap(cmap);
% Customize tick locations and labels
num_ticks = 6; % Number of ticks you want
tick_values = linspace(Min, Max, num_ticks); % Tick values
tick_positions = (tick_values - Min) / (Max - Min); % Normalized positions
% Update colorbar properties
cb.Ticks = tick_positions; % Set tick locations
cb.TickLabels = {'10^-^4', '10^-^3', '10^-^2', '10^-^1', '10^0', '10^1'}; % Set tick labels
cb.Label.String = 'Value'; % Add a label to the colorbar

%exportgraphics(f1,'FoodWebFig_clean.png','Resolution',300)