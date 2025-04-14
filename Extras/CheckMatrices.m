clear all
close all

%%This script, as the name implies, runs some initial statistics and
%%visualizations for model results mostly focused on error. For the
%%linegraph used in the paper to visualize these errors, see BasicPlots.m

%%Modify here to change cycle input
Cycle=5;

%%Loading input and output matrices
load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
if Cycle==1
    [num,txt,raw]  = xlsread('N15InverseModelNZC1.xlsx','NZ','B6:B127');
else
    [num,txt,raw]  = xlsread('N15InverseModelNZ.xlsx','NZ','B6:B145');
end

%%Basic summary stats
Columns = txt;
means = mean(MCMCmat)';
%start = centralval./wts;
l=size(MCMCmat);
final = MCMCmat(l(1),:)'./wts;
meanvals = mean(MCMCmat)'./wts;
%meanvalsplain = mean(MCMCmatplain)'./wts;
stdvals = std(MCMCmat)'./wts;
meantable = table(Columns,meanvals,stdvals);
MCMCcor = MCMCmat;
for i=1:length(MCMCmat(:,1))
    MCMCcor(i,:)=MCMCmat(i,:)'./wts;
end
meand15 = mean(MCMCmat)'./wts;
stdvalsd15 = std(MCMCmat)'./wts;
meantable = table(Columns,meanvals,stdvals);

%%Calculating approximate equality outputs and error
for i=1:length(MCMCmat(:,1))
    
    MCMC=MCMCmat(i,:);
    d15N=del15N(i,:);

    if Cycle==1
        [Aa] = ResetRN15C1(Aa0,MCMC,wts,d15N,d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
    else
        [Aa] = ResetRN15(Aa0,MCMC,wts,d15N,d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
    end

    AaOut(:,i)=Aa*MCMC';

    Error(:,i)=Aa*MCMC'-ba; %how far above or below Ae the solutions were
    ErrorNorm(:,i)=Error(:,i)./sdba; %how many STDs above or below Ae
    SSR(:,i)=sum(ErrorNorm(:,i).^2); %total error of the model
    InequalError(:,i)=G*MCMC'-h; %how far above or below the Gx>h constraints the solutions were
    
    
end
AaOutmean=mean(AaOut')';
%A%dding two columns to compare 50th percentiles minus the 20% burn in,
%%informative to convergence
len = length(AaOut);
temp = AaOut(:,round(len*0.2):len);
cut = length(temp)/2;
AaOutmean(:,2)=mean(temp(:,1:cut),2);
AaOutmean(:,3)=mean(temp(:,cut:end),2);
AaOutstd=std(AaOut')';
%%Same deal for the flow means to look at convergence more in depth
len = length(MCMCcor);
temp = MCMCcor(round(len*0.2):len,:)';
cut = length(temp)/2;
meanvals(:,2)=mean(temp(:,1:cut),2);
meanvals(:,3)=mean(temp(:,cut:end),2);
meanvals(:,4)=std(temp(:,:),0,2);

%%Label names for following figures
if length(Aa(:,1))==60
    %nonD15N Aa inputs
    Labels1={'NPP','STFlux','PicoNPP','FlagNPP','DtmNPP','PicoMort','FlagMort','DtmMort','MicGr','nvmMesoGr','nvmMacroGr','vmMacroGr',...
        'SalpAutoGr','SalpHetGr','Pico2Salp','Dtm2Salp','Other2Salp','SalpFPProd','BacProd','MesoTOMyct',...
        'MacroTOMyct','VMMacroTOMyct','AmphTOMyct','SalpTOMyct','GelTOMyct','MyctTOPisciv',...
        'VMMacroTOPisciv','AmphTOPisciv','SalpTOGelativ','GelTOGelativ','MesoTOHTL','MacroTOHTL',...
        'VMMacroTOHTL','AmphTOHTL','SalpTOHTL','GelTOHTL'};
    %Compartment names
    Labels2={'NO3','NH4','Pico','Dtm','Flag','HNF','Mic','nvmMeso','nvmMacro','vmMacro','Gel','Salp','Amph','PiscFish','GelatFish','Myct',...
        'Bac','sdet','mdet','ldet','SalpDet','Dom','SalpMort','Export'};
    %Some of the inequality constraints
    Labels3={'SalpFP','Exc_HNF','Exc_Mic','Exc_nvmMeso','Exc_nvmMacro','Exc_VMMacro','Exc_Gel','Exc_Salp','Exc_Amph','Exc_PiscFish',...
        'Exc_GelatFish','Exc_Myct',...
        'Resp_nvmMeso','Resp_nvmMacro','Resp_VMMacro','Resp_Amph','Resp_Gel','Resp_Myct','Resp_PiscFish','Resp_Salp',...
        'Resp_ULVMMacro','Resp_ULSalp',"Resp_ULMyct"};
elseif length(Aa(:,1))==64
    Labels1={'NPP','STFlux','PicoNPP','FlagNPP','DtmNPP','PicoMort','FlagMort','DtmMort','MicGr','nvmMesoGr','vmMesoGr','nvmMacroGr','vmMacroGr',...
        'SalpAutoGr','SalpHetGr','Pico2Salp','Dtm2Salp','Other2Salp','SalpFPProd','BacProd',...
        'MesoToMyct','VMMesoToMyct','MacroToMyct','VMMacroToMyct','AmphToMyct','SalpToMyct','GelToMyct',...
        'MyctToPisc','VMMacroToPisc','AmphToPisc','SalpToGelat','GelToGelat','MesoTOHTL',...
        'VMMesoTOHTL','MacroTOHTL','VMMacroTOHTL','AmphTOHTL','SalpTOHTL','GelTOHTL'};
    Labels2={'NO3','NH4','Pico','Dtm','Flag','HNF','Mic','nvmMeso','vmMeso','nvmMacro','vmMacro','Gel','Salp','Amph','PiscFish','GelatFish','Myct',...
        'Bac','sdet','mdet','ldet','SalpDet','Dom','SalpMort','Export'};
    Labels3={'SalpFP','Exc_HNF','Exc_Mic','Exc_nvmMeso','Exc_vmMeso','Exc_nvmMacro','Exc_VMMacro','Exc_Gel','Exc_Salp','Exc_Amph','Exc_PiscFish',...
        'Exc_GelatFish','Exc_Myct',...
        'Resp_nvmMeso','Resp_vmMeso','Resp_nvmMacro','Resp_VMMacro','Resp_Amph','Resp_Gel','Resp_Myct','Resp_PiscFish','Resp_Salp',...
        'Resp_ULVMMeso','Resp_ULVMMacro','Resp_ULSalp',"Resp_ULMyct"};
end

meanD15 = mean(del15N)';
stdvalsD15 = std(del15N)';
meantableD15 = table([Labels2,"Upwelling"]',meanD15,stdvalsD15);
if Cycle==1
    meantableD15{8,2:3}=d15NInputs(1,InputCol:InputCol+1); %nvmMeso
    meantableD15{9,2:3}=d15NInputs(2,InputCol:InputCol+1); %nvmMacro
    meantableD15{10,2:3}=d15NInputs(2,InputCol:InputCol+1); %vmMacro
    meantableD15{12,2:3}=d15NInputs(3,InputCol:InputCol+1); %Salp
elseif Cycle==2|Cycle==4
    meantableD15{8,2:3}=d15NInputs(1,InputCol:InputCol+1); %nvmMeso
    meantableD15{9,2:3}=d15NInputs(1,InputCol:InputCol+1); %nvmMeso
    meantableD15{10,2:3}=d15NInputs(2,InputCol:InputCol+1); %nvmMacro
    meantableD15{11,2:3}=d15NInputs(2,InputCol:InputCol+1); %vmMacro
    meantableD15{13,2:3}=d15NInputs(3,InputCol:InputCol+1); %Salp
elseif Cycle==3|Cycle==5
    meantableD15{8,2:3}=d15NInputs(1,InputCol:InputCol+1); %nvmMeso
    meantableD15{9,2:3}=d15NInputs(1,InputCol:InputCol+1); %nvmMeso
    meantableD15{10,2:3}=d15NInputs(2,InputCol:InputCol+1); %nvmMacro
    meantableD15{11,2:3}=d15NInputs(2,InputCol:InputCol+1); %vmMacro
end




%%Calculating and displaying measurements that are a particular source of
%%error (defined as >10% mean SSR)
ErrorNormMean=mean(ErrorNorm')';
SRmean=ErrorNormMean.^2;
l=sum(SRmean);
m=find(SRmean>.1*l);
[Labels1(m),' each represent >10% of the mean sum of squared residuals']
%%Root mean squared error, or the average number of standard errors the
%%model was from predictions
RMSEall=sqrt(mean((abs(AaOutmean-ba)).^2));
if Cycle==1
    RMSErate=sqrt(mean((abs(AaOutmean(1:19)-ba(1:19))).^2));
else
    RMSErate=sqrt(mean((abs(AaOutmean(1:20)-ba(1:20))).^2));
end




%%Normalized error of solutions as the run progressed
%Note only every 10000th sim is plotted
figure('Position',[100 100 1000 500])
plot(SSR)
%ylim([0 2000])
ylabel('Normalized Error Squared')

%Progression of approximate equalities
%Note only every 1,000th sim is saved and therefore plotted
AaPlot=AaOut';
for i=1:length(AaPlot)
    figure('Position',[100 100 1000 500])
    plot(AaPlot(:,i))
    yline(AaOutmean(i),'--r') %mean as red line
    yline(ba(i),'--k') %observation as black line
    ylabel(Labels1(i))
    ylim([0 max([max(AaPlot(:,i)) ba(i)])*1.1])
    xlim([0 length(AaPlot)*0.01])
end

%Progression of flows
close all
for i=80:100
    figure('Position',[100 100 1000 500])
    plot(MCMCmat(:,i)/wts(i))
    ylabel(Columns(i))
    yline(meanvals(i),'--r')
    ylim([0 max(MCMCmat(:,i)/wts(i))*1.1])
end

%Summed fish ingestions
if Cycle==1
    M=[51 57 64 72 86 96];%rows for ing flows
    P=[62 78 95];
    G=[73 93];
    %Myct
    figure('Position',[100 100 1000 500])
    plot(sum(MCMCcor(:,M),2))
    ylabel('Myct Total Ing')
    yline(sum(h([120 122 124 126 128 130])),'--b')%mins
    yline(mean(sum(MCMCcor(:,M),2)),'--k')
    yline(sum(abs(h([121 123 125 127 129 131]))),'--r')%max
    %Pisciv
    figure('Position',[100 100 1000 500])
    plot(sum(MCMCcor(:,P),2))
    ylabel('Pisciv Total Ing')
    yline(sum(h([132 134 136])),'--b')%mins
    yline(mean(sum(MCMCcor(:,P),2)),'--k')
    yline(sum(abs(h([133 135 137]))),'--r')%max
    %Gelativ
    figure('Position',[100 100 1000 500])
    plot(sum(MCMCcor(:,G),2))
    ylabel('Gelativ Total Ing')
    yline(sum(h([138 140])),'--b')%mins
    yline(mean(sum(MCMCcor(:,G),2)),'--k')
    yline(sum(abs(h([139 141]))),'--r')%max
else
    M=[55 63 71 78 86 100 110];
    P=[76 92 109];
    G=[87 107];
    %Myct
    figure('Position',[100 100 1000 500])
    plot(sum(MCMCcor(:,M),2))
    ylabel('Myct Total Ing')
    yline(sum(h([132 134 136 138 140 142 144])),'--b')%mins
    yline(mean(sum(MCMCcor(:,M),2)),'--k')
    yline(sum(abs(h([133 134 137 139 141 143 145]))),'--r')%max
    %Pisciv
    figure('Position',[100 100 1000 500])
    plot(sum(MCMCcor(:,P),2))
    ylabel('Pisciv Total Ing')
    yline(sum(h([146 148 150])),'--b')%mins
    yline(mean(sum(MCMCcor(:,P),2)),'--k')
    yline(sum(abs(h([147 149 151]))),'--r')%max
    %Gelativ
    figure('Position',[100 100 1000 500])
    plot(sum(MCMCcor(:,G),2))
    ylabel('Gelativ Total Ing')
    yline(sum(h([152 154])),'--b')%mins
    yline(mean(sum(MCMCcor(:,G),2)),'--k')
    yline(sum(abs(h([153 155]))),'--r')%max
end


%%Progression of d15N
del15Nplot=del15N';
for i=1:length(del15Nplot)
figure('Position',[100 100 1000 500])
plot(del15Nplot(i,:))
%ylim([min(min(AaOut(20:42,:))) max(max(AaOut(20:42,:)))])
ylabel(['d15N',Labels2(i)])
end


%%Absolute Error of the predicted approximate equalities (Measurements).
%%Note that the mean value represents accuracy relative to measurements 
%%whereas the spread is variation across simulations
figure('Position',[100 100 1000 500])
boxplot(Error(1:length(Labels1),:)')
yline(0,'--')
set(gca,'XTick',1:length(Labels1))
set(gca,'XTickLabel',Labels1)
xtickangle(45)
ylabel('Error (mmol N m-2 d-1) negative means a model underestimate')
title('Measurements')

%%Absolute Error of the modelled d15N, should all be 0 as it's a mass balance
figure('Position',[100 100 1000 500])
boxplot(Error(length(Labels1)+1:end,:)')
yline(0,'--')
set(gca,'XTick',1:length(Labels2))
set(gca,'XTickLabel',Labels2)
xtickangle(45)
ylabel('Error (mmol N m-2 d-1)')
title('d^1^5N')


%Error of the Ae relative to their measured uncertainty. I.e. the number
%of SDs above or below the measurement the model solutions fell.
figure('Position',[100 100 1000 500])
boxplot(ErrorNorm(1:length(Labels1),:)')
yline(0,'--')
set(gca,'XTick',1:length(Labels1))
set(gca,'XTickLabel',Labels1)
xtickangle(45)
ylim([-4 10])
ylabel('Normalized Error (mmol N m-2 d-1) negative means a model underestimate')
title('Measurements')

figure('Position',[100 100 1000 500])
boxplot(ErrorNorm(length(Labels1)+1:end,:)')
yline(0,'--')
set(gca,'XTick',1:length(Labels2))
set(gca,'XTickLabel',Labels2)
xtickangle(45)
ylabel('Normalized Error (mmol N m-2 d-1)')
title('d^1^5N')

%%How far the modeled constraints strayed from their min/max
figure('Position',[100 100 1000 500])
boxplot(InequalError(1:length(Labels3),:)')
yline(0,'--')
set(gca,'XTick',1:length(Labels3))
set(gca,'XTickLabel',Labels3)
xtickangle(45)
ylabel('Inequality Bounds (mmol N m-2 d-1)')

%%Plot of measurements vs inputs
figure('Position',[100 100 1000 500])
scatter(Inputs(:,InputCol),AaOutmean(1:length(Inputs)))
xlabel('Measurements (mmol N m-2 d-1)')
ylabel('Model (mmol N m-2 d-1)')
set(gca,'YScale','log')
set(gca,'XScale','log')
hold on
plot(xlim,ylim,'--r')
hold off
