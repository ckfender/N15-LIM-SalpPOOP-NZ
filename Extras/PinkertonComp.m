clear all
close all

%%Some figures comparing model outputs to those of the most comparable
%%prior study, the EcoPath model by Pinkerton 2011 for the Chatham Rise.
%%None of these made it into the paper and it's not as interesting of a 
%%comparison as it first seems since many of the inputs to the Pinkerton 
%%model were also used in this LIEM.

ProdMean = readtable('Prod.xlsx','Sheet','Mean');
ProdStd = readtable('Prod.xlsx','Sheet','Std');
TLMean = readtable('TL.xlsx','Sheet','Mean');
TLStd = readtable('TL.xlsx','Sheet','Std');
InputsNZ = readtable('InputsNZ.xlsx');
load('DietFluxMean.mat');
load('DietFluxStd.mat');
PinkComp = readtable('PinkComp.xlsx');
ProdCol = ProdMean.Properties.VariableNames;

f1=figure('Position',[100 100 1000 500]);
f2=figure('Position',[100 100 1000 500]);
f3=figure('Position',[100 100 1000 500]);
f4=figure('Position',[100 100 1000 500]);
f5=figure('Position',[100 100 1000 500]);
f6=figure('Position',[100 100 1000 500]);
for Cycle = 1:5
    load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
    if Cycle==1
        [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,Macro_col,...
        VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
        mdet_col,ldet_col,SalpDet_col,Dom_col,Export_col,upNO3_col] = GetColumnsC1(Aa0,MCMCmat,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
        VMMeso_col=NaN;
    else
        [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,...
        VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
        mdet_col,ldet_col,SalpDet_col,Dom_col,Export_col,upNO3_col] = GetColumns(Aa0,MCMCmat,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
    end
    %Calculating and populating observed biomasses. Note Pinkerton annual
    %average biomasses are assumed to be prepopulated.
    PicoBiom=InputsNZ{37,Cycle*2+2};
    DtmBiom=InputsNZ{38,Cycle*2+2};
    FlagBiom=InputsNZ{39,Cycle*2+2};
    HNFBiom=InputsNZ{40,Cycle*2+2};
    MicBiom=InputsNZ{41,Cycle*2+2};
    if Cycle == 1
        MesoBiom=mean([InputsNZ{42,Cycle*2+2} InputsNZ{43,Cycle*2+2}]);
        VMMesoBiom=0;
    else
        MesoBiom=InputsNZ{42,Cycle*2+2};
        VMMesoBiom=InputsNZ{43,Cycle*2+2}-InputsNZ{42,Cycle*2+2};
    end
    MacroBiom=InputsNZ{44,Cycle*2+2};
    VMMacroBiom=InputsNZ{45,Cycle*2+2}-InputsNZ{44,Cycle*2+2};
    GelBiom=InputsNZ{47,Cycle*2+2};
    SalpBiom=InputsNZ{49,Cycle*2+2};
    AmphBiom=InputsNZ{46,Cycle*2+2};
    PiscivBiom=InputsNZ{51,Cycle*2+2};
    GelativBiom=InputsNZ{52,Cycle*2+2};
    MyctBiom=InputsNZ{50,Cycle*2+2};
    PinkComp{1,Cycle+2}= PiscivBiom;%Pisciv
    PinkComp{2,Cycle+2}= GelativBiom;%Gelativ
    PinkComp{3,Cycle+2}= MyctBiom;%Myct
    PinkComp{4,Cycle+2}= SalpBiom+GelBiom;%Salp+Gel
    PinkComp{5,Cycle+2}= MacroBiom+VMMacroBiom+AmphBiom;%Macro+Amph
    PinkComp{6,Cycle+2}= MesoBiom+VMMesoBiom;%Meso
    PinkComp{7,Cycle+2}= MicBiom;%Mic
    PinkComp{8,Cycle+2}= HNFBiom;%HNF
    PinkComp{9,Cycle+2}= PicoBiom+DtmBiom+FlagBiom;%Phyto NPP
    PinkComp{10,Cycle+2}= NaN;%Bac
    %Populating Production
    PinkComp{1,Cycle+8}= ProdMean{Cycle,['Pisciv']};%Pisciv
    PinkComp{2,Cycle+8}= ProdMean{Cycle,['Gelativ']};%Gelativ
    PinkComp{3,Cycle+8}= ProdMean{Cycle,['Myct']};%Myct
    PinkComp{4,Cycle+8}= ProdMean{Cycle,['Salp']}+ProdMean{Cycle,['Gel']};%Salp+Gel
    PinkComp{5,Cycle+8}= ProdMean{Cycle,['Macro']}+ProdMean{Cycle,['VMMacro']}+ProdMean{Cycle,['Amph']};%Macro+Amph
    if Cycle==1
        PinkComp{6,Cycle+8}= ProdMean{Cycle,['Meso']};%Meso
    else
        PinkComp{6,Cycle+8}= ProdMean{Cycle,['Meso']}+ProdMean{Cycle,['VMMeso']};%Meso
    end
    PinkComp{7,Cycle+8}= ProdMean{Cycle,['Mic']};%Mic
    PinkComp{8,Cycle+8}= ProdMean{Cycle,['HNF']};%HNF
    PinkComp{9,Cycle+8}= ProdMean{Cycle,['Pico']}+ProdMean{Cycle,['Dtm']}+ProdMean{Cycle,['Flag']};%Phyto NPP
    PinkComp{10,Cycle+8}= ProdMean{Cycle,['Bac']};%Bac

    %Populating Consumption
    PinkComp{1,Cycle+14}= sum(DietFluxMean(:,PiscivFish_col,Cycle));%Pisciv
    PinkComp{2,Cycle+14}= sum(DietFluxMean(:,GelativFish_col,Cycle));%Gelativ
    PinkComp{3,Cycle+14}= sum(DietFluxMean(:,Myct_col,Cycle));%Myct
    PinkComp{4,Cycle+14}= sum(DietFluxMean(:,Salp_col,Cycle))+sum(DietFluxMean(:,Gel_col,Cycle));%Salp+Gel
    PinkComp{5,Cycle+14}= sum(DietFluxMean(:,Macro_col,Cycle))+sum(DietFluxMean(:,VMMacro_col,Cycle))+sum(DietFluxMean(:,Amph_col,Cycle));%Macro+Amph
    if Cycle==1
        PinkComp{6,Cycle+14}= sum(DietFluxMean(:,Meso_col,Cycle));%Meso
    else
        PinkComp{6,Cycle+14}= sum(DietFluxMean(:,Meso_col,Cycle))+sum(DietFluxMean(:,VMMeso_col,Cycle));%Meso
    end
    PinkComp{7,Cycle+14}= sum(DietFluxMean(:,Mic_col,Cycle));%Mic
    PinkComp{8,Cycle+14}= sum(DietFluxMean(:,HNF_col,Cycle));%HNF
    PinkComp{9,Cycle+14}= sum(DietFluxMean(:,Pico_col,Cycle))+sum(DietFluxMean(:,Dtm_col,Cycle))+sum(DietFluxMean(:,Flag_col,Cycle));%Phyto
    PinkComp{10,Cycle+14}= sum(DietFluxMean(:,bac_col,Cycle));%Bac

    %Populating Trophic Level
    %To be comparable with Pinkerton groups, groups that are a combination
    %of model groups are represented with biomass weighted mean TL
    PinkComp{1,Cycle+20}= TLMean{Cycle,['PiscFish']};%Pisciv
    PinkComp{2,Cycle+20}= TLMean{Cycle,['GelatFish']};%Gelativ
    PinkComp{3,Cycle+20}= TLMean{Cycle,['Myct']};%Myct
    PinkComp{4,Cycle+20}= ((TLMean{Cycle,['Salp']}*SalpBiom)+(TLMean{Cycle,['Gel']}*GelBiom))/(SalpBiom+GelBiom);%Salp+Gel
    PinkComp{5,Cycle+20}= ((TLMean{Cycle,['nvmMacro']}*MacroBiom)+(TLMean{Cycle,['vmMacro']}*VMMacroBiom)+(TLMean{Cycle,['Amph']}*AmphBiom))/(MacroBiom+VMMacroBiom+AmphBiom);%Macro+Amph
    if Cycle==1
        PinkComp{6,Cycle+20}= TLMean{Cycle,['nvmMeso']};%Meso
    else
        PinkComp{6,Cycle+20}= ((TLMean{Cycle,['nvmMeso']}*MesoBiom)+(TLMean{Cycle,['VMMeso']}*VMMesoBiom))/(MesoBiom+VMMesoBiom);%Meso
    end
    PinkComp{7,Cycle+20}= TLMean{Cycle,['Mic']};%Mic
    PinkComp{8,Cycle+20}= TLMean{Cycle,['HNF']};%HNF
    PinkComp{9,Cycle+20}= ((TLMean{Cycle,['Pico']}*PicoBiom)+(TLMean{Cycle,['Dtm']}*DtmBiom)+(TLMean{Cycle,['Flag']}*FlagBiom))/(PicoBiom+DtmBiom+FlagBiom);%Phyto NPP
    PinkComp{10,Cycle+20}= TLMean{Cycle,['Bac']};%Bac
end

%Plotting Biomass, Production normalized to Biomass, and Normalized Consumption vs Pinkerton's
for Cycle=[1:5]
    clearvars qw
    col = [0    0.4471    0.7412; %C1 dark blue
           0    0.6000    1.0000; %C2 blue
           0.6353    0.0784    0.1843; %C3 dark red
           0.3020    0.7451    0.9333; %C4 light blue
           1     0     0]; %C5 red
    mark = ['^'; '^'; '^'; '*'; 's'; 's'; 'o'; 'o'; 'd'; 'h'];
    mark2 = ['*'; 's'; 's'; 'o'; 'o'; 'd'; 'h'];
    %Biomass
    figure(f1)
    for i=4:9
        scatter(table2array(PinkComp(i,2)),table2array(PinkComp(i,Cycle+2)),'col',col(Cycle,:),...
            'Marker',mark2(i-3),'MarkerFaceColor',col(Cycle,:),'MarkerEdgeColor',col(Cycle,:))
        hold on
    end
    xlabel('Annual Mean Biomass (mmol N m^-^2 d^-^1)')
    ylabel('Observed Biomass (mmol N m^-^2 d^-^1)')
    xlim([10^-2 10^2])
    ylim([10^-2 10^2])
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        %Invisible dummy plots to act as handles for the legend
        qw{1}=scatter(nan,nan,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Jelly');
        qw{2}=scatter(nan,nan,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Zoo');
        qw{3}=scatter(nan,nan,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Protist');
        qw{4}=scatter(nan,nan,'Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Phyto');
        legend([qw{:}],{'Jelly','Zoo','Protist','Phyto'},'Location','southeast')
    end

    % Production
    figure(f2)
    for i=1:10
        scatter(table2array(PinkComp(i,8)),table2array(PinkComp(i,Cycle+8)),'col',col(Cycle,:),...
            'Marker',mark(i),'MarkerFaceColor',col(Cycle,:),'MarkerEdgeColor',col(Cycle,:))
        hold on
    end
    xlabel('Annual Mean Production (mmol N m^-^2 d^-^1)')
    ylabel('Modelled Production (mmol N m^-^2 d^-^1)')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10^-4 10^2])
    ylim([10^-4 10^2])
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        qw{1}=scatter(nan,nan,'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Fish');
        qw{2}=scatter(nan,nan,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Jelly');
        qw{3}=scatter(nan,nan,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Zoo');
        qw{4}=scatter(nan,nan,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Protist');
        qw{5}=scatter(nan,nan,'Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Phyto');
        qw{6}=scatter(nan,nan,'Marker','h','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Bac');
        legend([qw{:}],{'Fish','Jelly','Zoo','Protist','Phyto','Bac'},'Location','southeast')
    end

    %Consumption
    figure(f3)
    for i=1:10
        scatter(table2array(PinkComp(i,14)),table2array(PinkComp(i,Cycle+14)),'col',col(Cycle,:),...
            'Marker',mark(i),'MarkerFaceColor',col(Cycle,:),'MarkerEdgeColor',col(Cycle,:))
        hold on
    end
    xlabel('Annual Mean Consumption (mmol N m^-^2 d^-^1)')
    ylabel('Modelled Consumption (mmol N m^-^2 d^-^1)')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10^-4 10^2])
    ylim([10^-4 10^2])
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        qw{1}=scatter(nan,nan,'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Fish');
        qw{2}=scatter(nan,nan,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Jelly');
        qw{3}=scatter(nan,nan,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Zoo');
        qw{4}=scatter(nan,nan,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Protist');
        qw{5}=scatter(nan,nan,'Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Phyto');
        qw{6}=scatter(nan,nan,'Marker','h','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Bac');
        legend([qw{:}],{'Fish','Jelly','Zoo','Protist','Phyto','Bac'},'Location','southeast')
    end

    % P/B
    figure(f4)
    for i=1:10
        scatter(table2array(PinkComp(i,8))./table2array(PinkComp(i,2)),table2array(PinkComp(i,Cycle+8))./table2array(PinkComp(i,Cycle+2)),'col',col(Cycle,:),...
            'Marker',mark(i),'MarkerFaceColor',col(Cycle,:),'MarkerEdgeColor',col(Cycle,:))
        hold on
    end
    xlabel('Annual Mean Biomass Normalized Production (d^-^1)')
    ylabel('Modelled Biomass Normalized Production (d^-^1)')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10^-4 10^1])
    ylim([10^-4 10^1])
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        qw{1}=scatter(nan,nan,'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Fish');
        qw{2}=scatter(nan,nan,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Jelly');
        qw{3}=scatter(nan,nan,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Zoo');
        qw{4}=scatter(nan,nan,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Protist');
        qw{5}=scatter(nan,nan,'Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Phyto');
        qw{6}=scatter(nan,nan,'Marker','h','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Bac');
        legend([qw{:}],{'Fish','Jelly','Zoo','Protist','Phyto','Bac'},'Location','southeast')
    end

    %Norm Consumption
    figure(f5)
    for i=1:10
        scatter(table2array(PinkComp(i,14))./table2array(PinkComp(i,2)),table2array(PinkComp(i,Cycle+14))./table2array(PinkComp(i,Cycle+2)),'col',col(Cycle,:),...
            'Marker',mark(i),'MarkerFaceColor',col(Cycle,:),'MarkerEdgeColor',col(Cycle,:))
        hold on
    end
    xlabel('Annual Mean Biomass Normalized Consumption (d^-^1)')
    ylabel('Modelled Biomass Normalized Consumption (d^-^1)')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10^-4 10^1])
    ylim([10^-4 10^1])
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        qw{1}=scatter(nan,nan,'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Fish');
        qw{2}=scatter(nan,nan,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Jelly');
        qw{3}=scatter(nan,nan,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Zoo');
        qw{4}=scatter(nan,nan,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Protist');
        qw{5}=scatter(nan,nan,'Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Phyto');
        qw{6}=scatter(nan,nan,'Marker','h','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Bac');
        legend([qw{:}],{'Fish','Jelly','Zoo','Protist','Phyto','Bac'},'Location','southeast')
    end

    %Trophic Levels
    figure(f6)
    for i=1:10
        scatter(table2array(PinkComp(i,20)),table2array(PinkComp(i,Cycle+20)),'col',col(Cycle,:),...
            'Marker',mark(i),'MarkerFaceColor',col(Cycle,:),'MarkerEdgeColor',col(Cycle,:))
        hold on
    end
    xlabel('Annual Mean Trophic Level')
    ylabel('Modelled Tropic Level')
    xlim([0.5 5.5])
    ylim([0.5 5.5])
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        qw{1}=scatter(nan,nan,'Marker','^','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Fish');
        qw{2}=scatter(nan,nan,'Marker','*','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Jelly');
        qw{3}=scatter(nan,nan,'Marker','s','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Zoo');
        qw{4}=scatter(nan,nan,'Marker','o','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Protist');
        qw{5}=scatter(nan,nan,'Marker','d','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Phyto');
        qw{6}=scatter(nan,nan,'Marker','h','MarkerFaceColor','k','MarkerEdgeColor','k','DisplayName','Bac');
        legend([qw{:}],{'Fish','Jelly','Zoo','Protist','Phyto','Bac'},'Location','southeast')
    end
end