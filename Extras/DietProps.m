clear all
close all

%%Calculating proportion of each heterotrophs diet made up by each other
%%group. Requires an external Excel file called Feeding Matrix that
%%describes the consumptive portions of model flows, i.e. [column]
%%consumption of [row] is represented by flow [value] in the mass balance 
%%matrix

plot=1; %1 to make and save plots, 0 for just data output

if plot==1
    f1=figure('units','inches','Position',[0 0 6 7])%Salp and Fish Diet Prop
    f2=figure('units','inches','Position',[0 0 6 3])%Salp and Fish Diet Flux
    f3=figure('units','inches','Position',[0 0 1 10])%Legend
    f4=figure('units','inches','Position',[0 0 6 7])%Protistan Diet Prop
    f5=figure('units','inches','Position',[0 0 6 3])%Protistan Diet Flux
    f6=figure('units','inches','Position',[0 0 6 7])%Zoo Diet Prop
    f7=figure('units','inches','Position',[0 0 6 3])%Zoo Diet Flux
end

for Cycle=[1:5] 
    ['Cycle ' num2str(Cycle)]
    clearvars -except Cycle f1 f2 f3 f4 f5 f6 f7 DietPropMean DietPropStd DietPropLCI DietPropUCI ...
        DietFluxMean DietFluxStd DietFluxLCI DietFluxUCI HerbProp ProProp MesoProp MacroProp NPPProp plot
    load(['N15NZInverseCycle',num2str(Cycle),'.mat']);
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat']);
    
    %Loading Feeding Matrix, an excel relating Ae to flows (each column
    %gives the elements in flows that corresponds to the consumptive event
    %on the compartment in the respective row)
    %Note that for C2-5 a dummy column/row was added to prevent the
    %spreadsheet ending on a Z column because that somehow breaks ExcelCol
    sheet=['FeedingMatrix.xlsx'];
    if Cycle==1
        workbook='C1';
    else
        workbook='C2_5';
    end
    datsize=max(xlsread(sheet,workbook,'B2:B27'));
    rskip=2;    %rows to skip before data
    cskip=2;    %cols to skip before data
    ccl=cskip+1;	% col of upper left number to be read
    crl=rskip+1;    % row of upper left number to be read
    crr=rskip+datsize; %row of bottom right number to be read
    ccr=cskip+datsize;  %col of bottom right number to be read
    cread=[char(ExcelCol(ccl)),num2str(crl),':',char(ExcelCol(ccr)),num2str(crr)]; %Excel code for data to be read in
    [~,~,FeedingMat]=xlsread(sheet,workbook,cread); %The actual input matrix
    FeedingMat=cell2mat(FeedingMat);


    num=length(MCMCmat(:,1));
    MCMCmat(1:round(num*0.2),:)=[]; %Removing the first 20% of sims
    del15N(1:round(num*0.2),:)=[];
    len=length(MCMCmat(:,1));

    if Cycle==1
        [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,Macro_col,...
            VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
            mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumnsC1(Aa0,MCMCmat,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
        VMMeso_col=NaN;
    else
        [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,...
            VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
            mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumns(Aa0,MCMCmat,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
    end

    %Calcing Diet Proportion and total flux for every group in each MC sim
    FromNO3 = find(Ae0(NO3_col,:)==-1);
    FromNH4 = find(Ae0(NH4_col,:)==-1);
    ToPico = find(Ae0(Pico_col,:)==1);
    FromPico = find(Ae0(Pico_col,:)==-1);
    ToDtm = find(Ae0(Dtm_col,:)==1);
    FromDtm = find(Ae0(Dtm_col,:)==-1);
    ToFlag = find(Ae0(Flag_col,:)==1);
    FromFlag = find(Ae0(Flag_col,:)==-1);
    ToDom = find(Ae0(Dom_col,:)==1);
    NUTtoFLAG = intersect([FromNO3,FromNH4],ToFlag);
    NUTtoDTM = intersect([FromNO3,FromNH4],ToDtm);
    NUTtoPICO = intersect([FromNO3,FromNH4],ToPico);
    FLAGtoDON = intersect(FromFlag,ToDom);
    PICOtoDON = intersect(FromPico,ToDom);
    DTMtoDON = intersect(FromDtm,ToDom);
    for i=1:length(MCMCmat(:,1))
        flows=MCMCmat(i,:)./wts';
        [DietProp0,DietFlux0] = CalcDiet(Cycle,flows,del15N,d15NInputs,Ae0,InputCol,FeedingMat);
        DietProp_track(:,:,i)=DietProp0;
        DietFlux_track(:,:,i)=DietFlux0;
        %A few special summed measurements that have to be calculated
        %independently to get correct CIs. HerbProp is the proportion of
        %each groups diet made up by herbivory, and so on
        NPP(i) = sum(flows([NUTtoPICO NUTtoDTM NUTtoFLAG])) - sum(flows([PICOtoDON DTMtoDON FLAGtoDON]));
        HerbFlux_track(i,:) = sum(DietFlux0([Pico_col Dtm_col Flag_col],:),1);
        HerbProp_track(i,:) = HerbFlux_track(i,:)./sum(DietFlux0(:,:),1)*100;
        NPPProp_track(i,:) = HerbFlux_track(i,:)./NPP(i)*100;
        ProFlux_track(i,:) = sum(DietFlux0([HNF_col Mic_col],:),1);
        ProProp_track(i,:) = ProFlux_track(i,:)./sum(DietFlux0(:,:),1)*100;
        if Cycle == 1
            MesoFlux_track(i,:) = sum(DietFlux0(Meso_col,:),1);
        else
            MesoFlux_track(i,:) = sum(DietFlux0([Meso_col VMMeso_col],:),1);
        end
        MesoProp_track(i,:) = MesoFlux_track(i,:)./sum(DietFlux0(:,:),1)*100;
        if Cycle == 1
            MacroFlux_track(i,:) = sum(DietFlux0([Macro_col+1 VMMacro_col+1],:),1);
            MacroProp_track(i,:) = MacroFlux_track(i,:)./sum(DietFlux0(:,:),1)*100;
        else
            MacroFlux_track(i,:) = sum(DietFlux0([Macro_col VMMacro_col],:),1);
            MacroProp_track(i,:) = MacroFlux_track(i,:)./sum(DietFlux0(:,:),1)*100;
        end
    end
    DietPropSort = sort(DietProp_track,3);
    DietFluxSort = sort(DietFlux_track,3);
    DietPropMean(:,:,Cycle)=mean(DietProp_track,3);
    DietPropStd(:,:,Cycle)=std(DietProp_track,0,3);
    DietPropLCI(:,:,Cycle)=DietPropSort(:,:,round(0.025*len));
    DietPropUCI(:,:,Cycle)=DietPropSort(:,:,round(0.975*len));
    DietFluxMean(:,:,Cycle)=mean(DietFlux_track,3);
    DietFluxStd(:,:,Cycle)=std(DietFlux_track,0,3);
    DietFluxLCI(:,:,Cycle)=DietFluxSort(:,:,round(0.025*len));
    DietFluxUCI(:,:,Cycle)=DietFluxSort(:,:,round(0.975*len));
    %for special measurements these will be saved as smaller matrices
    HerbProp_track = sort(HerbProp_track); 
    HerbProp(1,:,Cycle) = mean(HerbProp_track);  %mean
    HerbProp(2,:,Cycle) = HerbProp_track(round(0.025*len),:); %LCI
    HerbProp(3,:,Cycle) = HerbProp_track(round(0.975*len),:); %UCI
    NPPProp_track = sort(NPPProp_track);
    NPPProp(1,:,Cycle) = mean(NPPProp_track);
    NPPProp(2,:,Cycle) = NPPProp_track(round(0.025*len),:);
    NPPProp(3,:,Cycle) = NPPProp_track(round(0.975*len),:);
    ProProp_track = sort(ProProp_track);
    ProProp(1,:,Cycle) = mean(ProProp_track);
    ProProp(2,:,Cycle) = ProProp_track(round(0.025*len),:);
    ProProp(3,:,Cycle) = ProProp_track(round(0.975*len),:);
    MesoProp_track = sort(MesoProp_track);
    MesoProp(1,:,Cycle) = mean(MesoProp_track);
    MesoProp(2,:,Cycle) = MesoProp_track(round(0.025*len),:);
    MesoProp(3,:,Cycle) = MesoProp_track(round(0.975*len),:);
    MacroProp_track = sort(MacroProp_track);
    MacroProp(1,:,Cycle) = mean(MacroProp_track);
    MacroProp(2,:,Cycle) = MacroProp_track(round(0.025*len),:);
    MacroProp(3,:,Cycle) = MacroProp_track(round(0.975*len),:);



    if plot == 1
        %Assigning colors for each compartment column
        %Initially in rgb triplets
        %To assign colors to plots the color matrix must include entries for
        %all compartments, even those I do not wish to plot
        col=[255 255 255; %NO3-white
            255 255 255; %NH4-white
            51 153 51;%Pico-dark green
            51 204 51;%Dtm-green
            102 255 153;%Flag-light green
            102 0 204;%Hnf-dark purple
            153 0 255;%Mic-purple
            255, 204, 0;%Meso-yellow
            255 153 51;%VMMeso-yellow orange
            204 102 0;%Macro-orange
            204 51 0;%VMMacro-orange red
            0 153 204;%Gel-blue
            0 204 255;%Salp-light blue
            5 16 108;%Amph-dark blue
            102 0 51;%PiscivFish-dark pink
            255 0 102;%GelativFish-neon pink
            204 102 153;%Myct-light pink
            120 120 120;%Bac-light grey
            153 153 102;%Sdet-light beige
            92 92 61;%Mdet-dark beige
            255 255 255;%Ldet-White
            255 255 255;%SalpDet-White
            255 255 255;%DOM-white
            0 0 0;%SalpMort-black
            156 12 12];%HTL-dark red
        col=col/255; %Converts to rgb 0-1

        %Plotting Salp/Fish Diet prop
        figure(f1)
        if Cycle==1
            text(50,0,'Salps')
            hold on
            text(50,6,'Myctophids')
            text(50,12,'Piscivorous Fish')
            text(50,18,'Gelativorous Fish')
            b=barh(Cycle,DietPropMean(:,Salp_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Salp_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            b=barh(Cycle+6,DietPropMean(:,Myct_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Myct_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+12,DietPropMean(:,PiscivFish_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,PiscivFish_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+18,DietPropMean(:,GelativFish_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,GelativFish_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
        else
            b=barh(Cycle,DietPropMean(:,Salp_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Salp_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            b=barh(Cycle+6,DietPropMean(:,Myct_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Myct_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+12,DietPropMean(:,PiscivFish_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,PiscivFish_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+18,DietPropMean(:,GelativFish_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,GelativFish_col,Cycle),1)
                b(k).CData = col(k,:);
            end
        end

        ylim([-1 24])
        xlim([0 100])
        set(gca,'YTick',[1:5 7:11 13:17 19:23])
        set(gca,'YTickLabel',{'Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5'})
        set(gca,'FontSize',10)
        set(gca,'YDir','reverse')
        set(gca,'XAxisLocation','top')
        %legend('on')

        %Absolute flux plot, split into 4 horizontal subplots
        figure(f2)
        if Cycle==1
            subplot(1,4,1)
            b=bar(Cycle,DietFluxMean(:,Salp_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Salp_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[3])
            set(gca,'XTickLabel',{'Salp'})
            set(gca,'FontSize',10)
            ylabel('Diet (mmol N m^-^2 d^-^1)')

            subplot(1,4,2)
            b=bar(Cycle+6,DietFluxMean(:,Myct_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Myct_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[9])
            set(gca,'XTickLabel',{'Myct'})
            set(gca,'FontSize',10)

            subplot(1,4,3)
            b=bar(Cycle+12,DietFluxMean(:,PiscivFish_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,PiscivFish_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[15])
            set(gca,'XTickLabel',{'Pisciv'})
            set(gca,'FontSize',10)

            subplot(1,4,4)
            b=bar(Cycle+18,DietFluxMean(:,GelativFish_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,GelativFish_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[21])
            set(gca,'XTickLabel',{'Gelativ'})
            set(gca,'FontSize',10)
        else
            subplot(1,4,1)
            b=bar(Cycle,DietFluxMean(:,Salp_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Salp_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[3])
            set(gca,'XTickLabel',{'Salp'})
            set(gca,'FontSize',10)
            ylabel('Diet (mmol N m^-^2 d^-^1)')

            subplot(1,4,2)
            b=bar(Cycle+6,DietFluxMean(:,Myct_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Myct_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[9])
            set(gca,'XTickLabel',{'Myct'})
            set(gca,'FontSize',10)

            subplot(1,4,3)
            b=bar(Cycle+12,DietFluxMean(:,PiscivFish_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,PiscivFish_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[15])
            set(gca,'XTickLabel',{'Pisciv'})
            set(gca,'FontSize',10)

            subplot(1,4,4)
            b=bar(Cycle+18,DietFluxMean(:,GelativFish_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,GelativFish_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[21])
            set(gca,'XTickLabel',{'Gelativ'})
            set(gca,'FontSize',10)
        end

        %Protistan Diet Prop
        figure(f4)
        if Cycle==1
            text(50,0,'Mixotrophic Flagellates')
            hold on
            text(50,6,'Heterotrophic Flagellates')
            text(50,12,'Microzoo')
            b=barh(Cycle,DietPropMean(:,Flag_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Flag_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            b=barh(Cycle+6,DietPropMean(:,HNF_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,HNF_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+12,DietPropMean(:,Mic_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Mic_col,Cycle),1)
                b(k).CData = col(k,:);
            end
        else
            b=barh(Cycle,DietPropMean(:,Flag_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Flag_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            b=barh(Cycle+6,DietPropMean(:,HNF_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,HNF_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+12,DietPropMean(:,Mic_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Mic_col,Cycle),1)
                b(k).CData = col(k,:);
            end
        end

        ylim([-1 18])
        xlim([0 100])
        set(gca,'YTick',[1:5 7:11 13:17])
        set(gca,'YTickLabel',{'Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5'})
        set(gca,'FontSize',10)
        set(gca,'YDir','reverse')
        set(gca,'XAxisLocation','top')
        %legend('on')

        %Absolute flux plot, split into 4 horizontal subplots
        figure(f5)
        subplot(1,3,1)
        b=bar(Cycle,DietFluxMean(:,Flag_col,Cycle),'stacked','FaceColor','flat')
        for k = 1:size(DietFluxMean(:,Flag_col,Cycle),1)
            b(k).CData = col(k,:);
        end
        hold on
        set(gca,'XTick',[3])
        set(gca,'XTickLabel',{'MixFlag'})
        set(gca,'FontSize',10)
        ylabel('Diet (mmol N m^-^2 d^-^1)')

        subplot(1,3,2)
        b=bar(Cycle+6,DietFluxMean(:,HNF_col,Cycle),'stacked','FaceColor','flat')
        for k = 1:size(DietFluxMean(:,HNF_col,Cycle),1)
            b(k).CData = col(k,:);
        end
        hold on
        set(gca,'XTick',[9])
        set(gca,'XTickLabel',{'HNF'})
        set(gca,'FontSize',10)

        subplot(1,3,3)
        b=bar(Cycle+12,DietFluxMean(:,Mic_col,Cycle),'stacked','FaceColor','flat')
        for k = 1:size(DietFluxMean(:,Mic_col,Cycle),1)
            b(k).CData = col(k,:);
        end
        hold on
        set(gca,'XTick',[15])
        set(gca,'XTickLabel',{'Micro'})
        set(gca,'FontSize',10)


        %Zoo Diet Prop
        figure(f6)
        if Cycle==1
            text(50,0,'Mesozoo')
            hold on
            text(50,6,'Macrozoo')
            text(50,12,'GelPred')
            text(50,18,'Amphipods')
            b=barh(Cycle,DietPropMean(:,Meso_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Meso_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            val=mean(DietPropMean(:,[Macro_col+1 VMMacro_col+1],Cycle),2);
            b=barh(Cycle+6,val,'stacked','FaceColor','flat')
            for k = 1:length(val)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+12,DietPropMean(:,Gel_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Gel_col+1),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+18,DietPropMean(:,Amph_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Amph_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
        else
            val=mean(DietPropMean(:,[Meso_col VMMeso_col],Cycle),2);
            b=barh(Cycle,val,'stacked','FaceColor','flat')
            for k = 1:length(val)
                b(k).CData = col(k,:);
            end
            hold on
            val=mean(DietPropMean(:,[Macro_col VMMacro_col],Cycle),2);
            b=barh(Cycle+6,val,'stacked','FaceColor','flat')
            for k = 1:length(val)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+12,DietPropMean(:,Gel_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Gel_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            b=barh(Cycle+18,DietPropMean(:,Amph_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietPropMean(:,Amph_col,Cycle),1)
                b(k).CData = col(k,:);
            end
        end

        ylim([-1 24])
        xlim([0 100])
        set(gca,'YTick',[1:5 7:11 13:17 19:23])
        set(gca,'YTickLabel',{'Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5','Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5'})
        set(gca,'FontSize',10)
        set(gca,'YDir','reverse')
        set(gca,'XAxisLocation','top')
        %legend('on')

        %Zoo Diet Flux
        figure(f7)
        if Cycle==1
            subplot(1,4,1)
            b=bar(Cycle,DietFluxMean(:,Meso_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Meso_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[3])
            set(gca,'XTickLabel',{'Meso'})
            set(gca,'FontSize',10)
            ylabel('Diet (mmol N m^-^2 d^-^1)')

            subplot(1,4,2)
            val=sum(DietFluxMean(:,[Macro_col+1 VMMacro_col+1],Cycle),2);
            b=bar(Cycle+6,val,'stacked','FaceColor','flat')
            for k = 1:length(val)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[9])
            set(gca,'XTickLabel',{'Macro'})
            set(gca,'FontSize',10)

            subplot(1,4,3)
            b=bar(Cycle+12,DietFluxMean(:,Gel_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Gel_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[15])
            set(gca,'XTickLabel',{'GelPred'})
            set(gca,'FontSize',10)

            subplot(1,4,4)
            b=bar(Cycle+18,DietFluxMean(:,Amph_col+1,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Amph_col+1,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[21])
            set(gca,'XTickLabel',{'Amph'})
            set(gca,'FontSize',10)
        else
            subplot(1,4,1)
            val=sum(DietFluxMean(:,[Meso_col VMMeso_col],Cycle),2);
            b=bar(Cycle,val,'stacked','FaceColor','flat')
            for k = 1:length(val)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[3])
            set(gca,'XTickLabel',{'Meso'})
            set(gca,'FontSize',10)
            ylabel('Diet (mmol N m^-^2 d^-^1)')

            subplot(1,4,2)
            val=sum(DietFluxMean(:,[Macro_col VMMacro_col],Cycle),2);
            b=bar(Cycle+6,val,'stacked','FaceColor','flat')
            for k = 1:length(val)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[9])
            set(gca,'XTickLabel',{'Macro'})
            set(gca,'FontSize',10)

            subplot(1,4,3)
            b=bar(Cycle+12,DietFluxMean(:,Gel_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Gel_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[15])
            set(gca,'XTickLabel',{'GelPred'})
            set(gca,'FontSize',10)

            subplot(1,4,4)
            b=bar(Cycle+18,DietFluxMean(:,Amph_col,Cycle),'stacked','FaceColor','flat')
            for k = 1:size(DietFluxMean(:,Amph_col,Cycle),1)
                b(k).CData = col(k,:);
            end
            hold on
            set(gca,'XTick',[21])
            set(gca,'XTickLabel',{'Amph'})
            set(gca,'FontSize',10)
        end

        %Legend as a separate figure
        col=[51 153 51;%Pico-dark green
            51 204 51;%Dtm-green
            102 255 153;%Flag-light green
            102 0 204;%Hnf-dark purple
            153 0 255;%Mic-purple
            255, 204, 0;%Meso-yellow
            255 153 51;%VMMeso-yellow orange
            204 102 0;%Macro-orange
            204 51 0;%VMMacro-orange red
            0 153 204;%Gel-blue
            0 204 255;%Salp-light blue
            5 16 108;%Amph-dark blue
            102 0 51;%PiscivFish-dark pink
            255 0 102;%GelativFish-neon pink
            204 102 153;%Myct-light pink
            120 120 120;%Bac-dark grey
            153 153 102;%Sdet-light beige
            92 92 61;%Mdet-dark beige
            0 0 0;%SalpMort-black
            156 12 12];%HTL-dark red
        col=col/255; %Converts to rgb 0-1
        if Cycle==1
            figure(f3)
            for i=1:(height(col))
                x=41-(2*i);
                r=rectangle('Position',[1 x 1 1])
                hold on
                r.FaceColor=col(i,:)
            end
            %ylim([0 37])
            xlim([0 3])
            axis off
        end
    end
end

%%Saving outputs as .mat files
save('DietFluxMean.mat','DietFluxMean');
save('DietFluxStd.mat','DietFluxStd');
save('DietFluxLCI.mat','DietFluxLCI');
save('DietFluxUCI.mat','DietFluxUCI');
save('DietPropMean.mat','DietPropMean');
save('DietPropStd.mat','DietPropStd');
save('DietPropLCI.mat','DietPropLCI');
save('DietPropUCI.mat','DietPropUCI');
save('HerbProp.mat','HerbProp');
save('NPPProp.mat','NPPProp');
save('ProProp.mat','ProProp');
save('MesoProp.mat','MesoProp')
save('MacroProp.mat','MacroProp')
if plot == 1
    exportgraphics(f1,'DietProp.png','Resolution',300)
    exportgraphics(f2,'DietFlux.png','Resolution',300)
    exportgraphics(f3,'DietLegend.png','Resolution',300)
    exportgraphics(f5,'ProDietProp.png','Resolution',300)
    exportgraphics(f7,'ZooDietProp.png','Resolution',300)
end