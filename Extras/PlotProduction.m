clear all
close all

%%A rather messy script that produces a variety of plots loosely related to
%%primary and secondary production. Requires PlotIndirect.m to be run first

f1=figure('Position',[100 100 1000 500]);
f2=figure('units','inches','Position',[0 0 6 3]);
f3=figure('units','inches','Position',[0 0 6 3]);
f4=figure('units','inches','Position',[0 0 6 3]);
f5=figure('Position',[100 100 1000 500]);
f6=figure('Position',[100 100 1000 500]);
f7=figure('Position',[100 100 1000 500]);
f8=figure('Position',[100 100 1000 500]);
f9=figure('Position',[100 100 1000 500]);
f10=figure('Position',[100 100 1000 500]);
f11=figure('Position',[100 100 1000 500]);
f12=figure('Position',[100 100 1000 500]);
load('GGEavg.mat'); %loading mean terminal pathway GGEs produced by PlotIndirect.m
load("prodmean.mat"); %Loading production matrix produced by PlotIndirect.m
for Cycle = 1:5 %[1:5]
    clearvars -except Cycle f1 f2 f3 f4 f5 f6 f7 f8 f9 f10 f11 f12 ProdMean ...
    ProdSort ProdLCI ProdUCI ProdStd SalpOutPropMean SalpOutPropStd ZooOutPropMean ZooOutPropStd ...
    GGEall GGEavg prodmean EEall EEpreyall TLEEall GGEpathall PicoPropall ProdMeanFront ProdMeanBack
    load(['N15NZInverseCycle',num2str(Cycle),'.mat']);
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat']);
    load(['TL_',num2str(Cycle),'.mat']); %Full array of Trophic Levels of each group
    TL_track=TL_track';

    num=length(MCMCmat(:,1));
    MCMCmat(1:round(num*0.2),:)=[]; %Removing the first 20% of sims
    %MCMCmat(1:round(length(MCMCmat(:,1))*0.5),:)=[]; %Front half of that
    %MCMCmat(round(length(MCMCmat(:,1))*0.5):round(length(MCMCmat(:,1))),:)=[]; %Back half
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

    %Denoting important flow locations
    ToNO3 = find(Ae0(NO3_col,:)==1);
    FromNO3 = find(Ae0(NO3_col,:)==-1);
    ToNH4 = find(Ae0(NH4_col,:)==1);
    FromNH4 = find(Ae0(NH4_col,:)==-1);
    ToPico = find(Ae0(Pico_col,:)==1);
    FromPico = find(Ae0(Pico_col,:)==-1);
    ToDtm = find(Ae0(Dtm_col,:)==1);
    FromDtm = find(Ae0(Dtm_col,:)==-1);
    ToFlag = find(Ae0(Flag_col,:)==1);
    FromFlag = find(Ae0(Flag_col,:)==-1);
    ToHNF = find(Ae0(HNF_col,:)==1);
    FromHNF = find(Ae0(HNF_col,:)==-1);
    ToMic = find(Ae0(Mic_col,:)==1);
    FromMic = find(Ae0(Mic_col,:)==-1);
    ToMeso = find(Ae0(Meso_col,:)==1);
    FromMeso = find(Ae0(Meso_col,:)==-1);
    if Cycle~=1
        ToVMMeso = find(Ae0(VMMeso_col,:)==1);
        FromVMMeso = find(Ae0(VMMeso_col,:)==-1);
    else
        ToVMMeso = NaN;
        FromVMMeso = NaN;
    end
    ToMacro = find(Ae0(Macro_col,:)==1);
    FromMacro = find(Ae0(Macro_col,:)==-1);
    ToVMMacro = find(Ae0(VMMacro_col,:)==1);
    FromVMMacro = find(Ae0(VMMacro_col,:)==-1);
    ToGel = find(Ae0(Gel_col,:)==1);
    FromGel = find(Ae0(Gel_col,:)==-1);
    ToSalp = find(Ae0(Salp_col,:)==1);
    FromSalp = find(Ae0(Salp_col,:)==-1);
    ToSalpMort = find(Ae0(SalpMort_col,:)==1);
    FromSalpMort = find(Ae0(SalpMort_col,:)==-1);
    ToAmph = find(Ae0(Amph_col,:)==1);
    FromAmph = find(Ae0(Amph_col,:)==-1);
    ToPisciv = find(Ae0(PiscivFish_col,:)==1);
    FromPisciv = find(Ae0(PiscivFish_col,:)==-1);
    ToGelativ = find(Ae0(GelativFish_col,:)==1);
    FromGelativ = find(Ae0(GelativFish_col,:)==-1);
    ToMyct = find(Ae0(Myct_col,:)==1);
    FromMyct = find(Ae0(Myct_col,:)==-1);
    ToBac = find(Ae0(bac_col,:)==1);
    FromBac = find(Ae0(bac_col,:)==-1);
    ToSdet = find(Ae0(sdet_col,:)==1);
    FromSdet = find(Ae0(sdet_col,:)==-1);
    ToMdet = find(Ae0(mdet_col,:)==1);
    FromMdet = find(Ae0(mdet_col,:)==-1);
    ToLdet = find(Ae0(ldet_col,:)==1);
    FromLdet = find(Ae0(ldet_col,:)==-1);
    ToSalpDet = find(Ae0(SalpDet_col,:)==1);
    FromSalpDet = find(Ae0(SalpDet_col,:)==-1);
    ToDom = find(Ae0(Dom_col,:)==1);
    FromDom = find(Ae0(Dom_col,:)==-1);
    FromPhy = sort([find(Ae0(Pico_col,:)==-1),find(Ae0(Dtm_col,:)==-1)]);
    FromFlag = find(Ae0(Flag_col,:)==-1);
    FromProtist = sort([find(Ae0(HNF_col,:)==-1),find(Ae0(Mic_col,:)==-1)]);
    ToProtist = sort([find(Ae0(HNF_col,:)==1),find(Ae0(Mic_col,:)==1)]);

    FLAGtoDON = intersect(FromFlag,ToDom);
    FLAGtoSdet = intersect(FromFlag,ToSdet);
    PICOtoDON = intersect(FromPico,ToDom);
    DTMtoDON = intersect(FromDtm,ToDom);
    NUTtoFLAG = intersect([FromNO3,FromNH4],ToFlag);
    NUTtoDTM = intersect([FromNO3,FromNH4],ToDtm);
    NUTtoPICO = intersect([FromNO3,FromNH4],ToPico);
    BACtoHTL = intersect(FromBac,[ToFlag,ToHNF,ToMic]);
    HNFtoHTL = intersect(FromHNF,[ToMic,ToMeso,ToVMMeso,ToMacro,ToVMMacro,ToSalp]);
    HNFtoSdet = intersect(FromHNF,ToSdet);
    MICtoHTL = intersect(FromMic,[ToMeso,ToVMMeso,ToMacro,ToVMMacro,ToSalp]);
    MICtoSdet = intersect(FromMic,ToSdet);
    MESOtoHTL = intersect(FromMeso,[ToMacro,ToVMMacro,ToGel,ToMyct]);
    MESOtoMdet = intersect(FromMeso,ToMdet);
    if Cycle~=1
        VMMESOtoHTL = intersect(FromVMMeso,[ToMacro,ToVMMacro,ToGel,ToMyct]);
        VMMESOtoMdet = intersect(FromVMMeso,ToMdet);
    end
    MACROtoHTL =  intersect(FromMacro,[ToPisciv,ToGel,ToMyct]);
    MACROtoLdet = intersect(FromMacro,ToLdet);
    VMMACROtoHTL =  intersect(FromVMMacro,[ToPisciv,ToGel,ToMyct]);
    VMMACROtoLdet = intersect(FromVMMacro,ToLdet);
    GELtoHTL = intersect(FromGel,[ToAmph,ToMyct,ToGelativ]);
    GELtoLdet = intersect(FromGel,ToLdet);
    MYCTtoHTL = intersect(FromMyct,ToPisciv);
    SALPtoHTL = intersect(FromSalp,[ToAmph,ToMyct]);
    SALPtoSalpDet = intersect(FromSalp,ToSalpDet);
    AMPHtoHTL = intersect(FromAmph,[ToPisciv ToMyct]);
    AMPHtoLdet = intersect(FromAmph,ToLdet);
    %Non-predatory secondary production as natural mortality and growth
    %Also egestion of fish as loss terms, Meso cannabalism, and predation
    %by unmodelled groups
    if Cycle==1
        SalpTOSalpMort = intersect(FromSalp,ToSalpMort);
        MesoTOMeso = 48;
        MesoTOPred = 53;
        MacroTOPred = 59;
        VMMacroTOPred = 66;
        gelTOPred = 75;
        MyctTOBiom=85;
        MyctTOPoop = 84;
        SalpTOPred = 92;
        AmphTOPred = 99;
        PiscivTOBiom = 106;
        PiscivTOPoop = 105;
        GelativTOBiom = 110;
        GelativTOPoop = 109;
    else
        SalpTOSalpMort = intersect(FromSalp,ToSalpMort);
        MesoTOMeso = 52;
        MesoTOVMMeso =53;
        MesoTOPred = 58;
        VMMesoTOMeso = 62;
        VMMesoTOVMMeso = 63;
        VMMesoTOPred = 68;
        MacroTOPred = 76;
        VMMacroTOPred = 83;
        gelTOPred = 92;
        MyctTOBiom = 102;
        MyctTOPoop = 101;
        SalpTOPred = 109;
        AmphTOPred = 116;
        PiscivTOBiom = 123;
        PiscivTOPoop = 122;
        GelativTOBiom = 127;
        GelativTOPoop = 126;
    end


    Prod = zeros(length(MCMCmat),16);
    GGE = zeros(length(MCMCmat),16);
    AE = zeros(length(MCMCmat),16);
    Columns = {'Pico','Dtm','Flag','Bac','HNF','Mic','Meso','VMMeso','Macro','VMMacro','Gel','Salp','Amph','Pisciv','Gelativ','Myct'};
    for i=1:len
        flows=MCMCmat(i,:)./wts';

        Prod(i,1) = sum(flows(NUTtoPICO)) - flows(PICOtoDON);
        Prod(i,2) = sum(flows(NUTtoDTM)) - flows(DTMtoDON);
        Prod(i,3) = sum(flows(NUTtoFLAG)) - flows(FLAGtoDON); %Note this is Flag NPP, does not include secondary production from heterotrophy

        Prod(i,4) = sum(flows(BACtoHTL));
        Prod(i,5) = sum(flows(HNFtoHTL));
        Prod(i,6) = sum(flows(MICtoHTL));
        if Cycle~=1
            Prod(i,7) = sum(flows([MesoTOMeso,MesoTOVMMeso,MESOtoHTL,MesoTOPred]));
            Prod(i,8) = sum(flows([VMMesoTOMeso,VMMesoTOVMMeso,VMMESOtoHTL,VMMesoTOPred]));
        else
            Prod(i,7) = sum(flows([MesoTOMeso,MESOtoHTL,MesoTOPred]));
            Prod(i,8) = NaN;
        end
        Prod(i,9) = sum(flows([MACROtoHTL,MacroTOPred]));
        Prod(i,10) = sum(flows([VMMACROtoHTL,VMMacroTOPred]));
        Prod(i,11) = sum(flows([GELtoHTL,gelTOPred]));
        Prod(i,12) = sum(flows([SALPtoHTL,SalpTOPred,SalpTOSalpMort]));
        Prod(i,13) = sum(flows([AMPHtoHTL,AmphTOPred]));
        Prod(i,14) = sum(flows(PiscivTOBiom));
        Prod(i,15) = sum(flows(GelativTOBiom));
        Prod(i,16) = sum(flows([MYCTtoHTL,MyctTOBiom]));
        Prod(i,17) = flows(MyctTOBiom); %Extra column for Myct biomass alone
        if Cycle ==1  %Extra column for summed consumption of unmodelled HTLs
            Prod(i,18) = sum(flows([MesoTOPred,MacroTOPred,VMMacroTOPred,gelTOPred,SalpTOPred,AmphTOPred]));
        else
            Prod(i,18) = sum(flows([MesoTOPred,VMMesoTOPred,MacroTOPred,VMMacroTOPred,gelTOPred,SalpTOPred,AmphTOPred]));
        end
        %Gross Growth Efficiency
        GGE(i,1:3) = NaN;
        GGE(i,4) = Prod(i,4)/sum(flows(ToBac));
        GGE(i,5) = Prod(i,5)/sum(flows(ToHNF));
        GGE(i,6) = Prod(i,6)/sum(flows(ToMic));
        GGE(i,7) = Prod(i,7)/sum(flows([ToMeso MesoTOMeso]));
        if Cycle ~= 1
            GGE(i,8) = Prod(i,8)/sum(flows([ToVMMeso VMMesoTOVMMeso]));
        else
            GGE(i,8) = NaN;
        end
        GGE(i,9) = Prod(i,9)/sum(flows(ToMacro));
        GGE(i,10) = Prod(i,10)/sum(flows(ToVMMacro));
        GGE(i,11) = Prod(i,11)/sum(flows(ToGel));
        GGE(i,12) = Prod(i,12)/sum(flows(ToSalp));
        GGE(i,13) = Prod(i,13)/sum(flows(ToAmph));
        GGE(i,14) = Prod(i,14)/sum(flows(ToPisciv));
        GGE(i,15) = Prod(i,15)/sum(flows(ToGelativ));
        GGE(i,16) = Prod(i,16)/sum(flows(ToMyct));

        %Assimilation Efficiency
        AE(i,1:4) = NaN;
        AE(i,5) = 1-(sum(flows(HNFtoSdet))/sum(flows(ToHNF)));
        AE(i,6) = 1-(sum(flows(MICtoSdet))/sum(flows(ToMic)));
        AE(i,7) = 1-(sum(flows(MESOtoMdet))/sum(flows([ToMeso MesoTOMeso])));
        if Cycle ~= 1
            AE(i,8) = 1-(sum(flows(VMMESOtoMdet))/sum(flows([ToVMMeso VMMesoTOVMMeso])));
        else
            AE(i,8) = NaN;
        end
        AE(i,9) = 1-(sum(flows(MACROtoLdet))/sum(flows(ToMacro)));
        AE(i,10) = 1-(sum(flows(VMMACROtoLdet))/sum(flows(ToVMMacro)));
        AE(i,11) = 1-(sum(flows(GELtoLdet))/sum(flows(ToGel)));
        AE(i,12) = 1-(sum(flows(SALPtoSalpDet))/sum(flows(ToSalp)));
        AE(i,13) = 1-(sum(flows(AMPHtoLdet))/sum(flows(ToAmph)));
        AE(i,14) = 1-(sum(flows(PiscivTOPoop))/sum(flows(ToPisciv)));
        AE(i,15) = 1-(sum(flows(GelativTOPoop))/sum(flows(ToGelativ)));
        AE(i,16) = 1-(sum(flows(MyctTOPoop))/sum(flows(ToMyct)));

        %Salp outputs
        %Amph,Myct,Gelativ,HTL,mort,nh4,dom,det
        SalpOut(i,1:2) = flows(SALPtoHTL);
        SalpOut(i,3) = flows(intersect(FromSalpMort,ToGelativ));
        SalpOut(i,4) = flows(SalpTOPred);
        SalpOut(i,5) = flows(FromSalpMort(2));
        SalpOut(i,6:7) = flows(FromSalp(3:4)) + flows(FromSalp(8:9));
        SalpOut(i,8) = flows(intersect(FromSalp,ToSalpDet));

        %Macrozoo outputs
        NVM = [0 flows(FromMacro)];
        VM = flows(FromVMMacro(1:7));
        VM(5:6) = VM(5:6) + flows(FromVMMacro(8:9)); %adding deep nh4/dom to shallow
        ZooOut(i,1:width(NVM)) = sum([NVM; VM]); %adding nvm and vm

        if Cycle == 1
            TL_htl(i) = ((flows(MesoTOPred)*TL_track(i,Meso_col) + ...
                     flows(MacroTOPred)*TL_track(i,Macro_col) + flows(VMMacroTOPred)*TL_track(i,VMMacro_col) + ...
                     flows(gelTOPred)*TL_track(i,Gel_col) + flows(SalpTOPred)*TL_track(i,Salp_col) + ...
                     flows(AmphTOPred)*TL_track(i,Amph_col))/Prod(i,18))+1;
        else
            TL_htl(i) = ((flows(MesoTOPred)*TL_track(i,Meso_col) + flows(VMMesoTOPred)*TL_track(i,VMMeso_col) + ...
                     flows(MacroTOPred)*TL_track(i,Macro_col) + flows(VMMacroTOPred)*TL_track(i,VMMacro_col) + ...
                     flows(gelTOPred)*TL_track(i,Gel_col) + flows(SalpTOPred)*TL_track(i,Salp_col) + ...
                     flows(AmphTOPred)*TL_track(i,Amph_col))/Prod(i,18))+1;
        end
    end

    %Saving mean and CIs
    ProdSort=sort(Prod(:,1:16));
    ProdMean(Cycle,:)=mean(Prod(:,1:16));
    ProdStd(Cycle,:)=std(Prod(:,1:16));
    ProdLCI(Cycle,:)=ProdSort(round(0.025*len),:);
    ProdUCI(Cycle,:)=ProdSort(round(0.975*len),:);
    GGESort=sort(GGE);
    GGEall(1,:,Cycle)=mean(GGESort);
    GGEall(2,:,Cycle)=GGESort(round(0.025*len),:);
    GGEall(3,:,Cycle)=GGESort(round(0.975*len),:);

    %Phytoplankton community composition by their contribution to total
    %phyto production (note this differs from NPP as it includes the
    %heterotrophy of flagellates)
    NPP = sum(Prod(:,1:3),2);
    PicoProp = sort(Prod(:,1)./NPP);
    DtmProp = sort(Prod(:,2)./NPP);
    FlagProp = sort(Prod(:,1)./NPP);


    %Ecosystem Efficiency as Production normalized to NPP for the 3 fish
    %and 1 'Terminal' group that includes the biomass accumulation of each
    %fish as well as the consumption by unmodelled higher trophic levels
    %(columns 1-4 respectively)
    Cols = [find(ismember(Columns,'Myct')) find(ismember(Columns,'Pisciv')) find(ismember(Columns,'Gelativ'))];
    for i=1:(length(Cols)+1)
        NPP = sum(Prod(:,1:3),2);
        if i<(length(Cols)+1)
            production = sum(Prod(:,Cols(i)),2);
        elseif i == (length(Cols)+1)
            production = sum(Prod(:,[14 15 17 18]),2);
        end
        EE(:,i)=production./NPP;
    end
    EESort=sort(EE);
    for i=1:4
        EEall(Cycle,i,1)=mean(EESort(:,i));
        EEall(Cycle,i,2)=EESort(round(0.025*len),i);
        EEall(Cycle,i,3)=EESort(round(0.975*len),i);
    end
    %Trophic levels for each of the EE groups. Fish are easy
    TL_EE(:,1)=TL_track(:,Myct_col);
    TL_EE(:,2)=TL_track(:,PiscivFish_col);
    TL_EE(:,3)=TL_track(:,GelativFish_col);
    %Terminal is harder as it's average TL weighted by each groups
    %contribution to production
    denom = sum(Prod(:,[14 15 17 18]),2);
    num=Prod(:,17).*TL_track(:,Myct_col)+Prod(:,14).*TL_track(:,PiscivFish_col)+...
        Prod(:,15).*TL_track(:,GelativFish_col)+Prod(:,18).*TL_htl';
    TL_EE(:,4)=num./denom;
    TLEESort=sort(TL_EE);
    TLEEall(1,Cycle)=mean(TLEESort(:,4));
    TLEEall(2,Cycle)=TLEESort(round(0.025*len),4);
    TLEEall(3,Cycle)=TLEESort(round(0.975*len),4);
    GGEpathSort=sort(GGEavg(1:len,Cycle));
    GGEpathall(1,Cycle)=mean(GGEpathSort);
    GGEpathall(2,Cycle)=GGEpathSort(round(0.025*len));
    GGEpathall(3,Cycle)=GGEpathSort(round(0.975*len));
    PicoProp=Prod(:,1)./sum(Prod(:,1:3),2);
    PicoProp=sort(PicoProp);
    PicoPropall(1,Cycle)=mean(PicoProp);
    PicoPropall(2,Cycle)=PicoProp(round(0.025*len));
    PicoPropall(3,Cycle)=PicoProp(round(0.975*len));

    %Summed Production Normalized to NPP for the prey of each fish
    for i=1:5
        NPP = sum(Prod(:,1:3),2);
        if i==1 %Myct
            if Cycle==1
                production = sum(Prod(:,[7 9:13]),2);
            else
                production = sum(Prod(:,7:13),2);
            end
        elseif i==2 %Pisciv
            production = sum(Prod(:,[10 13 16]),2);
        elseif i==3 %Gelativ
            production = sum(Prod(:,[11 12]),2);
        elseif i==4 %All
            if Cycle==1
                production = sum(Prod(:,[7 9:13 16]),2);
            else
                production = sum(Prod(:,[7:13 16]),2);
            end
        elseif i==5 %Just salp production for comparison
            production = sum(Prod(:,12),2);
        end
        EEprey(:,i)=production./NPP;
    end
    EESort=sort(EEprey);
    for i=1:5
        EEpreyall(Cycle,i,1)=mean(EESort(:,i));
        EEpreyall(Cycle,i,2)=EESort(round(0.025*len),i);
        EEpreyall(Cycle,i,3)=EESort(round(0.975*len),i);
    end


    %Plotting Production as one large series of boxplots
    %boxplot(TL_track(Columns,:)','Symbol','')
    figure(f1)
    for i=1:width(Prod)
        if Cycle==1 || Cycle==2 || Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        %try
        boxplot(Prod(:,i),'Positions',x,'Widths',.8)
        hold on
        %end
    end
    set(gca,'XTick',3:5:width(Prod)*5+2)
    set(gca,'XTickLabel',Columns)
    set(gca,'XTickLabelRotation',45)
    set(gca,'FontSize',20)
    set(gca, 'YScale', 'log')
    ylabel('Production (mmol N m^-^2 d^-^1)')
    xline(5.5:5:75.5,'--')
    xlim([0, width(Prod)*5+1])
    hold on

    %GGE as violin plots
    %Plotting TL
    figure(f9)
    ymax=45;
    ymin=0;
    %boxplot(TL_track(Columns,:)','Symbol','')
    for i=4:length(Columns)
        if Cycle==1 || Cycle==2 || Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*((i-3)-1);
        ViolinPlot(x,0.5,GGE(:,i)*100,200,col)
    end
    set(gca,'XTick',3:5:width(GGE)*5+2)
    set(gca,'XTickLabel',{'Bac','Nanoflag','Microzoo','Meso (NVM)','Meso (VM)','Macro (NVM)','Macro (VM)','Gel Pred','Salps','Amph','PiscivFish','GelativFish','Myct'})
    set(gca,'XTickLabelRotation',45)
    ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('GGE (%)')
    xline(5.5:5:65.5,'--')
    xlim([0, (width(GGE)-3)*5+1])
    hold on


    %AE as violin plots
    %Plotting TL
    figure(f10)
    ymax=95;
    ymin=45;
    %boxplot(TL_track(Columns,:)','Symbol','')
    for i=5:length(Columns)
        % if i==5 | i==6 | i==9 | i==10 %Pisciv and it's prey
        %     col = [0 0 1];
        % elseif i==7 | i==8 | i==11  %Gelativ and it's prey
        %     col = [1 0 0];
        % else
        %     col = [0.8 0.8 0];  %Myct and it's prey? Hard to do
        % end
        if Cycle==1 || Cycle==2 || Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*((i-4)-1);
        ViolinPlot(x,0.5,AE(:,i)*100,200,col)
    end
    set(gca,'XTick',3:5:width(AE)*5+2)
    set(gca,'XTickLabel',{'Nanoflag','Microzoo','Meso (NVM)','Meso (VM)','Macro (NVM)','Macro (VM)','Gel Pred','Salps','Amph','PiscivFish','GelativFish','Myct'})
    set(gca,'XTickLabelRotation',45)
    ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('AE (%)')
    xline(5.5:5:60.5,'--')
    xlim([0, (width(AE)-4)*5+1])
    hold on

    %Secondary Production as stacked barplots meant to match diet proportion
    %figures
    %Assigning colors for each compartment column
    %Initially in rgb triplets
    colors=[255 255 255; %NO3-white
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
         5 16 108;%Amph-maroon
         102 0 51;%PiscivFish-dark pink
         255 0 102;%GelativFish-neon pink
         204 102 153;%Myct-light pink
         120 120 120;%Bac-light grey
         0 0 0;%SalpMort-black
         156 12 12;]; %HTL
    colors=colors/255; %Converts to rgb 0-1

    figure(f2)
    prod=prodmean(:,:,Cycle)';
    prod([19:23 end-1],:)=[];%removing additional non-production components
    prod(:,[19:23 end-1])=[];
    if Cycle==1
        subplot(1,4,1)
        b=bar(Cycle,prod(:,Salp_col+1),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Salp_col+1),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[3])
        set(gca,'XTickLabel',{'Salp'})
        set(gca,'FontSize',10)
        ylabel('Production (mmol N m^-^2 d^-^1)')

        subplot(1,4,2)
        b=bar(Cycle+6,prod(:,Myct_col+1),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Myct_col+1),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[9])
        set(gca,'XTickLabel',{'Myct'})
        set(gca,'FontSize',10)

        subplot(1,4,3)
        b=bar(Cycle+12,prod(:,PiscivFish_col+1),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,PiscivFish_col+1),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[15])
        set(gca,'XTickLabel',{'Pisciv'})
        set(gca,'FontSize',10)

        subplot(1,4,4)
        b=bar(Cycle+18,prod(:,GelativFish_col+1),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,GelativFish_col+1),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[21])
        set(gca,'XTickLabel',{'Gelativ'})
        set(gca,'FontSize',10)
    else
        subplot(1,4,1)
        b=bar(Cycle,prod(:,Salp_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Salp_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[3])
        set(gca,'XTickLabel',{'Salp'})
        set(gca,'FontSize',10)
        ylabel('Production (mmol N m^-^2 d^-^1)')

        subplot(1,4,2)
        b=bar(Cycle+6,prod(:,Myct_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Myct_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[9])
        set(gca,'XTickLabel',{'Myct'})
        set(gca,'FontSize',10)

        subplot(1,4,3)
        b=bar(Cycle+12,prod(:,PiscivFish_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,PiscivFish_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[15])
        set(gca,'XTickLabel',{'Pisciv'})
        set(gca,'FontSize',10)

        subplot(1,4,4)
        b=bar(Cycle+18,prod(:,GelativFish_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,GelativFish_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[21])
        set(gca,'XTickLabel',{'Gelativ'})
        set(gca,'FontSize',10)
    end

    %Protists
    figure(f3)
    subplot(1,3,1)
    b=bar(Cycle,prod(:,Flag_col),'stacked','FaceColor','flat')
    for k = 1:size(prod(:,Flag_col),1)
          b(k).CData = colors(k,:);
    end
    hold on
    set(gca,'XTick',[3])
    set(gca,'XTickLabel',{'MixFlag'})
    set(gca,'FontSize',10)
    ylabel('Production (mmol N m^-^2 d^-^1)')

    subplot(1,3,2)
    b=bar(Cycle+6,prod(:,HNF_col),'stacked','FaceColor','flat')
    for k = 1:size(prod(:,HNF_col),1)
          b(k).CData = colors(k,:);
    end
    hold on
    set(gca,'XTick',[9])
    set(gca,'XTickLabel',{'HNF'})
    set(gca,'FontSize',10)

    subplot(1,3,3)
    b=bar(Cycle+12,prod(:,Mic_col),'stacked','FaceColor','flat')
    for k = 1:size(prod(:,Mic_col),1)
          b(k).CData = colors(k,:);
    end
    hold on
    set(gca,'XTick',[15])
    set(gca,'XTickLabel',{'Micro'})
    set(gca,'FontSize',10)

    %Zoos
    figure(f4)
    if Cycle==1
        subplot(1,4,1)
        b=bar(Cycle,prod(:,Meso_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Meso_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[3])
        set(gca,'XTickLabel',{'Meso'})
        set(gca,'FontSize',10)
        ylabel('Production (mmol N m^-^2 d^-^1)')

        subplot(1,4,2)
        val=prod(:,Macro_col+1)+prod(:,VMMacro_col+1);
        b=bar(Cycle+6,val,'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Macro_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[9])
        set(gca,'XTickLabel',{'Macro'})
        set(gca,'FontSize',10)

        subplot(1,4,3)
        b=bar(Cycle+12,prod(:,Gel_col+1),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Gel_col+1),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[15])
        set(gca,'XTickLabel',{'GelPred'})
        set(gca,'FontSize',10)

        subplot(1,4,4)
        b=bar(Cycle+18,prod(:,Amph_col+1),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Amph_col+1),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[21])
        set(gca,'XTickLabel',{'Amph'})
        set(gca,'FontSize',10)
    else
        subplot(1,4,1)
        val=prod(:,Meso_col)+prod(:,VMMeso_col);
        b=bar(Cycle,val,'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Meso_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[3])
        set(gca,'XTickLabel',{'Meso'})
        set(gca,'FontSize',10)
        ylabel('Production (mmol N m^-^2 d^-^1)')

        subplot(1,4,2)
        val=prod(:,Macro_col)+prod(:,VMMacro_col);
        b=bar(Cycle+6,val,'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Macro_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[9])
        set(gca,'XTickLabel',{'Macro'})
        set(gca,'FontSize',10)

        subplot(1,4,3)
        b=bar(Cycle+12,prod(:,Gel_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Gel_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[15])
        set(gca,'XTickLabel',{'GelPred'})
        set(gca,'FontSize',10)

        subplot(1,4,4)
        b=bar(Cycle+18,prod(:,Amph_col),'stacked','FaceColor','flat')
        for k = 1:size(prod(:,Amph_col),1)
            b(k).CData = colors(k,:);
        end
        hold on
        set(gca,'XTick',[21])
        set(gca,'XTickLabel',{'Amph'})
        set(gca,'FontSize',10)
    end

    

    %Bargraph of relative salp outputs
    SalpOutProp=SalpOut;
    for j=1:height(SalpOutProp)
        total=sum(SalpOut(j,:));
        for k=1:width(SalpOut)
            temp=SalpOut(j,k);
            SalpOutProp(j,k)=temp/total*100;
        end
    end
    SalpOutPropMean(Cycle,:)=mean(SalpOutProp);
    SalpOutPropStd(Cycle,:)=std(SalpOutProp);
    figure(f5)
    b=bar(Cycle,mean(SalpOutProp),'stacked','FaceColor','flat')
    col=jet(width(SalpOutProp));
    for k=1:width(SalpOutProp)
        b(k).CData = col(k,:);
    end
    legend({'Amph','Myct','Gelativ','HTL','Mortality','NH4','DOM','Det'})
    set(gca,'XTick',1:5)
    set(gca,'XTickLabel',{'1','2','3','4','5'})
    set(gca,'FontSize',20)
    ylabel('Relative Output (%)')
    xlabel('Cycle')
    xlim([0, 6])
    ylim([0,100])
    hold on
    %Boxplot of salp outputs
    figure(f6)
    subplot(1,2,1)
    for i=1:5
        if Cycle==1 || Cycle==2 || Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        boxplot(SalpOut(:,i),'Positions',x,'Widths',0.8,'Colors',col)
        hold on
    end
    set(gca,'XTick',3:5:5*5+2)
    set(gca,'XTickLabel',{'Amph','Myct','Gelativ','HTL','Mort'})
    set(gca,'XTickLabelRotation',45)
    %ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('Flux (mmol N m^-^2 d^-^1)')
    xline([5.5 10.5 15.5 20.5],'--')
    xlim([0, 26])
    hold on

    subplot(1,2,2)
    for i=6:8
        if Cycle==1 || Cycle==2 || Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        boxplot(SalpOut(:,i),'Positions',x,'Widths',0.8,'Colors',col)
        hold on
    end
    set(gca,'XTick',23:5:43)
    set(gca,'XTickLabel',{'NH4','DOM','Det',})
    set(gca,'XTickLabelRotation',45)
    %ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('Flux (mmol N m^-^2 d^-^1)')
    xline(25.5:5:30.5,'--')
    xlim([20, 36])
    hold on

    %Bargraph of relative Macrozoo outputs
    ZooOutProp=ZooOut;
    for j=1:height(ZooOutProp)
        total=sum(ZooOut(j,:));
        for k=1:width(ZooOut)
            temp=ZooOut(j,k);
            ZooOutProp(j,k)=temp/total*100;
        end
    end
    ZooOutPropMean(Cycle,:)=mean(ZooOutProp);
    ZooOutPropStd(Cycle,:)=std(ZooOutProp);
    figure(f7)
    b=bar(Cycle,mean(ZooOutProp),'stacked','FaceColor','flat')
    col=jet(width(ZooOutProp));
    for k=1:width(ZooOutProp)
        b(k).CData = col(k,:);
    end
    legend({'Pisciv','Gel','Myct','HTL','NH4','DOM','Det'})
    set(gca,'XTick',1:5)
    set(gca,'XTickLabel',{'1','2','3','4','5'})
    xlim([0, 6])
    ylim([0,100])
    hold on
    %Boxplot of Zoo outputs
    figure(f8)
    subplot(1,2,1)
    for i=1:3
        if Cycle==1 || Cycle==2 || Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        boxplot(ZooOut(:,i),'Positions',x,'Widths',0.8,'Colors',col)
        hold on
    end
    set(gca,'XTick',3:5:3*5+2)
    set(gca,'XTickLabel',{'Pisciv','Gel','Myct'})
    set(gca,'XTickLabelRotation',45)
    %ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('Flux (mmol N m^-^2 d^-^1)')
    xline([5.5 10.5],'--')
    xlim([0, 16])
    hold on

    subplot(1,2,2)
    for i=4:6
        if Cycle==1 | Cycle==2 | Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        boxplot(ZooOut(:,i),'Positions',x,'Widths',0.8,'Colors',col)
        hold on
    end
    set(gca,'XTick',18:5:28)
    set(gca,'XTickLabel',{'NH4','DOM','Ldet'})
    set(gca,'XTickLabelRotation',45)
    %ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('Flux (mmol N m^-^2 d^-^1)')
    xline([20.5:5:35.5],'--')
    xlim([15, 31])
    hold on

    %Ecosystem Efficiency as Fish Production normalized to NPP
    ymax=10;
    ymin=10^-4;
    figure(f11)
    for i=1:width(EE)
        if Cycle==1 | Cycle==2 | Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        ViolinPlot(x,0.5,EE(:,i)*100,200,col)
    end
    x=Cycle+15;
    set(gca,'XTick',3:5:18)
    set(gca,'XTickLabel',{'Myct','PiscivFish','GelativFish','Terminal'})
    %set(gca,'XTickLabelRotation',45)
    set(gca, 'YScale', 'log')
    ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('Ecosystem Efficiency (%)')
    xline([5.5 10.5 15.5],'--')
    xlim([0, (length(Cols)+1)*5+1])
    hold on

    %Ecosystem Efficiency of terminal trophic groups
    figure(f12)
    subplot(1,4,1)
    col=[117 197 250; %Cycle 1 - salp - light blue
        3 154 255;  %Cycle 2 - salp - blue
        252 111 118;%Cycle 3 - nonsalp - light red
        3 79 130;  %Cycle 4 - salp - dark blue
        252 3 15];  %Cycle 5 - nonsalp - red
    col=col/255;
    ymax=3.5;
    ymin=0;
    ViolinPlot(Cycle,0.5,EE(:,4)*100,200,col(Cycle,:))
    set(gca,'XTick',1:5)
    set(gca,'XTickLabel',{'Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5'})
    set(gca,'XTickLabelRotation',45)
    set(gca, 'YScale', 'linear')
    ylim([ymin ymax]);
    set(gca,'FontSize',10)
    ylabel('Ecosystem Efficiency (%)')
    hold on

    subplot(1,4,2)
    ind=find(GGEavg(:,Cycle)>0);
    x=GGEavg(ind,Cycle)*100;
    y=EE(:,4)*100;
    scatter(mean(x),mean(y),'MarkerEdgeColor',col(Cycle,:),'MarkerFaceColor',col(Cycle,:))
    hold on
    handle_ellipse= plot_ellipse(x,y); %as an ellipse around 2sd (note loses asymmetry of uncertainty)
    handle_ellipse.Color=col(Cycle,:);
    xlabel('Mean GGE (%)')
    ylabel('Ecosystem Efficiency (%)')
    ylim([ymin ymax]);
    xlim([15 45])
    set(gca,'FontSize',10)
    set(gca,'box','on')
    hold on

    subplot(1,4,3)
    x=TL_EE(:,4);
    y=EE(:,4)*100;
    scatter(mean(x),mean(y),'MarkerEdgeColor',col(Cycle,:),'MarkerFaceColor',col(Cycle,:))
    hold on
    handle_ellipse= plot_ellipse(x,y);
    handle_ellipse.Color=col(Cycle,:);
    xlabel('Mean Terminal Food Chain Length')
    ylabel('Ecosystem Efficiency (%)')
    ylim([ymin ymax]);
    xlim([3 5.5])
    set(gca,'FontSize',10)
    set(gca,'box','on')
    hold on

    subplot(1,4,4)
    x=Prod(:,1)./sum(Prod(:,1:3),2)*100;
    y=EE(:,4)*100;
    p=scatter(mean(x),mean(y),'MarkerEdgeColor',col(Cycle,:),'MarkerFaceColor',col(Cycle,:));
    hold on
    handle_ellipse= plot_ellipse(x,y);
    handle_ellipse.Color=col(Cycle,:);
    xlabel('Pico Proportion (%)')
    ylabel('Ecosystem Efficiency (%)')
    ylim([ymin ymax]);
    xlim([10 100])
    set(gca,'FontSize',10)
    set(gca,'box','on')
    hold on

end

%%Exporting outputs into various Excel tables and graphics
writetable(array2table(ProdMean,'VariableNames',Columns),'Prod.xlsx','Sheet','Mean');
writetable(array2table(ProdStd,'VariableNames',Columns),'Prod.xlsx','Sheet','Std');
writetable(array2table(ProdLCI,'VariableNames',Columns),'Prod.xlsx','Sheet','LCI');
writetable(array2table(ProdUCI,'VariableNames',Columns),'Prod.xlsx','Sheet','UCI');
writetable(array2table(EEall(:,:,1),'VariableNames',{'Myct','Pisciv','Gelativ','Terminal'}),'EEall.xlsx','Sheet','Mean');
writetable(array2table(EEall(:,:,2),'VariableNames',{'Myct','Pisciv','Gelativ','Terminal'}),'EEall.xlsx','Sheet','LCI');
writetable(array2table(EEall(:,:,3),'VariableNames',{'Myct','Pisciv','Gelativ','Terminal'}),'EEall.xlsx','Sheet','UCI');
writetable(array2table(EEpreyall(:,:,1),'VariableNames',{'Myct','Pisciv','Gelativ','Terminal','Salp'}),'EEprey.xlsx','Sheet','Mean');
writetable(array2table(EEpreyall(:,:,2),'VariableNames',{'Myct','Pisciv','Gelativ','Terminal','Salp'}),'EEprey.xlsx','Sheet','LCI');
writetable(array2table(EEpreyall(:,:,3),'VariableNames',{'Myct','Pisciv','Gelativ','Terminal','Salp'}),'EEprey.xlsx','Sheet','UCI');
writematrix(EEall,'EE.xlsx');
writematrix(TLEEall,'TLEE.xlsx');
writematrix(GGEpathall,'GGEpath.xlsx');
writematrix(PicoPropall,'PicoProp.xlsx');
exportgraphics(f2,'ProdSF.png','Resolution',300)
exportgraphics(f3,'ProdProtist.png','Resolution',300)
exportgraphics(f4,'ProdMeso.png','Resolution',300)
exportgraphics(f5,'SalpOut.png','Resolution',300)
exportgraphics(f7,'MacroOut.png','Resolution',300)
exportgraphics(f9,'GGE.png','Resolution',300)
exportgraphics(f10,'AE.png','Resolution',300)
exportgraphics(f11,'FishEff.png','Resolution',300)
exportgraphics(f12,'TTE.png','Resolution',300)