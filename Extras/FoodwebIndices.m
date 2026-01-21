clear all
clearvars
close all

%This script calculates some ancillary food web metrics like nutrient
%support to phytos in the model and relative grazing by different
%heterotrophs

for Cycle=5
    clearvars -except Cycle

    load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
    clearvars -except Ae0 MCMCmat Aa Aa0 d15NInputs del15N InputCol wts Cycle
    num=length(MCMCmat(:,1));
    MCMCmat(1:round(num*0.2),:)=[];
    del15N(1:round(num*0.2),:)=[];
    len=length(MCMCmat(:,1));
    if Cycle==1
        [num,txt,raw]  = xlsread('N15InverseModelNZC1.xlsx','NZ','B6:B127');
    else
        [num,txt,raw]  = xlsread('N15InverseModelNZ.xlsx','NZ','B6:B145');
    end
    Columns = txt;
    meanvals = mean(MCMCmat)'./wts;
    stdvals = std(MCMCmat)'./wts;
    meantable = table(Columns,meanvals,stdvals);


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
    if Cycle==1
        ToMesozoo = sort([find(Ae0(Meso_col,:)==1),find(Ae0(Macro_col,:)==1),find(Ae0(VMMacro_col,:)==1),find(Ae0(Gel_col,:)==1),find(Ae0(Amph_col,:)==1)]);
        FromMesozoo = sort([find(Ae0(Meso_col,:)==-1),find(Ae0(Macro_col,:)==-1),find(Ae0(VMMacro_col,:)==-1),find(Ae0(Gel_col,:)==-1),find(Ae0(Amph_col,:)==-1)]);
    else
        ToMesozoo = sort([find(Ae0(Meso_col,:)==1),find(Ae0(VMMeso_col,:)==1),find(Ae0(Macro_col,:)==1),find(Ae0(VMMacro_col,:)==1),find(Ae0(Gel_col,:)==1),find(Ae0(Amph_col,:)==1)]);
        FromMesozoo = sort([find(Ae0(Meso_col,:)==-1),find(Ae0(VMMeso_col,:)==-1),find(Ae0(Macro_col,:)==-1),find(Ae0(VMMacro_col,:)==-1),find(Ae0(Gel_col,:)==-1),find(Ae0(Amph_col,:)==-1)]);
    end
    ToSalp = find(Ae0(Salp_col,:)==1);
    FromSalp = find(Ae0(Salp_col,:)==-1);
    FromHTL = sort([find(Ae0(PiscivFish_col,:)==-1),find(Ae0(GelativFish_col,:)==-1),find(Ae0(Myct_col,:)==-1)]);
    FromPhy = sort([find(Ae0(Pico_col,:)==-1),find(Ae0(Dtm_col,:)==-1)]);
    FromMix = find(Ae0(Flag_col,:)==-1);
    ToMix = find(Ae0(Flag_col,:)==1);
    FromPico = find(Ae0(Pico_col,:)==-1);
    ToPico = find(Ae0(Pico_col,:)==1);
    FromNO3 = find(Ae0(NO3_col,:)==-1);
    FromNH4 = find(Ae0(NH4_col,:)==-1);
    ToNH4 = find(Ae0(NH4_col,:)==1);
    ToDon = find(Ae0(Dom_col,:)==1);
    FromBac = find(Ae0(bac_col,:)==-1);
    FromProtist = sort([find(Ae0(HNF_col,:)==-1),find(Ae0(Mic_col,:)==-1)]);
    ToProtist = sort([find(Ae0(HNF_col,:)==1),find(Ae0(Mic_col,:)==1)]);
    FromDtm = find(Ae0(Dtm_col,:)==-1);
    ToDtm = find(Ae0(Dtm_col,:)==1);
    ToSdet = find(Ae0(sdet_col,:)==1);
    ToMdet = find(Ae0(mdet_col,:)==1);

    PHYtoMESOZOO = intersect(ToMesozoo,FromPhy);
    PHYtoSalp = intersect(ToSalp,FromPhy);
    FLAGtoMESOZOO = intersect(ToMesozoo,FromMix);
    FLAGtoSalp = intersect(ToSalp,FromMix);
    BACtoFLAG = intersect(FromBac,ToMix);
    PICOtoFLAG = intersect(FromPico,ToMix);
    FLAGtoDON = intersect(FromMix,ToDon);
    PICOtoDON = intersect(FromPico,ToDon);
    DTMtoDON = intersect(FromDtm,ToDon);
    NUTtoFLAG = intersect([FromNO3,FromNH4],ToMix);
    NUTtoDTM = intersect([FromNO3,FromNH4],ToDtm);
    NUTtoPICO = intersect([FromNO3,FromNH4],ToPico);
    PICOtoDET = intersect(FromPico,ToSdet);
    FLAGtoDET = intersect(FromMix,ToSdet);
    DTMtoDET = intersect(FromDtm,ToMdet);
    PROTISTtoMESOZOO = intersect(ToMesozoo,FromProtist);
    BACtoNH4 = intersect(FromBac,ToNH4);
    BACtoPROTIST = intersect(FromBac,ToProtist);
    PROTISTtoDON = intersect(FromProtist,ToDon);
    PROTISTtoNH4 = intersect(FromProtist,ToNH4);
    NO3toPHY = intersect(FromNO3,[ToDtm,ToPico,ToMix]);
    NO3toPico = intersect(FromNO3,ToPico);
    NO3toDtm = intersect(FromNO3,ToDtm);
    NO3toFlag = intersect(FromNO3,ToMix);
    NH4toPHY = intersect(FromNH4,[ToDtm,ToPico,ToMix]);
    NH4toPico = intersect(FromNH4,ToPico);
    NH4toDtm = intersect(FromNH4,ToDtm);
    NH4toFlag = intersect(FromNH4,ToMix);
    PHYtoPROTIST = intersect([FromPhy,FromMix],ToProtist);
    PICOtoPROTIST = intersect(FromPico,ToProtist);
    DTMtoPROTIST = intersect(FromDtm,ToProtist);
    FLAGtoPROTIST = intersect(FromMix,ToProtist);
    PHYtoDON = intersect(FromPhy,ToDon);
    MIXtoDON = intersect(FromMix,ToDon);
    MESOZOOtoDON = intersect(FromMesozoo,ToDon);
    HTLtoDON = intersect(FromHTL,ToDon);
    BACtoNH4 = intersect(FromBac,ToNH4);
    MIXtoNH4 = intersect(FromMix,ToNH4);
    MESOZOOtoNH4 = intersect(FromMesozoo,ToNH4);
    HTLtoNH4 = intersect(FromHTL,ToNH4);



    for i=1:len
        flows=MCMCmat(i,:)./wts';

        PICOnpp(i) = sum(flows(NUTtoPICO)) - flows(PICOtoDON);
        DTMnpp(i) = sum(flows(NUTtoDTM)) - flows(DTMtoDON);
        FLAGnpp(i) = sum(flows(NUTtoFLAG)) - flows(FLAGtoDON);

        FLAGmix = flows(BACtoFLAG)+flows(PICOtoFLAG);
        %Per Stukel et al 2012 "the sum of direct nitrogen flux from phytoplankton to metazoan zooplankton"
        MesoHerbFoodChain(i) = sum(flows(PHYtoMESOZOO)) + sum(flows(FLAGtoMESOZOO))*FLAGnpp(i)/(FLAGnpp(i)+FLAGmix);
        SalpHerbFoodChain(i) = sum(flows(PHYtoSalp)) + sum(flows(FLAGtoSalp))*FLAGnpp(i)/(FLAGnpp(i)+FLAGmix);
        HerbFoodChain(i) = sum(MesoHerbFoodChain) + sum(SalpHerbFoodChain);
        %the sum of nitrogen flux that reaches metazoan zooplankton after passing through protistan grazers"
        MultivFoodChain(i) = sum(flows(PROTISTtoMESOZOO)) + sum(flows(FLAGtoMESOZOO))*sum(flows(PICOtoFLAG))/(FLAGnpp(i)+FLAGmix);

        ProtistRespirationPercent = sum(flows(PROTISTtoNH4))/sum(flows(ToProtist));
        %the sum of bacterial respiration and the fraction of protistan respiration that was supported by bacterial production
        MicrobialLoop(i) = flows(BACtoNH4) + sum(flows(BACtoPROTIST))*ProtistRespirationPercent;

        NPP(i) = PICOnpp(i) + DTMnpp(i) + FLAGnpp(i);

        NO3UPTAKE(i) = sum(flows(NO3toPHY));
        NO3UPTAKEPico(i) = sum(flows(NO3toPico));
        NO3UPTAKEDtm(i) = sum(flows(NO3toDtm));
        NO3UPTAKEFlag(i) = sum(flows(NO3toFlag));
        NH4UPTAKE(i) = sum(flows(NH4toPHY));
        NH4UPTAKEPico(i) = sum(flows(NH4toPico));
        NH4UPTAKEDtm(i) = sum(flows(NH4toDtm));
        NH4UPTAKEFlag(i) = sum(flows(NH4toFlag));

        PicoMort(i) = sum(flows(PICOtoDET));
        DtmMort(i) = sum(flows(DTMtoDET));
        FlagMort(i) = sum(flows(FLAGtoDET));
        PhytoMort(i) = PicoMort(i) + DtmMort(i) + FlagMort(i);

        PROTISTGRAZING(i) = sum(flows(PHYtoPROTIST)) + sum(flows(PICOtoFLAG));
        PROTISTGRAZINGPico(i) = sum(flows(PICOtoPROTIST)) + sum(flows(PICOtoFLAG));
        PROTISTGRAZINGDtm(i) = sum(flows(DTMtoPROTIST));
        PROTISTGRAZINGFlag(i) = sum(flows(FLAGtoPROTIST));
        MESOZOOGRAZING(i) = sum(flows(PHYtoMESOZOO)) + sum(flows(FLAGtoMESOZOO));

        PhytoplanktonDONproduction(i)=sum(flows(PHYtoDON)) + sum(flows(MIXtoDON))*FLAGnpp(i)/(FLAGnpp(i)+FLAGmix);
        ProtistDONproduction(i)=sum(flows(PROTISTtoDON)) + sum(flows(MIXtoDON))*FLAGmix/(FLAGnpp(i)+FLAGmix);
        MesozooDONproduction(i)=sum(flows(MESOZOOtoDON));
        HTLDONproduction(i)=sum(flows(HTLtoDON));
        TotalDONproduction(i)=sum(flows(ToDon));

        BacNH4production(i)=sum(flows(BACtoNH4));
        ProtistNH4production(i)=sum(flows(PROTISTtoNH4)) + sum(flows(MIXtoNH4))*FLAGmix/(FLAGnpp(i)+FLAGmix);
        MesozooNH4production(i)=sum(flows(MESOZOOtoNH4));
        HTLNH4production(i)=sum(flows(HTLtoNH4));
        TotalNH4production(i)=sum(flows(ToNH4));
    end

    % col=[0 0.7 0;
    %     0 0 1;
    %     1 0.5 0;];
    % figure('Position',[50 50 600 400])
    % hold on
    % ViolinPlot(1,0.5,HerbFoodChain,200,col(1,:))
    % ViolinPlot(2,0.5,MultivFoodChain,200,col(2,:))
    % ViolinPlot(3,0.5,MicrobialLoop,200,col(3,:))
    % plot([4,4],[0,3],'-k')
    % fill([0.1, 0.1, 1.9, 1.9, 0.1],[15 11 11 15 15],'w')
    % plot(0.2,14,'dk','MarkerFaceColor',col(1,:),'MarkerSize',12)
    % text(0.3,14,'Herbiv Food Chain','FontSize',16)
    % plot(0.2,13,'dk','MarkerFaceColor',col(2,:),'MarkerSize',12)
    % text(0.3,13,'Multiv Food Chain','FontSize',16)
    % plot(0.2,12,'dk','MarkerFaceColor',col(3,:),'MarkerSize',12)
    % text(0.3,12,'Microbial Loop','FontSize',16);
    % set(gca,'box','on')
    % set(gca,'XTick',[2,6])
    % set(gca,'XTickLabel',{'Shallow Eup Zone'})
    % %set(gca,'XTickLabelRotation',45)
    % ylabel('Flux (mmol N m^-^2 d^-^1)')
    % set(gca,'FontSize',16)
    % title(['Cycle ',num2str(Cycle)])

    % figure('Position',[50 50 600 400])
    % hold on
    % ViolinPlot(1,0.5,HerbFoodChain./NPP,200,col(1,:))
    % ViolinPlot(2,0.5,MultivFoodChain./NPP,200,col(2,:))
    % ViolinPlot(3,0.5,MicrobialLoop./NPP,200,col(3,:))
    % plot([4,4],[0,3],'-k')
    % set(gca,'box','on')
    % set(gca,'XTick',[2,6])
    % set(gca,'XTickLabel',{'Shallow Eup Zone'})
    % %set(gca,'XTickLabelRotation',45)
    % ylabel('Flux / NPP')
    % set(gca,'FontSize',16)
    % ylim([0 1.5])

    NH4UPTAKEpercent = sort(NH4UPTAKE./(NH4UPTAKE+NO3UPTAKE))*100;
    NO3UPTAKEpercent = sort(NO3UPTAKE./(NH4UPTAKE+NO3UPTAKE))*100;
    NH4UPTAKEpercentPico = sort(NH4UPTAKEPico./(NH4UPTAKEPico+NO3UPTAKEPico))*100;
    NO3UPTAKEpercentPico = sort(NO3UPTAKEPico./(NH4UPTAKEPico+NO3UPTAKEPico))*100;
    NH4UPTAKEpercentDtm = sort(NH4UPTAKEDtm./(NH4UPTAKEDtm+NO3UPTAKEDtm))*100;
    NO3UPTAKEpercentDtm = sort(NO3UPTAKEDtm./(NH4UPTAKEDtm+NO3UPTAKEDtm))*100;
    NH4UPTAKEpercentFlag = sort(NH4UPTAKEFlag./(NH4UPTAKEFlag+NO3UPTAKEFlag))*100;
    NO3UPTAKEpercentFlag = sort(NO3UPTAKEFlag./(NH4UPTAKEFlag+NO3UPTAKEFlag))*100;
    ['NH4 uptake by all phytoplankton was ','(mean = ',num2str(mean(NH4UPTAKEpercent),2),...
        '%; 95% C.I. = ',num2str(NH4UPTAKEpercent(round(0.025*len)),2),' - ',num2str(NH4UPTAKEpercent(round(0.975*len)),2),...
        '%)',' whereas NO3 uptake was (',num2str(mean(NO3UPTAKEpercent),2),...
        '%; ',num2str(NO3UPTAKEpercent(round(0.025*len)),2),' - ',num2str(NO3UPTAKEpercent(round(0.975*len)),2),...
        '%).']

    ['NH4 uptake by Pico was ','(mean = ',num2str(mean(NH4UPTAKEpercentPico),2),...
        '%; 95% C.I. = ',num2str(NH4UPTAKEpercentPico(round(0.025*len)),2),' - ',num2str(NH4UPTAKEpercentPico(round(0.975*len)),2),...
        '%)',' whereas NO3 uptake was (',num2str(mean(NO3UPTAKEpercentPico),2),'%; ',num2str(NO3UPTAKEpercentPico(round(0.025*len)),2),' - ',...
        num2str(NO3UPTAKEpercentPico(round(0.975*len)),2),'%).']

    ['NH4 uptake by Dtm was ','(mean = ',num2str(mean(NH4UPTAKEpercentDtm),2),...
        '%; 95% C.I. = ',num2str(NH4UPTAKEpercentDtm(round(0.025*len)),2),' - ',num2str(NH4UPTAKEpercentDtm(round(0.975*len)),2),...
        '%)',' whereas NO3 uptake was (',num2str(mean(NO3UPTAKEpercentDtm),2),'%; ',num2str(NO3UPTAKEpercentDtm(round(0.025*len)),2),' - ',...
        num2str(NO3UPTAKEpercentDtm(round(0.975*len)),2),'%).']

    ['NH4 uptake by Flag was ','(mean = ',num2str(mean(NH4UPTAKEpercentFlag),2),...
        '%; 95% C.I. = ',num2str(NH4UPTAKEpercentFlag(round(0.025*len)),2),' - ',num2str(NH4UPTAKEpercentFlag(round(0.975*len)),2),...
        '%)',' whereas NO3 uptake was (',num2str(mean(NO3UPTAKEpercentFlag),2),'%; ',num2str(NO3UPTAKEpercentFlag(round(0.025*len)),2),' - ',...
        num2str(NO3UPTAKEpercentFlag(round(0.975*len)),2),'%).']


    NPPsort = sort(NPP);
    ['Mean NPP was (',num2str(mean(NPPsort),2),'; ',num2str(NPPsort(round(0.025*len)),2),...
        ' - ',num2str(NPPsort(round(0.975*len)),2),')']

    DTMNPPpercent = sort(DTMnpp./NPP)*100;
    FLAGNPPpercent = sort(FLAGnpp./NPP)*100;
    PICONPPpercent = sort(PICOnpp./NPP)*100;
    ['Picophytoplankton NPP was ',num2str(mean(PICONPPpercent),2),'% (',num2str(PICONPPpercent(0.025*len),2),' - ',num2str(PICONPPpercent(0.975*len),2),'%)']
    ['Flagellate NPP was ',num2str(mean(FLAGNPPpercent),2),'% (',num2str(FLAGNPPpercent(0.025*len),2),' - ',num2str(FLAGNPPpercent(0.975*len),2),'%)']
    ['Diatom NPP was ',num2str(mean(DTMNPPpercent),2),'% (',num2str(DTMNPPpercent(0.025*len),2),' - ',num2str(DTMNPPpercent(0.975*len),2),'%)']

    ProtistGrazingPercent = sort(PROTISTGRAZING./NPP)*100;
    ProtistPicoGrazingPercent = sort(PROTISTGRAZINGPico./PROTISTGRAZING)*100;
    ProtistDtmGrazingPercent = sort(PROTISTGRAZINGDtm./PROTISTGRAZING)*100;
    ProtistFlagGrazingPercent = sort(PROTISTGRAZINGFlag./PROTISTGRAZING)*100;
    ['Protists (including mixotrophic flagellates) consumed ',num2str(mean(ProtistGrazingPercent),2),'% (',num2str(ProtistGrazingPercent(round(0.025*len)),2),' - ',num2str(ProtistGrazingPercent(round(0.975*len)),2),'%) of phytoplankton production.']
    ['Of which ',num2str(mean(ProtistPicoGrazingPercent),2),'% (',num2str(ProtistPicoGrazingPercent(round(0.025*len)),2),' - ',num2str(ProtistPicoGrazingPercent(round(0.975*len)),2),'%) was on pico,']
    [num2str(mean(ProtistDtmGrazingPercent),2),'% (',num2str(ProtistDtmGrazingPercent(round(0.025*len)),2),' - ',num2str(ProtistDtmGrazingPercent(round(0.975*len)),2),'%) was on diatoms,']
    ['and ',num2str(mean(ProtistFlagGrazingPercent),2),'% (',num2str(ProtistFlagGrazingPercent(round(0.025*len)),2),' - ',num2str(ProtistFlagGrazingPercent(round(0.975*len)),2),'%) was on flagellates.']

    MesozooGrazingPercent = sort((MesoHerbFoodChain)./(NPP))*100;
    ['Metazoan zooplankton (excluding salps) consumed ',num2str(mean(MesozooGrazingPercent),2),'% (',num2str(MesozooGrazingPercent(round(0.025*len)),2),' - ',num2str(MesozooGrazingPercent(round(0.975*len)),2),'%) of phytoplankton production.']
    SalpGrazingPercent = sort((SalpHerbFoodChain)./(NPP))*100;
    ['Salps consumed ',num2str(mean(SalpGrazingPercent),2),'% (',num2str(SalpGrazingPercent(round(0.025*len)),2),' - ',num2str(SalpGrazingPercent(round(0.975*len)),2),'%) of phytoplankton production.']

    PhytoMortPercent = sort(PhytoMort./(NPP))*100;
    DTMMortPercent = sort(DtmMort./NPP)*100;
    FLAGMortPercent = sort(FlagMort./NPP)*100;
    PICOMortPercent = sort(PicoMort./NPP)*100;
    ['Natural mortality of phytos was ',num2str(mean(PhytoMortPercent),2),'% (',num2str(PhytoMortPercent(round(0.025*len)),2),' - ',num2str(PhytoMortPercent(round(0.975*len)),2),'%) of phytoplankton production.']
    ['Natural mortality of picos was ',num2str(mean(PICOMortPercent),2),'% (',num2str(PICOMortPercent(round(0.025*len)),2),' - ',num2str(PICOMortPercent(round(0.975*len)),2),'%) of total production.']
    ['Natural mortality of diatoms was ',num2str(mean(DTMMortPercent),2),'% (',num2str(DTMMortPercent(round(0.025*len)),2),' - ',num2str(DTMMortPercent(round(0.975*len)),2),'%) of total production.']
    ['Natural mortality of flags was ',num2str(mean(FLAGMortPercent),2),'% (',num2str(FLAGMortPercent(round(0.025*len)),2),' - ',num2str(FLAGMortPercent(round(0.975*len)),2),'%) of total production.']

    HerbFoodChainPercent=sort(HerbFoodChain./NPP)*100;
    MultivFoodChainPercent=sort(MultivFoodChain./NPP)*100;
    MicrobialLoopPercent=sort(MicrobialLoop./NPP)*100;
end