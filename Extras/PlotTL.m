clear all
close all

%%Calculates Trophic Level for heterotrophic group and plots distribution
%%of outputs as violin plots

f1=figure('Position',[100 100 1000 500])
for Cycle=[1:5] 
    clearvars -except Cycle TLMean TLStd TLLCI TLUCI TL_all f1
    load(['N15NZInverseCycle',num2str(Cycle),'.mat']);
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat']);
    
    num=length(MCMCmat(:,1));
    MCMCmat(1:round(num*0.2),:)=[]; %Removing the first 20% of sims as burnin
    del15N(1:round(num*0.2),:)=[];
    len=length(MCMCmat(:,1));

    Columns={'NO3','NH4','Pico','Dtm','Flag','HNF','Mic','nvmMeso','VMMeso','nvmMacro','vmMacro','Gel','Salp','Amph','PiscFish','GelatFish','Myct',...
        'Bac','sdet','mdet','ldet','SalpDet','Dom','SalpMort'};
    
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

    %Calcing TL for every group in each MC sim
    for i=1:length(MCMCmat(:,1))
        flows=MCMCmat(i,:)./wts';
        TL0 = CalcTL(Cycle,flows,del15N,d15NInputs,Ae0,InputCol);
        TL_track(:,i)=TL0;
    end
    
    PlotCols = [HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,...
            PiscivFish_col,GelativFish_col,Myct_col];


    %Plotting TL
    ymax=6.3;
    ymin=1.8;
    %boxplot(TL_track(PlotCols,:)','Symbol','')
    for i=1:length(PlotCols)
        % if i==5 | i==6 | i==9 | i==10 %Pisciv and it's prey
        %     col = [0 0 1];
        % elseif i==7 | i==8 | i==11  %Gelativ and it's prey
        %     col = [1 0 0];
        % else
        %     col = [0.8 0.8 0];  %Myct and it's prey? Hard to do
        % end
        if Cycle==1 | Cycle==2 | Cycle==4
            col = [0 0 1];
        else
            col = [1 0 0];
        end
        x=Cycle+5*(i-1);
        try
            ViolinPlot(x,0.5,TL_track(PlotCols(i),:),200,col)
        end
    end
    set(gca,'XTick',3:5:length(PlotCols)*5+2)
    set(gca,'XTickLabel',{'Nanoflag','Microzoo','Meso (NVM)','Meso (VM)','Macro (NVM)','Macro (VM)','Gel Pred','Salps','Amph','PiscivFish','GelativFish','Myct'})
    set(gca,'XTickLabelRotation',45)
    ylim([ymin ymax]);
    set(gca,'FontSize',20)
    ylabel('Trophic Level')
    xline([5.5:5:55.5],'--')
    xlim([0, length(PlotCols)*5+1])
    hold on
 
    
    %Saving mean and CIs
    TLSort=sort(TL_track');
    if Cycle==1
        TLMean(Cycle,1:7)=mean(TLSort(:,1:7));
        TLStd(Cycle,1:7)=std(TLSort(:,1:7));
        TLLCI(Cycle,1:7)=TLSort(round(0.025*len),1:7);
        TLUCI(Cycle,1:7)=TLSort(round(0.975*len),1:7);
        TLMean(Cycle,8)=NaN;
        TLStd(Cycle,8)=NaN;
        TLLCI(Cycle,8)=NaN;
        TLUCI(Cycle,8)=NaN;
        TLMean(Cycle,9:24)=mean(TLSort(:,8:23));
        TLStd(Cycle,9:24)=std(TLSort(:,8:23));
        TLLCI(Cycle,9:24)=TLSort(round(0.025*len),8:23);
        TLUCI(Cycle,9:24)=TLSort(round(0.975*len),8:23);
    else
        TLMean(Cycle,:)=mean(TLSort(:,:));
        TLStd(Cycle,:)=std(TLSort(:,:));
        TLLCI(Cycle,:)=TLSort(round(0.025*len),:);
        TLUCI(Cycle,:)=TLSort(round(0.975*len),:);
    end
    save(['TL_',num2str(Cycle),'.mat'],'TL_track');


    %Text output of summary stats
    TLNano = sort(TL_track(HNF_col,:));
    TLMicro = sort(TL_track(Mic_col,:));
    if Cycle==1
        TLMeso = sort(TL_track(Meso_col,:));
    else
        TLMeso = sort(mean(TL_track([Meso_col,VMMeso_col],:)));
    end
    TLMacro = sort(mean(TL_track([Macro_col,VMMacro_col],:)));
    TLGel = sort(TL_track(Gel_col,:));
    TLSalp = sort(TL_track(Salp_col,:));
    TLAmph = sort(TL_track(Amph_col,:));
    TLGelat = sort(TL_track(GelativFish_col,:));
    TLMyct = sort(TL_track(Myct_col,:));
    TLPisc = sort(TL_track(PiscivFish_col,:));

    fn = 'TrophicLevels.txt'
    if Cycle==1
        fileID=fopen(fn,'w');
        tmp=datetime;
        formatSpec='Last run at %s \n \n';
        fprintf(fileID,formatSpec,tmp)
        fclose(fileID)
    end
    text=['In Cycle ',num2str(Cycle),', Nanoflagellates had a TL of ',num2str(mean(TLNano),2),...
        ' (95% C.I. = ',num2str(TLNano(round(0.025*len)),2),' - ',num2str(TLNano(round(0.975*len)),2),...
        '), while',char(10),'Protists had a TL of ',num2str(mean(TLMicro),2),'(',...
        num2str(TLMicro(round(0.025*len)),2),' - ',num2str(TLMicro(round(0.975*len)),2),').',...
        ' Mesozoos had a TL of ',num2str(mean(TLMeso),2),...
        ' (95% C.I. = ',num2str(TLMeso(round(0.025*len)),2),' - ',num2str(TLMeso(round(0.975*len)),2),...
        '), while',char(10),'Macrozoos had a TL of ',num2str(mean(TLMacro),2),'(',...
        num2str(TLMacro(round(0.025*len)),2),' - ',num2str(TLMyct(round(0.975*len)),2),').',...
        ' Gel Preds had a TL of ',num2str(mean(TLGel),2),...
        ' (95% C.I. = ',num2str(TLGel(round(0.025*len)),2),' - ',num2str(TLGel(round(0.975*len)),2),...
        '), while',char(10),'Salps had a TL of ',num2str(mean(TLSalp),2),'(',...
        num2str(TLSalp(round(0.025*len)),2),' - ',num2str(TLSalp(round(0.975*len)),2),')',...
        ' and Amphipods had a TL of ',num2str(mean(TLAmph),2),'(',num2str(TLAmph(round(0.025*len)),2),...
        ' - ',num2str(TLAmph(round(0.975*len)),2),').',char(10),...
        ' Gelativorous Fish had a TL of ',num2str(mean(TLGelat),2),...
        ' (95% C.I. = ',num2str(TLGelat(round(0.025*len)),2),' - ',num2str(TLGelat(round(0.975*len)),2),...
        '), while',char(10),'Myctophids had a TL of ',num2str(mean(TLMyct),2),'(',...
        num2str(TLMyct(round(0.025*len)),2),' - ',num2str(TLMyct(round(0.975*len)),2),')',...
        ' and Piscivorous Fish had a TL of ',num2str(mean(TLPisc),2),'(',num2str(TLPisc(round(0.025*len)),2),...
        ' - ',num2str(TLPisc(round(0.975*len)),2),').']
    writelines(text,fn,WriteMode="append")

end
%writetable(array2table(TLMean,'VariableNames',Columns),'TL.xlsx','Sheet','Mean');
%writetable(array2table(TLStd,'VariableNames',Columns),'TL.xlsx','Sheet','Std');
%writetable(array2table(TLLCI,'VariableNames',Columns),'TL.xlsx','Sheet','LCI');
%writetable(array2table(TLUCI,'VariableNames',Columns),'TL.xlsx','Sheet','UCI');
%exportgraphics(f1,'Trophic Level.png','Resolution',300)