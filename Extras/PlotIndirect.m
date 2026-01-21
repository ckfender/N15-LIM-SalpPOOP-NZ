clear all
close all

%%Calculations and visualizations of indirect support per Hannon 1973 "The
%%Structure of Ecosystems." Actual indirect support didn't make it into the
%%paper but some of the matrices produced by this script are used elsewhere
%%as it is a useful way to conceptualize production matrices

%Preallocating some result matrices
prodmean = zeros(26,26,5);
prodLCI = zeros(26,26,5);
prodUCI = zeros(26,26,5);
SuppMean = zeros(25,25,5);
SuppLCI = zeros(25,25,5);
SuppUCI = zeros(25,25,5);
NormSuppMean = zeros(25,25,5);
NormSuppLCI = zeros(25,25,5);
NormSuppUCI = zeros(25,25,5);
%f6=figure('Position',[50 50 1000 500]);
for Cycle = 1:5
    load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
    clearvars -except Ae0 MCMCmat Aa Aa0 d15NInputs del15N InputCol wts Cycle f6...
        prodmean prodLCI prodUCI SuppMean SuppLCI SuppUCI NormSuppMean NormSuppLCI ...
        NormSuppUCI GGEavg GGESalp GGENonSalp
    %MCMCmat=MCMCmat(1:10000,:);
    num=length(MCMCmat(:,1));
    MCMCmat(1:round(num*0.2),:)=[]; %Removing first 20% of sims as burnin
    del15N(1:round(num*0.2),:)=[];
    MCMCmat = MCMCmat(1:10:end,:); %If necessary, removing all but every 10th result
    del15N = del15N(1:10:end,:);
    len=length(MCMCmat(:,1));
    

    % f1=figure('Position',[50 50 300 600]);
    % f2=figure('Position',[150 50 300 600]);
    % f3=figure('Position',[250 50 300 600]);
    % f4=figure('Position',[450 50 500 600]);

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

    for i=1:len
        %['Cycle ' num2str(Cycle) ' simulation ' num2str(i) ' out of ' num2str(len)]
        flows=MCMCmat(i,:)./wts';

        [P,struct,ext,r,e,p,Pspec,structspec] = IndirectFlowsRespiration(Cycle,Ae0,flows);

        %Per Hannon 1973, multiplying struct (the inverse structure matrix) by
        %the respiration matrix produces the special energy matrix (Table 6)
        %where each entry is the total direct (through consumption of production)
        %and indirect energy consumed by the
        %row components to allow the respiration of the column component.
        %%Relative to the previous config these are read as "row element
        %%production respired by column element" normalized to total
        %%production of the row element. So what was Salp2Gelativ is now
        %%found in NormSupport(Salp_col,GelativFish_col)
        Support(:,:,i)=struct.*r';
        SuppSpec(:,:) = structspec.*[r;sum(Pspec(:,end))]';
        NormSupport(:,:,i)=Support(:,:,i)./p';
        P_all(:,:,i)=Pspec;


        %Calculating GGE from the P matrix (the output production i.e. row sum
        %divided by the input ingestion i.e. column sum)
        temp=P;
        temp(:,end+1)=NaN;
        temp(end+1,:)=NaN;
        temp(:,end)=Pspec(:,end);
        temp(end,:)=Pspec(end,:);
        GGE = sum(Pspec(:,[1:(sdet_col-1) SalpMort_col end]),2)./sum(temp(Pico_col:end,:),1)';
        GGE([1:Flag_col sdet_col:end-1]) = 0; %GGE for anything but heterotrophs is meaningless
        %Estimating the mean GGE of the trophic pathway leading to Myct,
        %Pisciv, Gelativ, and HTL as the mean of GGEs of the pathway's
        %constituents weighted by their direct and indirect support of the
        %terminal group
        M_path = sum(GGE.*SuppSpec(:,Myct_col))/sum(SuppSpec(HNF_col:bac_col,Myct_col));
        P_path = sum(GGE.*SuppSpec(:,PiscivFish_col))/sum(SuppSpec(HNF_col:bac_col,PiscivFish_col));
        G_path = sum(GGE.*SuppSpec(:,GelativFish_col))/sum(SuppSpec(HNF_col:bac_col,GelativFish_col));
        H_path = sum(GGE.*SuppSpec(:,end))/sum(SuppSpec(HNF_col:bac_col,end));
        Mprod = sum(Pspec(Myct_col,Myct_col));
        Pprod = sum(Pspec(PiscivFish_col,PiscivFish_col));
        Gprod = sum(Pspec(GelativFish_col,GelativFish_col));
        Hprod = sum(Pspec(:,end));
        %The mean GGE of the ecosystem as the mean of the above 4
        %non-steady state pathways through which energy may exit, weighted
        %by the contribution of each pathway to outgoing production.
        GGEavg(i,Cycle) = (M_path*Mprod + P_path*Pprod + G_path*Gprod + H_path*Hprod)/sum([Mprod Pprod Gprod Hprod]);

        %Mean GGEs of the salp and non-salp pathways for GelativFish
        if Cycle == 1
            path = [HNF_col Mic_col Salp_col+1];
        else
            path = [HNF_col Mic_col Salp_col];
        end
        GGE_pathSalp(i) = sum(GGE(path).*SuppSpec(path,GelativFish_col))/sum(SuppSpec(path,GelativFish_col));
        if Cycle == 1
            path = [HNF_col Mic_col Meso_col Macro_col+1 VMMacro_col+1 Gel_col+1];
        else
            path = [HNF_col Mic_col Meso_col VMMeso_col Macro_col VMMacro_col Gel_col];
        end
        GGE_pathNonSalp(i) = sum(GGE(path).*SuppSpec(path,GelativFish_col))/sum(SuppSpec(path,GelativFish_col));
    end

    %%Saving output summary metrics
    %Absolute support
    if Cycle==1 %inserting row and column for empty VMMeso compartment
        temp = mean(Support,3);
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        SuppMean(:,:,Cycle) = temp2;

        suppsort=sort(Support,3);
        temp = suppsort(:,:,round(0.025*len));
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        SuppLCI(:,:,Cycle) = temp2;
        temp = suppsort(:,:,round(0.975*len));
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        SuppUCI(:,:,Cycle) = temp2;
    else
        suppsort = sort(Support,3);
        SuppMean(:,:,Cycle) = mean(Support,3);
        SuppLCI(:,:,Cycle)=suppsort(:,:,round(0.025*len));
        SuppUCI(:,:,Cycle)=suppsort(:,:,round(0.975*len));
    end
    %Normalized support
    if Cycle==1 %inserting row and column for empty VMMeso compartment
        temp = mean(NormSupport,3);
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        NormSuppMean(:,:,Cycle) = temp2;

        normsuppsort=sort(NormSupport,3);
        temp = normsuppsort(:,:,round(0.025*len));
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        NormSuppLCI(:,:,Cycle) = temp2;
        temp = normsuppsort(:,:,round(0.975*len));
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        NormSuppUCI(:,:,Cycle) = temp2;
    else
        normsuppsort = sort(NormSupport,3);
        NormSuppMean(:,:,Cycle) = mean(NormSupport,3);
        NormSuppLCI(:,:,Cycle)=normsuppsort(:,:,round(0.025*len));
        NormSuppUCI(:,:,Cycle)=normsuppsort(:,:,round(0.975*len));
    end

    %p matrix for the given cycle into an array for use in other scripts
    if Cycle==1 %inserting row and column for empty VMMeso compartment
        temp = mean(P_all,3);
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        prodmean(:,:,Cycle) = temp2;

        prodsort=sort(P_all,3);
        temp = prodsort(:,:,round(0.025*len));
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        prodLCI(:,:,Cycle) = temp2;
        temp = prodsort(:,:,round(0.975*len));
        newrow = zeros(1,size(temp,2));
        temp2 = [temp(1:8,:); newrow; temp(9:end,:)];
        newcol = zeros(size(temp2,1),1);
        temp2 = [temp2(:,1:8) newcol temp2(:,9:end)];
        prodUCI(:,:,Cycle) = temp2;
    else
        prodsort = sort(P_all,3);
        prodmean(:,:,Cycle) = mean(P_all,3);
        prodLCI(:,:,Cycle)=prodsort(:,:,round(0.025*len));
        prodUCI(:,:,Cycle)=prodsort(:,:,round(0.975*len));
    end

    %GelativFish pathway mean GGEs
    GGEsort = sort(GGE_pathSalp);
    GGESalp(1,Cycle) = mean(GGEsort);
    GGESalp(2,Cycle) = GGEsort(round(0.025*len));
    GGESalp(3,Cycle) = GGEsort(round(0.975*len));
    GGEsort = sort(GGE_pathNonSalp);
    GGENonSalp(1,Cycle) = mean(GGEsort);
    GGENonSalp(2,Cycle) = GGEsort(round(0.025*len));
    GGENonSalp(3,Cycle) = GGEsort(round(0.975*len));

    %
    % 
    % %Sub Figure 1 - Absolute support of fish from phytos
    % figure(f1)
    % hold on
    % ViolinPlot(1,0.5,Support(Pico_col,Myct_col,:),200,[0.9 0.9 0])
    % ViolinPlot(2,0.5,Support(Dtm_col,Myct_col,:),200,[0.2 1 0])
    % ViolinPlot(3,0.5,Support(Flag_col,Myct_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(6,0.5,Support(Pico_col,PiscivFish_col,:),200,[0.9 0.9 0])
    % ViolinPlot(7,0.5,Support(Dtm_col,PiscivFish_col,:),200,[0.2 1 0])
    % ViolinPlot(8,0.5,Support(Flag_col,PiscivFish_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(11,0.5,Support(Pico_col,GelativFish_col,:),200,[0.9 0.9 0])
    % ViolinPlot(12,0.5,Support(Dtm_col,GelativFish_col,:),200,[0.2 1 0])
    % ViolinPlot(13,0.5,Support(Flag_col,GelativFish_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(16,0.5,Support(Pico_col,Myct_col,:)+Support(Pico_col,PiscivFish_col,:)+Support(Pico_col,GelativFish_col,:),200,[0.9 0.9 0])
    % ViolinPlot(17,0.5,Support(Dtm_col,Myct_col,:)+Support(Dtm_col,PiscivFish_col,:)+Support(Dtm_col,GelativFish_col,:),200,[0.2 1 0])
    % ViolinPlot(18,0.5,Support(Flag_col,Myct_col,:)+Support(Flag_col,PiscivFish_col,:)+Support(Flag_col,GelativFish_col,:),200,[0 0.6 0])
    % 
    % set(gca,'XTick',[2.5 7.5 12.5 17.5])
    % set(gca,'XTickLabel',{'Myct','Pisciv','Gelativ','All'})
    % set(gca,'box','on')
    % ylabel(['Indirect support of fish',char(10),'(mmol N m^-^2 d^-^1)'])
    % set(gca,'XTickLabelRotation',45)
    % set(gca,'FontSize',10)
    % 
    % %Sub Figure 2 - Normalized indirect support of fish from phytos
    % figure(f2)
    % hold on
    % ViolinPlot(1,0.5,NormSupport(Pico_col,Myct_col,:),200,[0.9 0.9 0])
    % ViolinPlot(2,0.5,NormSupport(Dtm_col,Myct_col,:),200,[0.2 1 0])
    % ViolinPlot(3,0.5,NormSupport(Flag_col,Myct_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(6,0.5,NormSupport(Pico_col,PiscivFish_col,:),200,[0.9 0.9 0])
    % ViolinPlot(7,0.5,NormSupport(Dtm_col,PiscivFish_col,:),200,[0.2 1 0])
    % ViolinPlot(8,0.5,NormSupport(Flag_col,PiscivFish_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(11,0.5,NormSupport(Pico_col,GelativFish_col,:),200,[0.9 0.9 0])
    % ViolinPlot(12,0.5,NormSupport(Dtm_col,GelativFish_col,:),200,[0.2 1 0])
    % ViolinPlot(13,0.5,NormSupport(Flag_col,GelativFish_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(16,0.5,NormSupport(Pico_col,Myct_col,:)+NormSupport(Pico_col,PiscivFish_col,:)+NormSupport(Pico_col,GelativFish_col,:),200,[0.9 0.9 0])
    % ViolinPlot(17,0.5,NormSupport(Dtm_col,Myct_col,:)+NormSupport(Dtm_col,PiscivFish_col,:)+NormSupport(Dtm_col,GelativFish_col,:),200,[0.2 1 0])
    % ViolinPlot(18,0.5,NormSupport(Flag_col,Myct_col,:)+NormSupport(Flag_col,PiscivFish_col,:)+NormSupport(Flag_col,GelativFish_col,:),200,[0 0.6 0])
    % 
    % set(gca,'XTick',[2.5 8.5 12.5 17.5])
    % set(gca,'XTickLabel',{'Myct','Pisciv','Gelativ','All'})
    % set(gca,'box','on')
    % ylabel(['Indirect support of fish',char(10),'(Normalized to phytoplankton production)'])
    % set(gca,'XTickLabelRotation',45)
    % set(gca,'FontSize',10)
    % 
    % 
    % %Sub Figure 3 - Indirect support of protists
    % figure(f3)
    % hold on
    % ViolinPlot(1,0.5,NormSupport(Pico_col,HNF_col,:),200,[0.9 0.9 0])
    % ViolinPlot(2,0.5,NormSupport(Dtm_col,HNF_col,:),200,[0.2 1 0])
    % ViolinPlot(3,0.5,NormSupport(Flag_col,HNF_col,:),200,[0 0.6 0])
    % 
    % ViolinPlot(6,0.5,NormSupport(Pico_col,Mic_col,:),200,[0.9 0.9 0])
    % ViolinPlot(7,0.5,NormSupport(Dtm_col,Mic_col,:),200,[0.2 1 0])
    % ViolinPlot(8,0.5,NormSupport(Flag_col,Mic_col,:),200,[0 0.6 0])
    % 
    % set(gca,'XTick',[2.5 7.5])
    % set(gca,'XTickLabel',{'Nanoflag','Microzoo'})
    % set(gca,'box','on')
    % ylabel(['Indirect support of heterotrophic protists',char(10),'(Normalized to phytoplankton production)'])
    % set(gca,'XTickLabelRotation',45)
    % set(gca,'FontSize',10)
    % 
    % %Sub Figure 4 - Indirect Support of mesozoos
    % figure(f4)
    % if Cycle ~= 1
    %     hold on
    %     ViolinPlot(13,0.5,NormSupport(Pico_col,Meso_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(14,0.5,NormSupport(Dtm_col,Meso_col,:),200,[0.2 1 0])
    %     ViolinPlot(15,0.5,NormSupport(Flag_col,Meso_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(18,0.5,NormSupport(Pico_col,VMMeso_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(19,0.5,NormSupport(Dtm_col,VMMeso_col,:),200,[0.2 1 0])
    %     ViolinPlot(20,0.5,NormSupport(Flag_col,VMMeso_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(23,0.5,NormSupport(Pico_col,Macro_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(24,0.5,NormSupport(Dtm_col,Macro_col,:),200,[0.2 1 0])
    %     ViolinPlot(25,0.5,NormSupport(Flag_col,Macro_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(28,0.5,NormSupport(Pico_col,VMMacro_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(29,0.5,NormSupport(Dtm_col,VMMacro_col,:),200,[0.2 1 0])
    %     ViolinPlot(30,0.5,NormSupport(Flag_col,VMMacro_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(33,0.5,NormSupport(Pico_col,Gel_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(34,0.5,NormSupport(Dtm_col,Gel_col,:),200,[0.2 1 0])
    %     ViolinPlot(35,0.5,NormSupport(Flag_col,Gel_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(38,0.5,NormSupport(Pico_col,Salp_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(39,0.5,NormSupport(Dtm_col,Salp_col,:),200,[0.2 1 0])
    %     ViolinPlot(40,0.5,NormSupport(Flag_col,Salp_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(43,0.5,NormSupport(Pico_col,Amph_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(44,0.5,NormSupport(Dtm_col,Amph_col,:),200,[0.2 1 0])
    %     ViolinPlot(45,0.5,NormSupport(Flag_col,Amph_col,:),200,[0 0.6 0])
    % 
    %     set(gca,'XTick',[14.5 19.5 24.5 29.5 34.5 39.5 44.5])
    %     set(gca,'XTickLabel',{'Meso (NVM)','Meso (VM)','Macro (NVM)','Macro (VM)','GelPred','Salps','Amphipods'})
    %     set(gca,'box','on')
    %     ylabel(['Indirect support of mesozooplankton',char(10),'(Normalized to phytoplankton production)'])
    %     set(gca,'XTickLabelRotation',45)
    %     set(gca,'FontSize',10)
    %     %ylim([0 0.08])
    %     xlim([12.5 45.5])
    % 
    % elseif Cycle==1
    %     hold on
    %     ViolinPlot(13,0.5,NormSupport(Pico_col,Meso_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(14,0.5,NormSupport(Dtm_col,Meso_col,:),200,[0.2 1 0])
    %     ViolinPlot(15,0.5,NormSupport(Flag_col,Meso_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(18,0.5,NormSupport(Pico_col,Macro_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(19,0.5,NormSupport(Dtm_col,Macro_col,:),200,[0.2 1 0])
    %     ViolinPlot(20,0.5,NormSupport(Flag_col,Macro_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(23,0.5,NormSupport(Pico_col,VMMacro_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(24,0.5,NormSupport(Dtm_col,VMMacro_col,:),200,[0.2 1 0])
    %     ViolinPlot(25,0.5,NormSupport(Flag_col,VMMacro_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(28,0.5,NormSupport(Pico_col,Gel_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(29,0.5,NormSupport(Dtm_col,Gel_col,:),200,[0.2 1 0])
    %     ViolinPlot(30,0.5,NormSupport(Flag_col,Gel_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(33,0.5,NormSupport(Pico_col,Salp_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(34,0.5,NormSupport(Dtm_col,Salp_col,:),200,[0.2 1 0])
    %     ViolinPlot(35,0.5,NormSupport(Flag_col,Salp_col,:),200,[0 0.6 0])
    % 
    %     ViolinPlot(38,0.5,NormSupport(Pico_col,Amph_col,:),200,[0.9 0.9 0])
    %     ViolinPlot(39,0.5,NormSupport(Dtm_col,Amph_col,:),200,[0.2 1 0])
    %     ViolinPlot(40,0.5,NormSupport(Flag_col,Amph_col,:),200,[0 0.6 0])
    % 
    %     set(gca,'XTick',[14.5 19.5 24.5 29.5 34.5 39.5])
    %     set(gca,'XTickLabel',{'Meso (NVM)','Macro (NVM)','Macro (VM)','GelPred','Salps','Amphipods'})
    %     set(gca,'box','on')
    %     ylabel(['Indirect support of mesozooplankton',char(10),'(Normalized to phytoplankton production)'])
    %     set(gca,'XTickLabelRotation',45)
    %     set(gca,'FontSize',10)
    %     %ylim([0 0.08])
    %     xlim([12.5 40.5])
    % 
    %     upper=0.39;
    %     lower=0.33;
    %     fill([16 32.5 32.5 16 16],[lower lower upper upper lower],'w')
    %     fill([16.4 17.3 17.3 16.4 16.4],[upper-(upper-lower)*1/7 upper-(upper-lower)*1/7 upper-(upper-lower)*2/7 upper-(upper-lower)*2/7  upper-(upper-lower)*1/7],'r','FaceColor',[0.9 0.9 0])
    %     text(17.5,upper-(upper-lower)*1.5/7,'Picophytoplankton','FontSize',6)
    %     fill([16.4 17.3 17.3 16.4 16.4],[upper-(upper-lower)*3/7 upper-(upper-lower)*3/7 upper-(upper-lower)*4/7 upper-(upper-lower)*4/7  upper-(upper-lower)*3/7],'r','FaceColor',[0.2 1 0])
    %     text(17.5,upper-(upper-lower)*3.5/7,'Diatoms','FontSize',6)
    %     fill([16.4 17.3 17.3 16.4 16.4],[upper-(upper-lower)*5/7 upper-(upper-lower)*5/7 upper-(upper-lower)*6/7 upper-(upper-lower)*6/7  upper-(upper-lower)*5/7],'r','FaceColor',[0 0.6 0])
    %     text(17.5,upper-(upper-lower)*5.5/7,'Mix. Flagellates','FontSize',6)
    % end
    % 
    % %Combining all open figures into 1 new figure
    % figs=[f1 f2 f3 f4];
    % f5=figure('Position',[50 50 1000 500]);
    % figure(f5)
    % tcl=tiledlayout(f5,'horizontal');
    % for n=1:length(figs)
    %     figure(figs(n));
    %     ax=gca;
    %     ax.Parent=tcl;
    %     ax.Layout.Tile=n;
    % end
    % 
    % %Figure 2 - Normalized indirect support of fish from salps
    % figure(f6)
    % hold on
    % if Cycle==1 || Cycle==2 || Cycle==4
    %    col = [0 0 1];
    % else
    %    col = [1 0 0];
    % end
    % ViolinPlot(Cycle,0.5,NormSupport(Salp_col,Myct_col,:),200,col)
    % ViolinPlot(Cycle+5,0.5,NormSupport(Salp_col,PiscivFish_col,:),200,col)
    % ViolinPlot(Cycle+10,0.5,NormSupport(Salp_col,GelativFish_col,:),200,col)
    % ViolinPlot(Cycle+15,0.5,NormSupport(Salp_col,Myct_col,:)+NormSupport(Salp_col,PiscivFish_col,:)+NormSupport(Salp_col,GelativFish_col,:),200,col)
    % set(gca,'XTick',[3 8 13 17.5])
    % set(gca,'XTickLabel',{'Myct','Pisciv','Gelativ','All'})
    % set(gca,'box','on')
    % ylabel(['Indirect support of fish',char(10),'(Normalized to salp production)'])
    % xline([5.5:5:15.5],'--')
    % ylim([0 0.01])
    % set(gca,'XTickLabelRotation',45)
    % set(gca,'FontSize',10)
    % 
    % filename=['IndirectC',num2str(Cycle),'.png'];
    % exportgraphics(f5,filename,'Resolution',300)
    % close Figure 2 Figure 3 Figure 4 Figure 5 Figure 6
end
save('prodmean.mat','prodmean');
save('prodLCI.mat','prodLCI');
save('prodUCI.mat','prodUCI');
save('GGEavg.mat','GGEavg');
writematrix(GGESalp,'GGEgelpath.xlsx','Sheet','Salp');
writematrix(GGENonSalp,'GGEgelpath.xlsx','Sheet','NonSalp');
%exportgraphics(f6,'IndirectSF','Resolution',300)