function [DietProp,DietFlux] = CalcDiet(Cycle,flows,del15N,d15NInputs,Ae,InputCol,FeedingMat)

% load('N15GoMInverseCycle1.mat')
% load('N15GoMInverseCycle1Routputs.mat')
% data = MCMCmat;
% flows=data(end,:)./wts';
% Ae=Ae0;

if Cycle==1
    [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,Macro_col,...
        VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
        mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumnsC1(Ae,flows,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));

    Consumers = [Flag_col,HNF_col,Mic_col,Meso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,...
        PiscivFish_col,GelativFish_col,Myct_col];
else
    [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,...
        VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
        mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumns(Ae,flows,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));

    Consumers = [Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,...
        PiscivFish_col,GelativFish_col,Myct_col];
end

DietProp=zeros(height(Ae),height(Ae));
DietFlux=zeros(height(Ae),height(Ae));
for j=1:length(Consumers)
    ToConsumer = find(~isnan(FeedingMat(:,Consumers(j))));
    total=sum(flows(FeedingMat(ToConsumer,Consumers(j))));
    for k=1:length(ToConsumer)
        temp=flows(FeedingMat(ToConsumer(k),Consumers(j)));
        DietFlux(ToConsumer(k),Consumers(j))=temp;
        DietProp(ToConsumer(k),Consumers(j))=temp/total*100;
    end
end 

if Cycle==1
    %Adding VMMeso Column
    T=zeros(23,1);
    DietProp=[DietProp(:,1:8),T,DietProp(:,9:23)];
    DietFlux=[DietFlux(:,1:8),T,DietFlux(:,9:23)];
    %Adding VMMeso Row
    T=zeros(1,24);
    DietProp=[DietProp(1:8,:);T;DietProp(9:23,:)];
    DietFlux=[DietFlux(1:8,:);T;DietFlux(9:23,:)];
end