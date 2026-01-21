function [TL] = CalcTL(Cycle,flows,del15N,d15NInputs,Ae,InputCol)

% load('N15GoMInverseCycle1.mat')
% load('N15GoMInverseCycle1Routputs.mat')
% data = MCMCmat;
% flows=data(end,:)./wts';
% Ae=Ae0;

if Cycle==1
    [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,Macro_col,...
        VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
        mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumnsC1(Ae,flows,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));

    % Consumers = [HNF_col,Mic_col,Meso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,...
    %     PiscivFish_col,GelativFish_col,Myct_col];
    Consumers = [HNF_col,Mic_col,Meso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,...
        Myct_col,GelativFish_col,PiscivFish_col];
else
    [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,...
        VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,...
        mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumns(Ae,flows,del15N(1,:),d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));

    Consumers = [HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,...
        PiscivFish_col,GelativFish_col,Myct_col];
end
Dets = [sdet_col,mdet_col,ldet_col,SalpDet_col];
Dons = Dom_col;
Flags = Flag_col;
HBacs = bac_col;
Phys=[Pico_col,Dtm_col];
Nuts = [NO3_col,NH4_col];

%TL_final = sum(TL_source*input_source)/input_total
%for i=1:length(data(:,1))
    %flows=data(i,:);
    TL=ones(length(Ae(:,1)),1)*2;
    for iter=1:10 %The whole process is done a number of times so the TL of the last calc affects the next
        TL(Phys)=1;
        TL(Nuts)=0;
        for j=1:length(Flags)
            ToFlags = find(Ae(Flags(j),:)==1);
            total=0;
            temp=0;
            for k=1:length(ToFlags)
                ind=find(Ae(:,ToFlags(k))==-1); %the column of the ingestion source
                if ind~=Flags(1) %if the outflow is not...intraguild predation?
                    total=total+flows(ToFlags(k)); %A running sum of ingestion
                    temp=temp+flows(ToFlags(k))*TL(ind); %Running sum of TL*ingestion
                end
            end
            TL(Flags(j),1)=temp./total+1; %Final TL
        end
        for j=1:length(Consumers)
            if Consumers(j) == Meso_col
                if Cycle==1
                    can=48;
                    ToConsumers = sort([find(Ae(Consumers(j),:)==1) can]);
                elseif Cycle ~= 1
                    can=52;
                    ToConsumers = sort([find(Ae(Consumers(j),:)==1) can]);
                end
                total=0;
                temp=0;
                for k=1:length(ToConsumers)
                    ind=find(Ae(:,ToConsumers(k))==-1);
                    if isempty(ind)
                        ind=Meso_col;
                    end
                    total=total+flows(ToConsumers(k));
                    temp=temp+flows(ToConsumers(k))*TL(ind);
                end
            elseif Cycle~=1 & Consumers(j) == VMMeso_col
                can=63;
                ToConsumers = sort([find(Ae(Consumers(j),:)==1) can]);
                total=0;
                temp=0;
                for k=1:length(ToConsumers)
                    ind=find(Ae(:,ToConsumers(k))==-1);
                    if isempty(ind)
                        ind=VMMeso_col;
                    end
                    total=total+flows(ToConsumers(k));
                    temp=temp+flows(ToConsumers(k))*TL(ind);
                end
            else
                ToConsumers = find(Ae(Consumers(j),:)==1);
                total=0;
                temp=0;
                for k=1:length(ToConsumers)
                    ind=find(Ae(:,ToConsumers(k))==-1);
                    total=total+flows(ToConsumers(k));
                    temp=temp+flows(ToConsumers(k))*TL(ind);
                end
            end
            TL(Consumers(j),1)=temp./total+1;
        end
        for j=1:length(Dons)
            ToDons = find(Ae(Dons(j),:)==1);
            total=0;
            temp=0;
            for k=1:length(ToDons)
                ind=find(Ae(:,ToDons(k))==-1);
                if length(ind)==1  %Note this is necessary to avoid issues associated with external DOM not having a TL.
                    total=total+flows(ToDons(k));
                    temp=temp+flows(ToDons(k))*TL(ind);
                end
            end
            TL(Dons(j),1)=temp./total;
        end
        for j=1:length(HBacs)
            ToHBacs = find(Ae(HBacs(j),:)==1);
            total=0;
            temp=0;
            for k=1:length(ToHBacs)
                ind=find(Ae(:,ToHBacs(k))==-1);
                total=total+flows(ToHBacs(k));
                temp=temp+flows(ToHBacs(k))*TL(ind);
            end
            TL(HBacs(j),1)=temp./total+1;
        end
        for j=1:length(Dets)
            ToDets = find(Ae(Dets(j),:)==1);
            total=0;
            temp=0;
            for k=1:length(ToDets)
                ind=find(Ae(:,ToDets(k))==-1);
                if length(ind)==1  %Note this is necessary to avoid issues associated with external Det not having a TL.
                    total=total+flows(ToDets(k));
                    if length(intersect(ind,[Phys,Flags]))==1
                        temp=temp+flows(ToDets(k))*TL(ind);
                    else
                        temp=temp+flows(ToDets(k))*(TL(ind)-1);   %Note assuming that egesta has a TL equal to your food
                    end
                end
            end
            TL(Dets(j),1)=temp./total;
        end
        TL(length(TL),1)=TL(Salp_col);
    end
%end
        