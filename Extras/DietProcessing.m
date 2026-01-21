load("DietPropMean.mat") %Array of full diet constituent proportion means
load("DietPropLCI.mat")  %Full diet LCIS 
load("DietPropUCI.mat")  %Full diet UCIs
load("HerbProp.mat")     %Mean, LCI, and UCI of just herbivory in each group
load("NPPProp.mat")      %Mean, LCI, and UCI of proportion of NPP consumed by each group
load("ProProp.mat")      %Mean, LCI, and UCI of feeding on HNF+Mic
load("MesoProp.mat")     %Meso+VMMeso
load("MacroProp.mat")    %Macro+VMMacro

for i = 1:5
    %Aiding visualization by compressing arrays into matrices of each group
    %Noting these are feeding BY these groups, not feeding ON them
    SalpDietMean(:,i) = DietPropMean(:,13,i);
    SalpDietLCI(:,i) = DietPropLCI(:,13,i);
    SalpDietUCI(:,i) = DietPropUCI(:,13,i);

    MyctDietMean(:,i) = DietPropMean(:,17,i);
    MyctDietLCI(:,i) = DietPropLCI(:,17,i);
    MyctDietUCI(:,i) = DietPropUCI(:,17,i);

    PiscivDietMean(:,i) = DietPropMean(:,15,i);
    PiscivDietLCI(:,i) = DietPropLCI(:,15,i);
    PiscivDietUCI(:,i) = DietPropUCI(:,15,i);

    GelatDietMean(:,i) = DietPropMean(:,16,i);
    GelatDietLCI(:,i) = DietPropLCI(:,16,i);
    GelatDietUCI(:,i) = DietPropUCI(:,16,i);

    if i~=1

        MesoDietPropMean(:,i) = mean(DietPropMean(:,[Meso_col VMMeso_col],i),2);
        MesoDietPropLCI(:,i) = mean(DietPropLCI(:,[Meso_col VMMeso_col],i),2);
        MesoDietPropUCI(:,i) = mean(DietPropUCI(:,[Meso_col VMMeso_col],i),2);

        MacroDietPropMean(:,i) = mean(DietPropMean(:,[Macro_col VMMacro_col],i),2);
        MacroDietPropLCI(:,i) = mean(DietPropLCI(:,[Macro_col VMMacro_col],i),2);
        MacroDietPropUCI(:,i) = mean(DietPropUCI(:,[Macro_col VMMacro_col],i),2);

        MesoHerbProp(:,i) = mean(HerbProp(:,[Meso_col VMMeso_col],i),2);
        MesoProProp(:,i) = mean(ProProp(:,[Meso_col VMMeso_col],i),2);

        MacroHerbProp(:,i) = mean(HerbProp(:,[Macro_col VMMacro_col],i),2);
        MacroProProp(:,i) = mean(ProProp(:,[Macro_col VMMacro_col],i),2);

    else
        MesoDietPropMean(:,i) = mean(DietPropMean(:,[Meso_col],i),2);
        MesoDietPropLCI(:,i) = mean(DietPropLCI(:,[Meso_col],i),2);
        MesoDietPropUCI(:,i) = mean(DietPropUCI(:,[Meso_col],i),2);

        MacroDietPropMean(:,i) = mean(DietPropMean(:,[Macro_col VMMacro_col],i),2);
        MacroDietPropLCI(:,i) = mean(DietPropLCI(:,[Macro_col VMMacro_col],i),2);
        MacroDietPropUCI(:,i) = mean(DietPropUCI(:,[Macro_col VMMacro_col],i),2);

        MesoHerbProp(:,i) = mean(HerbProp(:,[Meso_col],i),2);
        MesoProProp(:,i) = mean(ProProp(:,[Meso_col],i),2);

        MacroHerbProp(:,i) = mean(HerbProp(:,[Macro_col VMMacro_col],i),2);
        MacroProProp(:,i) = mean(ProProp(:,[Macro_col VMMacro_col],i),2);
    end
end
