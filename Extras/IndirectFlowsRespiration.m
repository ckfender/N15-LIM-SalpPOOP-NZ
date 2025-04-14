function [P,struct,ext,r,e,p,Pspec,structspec] = IndirectFlowsRespiration(Cycle,Ae_keep,flows)

%Note that this version of the code treats NH4 excretion as a loss from the ecosystem.
%It also differs from Stukel's previous versions in that it includes deep
%respiration in the r vector, allowing it to represent indirect support of
%migrators and benthic groups

% clear all
% close all
% 
% Cycle=1
% load(['N15GoMInverseCycle',num2str(Cycle),'.mat'])
% load(['N15GoMInverseCycle',num2str(Cycle),'Routputs.mat'])
% Ae_keep = Ae;
% flows=mean(MCMCmat);
% clearvars -except Ae_keep flows


load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
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
%clearvars -except NH4_col flows Ae_keep
Ae=Ae_keep;
clear Ae_keep p


%This code is from methods explained in Hannon (1973)



%Calculating the number of different inputs to the ecosystem.
%tmp represents each flow, with values of 1 and -1 representing flows that
%are net inputs of energy (lateral input of DON, deep upwelling, etc) or loss 
%terms (ingestion by unmodelled HTL, export out of the deeper depth) respectively
tmp=sum(Ae);
numInputs = sum(tmp==1);
clear tmp
numCompartments = length(Ae(:,1));

%Lookup will be a tool to map the solution vector into a matrix.
%The first [numCompartments] rows/columns correspond to the compartments, with 
%columns being inputs and rows being outputs. So if the first entry in the 
%first column is a value of 4 at the second row, that reads as the first 
%balanced input to the first compartment occurs in the 4th flow. By the same 
%token if the first entry in the 25th row is a 1, it represents that the first 
%output for 25th compartment is coming from flow 1. The
%ones after that (numCompartments+1) represent each sequential input
%Thus an entry of 2 in the 22nd column of the [numCompartments+1]th row
%represents the first input occurs in the second flow feeding to the 22nd
%compartment.
%TLDR the end result is a matrix where the elements in a given column
%represent the flows from the column compartment to the row compartment
%(excluding instances or net input of output including respiration or sinks)
Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
%r_Lookup records the flows to ammonia (shallow respiration). Both it and ext_Lookup seem to be
%formatted differently with sequential respiration or output events as each
%column. Thus an ext_Lookup first column value of 107 in row 15 means the
%first output event happens in flow 107 out of compartment 15
%%This structure only allows for the calculation of indirect support for 
%%exclusively shallow dwelling compartments. 
r_Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
%%An additional r_Lookup just for deep respiration, done semi-manually and
%%to be recombined with r later
rdeep_Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
if Cycle==1
    deepresp = [70 82 93 103 107];
else
    deepresp = [72 87 99 110 120 124];
end
%ext_Lookup records loss terms out of the system (excluding deep
%respiration)
ext_Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
%An additional Lookup for losses to HTL or biomass that should be counted
%as production
extprod_Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
extprodHTL_Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
extprodBiom_Lookup = zeros(numInputs+numCompartments,numInputs+numCompartments);
if Cycle==1
    extprod = [53 59 66 75 85 92 99 106 110];
    extprodHTL = [53 59 66 75 92 99];
    extprodBiom = [85 106 110];
else
    extprod = [58 68 76 83 92 102 109 116 123 127];
    extprodHTL = [58 68 76 83 92 109 116];
    extprodBiom = [102 123 127];
end

%%Hannon's theory assumes a perfectly balanced system in which the
%%difference between inputs and outputs is due solely to respiration, but
%%this is not the case in our system as we have multiple things that might
%%contribute to that difference. Thus we must construct the r vector
%%manually.
%Initializing stepping variables to keep track of special cases, being
%flows that are inputs, outputs (loss terms), or respiration to NH4
whichinput=1;  %Need to keep track of which input we are on
whichoutput=1;  %This is obnoxious, but necessary because we have multiple outputs from some compartments
whichr=1; %This is obnoxious, but necessary because some organisms excrete to surface and others excrete to depth
whichrdeep=1;
whichprod=1;
whichprodHTL=1;
whichprodBiom=1;
if Cycle == 1
    Mcan = 48;
else
    Mcan = 52;
    VMcan = 63;
end

%Constructing the Lookup matrices
%for each flow in the model, 
for i=1:length(Ae(1,:))
    if sum(Ae(:,i))==0 %if the flow is balanced
        tmp1=find(Ae(:,i)==1); %record the input columns
        tmp2=find(Ae(:,i)==-1); %record the output columns
        if Lookup(tmp2,tmp1)==0 & tmp1~=NH4_col
            Lookup(tmp2,tmp1)=i; % if not shallow respiration, record the input and output columns in Lookup
        elseif tmp1==NH4_col
            r_Lookup(tmp2,whichr)=i; %if it is shallow respiration, record output column in r_Lookup instead
            whichr=whichr+1; %we then increment a sequential step per respiratory flow
        elseif i==Mcan %Manually adding meso canibalism
            Lookup(Meso_col,Meso_col)=i;
        elseif i==VMcan
            Lookup(VMMeso_col,VMMeso_col)=i;
        else
            ['Error at ',num2str(i)]
            Lookup(tmp2,tmp1)=i;
        end
    elseif sum(Ae(:,i))==1 %if the flow is a net source of energy (i.e. upwelling)
        tmp1=find(Ae(:,i)==1); %record the input column
        if Lookup(numCompartments+whichinput,tmp1)==0
            Lookup(numCompartments+whichinput,tmp1)=i; %record input column at row corresponding to whichinput
        else
            Lookup(numCompartments+whichinput,tmp1)=i;
            ['Error at ',num2str(i)]
        end
        whichinput=whichinput+1;
    elseif sum(Ae(:,i))==-1 & ismember(i,deepresp) %if the flow is deep resp
        tmp2=find(Ae(:,i)==-1); %record output compartment
        rdeep_Lookup(tmp2,whichrdeep)=i; %record compartment as entry with the flow at the output compartment's row
        whichrdeep=whichrdeep+1; %increment output counter
    elseif sum(Ae(:,i))==-1 & ismember(i,extprod) %if the flow is loss to extra production (biom,HTL)
        tmp2=find(Ae(:,i)==-1); 
        extprod_Lookup(tmp2,whichprod)=i; 
        whichprod=whichprod+1; 
    elseif sum(Ae(:,i))==-1 %if the flow is any other output (deep dom, sinking, egestion)
        tmp2=find(Ae(:,i)==-1); 
        ext_Lookup(tmp2,whichoutput)=i; 
        whichoutput=whichoutput+1; 
    end
end
%secondary look for special lookups
for i=1:length(Ae(1,:))
    if sum(Ae(:,i))==-1 & ismember(i,extprodHTL) %if the flow is loss to HTL consumption
        tmp2=find(Ae(:,i)==-1); 
        extprodHTL_Lookup(tmp2,whichprodHTL)=i; 
        whichprodHTL=whichprodHTL+1;
    elseif sum(Ae(:,i))==-1 & ismember(i,extprodBiom) %if the flow is loss to biomass production
        tmp2=find(Ae(:,i)==-1); 
        extprodBiom_Lookup(tmp2,whichprodBiom)=i; 
        whichprodBiom=whichprodBiom+1;
    end
end

%P is the input-output matrix of production energy flows (Table 2 in Hannon), with the inputs
%to a given compartment listed in the column of that compartment and the
%outputs as the rows. Thus each element of the matrix represents the energy
%output of the row compartment used as input by the column compartment 
%to be ultimately used for production and respiration.
%Practically it is the Lookup matrix but replacing the flow indices in said
%matrix with the actual MCMC flow values. Note that because our model
%structure includes non-respiratory loss terms its interpretation is a
%little different The difference between column sum and row sum is not
%respiration as it was for Hannon, for example. The row sum is also not
%necessarily production for groups which have production losses as sink
%terms (these will be added back in manually later)
P = zeros(numInputs+numCompartments,numInputs+numCompartments);
for i=1:numInputs+numCompartments
    for j=1:numInputs+numCompartments
        if Lookup(i,j)~=0
            P(i,j)=flows(Lookup(i,j));
        end
    end
end
%A quick tweak to the P matrix as the flow of salps to Gelativ is actually
%represented by a separate compartment SalpMort
P(Salp_col,GelativFish_col)=P(SalpMort_col,GelativFish_col);
P(Salp_col,SalpMort_col)=P(Salp_col,SalpMort_col)-P(Salp_col,GelativFish_col);
P(SalpMort_col,GelativFish_col)=0;

%As above but with the extprod flows
prod=zeros(size(extprod_Lookup));
for i=1:length(extprod_Lookup(:,1))
    for j=1:length(extprod_Lookup(1,:))
        if extprod_Lookup(i,j)~=0
            prod(i,j)=flows(extprod_Lookup(i,j));
        end
    end
end
prod=sum(prod')'; %Summing output flows across each compartment and saving as a vector
%Special P matrix with biomass flows as identity and an extra compartment
%for HTL flows. This is not to be used with any Hanson calculations and is
%for other purposes requiring a more complete production matrix.
Pspec=P;
Pspec(:,end+1)=0;
Pspec(end+1,:)=0;
for i=1:length(extprodHTL_Lookup(:,1))
    for j=1:length(extprodHTL_Lookup(1,:))
        if extprodHTL_Lookup(i,j)~=0
            Pspec(i,end)=flows(extprodHTL_Lookup(i,j));
        end
    end
end
for i=1:length(extprodBiom_Lookup(:,1))
    for j=1:length(extprodBiom_Lookup(1,:))
        if extprodBiom_Lookup(i,j)~=0
            Pspec(i,i)=flows(extprodBiom_Lookup(i,j));
        end
    end
end
%Regrettably yet another version of the P matrix that DOES include HTL but
%does NOT include the biomass production of fish
Pspec2 = P;
Pspec2(:,end+1)=0;
Pspec2(end+1,:)=0;
Pspec2(:,end)=Pspec(:,end);


%As above but with the output flows
ext=zeros(size(ext_Lookup));
for i=1:length(ext_Lookup(:,1))
    for j=1:length(ext_Lookup(1,:))
        if ext_Lookup(i,j)~=0
            ext(i,j)=flows(ext_Lookup(i,j));
        end
    end
end
ext=sum(ext')'; %Summing output flows across each compartment and saving as a vector

%As above for the respiration flows
r=zeros(size(r_Lookup));
for i=1:length(r_Lookup(:,1))
    for j=1:length(r_Lookup(1,:))
        if r_Lookup(i,j)~=0
            r(i,j)=flows(r_Lookup(i,j));
        end
    end
end
for i=1:length(rdeep_Lookup(:,1))  %adding deep respiration to the r matrix
    for j=1:length(rdeep_Lookup(1,:))
        if rdeep_Lookup(i,j)~=0
            r(i,j)=r(i,j)+flows(rdeep_Lookup(i,j));
        end
    end
end
r=sum(r')'; %Summing respiration flows across each compartment (rows), transposing, and saving as a vector

%The normalized production matrix G (Table 3), i.e. the relative proportions of total 
%input made up by each flow to a given compartment. If the element in the
%second row of the first column is .95, this represents NO3 getting 95% of
%it's inputs from NH4 (through ammonification).
%Note that as a result of segregating the respiration and export flows,
%they are not part of the P matrix that we are doing the math on.
G = zeros(numInputs+numCompartments,numInputs+numCompartments);
for i=1:numInputs+numCompartments %for each column of the matrix (each compartment)
    if sum(P(:,i))>0 %If the sum of balanced flows is positive (i.e. it's not an NH4 column)
        G(:,i)=P(:,i)/sum(P(:,i)); %record the proportion of the total inflow that flow represents
    end
    e(i)=sum(P(:,i)); %e is the direct energy vector, the summed inflows (columns) for a given compartment
    p(i)=sum(P(i,:)); %p is production vector, the summed outflows (rows) of a compartment ONLY into the shallows
    %entries for which e and p are not the same are those that contain net
    %inflows or outflows to the model (NH4 or any heterotroph column due to 
    %all respiration being treated as a loss)
end
%Adding extra production flows to production vector
p=p+prod';
Gspec = zeros(numInputs+numCompartments+1,numInputs+numCompartments+1);
for i=1:(numInputs+numCompartments+1) %for each column of the matrix (each compartment)
    if sum(Pspec2(:,i))>0 %If the sum of balanced flows is positive (i.e. it's not an NH4 column)
        Gspec(:,i)=Pspec2(:,i)/sum(Pspec2(:,i)); %record the proportion of the total inflow that flow represents
    end
end

%An identity matrix (diagonal line of 1s) with the same size as P/G. 
I = eye(numInputs+numCompartments);
Ispec = eye(numInputs+numCompartments+1);

struct=inv(I-G);
structspec=inv(Ispec-Gspec);

% (I-G)^-1 codifies the structural elements of the ecosysten since r
% can be varied to produce changes in e without changing (I-G)-‘. Each
% element represents the total energy flow both directly and indirectly from the
% ith column to the row component per unit of respiration of the row component. These
% elements are the key to understanding the interdependence of each
% component on the other.  (page 537 of Hannon (1973) )