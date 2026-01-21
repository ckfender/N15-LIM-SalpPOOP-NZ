clear all
close all

%%Outputting the modelled flows and their error into an excel file

for Cycle=1:5
    load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
    if Cycle==1
        [num,txt,raw]  = xlsread('N15InverseModelNZC1.xlsx','NZ','B6:B127');
    else
        [num,txt,raw]  = xlsread('N15InverseModelNZ.xlsx','NZ','B6:B145');
    end

    lng=length(MCMCmat(:,1));
    MCMCmat(1:round(lng*0.2),:)=[];

    Columns = txt;
    meanvals = mean(MCMCmat)'./wts;
    stdvals = std(MCMCmat)'./wts;
    if Cycle==1
        meantable = table('Size',[140 3],'VariableTypes',{'string','double','double'},'VariableNames',{'Column','Mean','Std'});
        meantable.Column(1:length(Columns)) = Columns;
        meantable.Mean(1:length(meanvals)) = meanvals;
        meantable.Std(1:length(stdvals)) = stdvals;
    else
        meantable = addvars(meantable,Columns,meanvals,stdvals);
    end
end
%writetable(meantable,'ModelOut.xlsx');