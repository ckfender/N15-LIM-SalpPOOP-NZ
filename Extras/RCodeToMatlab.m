clearvars
close all

%%A script written by Mike Stukel and tweaked by Christian Fender to help
%%in translating the portions of the model written in R into Matlab as is
%%required for certain functions like GetColumns.m to avoid a lot of
%%copy-pasting followed by manual grammatical changes. Unfortunately it
%%doesn't work very well and only catches half of what I would want it to.

fid = fopen('ExternalFunctionsNZC1.R');
C = textscan(fid,'%s','Delimiter','@@');

lines=length(C{1,1});

for i=1:lines
    temp=char(C{1,1}(i));
    for j=1:length(temp)-3
        if strcmp(temp(j:j+3),' <- ')
            temp(j:j+3) = ' =  ';
        end
    end
    for j=1:length(temp)-1
        if strcmp(temp(j:j+1),',]')
            temp(j+2:end+1) = temp(j+1:end);
            temp(j:j+2)=',:]';
        end
    end
    for j=1:length(temp)
        if strcmp(temp(j),'[')
            temp(j) = '(';
        end
        if strcmp(temp(j),']')
            temp(j) = ')';
        end
        if strcmp(temp(j),'#')
            temp(j) = '%';
        end
    end
    if strcmp(temp(length(temp)),'+')
        temp(end+1:end+3)='...';
    elseif length(temp)>1
        if strcmp(temp(end-1:end),'+ ')
            temp(end+1:end+3)='...';
        end
    end
    if length(temp)>22
        if strcmp(temp(1:22),'ResetRN15 =  function(')
            temp(31:end+8)=temp(23:end);
            temp(1:30)='function [output] = ResetRN15(';
            temp(end)=[];
        else
            temp(length(temp)+1)=';';
        end
    elseif length(temp)>8
        if strcmp(temp(1:7),'return(')
            temp(1:7)='output=';
            temp(end)=';';
            endline=i;
        else
            temp(length(temp)+1)=';';
        end
    else
        temp(length(temp)+1)=';';
    end
    C{1,1}(i)={temp};
end
fclose(fid)

filenameout='ResetRN15test.m';
%fid=fopen(filenameout);
i=1;
dlmwrite(filenameout,C{1,1}{i,1},'delimiter','');
dlmwrite(filenameout,'','-append','delimiter','','newline','pc');
for i=2:endline
    dlmwrite(filenameout,C{1,1}{i,1},'-append','delimiter','');
    dlmwrite(filenameout,'','-append','delimiter','','newline','pc');
end


filenameout='GetColumnstest.m';
%fid=fopen(filenameout);
i=1;
inputs=C{1,1}{1,1}(30:end);
functionline='function [NO3_col,NH4_col,Pico_col,Dtm_col,Flag_col,HNF_col,Mic_col,Meso_col,VMMeso_col,Macro_col,VMMacro_col,Gel_col,Salp_col,Amph_col,PiscivFish_col,GelativFish_col,Myct_col,bac_col,sdet_col,mdet_col,ldet_col,SalpDet_col,Dom_col,SalpMort_col,Export_col,upNO3_col] = GetColumns(input)';
functionline(484:490)=[];
functionline=[functionline,inputs]
dlmwrite(filenameout,functionline,'delimiter','');
dlmwrite(filenameout,'','-append','delimiter','','newline','pc');
for i=2:endline
    dlmwrite(filenameout,C{1,1}{i,1},'-append','delimiter','');
    dlmwrite(filenameout,'','-append','delimiter','','newline','pc');
end