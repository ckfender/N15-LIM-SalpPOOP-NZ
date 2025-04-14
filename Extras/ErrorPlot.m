clear all
close all

f1=figure('Position',[100 100 1000 500]);
for Cycle=1:5
    clearvars -except Cycle f1

    load(['N15NZInverseCycle',num2str(Cycle),'.mat'])
    load(['N15NZInverseCycle',num2str(Cycle),'Routputs.mat'])
    if Cycle==1
        [num,txt,raw]  = xlsread('N15InverseModelNZC1.xlsx','NZ','B6:B127');
    else
        [num,txt,raw]  = xlsread('N15InverseModelNZ.xlsx','NZ','B6:B145');
    end

    Columns = txt;
    means(:,Cycle) = mean(MCMCmat)';
    l=size(MCMCmat);
    final(:,Cycle) = MCMCmat(l(1),:)'./wts;
    meanvals(:,Cycle) = mean(MCMCmat)'./wts;
    meanvalsplain(:,Cycle) = mean(MCMCmatplain)'./wts;
    stdvals(:,Cycle) = std(MCMCmat)'./wts;
    InputCols(Cycle)=InputCol;

    for i=1:length(MCMCmat(:,1))

        MCMC=MCMCmat(i,:);
        d15N=del15N(i,:);

        if Cycle==1
            [Aa] = ResetRN15C1(Aa0,MCMC,wts,d15N,d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
        else
            [Aa] = ResetRN15(Aa0,MCMC,wts,d15N,d15NInputs(1,InputCol),d15NInputs(2,InputCol),d15NInputs(3,InputCol));
        end

        AaOut(:,i)=Aa*MCMC';
    end

    AaOutmean(:,Cycle)=mean(AaOut')';
    AaOutstd(:,Cycle)=std(AaOut')';

    col=[117 197 250; %Cycle 1 - salp - light blue
        3 154 255;  %Cycle 2 - salp - blue
        252 111 118;%Cycle 3 - nonsalp - light red
        3 79 130;  %Cycle 4 - salp - dark blue
        252 3 15];  %Cycle 5 - nonsalp - red
    col=col/255;

    %Plot of measurements vs inputs
    %Circles will be in situ measurments, squares will be Pinkerton fish ing
    %assumptions
    sz=20; %font size for markers
    figure(f1)
    hold on
    if Cycle ==1
        xmean=Inputs(:,Cycle*2-1);
        xmin=Inputs(:,Cycle*2);
        xminmod=xmin;
        isneg=Inputs(:,Cycle*2-1)-xminmod<=0;
        xminmod(isneg)=Inputs(isneg,Cycle*2-1)-1E-9;
        ymean=AaOutmean(:,Cycle);
        ymin=AaOutstd(:,Cycle);
        yminmod=ymin;
        isneg=AaOutmean(:,Cycle)-yminmod<=0;
        yminmod(isneg)=AaOutmean(isneg,Cycle)-1E-9;
        errorbar(xmean(1:19),ymean(1:19),yminmod(1:19),AaOutstd(1:19,Cycle),"vertical","LineStyle","none",'CapSize',0,'col','k')
        errorbar(xmean(1:19),ymean(1:19),xminmod(1:19),Inputs(1:19,Cycle*2),"horizontal","LineStyle","none",'CapSize',0,'col','k')
        errorbar(xmean(20:36),ymean(20:36),yminmod(20:36),AaOutstd(20:36,Cycle),"vertical","LineStyle","none",'CapSize',0,'col','k')
        errorbar(xmean(20:36),ymean(20:36),xminmod(20:36),Inputs(20:36,Cycle*2),"horizontal","LineStyle","none",'CapSize',0,'col','k')
        xmean=Inputs(:,Cycle*2-1);
        ymean=AaOutmean(:,Cycle);
        scatter(xmean(1:19),ymean(1:19),sz,col(Cycle,:),'filled') %rate meas as circles
        scatter(xmean(20:36),ymean(20:36),sz,col(Cycle,:),'filled','s') %fish ing as squares
    else
        xmean=Inputs(:,Cycle*2-1);
        xmin=Inputs(:,Cycle*2);
        xminmod=xmin;
        isneg=Inputs(:,Cycle*2-1)-xminmod<=0;
        xminmod(isneg)=Inputs(isneg,Cycle*2-1)-1E-9;
        ymean=AaOutmean(:,Cycle);
        ymin=AaOutstd(:,Cycle);
        yminmod=ymin;
        isneg=AaOutmean(:,Cycle)-yminmod<=0;
        yminmod(isneg)=AaOutmean(isneg,Cycle)-1E-9;
        errorbar(xmean(1:20),ymean(1:20),yminmod(1:20),AaOutstd(1:20,Cycle),"vertical","LineStyle","none",'CapSize',0,'col','k')
        errorbar(xmean(1:20),ymean(1:20),xminmod(1:20),Inputs(1:20,Cycle*2),"horizontal","LineStyle","none",'CapSize',0,'col','k')
        errorbar(xmean(21:39),ymean(21:39),yminmod(21:39),AaOutstd(21:39,Cycle),"vertical","LineStyle","none",'CapSize',0,'col','k')
        errorbar(xmean(21:39),ymean(21:39),xminmod(21:39),Inputs(21:39,Cycle*2),"horizontal","LineStyle","none",'CapSize',0,'col','k')
        xmean=Inputs(:,Cycle*2-1);
        ymean=AaOutmean(:,Cycle);
        scatter(xmean(1:20),ymean(1:20),sz,col(Cycle,:),'filled')
        scatter(xmean(21:39),ymean(21:39),sz,col(Cycle,:),"filled",'s')
    end


    xlabel('Measurements (mmol N m^{-2} d^{-1})')
    ylabel('Model (mmol N m^{-2} d^{-1})')
    set(gca,'YScale','log')
    set(gca,'XScale','log')
    xlim([10^-8 10^2])
    ylim([10^-8 10^2])
    plot(xlim,ylim,'--k')
    hold on
    if Cycle==5
        plot(xlim,ylim,'--k')
        qw{1}=scatter(nan,nan,'MarkerFaceColor',col(1,:),'MarkerEdgeColor',col(1,:),'DisplayName','Cycle 1');
        qw{2}=scatter(nan,nan,'MarkerFaceColor',col(2,:),'MarkerEdgeColor',col(2,:),'DisplayName','Cycle 2');
        qw{3}=scatter(nan,nan,'MarkerFaceColor',col(3,:),'MarkerEdgeColor',col(3,:),'DisplayName','Cycle 3');
        qw{4}=scatter(nan,nan,'MarkerFaceColor',col(4,:),'MarkerEdgeColor',col(4,:),'DisplayName','Cycle 4');
        qw{5}=scatter(nan,nan,'MarkerFaceColor',col(5,:),'MarkerEdgeColor',col(5,:),'DisplayName','Cycle 5');
        legend([qw{:}],{'Cycle 1','Cycle 2','Cycle 3','Cycle 4','Cycle 5'},'Location','southeast')
    end
end
%exportgraphics(f1,'Error.png','Resolution',300)

