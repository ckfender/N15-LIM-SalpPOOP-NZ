
% Creating data
%close all
clearvars -except TL_EE EE Prod GGE_EE
clc

figure
for i=1:3
%     x0=normrnd(rand,rand,[60,1]);
%     y0=normrnd(rand, rand,[60,1]);
    x0=randn([60,1])*rand+rand;
    y0=randn([60,1])*rand+rand;
    orientation_rad =rand*pi;
    cos_phi = cos( orientation_rad );
    sin_phi = sin( orientation_rad );
    R = [ cos_phi sin_phi; -sin_phi cos_phi ];
    xy=[x0,y0]*R;
    x{i}=xy(:,1);
    y{i}=xy(:,2);
end


% Plot the data with ellipse

for i=1:length(x)
    %p(i)=plot(x{i},y{i},'o');
    p(i)=plot(mean(x{i}),mean(y{i}),'o');
    hold on
    handle_ellipse= plot_ellipse(x{i},y{i});
    handle_ellipse.Color=p(i).Color;
end
