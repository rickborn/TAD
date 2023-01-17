function [] = bias_plot(bias,pow)
% bias_plot.m: replicate Figure 1 in Ioannidis 2005

% Values used by Ioanidis:
%pow = [0.8,0.5,0.2];
%bias = [0.05,0.2,0.5,0.8];

R = 0:0.01:1;
alpha = 0.05;
colors = lines(6); % Generate color values
my_colors = [colors(1,:);colors(2,:);colors(5,:);colors(3,:)];

figure('Position',[50 10 600 900],'Name','Fig. 1 of Ioannidis 2005');
for k = 1:length(pow)
    subplot(length(pow),1,k);
    leg_str = cell(1,length(bias));
    for m = 1:length(bias)
        y = ppv2(R,alpha,pow(k),bias(m));
        h = plot(R,y,'Color',my_colors(m,:));
        set(h,'LineWidth',2);
        leg_str(1,m) = {num2str(bias(m))};
        hold on
    end
    axis([-0.05,1,-0.05,1.05])
    xlabel('Pre-study odds, R')
    ylabel('Post-study probability, PPV');
    title(['Power = ' num2str(pow(k))]);
    lgd = legend(leg_str,'Location','northwest');
    title(lgd,'Bias');
    grid on
end