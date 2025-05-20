% JSdemo.m: Demo of the James-Stein estimator with batting averages
%
% This is a good introduction to the important concept of "shrinkage" that
% has grown of out Empirical Bayesian approaches to statistics. Brad Efron
% refers to this as "learning from the experience of others." The deep idea
% is that there is information about Player A's batting ability in the
% observed batting averages of players B, C, D, etc.
%
% RTB wrote it 05 Oct. 2016, procrastinating from reading T32s

%% Read in batting average data
ds = readtable('BattingAverages.xlsx');

% Data consists of data for 18 professional baseball players (rows):
% column 1: 'Player', player's name
% column 2: # hits in first 45 at-bats ("sample")
% column 3: season batting average ("ground truth" mu)

%% Calculate and plot the maximum likelihood estimator

% We first need to convert the text batting averages ('18/45') to numbers
muMLE = zeros(height(ds),1);
for k = 1:height(ds)
    muMLE(k) = eval(ds.hitAB{k});
end

x = [1:height(ds)];     % x values for plots
dkGold = [0.9290    0.6940    0.1250];
%p1 = plot(x,muMLE,'k*');
p1 = plot(x,muMLE,'^','Color',dkGold);
set(p1,'MarkerFaceColor',dkGold);
hold on
set(gca,'XTick',[1:18],'XTickLabel',ds.Player,'XTickLabelRotation',60);
ylabel('Batting Average');

% line for the grand mean across players
hl=line([1,height(ds)],[mean(muMLE),mean(muMLE)]);
set(hl,'Color','k','LineStyle','-');

%% Calculate and plot the James-Stein estimator
dkRed = [0.6350    0.0780    0.1840];
muJS = JSestimator(muMLE,45);
p2 = plot(x,muJS,'o','Color',dkRed);
set(p2,'MarkerFaceColor',dkRed);

%% Plot the season batting average--this is the ground truth value
p3 = plot(x,ds.SeasonBA,'bs');
set(p3,'MarkerFaceColor','b');

% add a legend
legend([p1,p2,p3],'mu^M^L^E','mu^J^S','mu');

%% Compute the ratio of prediction errors for the two estimators

RPE = sum((muJS-ds.SeasonBA).^2) / sum((muMLE - ds.SeasonBA).^2);
tStr = sprintf('Ratio of prediction errors: mu^J^S/mu^M^L^E = %0.2f', RPE);
title(tStr);