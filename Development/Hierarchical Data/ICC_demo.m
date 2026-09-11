function [allModels] = ICC_demo(panel,nObsPerGroup)
%
% reproducing figures 2 & 3 of Aarts et al. 2014

%% Defaults
if nargin < 2, nObsPerGroup = 10; end
if nargin < 1, panel = 'c'; end

nGroups = 4;
switch lower(panel)
    case 'a'
        grpMeans = [10,10,10,10];
    case 'b'
        grpMeans = [8.5,8.5,11.5,11.5];
    case 'c'
        grpMeans = [9.5,8.5,11.5,10.5];
    otherwise
        grpMeans = [9.5,8.5,11.5,10.5];
end

%% Set up data table
T = table();

% Generate the data:
y = randn(nObsPerGroup,nGroups);
for k = 1:nGroups
    y(:,k) = y(:,k) + grpMeans(k);
end
T.y = y(:);

% Assign group IDs
T.grpID = [ones(nObsPerGroup,1);ones(nObsPerGroup,1).*2;ones(nObsPerGroup,1).*3;ones(nObsPerGroup,1).*4];

% Assign unique observation IDs
T.obsID = repmat((1:nObsPerGroup)',nGroups,1);

% Assign treatment:
T.rx = [zeros(nObsPerGroup*2,1);ones(nObsPerGroup*2,1)];

%% Calculate the intracluster correlations (ICC)
allGroups = unique(T.grpID);
varSum = 0;
gMeans = zeros(nGroups,1);
for k = 1:nGroups
    thisGroup = allGroups(k);
    thisGroupIdx = T.grpID == thisGroup;
    % calculate the variance for this group & add to varSum:
    varSum = varSum + var(T.y(thisGroupIdx),1);
    % calculate the mean for this well:
    gMeans(k) = mean(T.y(thisGroupIdx));
end

% I'm a bit puzzled here as to whether the 'within group variance' (varW)
% should be the sum of the variance for the 4 groups (i.e. varSum) or
% should be the average of the within-group variances. The former seems
% more correct, but then the two variances (varB + varW) don't add up to
% the total variance in y: var(y(:),1)
varW = varSum / nGroups;
varB = var(gMeans,1);
ICCgroup = varB / (varB + varW);
fprintf(1,'ICC: %.3f\n', ICCgroup);

%% Plot the raw data:
x = jitter(T.grpID,0.2);
f1 = figure('Name','Aarts et al. 2014, Figs. 2 & 3','Position',[50 600 1200 400]);
subplot(1,3,1);
gscatter(x,T.y,T.grpID);
hold on
% add dashed lines for mean of each group
for k = 1:4
    hl = line([k-0.5,k+0.5],[gMeans(k),gMeans(k)]);
    set(hl,'LineWidth',1.5,'LineStyle','--', 'Color','k');
end
ax = axis;
legend off
xlabel('Group')
ylabel('y')
title(['\rm Variance \bf within \rm groups: ' num2str(varW,2)]);

%% Plot the means
subplot(1,3,2);
myMap = colormap('lines');
hold on
for k = 1:4
    hp = plot(k,gMeans(k),'o','Color',myMap(k,:));
    set(hp,'MarkerSize',8,'MarkerFaceColor',myMap(k,:));
end
hl = line([1,nGroups],[mean(gMeans),mean(gMeans)]);
set(hl,'LineWidth',1.5,'LineStyle','--', 'Color','k');
axis(ax);
xlabel('Group')
ylabel('y')
title(['\rm Variance \bf between \rm groups: ' num2str(varB,2)]);
text(2,ax(4)-1,['\bf ICC \rm= ' num2str(ICCgroup,2)]);

%% Create a new table of just the means and the treatment
Tm = table();
Tm.y = gMeans;
Tm.rx = [0;0;1;1];

%% Fit models

% Fit 3 models:
% 1: wrong one, treating all observations as independent:
lme1 = fitlm(T,'y ~ 1 + rx');
% 2: wrong one, fitting only the means:
lme2 = fitlm(Tm,'y ~ 1 + rx');
% right one, group-specific y-intercept (see fig. 3 of Aarts et al. 2014)
lme3 = fitlme(T,'y ~ 1 + rx + (1|grpID)');

%% Bundle models for reporting back
allModels = {lme1,lme2,lme3};
allModelNames = {'lme1','lme2','lme3'};
nModels = length(allModels);

% Which model to plot and do regression diagnostics on:
model = lme3;
mText = '3';

%% Read out the random effects:
% Now we need to read out the individual intercepts from the model
betaFit = fixedEffects(lme3);          % give us the fixed effects (slope & intercept)
[~,~,STATS] = randomEffects(lme3);      % Compute the random-effects statistics
STATS.Level = nominal(STATS.Level);     % declare as a nominal variable
% Add the random effects to fixed effect to get the y-intercept for each
% group
beta0Grps = zeros(nGroups,1);
beta0Grps(1) = betaFit(1) + STATS.Estimate(STATS.Level=='1');
beta0Grps(2) = betaFit(1) + STATS.Estimate(STATS.Level=='2');
beta0Grps(3) = betaFit(1) + STATS.Estimate(STATS.Level=='3');
beta0Grps(4) = betaFit(1) + STATS.Estimate(STATS.Level=='4');

% fit coefficients for model 1 (same as fixed effects for lme2):
%betaFit1 = lme1.Coefficients.Estimate;
%% Plot models
% plot replicating figure 3 of Aarts et al. 2014
subplot(1,3,3);
x = jitter(T.rx,0.1);
gscatter(x,T.y,T.grpID);
hold on
xlabel('Treatment');
ylabel('y')

% plot the fits:
% 1. fixed effects:
h0 = line([0,1],[betaFit(1), betaFit(1) + betaFit(2)]);
set(h0,'LineWidth',2,'Color','k');
% 2. random effects
for k = 1:nGroups
    h1 = line([0,1],[beta0Grps(k),beta0Grps(k) + betaFit(2)]);
    set(h1,'LineWidth',2,'LineStyle','--','Color',myMap(k,:));
end

lgd = legend('1','2','3','4');
legend('Location','north');
lgd.Title.String = 'Group';
title('Model fit')

%% Print out t-stat & p-value for the Rx param for each model:
for k = 1:nModels
    fprintf(1,'Model: %d, R^2: %.2f, Beta1: %.2f, tStat: %.2f, pVal: %.4f\n',...
        k, allModels{k}.Rsquared.Adjusted,...
        allModels{k}.Coefficients.Estimate(2),...
        allModels{k}.Coefficients.tStat(2),...
        allModels{k}.Coefficients.pValue(2));
end