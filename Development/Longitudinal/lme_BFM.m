function [p,lme1,lme2,cv1,cv2] = lme_BFM(fileName)

% lme_BFM.m: linear mixed-effects model of bodyfat/menarche data
%
% [p,lme1,lme2] = lme_BFM(fileName)
%
% e.g. [p,M1,M2,cv] = lme_BFM('bodyfat.xlsx')
%
% Inputs:
% - fileName: Excel file containing data
%
% Outputs:
% - p: p-value for the critical test of our H0
% - lme1: mixed-effects model fit by Garrett Fitzmaurice
% - lme2: mixed-effects model fit in Naumova et al. 2001
% - covar1: the covariance matrix for the random effects for model #1
% - covar2: the covariance matrix for the random effects for model #1
% 
% This example is from:
%
% Naumova, E.N., Must, A. and Laird, N.M., 2001. Tutorial in biostatistics:
% evaluating the impact of ‘critical periods’ in longitudinal studies of
% growth using piecewise mixed effects models. International journal of
% epidemiology, 30(6), pp.1332-1341.
%
% It was also covered by Garrett Fitzmaurice in lecture #25 of Brian
% Healy's "Biostatistics Certificate" course. See my:
%
% 'Biostatistics-Lecture_Notes-RTB.docx', sections 25.7 to 25.10 for
% details
%
% RTB wrote it, June 3, 2025; bending sails (jib + mainsail) with David Cardozo

% Scientific question: Does the rate of body fat accumulation in girls
% change at menarche?
%
% Experiment:
% 1. Prospective study on body fat accretion in a cohort of 162 girls from the
% MIT Growth and Development Study
% 2. At start of study, all the girls were pre-menarchal and non-obese.
% 3. All girls were followed over time according to a schedule of annual
% measurements until four years afer menarche.
% 4. The final measurement was scheduled on the fourth anniversary of their
% reported date of menarche.
% 5. At each examination, a measure of body fatness was obtained based on
% bioelectric impedance analysis.
%
% Data: (http://alecri.github.io/downloads/data/bodyfat.RData)
% I converted it to an Excel file using my 'ds2xls.R' script.
% A total of 1049 measurements from 162 girls. Each row corresponds to one
% measurement of % body fat:
%
% Variables:
% 1. 'id' - individual subject id number (1 to 162)
% 2. 'agec' - girl's age in years at the time of the measurement
% 3. 'agem' - girl's age at menarche
% 4. 'time' - girl's age relative to menarche (col. 2 - col. 3)
% 5. 'pbf' - percent body fat
% 6. 'occasion' - sequential # of the visit for that girl

% Analysis plan:
% We will fit a piece-wise linear model to the data, with the 'knot' or
% 'hinge' at the time of menarche ('time' = 0). This means that there will
% be a slope before menarche and a slope after menarche. We'll also allow
% for random effects for both slopes and for the y-intercept.

%% Defaults
if nargin < 1, fileName = 'bodyfat.xlsx'; end

%% Load the data
ds = readtable(fileName);
%varNames = ds.Properties.VariableNames;
%disp(varNames);

%% Gather some information
allIDs = unique(ds.id);
nGirls = length(allIDs);
% The 2 examples used in lecture 25 by GF were #16 and #4
exampleGirls = [16,4];

%% Distribution of numbers of measurments for different girls
% illustrates that the data are highly unbalanced
allMs = zeros(nGirls,1);
for k = 1:nGirls
    allMs(k) = max(ds.occasion(ds.id == k));
end
figure('Name','Unbalanced data');
histogram(allMs);
xlabel('Total number of measurements');
ylabel('Number of subjects');
title('Compare with Table 1 of Naumova et al. 2001');

%% Spaghetti plot of all the data

% Label the figure with the name of the Excel file:
figName = fileName(1:strfind(fileName,'.')-1);
figName = [figName, ' data analysis'];
f1 = figure('Name',figName,'Position',[50 450 1200 400]);
myMap = colormap('lines');
subplot(1,2,1)
% plot each girl as a connected line
% It doesn't look like 'gscatter' will connect with lines, so we'll have to
% do it with a 'for' loop
for k = 1:nGirls
    h = plot(ds.time(ds.id == k),ds.pbf(ds.id == k),'Color',myMap(k,:),...
        'Marker','.','LineStyle','-');
    hold on
    if any(k == exampleGirls)
        set(h,'Marker','^','MarkerFaceColor',myMap(k,:),'LineWidth',1.5);
    end
end

%% Fit data with a smoothing spline
% Set up fittype and options.
ft = fittype( 'smoothingspline' );
opts = fitoptions( 'Method', 'SmoothingSpline' );
opts.SmoothingParam = 0.0302091785132346;

% Fit model to data.
[f,~] = fit(ds.time,ds.pbf,ft,opts);

% Plot fit.
hSpline = plot(f);
set(hSpline,'Color',myMap(3,:),'LineWidth',2);
legend off
ax = axis;
axis([min(ds.time),max(ds.time),ax(3),ax(4)]);
xlabel('Time relative to menarche (yrs)');
ylabel('Percent body fat');
title('Spaghetti plot');

% draw a dashed line at 'time' = 0
line([0,0],[ax(3),ax(4)],'Color','k','LineStyle','--');

%% Set up the piece-wise model

% 1. Let t_ij denote the time of the jth measurement on the ith subject
% (t_ij = 0 at menarche)
%
% 2. To create the “broken stick” model (i.e. piece-wise linear) with a
% “knot” at menarche, we create a new variable: (tij)+ such that:
%
%   (t_ij)+ =  t_ij  if  t_ij  > 0
% 	(t_ij)+ =  0     if  t_ij  <= 0
%
% 3. Think of this new variable as “the positive part of time” (since menarche)
%
% 4. Then we fit the following model:
%     Y_ij = B1 + B2t_ij + B3(t_ij)+
%
% B2 is the slope prior to menarche
% (B2 + B3) is the slope after menarche
%
% Our scientific question is answered by testing B3 = 0

% Create a new variable call 'posTime':
ds.posTime = ds.time;
ds.posTime(ds.posTime <= 0) = 0;

%% Fit the linear mixed-effects model

% The first way I tried it was wrong. When the random effects are specified
% as separate parentheticals separated by '+' signs, the model assumes the
% random effects are uncorrelated.
% Wrong way:
% lme1 = fitlme(ds,'pbf ~ time + posTime + (1|id) + (time-1|id) + (posTime-1|id)');
% Right way:
lme1 = fitlme(ds,'pbf ~ time + posTime + (1 + time + posTime|id)');
% The main p-value of interest is for the fixed of 'posTime'
p = lme1.Coefficients.pValue(3);

%% Extract the covariance matrix for the random effects:
psi = covarianceParameters(lme1);
cv1 = psi{1};
% The main diagonal gives the variances of the 3 random effects
% Off-diagonals are covariances between pairs of random effects

% NOTE: It's easy to convert the covariance matrix to give 1) the standard
% deviations (sqrt of the variances on the main diagonal and 2) the
% correlation matrix:
% [sd,corr] = cov2corr(covar);

%% Second model, as instantiated in Naumova et al. 2001
% See p. 1336. They used a slightly different way of coding time w/r/t
% menarche. Essentially they created a 'pre-time' and a 'post-time' so that
% the 2 slopes would correspond to the pre- and post-menarchal slopes
% directly. In our first approach, the post-menarchal slope was (B2+B3).

% indicator variable that is '1' for ds.time < 0:
indPre = ds.time < 0;
ds.preTime = ds.time .* indPre;
ds.postTime = ds.time .* ~indPre;

lme2 = fitlme(ds,'pbf ~ preTime + postTime + (1+preTime+postTime|id)');

% extract covariance matrix for the random effects:
psi = covarianceParameters(lme2);
cv2 = psi{1};
% results are identical:
% For model 1, b2 + b3 = 2.4640 (post-menarchal slope)
% For model 2, b3 = 2.4640 (post-menarchal slope)

% NOTE: I think model #1 is more elegant, since its 'b3' represents the
% *difference* in the pre-post slopes, and therefore allows us to directly
% test our H0 that there is no difference after menarche. We can
% definitively reject H0 (p = 9 x 10^-27).

%% Plot the fixed effects fit

% x data for 2 halves of fit
xNeg = -6:1:0;
xPos = 0:1:4;

bFix = fixedEffects(lme1); % fixed effect coefficients

yNeg = bFix(1) + (xNeg .* bFix(2));
yPos = bFix(1) + (xPos .* (bFix(2)+bFix(3)));
hFit = plot(xNeg,yNeg,'Color','k','LineStyle','-','LineWidth',2);
plot(xPos,yPos,'Color','k','LineStyle','-','LineWidth',2);

legend([hSpline,hFit],{'Smoothing spline','Model fit'},'Location','Northwest');

%% Plots for example subjects

% first figure
figure(f1); subplot(1,2,2);

% re-plot the fit to the population mean:
plot(xNeg,yNeg,'Color','k','LineStyle','-','LineWidth',2);
hold on
plot(xPos,yPos,'Color','k','LineStyle','-','LineWidth',2);

for m = 1:length(exampleGirls)
    thisGirl = exampleGirls(m);
    T = ds(ds.id == thisGirl,:);
    
    % find appropriate range of times for this girl:
    xNeg = floor(min(T.time)):1:0;
    xPos = 0:1:ceil(max(T.time));
    
    % plot her raw data with larger symbols
    plot(T.time,T.pbf,...
        'Color',myMap(thisGirl,:),'Marker','o','LineStyle','none',...
        'MarkerFaceColor',myMap(thisGirl,:));
    
    % fit ML line to just her data:
    bML1 = glmfit([T.time,T.posTime],T.pbf);
    % calculate and plot the lines:
    yNeg = bML1(1) + (xNeg .* bML1(2));
    yPos = bML1(1) + (xPos .* (bML1(2)+bML1(3)));
    h1 = plot(xNeg,yNeg,'Color',myMap(thisGirl,:),'LineStyle','-','LineWidth',1.5);
    plot(xPos,yPos,'Color',myMap(thisGirl,:),'LineStyle','-','LineWidth',1.5);
    
    % Plot the BLUP for this girl:
    % First find the random effects for this girl
    [~,~,stats] = randomEffects(lme1);
    % Find appropriate rows in 'stats':
    L = strcmp(num2str(thisGirl),stats.Level);
    bRE = stats.Estimate(L);
    % add to the fixed effects:
    bBLUP = bFix + bRE;
    yNeg = bBLUP(1) + (xNeg .* bBLUP(2));
    yPos = bBLUP(1) + (xPos .* (bBLUP(2)+bBLUP(3)));
    h2 = plot(xNeg,yNeg,'Color',myMap(thisGirl,:),'LineStyle','--','LineWidth',1.5);
    plot(xPos,yPos,'Color',myMap(thisGirl,:),'LineStyle','--','LineWidth',1.5);
    
end

% Make it pretty
axis([min(ds.time),max(ds.time),ax(3),ax(4)]);
% draw a dashed line at 'time' = 0
line([0,0],[ax(3),ax(4)],'Color','k','LineStyle','--');
xlabel('Time relative to menarche (yrs)');
ylabel('Percent body fat');
title('Shrinkage is real!');
legend([h1,h2],{'ML fit','BLUP'});

%% Regression diagnostics

model = lme1;
mText = '1';

figure('Name','Regression Diagnostics','Position',[50 50 1200 400]);
subplot(1,2,1);
plot(model.Fitted,model.Residuals.Raw,'ko');
hold on;
ax = axis;
line([ax(1),ax(2)],[0,0]);
xlabel('Linear Predictor'); ylabel('Residual');
title(['Residuals vs. Fitted for Model #' mText]);

subplot(1,2,2);
qqplot(model.Residuals.Raw);
title(['Q-Q Plot of Model #' mText ' Residuals']);