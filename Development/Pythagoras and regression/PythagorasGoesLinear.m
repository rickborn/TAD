% PythagorasGoesLinear.m
%
% from p. 160-1 of Gelman & Nolan's "Teaching statistics: A bag of tricks"
%
% From p. 161: "What is the point of this example? At one level, it shows
% the power of multiple regression--even when the data come from an
% entirely different model, the regression can fit well. There is also a
% cautionary message, showing the limitations of any purely data-analytic
% method for finding true underlying relations. As we tell the students, if
% Pythagoras knew about multiple regression, he might never have discovered
% his famous theorem!"
%
% IMHO this is also a good lesson on the importance of regression
% diagnostics: the residuals-vs-fitted and qqplots clearly reveal that our
% assumptions for linear regression are not met.
%
% RTB wrote it, 26 Sept. 2017

%% Generate a Pythagorean data set

nSamp = 50;
x1 = unidrnd(10,nSamp,1);
x2 = unidrnd(20,nSamp,1);
% x1 = rand(nSamp,1) .* 10;
% x2 = rand(nSamp,1) .* 20;
y = sqrt(x1.^2 + x2.^2);

% consider adding noise

% or load in existing file:
% load pyData.mat

%% Model it with linear regression
const = ones(length(y),1);
[betaFit,betaCI,resid,residInt,stats] = regress(y,[const,x1,x2]);

% 'stats' contains the R2 statistic, the F-statistic and its p-value,
% and an estimate of the error variance.

%% Plot residuals vs. fitted

% model predictions:
yPred = betaFit(1) + betaFit(2).*x1 + betaFit(3).*x2;

main = figure('position',[50 50 500 850]);
subplot(3,1,1)
plot(yPred,resid,'ko');
xlabel('Predicted values');
ylabel('Residuals');
tStr = sprintf('Linear prediction R^2 = %.3f', stats(1));
title(tStr);

%% Q-Q Plot of residuals

subplot(3,1,2)
qqplot(resid);
title('QQ Plot of Residuals vs. Standard Normal');

%% Finally, plot the data and the fit:

subplot(3,1,3)
hp=plot3(x1,x2,y,'ro');
set(hp,'MarkerFaceColor','r');
hold on
ax = axis;

xx = linspace(ax(1),ax(2),length(x1));
yy = linspace(ax(3),ax(4),length(x2));
[X,Y] = meshgrid(xx,yy);
Z = betaFit(1) + betaFit(2).*X + betaFit(3).*Y;
surf(X,Y,Z);
xlabel('x1');
ylabel('x2');
zlabel('y');
title(tStr);

%% But how does our regression model do beyond the original data?

% create an anonymous function that calculates the hypotenuse:
hypot = @(a,b) sqrt(a.^2 + b.^2);

% large triangle sides:
%newX1 = 500; newX2 = 750;

% tiny triangle sides:
newX1 = 0.1; newX2 = 0.3;

% regression model prediction for the hypotenuse:
regY = betaFit(1) + betaFit(2).*newX1 + betaFit(3).*newX2;

% prediction by the correct model:
pythagY = hypot(newX1,newX2);

perCentError = (abs(pythagY - regY) / pythagY) * 100