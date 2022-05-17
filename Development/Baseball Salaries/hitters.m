% hitters.m
%
% Builds intuition on regression trees for machine learning
%
% RTB wrote it, summer 2021 for use in code4fun ML course with Andrei
% Grigoriev.

%% Data is from "Introduction to Statistical Learning"
% see "Hitters" in the file ISLR.pdf
% 
% This data was originally part of the 'R' data base for ISLR:
% see my 'ch6.R' for example of model selection:
% ### Model selection with 'R'
% install.packages("ISLR")
% library(ISLR)
% summary(Hitters)
% 
% # There are some missing values, so we will remove them:
% Hitters = na.omit(Hitters)
% with(Hitters, sum(is.na(Salary)))
% 
% # Write out data as a CSV file for importing into MATLAB
% write.csv(Hitters,"Hitters.csv",row.names=FALSE)

ds = readtable("Hitters.csv");

% There are a couple of players who had very few at-bats; we'll eliminate
% them:
ds = ds(ds.AtBat >= 100,:);

figure
% Since "Years in league" is discrete, might help to jitter a bit:
jX = jitter(ds.Years,0.3);
scatter(jX,ds.Hits,25,ds.Salary,'filled');
colormap(gca,'jet')
c = colorbar;
c.Label.String = 'Salary (thousands of $)';
xlabel('Years in league');
ylabel('Hits');

%% Calculate RSS as a function of different values for Years:

xCrit = 1.5:15.5;
yRSS = zeros(size(xCrit));

for k = 1:length(xCrit)
    % partition the data: x < xCrit vs. x >= xCrit
    xLTsel = ds.Years < xCrit(k);
    
    yRSS(k) = sum((ds.Salary(xLTsel) - mean(ds.Salary(xLTsel))) .^ 2) + ...
        sum((ds.Salary(~xLTsel) - mean(ds.Salary(~xLTsel))) .^ 2);
end

% find the criterion that minimizes RSS:
myCrit = xCrit(yRSS == min(yRSS));

figure
plot(xCrit,yRSS,'ko-');
hold on
plot(myCrit, min(yRSS),'ro','MarkerFaceColor','r');

xlabel('Criterion (Years in league)');
ylabel('Residual sum of squares');

%% Plot this on figure 1

figure(1);
ax = axis;

hl = line([myCrit,myCrit], [ax(3),ax(4)], 'Color','k','LineStyle','--');

%% Make a new plot of just salary vs. years in league

figure
%plot(ds.Years, ds.Salary, 'bo');
plot(jX, ds.Salary, 'bo');
xlabel('Years in league');
ylabel('Salary (thousands of $)');

ax = axis;
line([myCrit,myCrit], [ax(3),ax(4)], 'Color','k','LineStyle','--');

% now plot a line for the mean salary on each side:
xLTsel = ds.Years < myCrit;
y1 = mean(ds.Salary(xLTsel));
y2 = mean(ds.Salary(~xLTsel));

line([ax(1),myCrit], [y1,y1], 'Color','r','LineStyle','-');
line([myCrit,ax(2)], [y2,y2], 'Color','r','LineStyle','-');

