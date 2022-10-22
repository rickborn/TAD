% incomeByHeightDemo.m
%
% RTB trying to emulate Gelman & Nolan demo (pp. 49-51)
% "Teaching Statistics: A Bag of Tricks"

%% Load and plot data

% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\TAD Fall 2017\TAD2017\Last class'
fileName = 'IncomeByHgtData.xlsx';
ds = readtable(fileName);

figure;
% jitter x-values for better viewing:
jFactor = 0.5;
jitterX = (rand(length(ds.Hgt),1) .* jFactor) - (jFactor/2);

% We want to use actual height for regression, but jittered height for
% plotting:
jHgt = ds.Hgt+jitterX;
hp = plot(jHgt, ds.Income, 'b.');
set(hp,'MarkerSize',10);
hold on
% draw the least-squares regression line:
hl = lsline;
set(hl,'Color','k','LineWidth',2);
xlabel('Height (inches)');
ylabel('Income (thousands of $)');
x = xlim;

set(gca,'FontSize',12);
r = corrcoef(ds.Hgt,ds.Income);
text(57,110,['r = ' num2str(r(2,1),2)],'FontSize',14);

%% Do simple regression
modelspec1 = 'Income ~ Hgt';
mdl1 = fitglm(ds,modelspec1,'Distribution','normal');
[b1,dev1,stats1] = glmfit(ds.Hgt,ds.Income);

%% Suppose we didn't know about regression.
% How might we develop an intuition about whether the slope is not 0?

% resample pairs from data set (with replacement) and replot a line:
nBoot = 1000;
nPairs = length(ds.Hgt);
allSlopes = zeros(nBoot,1);

for k = 1:nBoot
    % random sample with replacement
    randNdx = unidrnd(nPairs,nPairs,1);
    % calculate regression parameters:
    [allSlopes(k),b] = l_regression(ds.Hgt(randNdx),ds.Income(randNdx));
    % calculate the regression line:
    y = b + (allSlopes(k) .* x);
    hl = line(x,y,'Color','r','LineWidth',0.3);
end
    
figure
histogram(allSlopes);

% Motivated by Gelman & Carlin 2014:
% Beyond Power Calculations: Assessing Type S (Sign) and Type M (Magnitude)Errors
%
% Abstract: Statistical power analysis provides the conventional approach to
% assess error rates when designing a research study. However, power
% analysis is flawed in that a narrow emphasis on statistical significance
% is placed as the primary focus of study design. In noisy, small-sample
% settings, statistically significant results can often be misleading. To
% help researchers address this problem in the context of their own
% studies, we recommend design calculations in which (a) the probability of
% an estimate being in the wrong direction (Type S [sign] error) and (b)
% the factor by which the magnitude of an effect might be overestimated
% (Type M [magnitude] error or exaggeration ratio) are estimated. We
% illustrate with examples from recent published research and discuss the
% largest challenge in a design calculation: coming up with reasonable
% estimates of plausible effect sizes based on external information.
%
% This would be the probability of a Type S error:
pWrongDirection = sum(allSlopes < 0) / nBoot;

%% Lurking variable: sex
% men taller than women and make more $$$

figure;
ds.Male = logical(ds.Male);
h1 = plot(jHgt(ds.Male),ds.Income(ds.Male),'k+');
h2 = plot(jHgt(~ds.Male),ds.Income(~ds.Male),'ro');
legend([h1,h2],{'Male','Female'},'Location','NorthWest');

%% Add sex to regression
modelspec2 = 'Income ~ Hgt + Male';
mdl2 = fitglm(ds,modelspec2,'Distribution','normal');
[b2,dev2,stats2] = glmfit([ds.Hgt,ds.Male],ds.Income);

%% Show separate regression lines for men and women

% re-plot original data:
figure
hp = plot(jHgt, ds.Income, 'b.');
set(hp,'MarkerSize',10);
hold on
h1 = plot(jHgt(ds.Male),ds.Income(ds.Male),'ks');
h2 = plot(jHgt(~ds.Male),ds.Income(~ds.Male),'ro');
%legend([h1,h2],{'Male','Female'},'Location','NorthWest');
xlabel('Height (inches)');
ylabel('Income (thousands of $)');

set(gca,'FontSize',12);
%r = corrcoef(ds.Hgt(~ds.Male),ds.Income(~ds.Male));
%text(57,110,['r = ' num2str(r(2,1),2)],'FontSize',14);

%% Get predictions and CIs for men:

xVals = (min(ds.Hgt):max(ds.Hgt))';
men = ones(size(xVals));

[yMen,yMenCI] = predict(mdl2,[xVals,men]);
% Plot 'em
plot(xVals,yMen,'k-','LineWidth',2);
plot(xVals,yMenCI(:,1),'k--');
plot(xVals,yMenCI(:,2),'k--');

% Get predictions and CIs for women:
[yWomen,yWomenCI] = predict(mdl2,[xVals,~men]);
% Plot 'em
plot(xVals,yWomen,'r-','LineWidth',2);
plot(xVals,yWomenCI(:,1),'r--');
plot(xVals,yWomenCI(:,2),'r--');

legend([h1,h2],{'Male','Female'},'Location','NorthWest');

%% Bonus: Included-variable bias:

% Inspired by: 
% Jung, J., Corbett-Davies, S., Shroff, R. and Goel, S., 2018. Omitted and
% included variable bias in tests for disparate impact. arXiv preprint
% arXiv:1809.05651.

% If we include a "proxy" variable--in this case, "voice register"--what
% happens to our coefficient for sex?

fileName = 'IncomeByHgtData2.xlsx';
ds2 = readtable(fileName);

% Recode to make 'S' and 'A' into 'hi' and 'T' and 'B' into 'lo'
ds2.Voice = categorical(ds2.Voice);
ds2.Voice(ds2.Voice == 'S' | ds2.Voice == 'A') = 'high';
ds2.Voice(ds2.Voice == 'T' | ds2.Voice == 'B') = 'low';
ds2.Voice = categorical(ds2.Voice);

% Now fit the model:
modelspec3 = 'Income ~ Hgt + Male + Voice';
mdl3 = fitglm(ds2,modelspec3,'Distribution','normal');
%[b3,dev3,stats3] = glmfit([ds2.Hgt,ds2.Male,ds2.Voice],ds2.Income);

%% Create a data set with vocal ranges of 'S,A,T,B'

% Add another column that correlates with sex: voice = S,A,T,B
t1 = randi([3,4],100,1);    % men are 'tenor' or 'bass'
t2 = randi([1,2],100,1);    % women are 'soprano' or 'alto'
ds.Voice(ds.Male == 1) = t1;
ds.Voice(ds.Male == 0) = t2;
ds.Voice = categorical(ds.Voice,[1,2,3,4],{'S','A','T','B'});
writetable(ds,'incomeByHgtData2.xlsx');