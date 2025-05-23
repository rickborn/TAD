% WRSTtutorial.m: A brief tutorial on computing ranks with MATLAB
%
% https://www.stat.auckland.ac.nz/~wild/ChanceEnc/Ch10.wilcoxon.pdf
%
% RTB wrote it, 30 October 2020, unexpected blizzard!

% cd 'C:\usr\rick\doc\Committees\PIN\PIN Director\Courses\Stats\TAD\TAD Code\Development\WRST Tutorial'
%% Example 1

% Useful for labeling later:
NA = 1;
C = 0;

% Two different measures to illustrate basics and then ranking when there
% are ties in the data:
% 1. MSCE (mean sister chromatid exchage)
% 2. Dispersion

msceFlag = 0;
if msceFlag
    % We shall compare the groups labeled “Native American” and “Caucasian”
    % with respect to the variable MSCE (mean sister chromatid exchange).
    natAm = [8.5,9.48,8.65,8.16,8.83,7.76,8.63];
    cauc = [8.27,8.2,8.25,8.14,9,8.1,7.2,8.32,7.7];
    xStr = 'MSCE';
else
    natAm = [0.61,0.63,0.45,0.50,0.75,0.93,0.85];
    cauc = [0.56,0.58,0.45,0.44,0.58,0.79,0.65,0.52,0.53];
    xStr = 'Dispersion';
end

% Create an indicator variable: 0 = Caucasian, 1 = Native American
iNatAm = ones(size(natAm));
iCauc = zeros(size(cauc));

% Create a single data vector for the values and the indicators:
allData = [natAm,cauc];
allInd = [iNatAm,iCauc];

%% Plot the data:
figure
plot(allData,allInd,'bo');

% fix axes
ax = axis;
axis([ax(1),ax(2),-1,2]);

% make it pretty:
yVals = get(gca,'YTick');
yLabels = get(gca,'YTickLabel');

for k = 1:length(yVals)
    if yVals(k) == NA
        yLabels{k} = 'Native Am.';
    elseif yVals(k) == C
        yLabels{k} = 'Caucasian';
    else
        yLabels{k} = ' ';
    end
end
set(gca,'YTickLabel',yLabels);
xlabel(xStr);
title(['Comparing ',xStr, ' measurements']);

%% Rank the pooled data:

[rnk,tiedJ] = tiedrank(allData);

% compute the rank sums for each group:
wNA = sum(rnk(allInd == NA));
wC = sum(rnk(allInd == C));

%% Significance testing

% We can either compare of test statistic (wNA, or the 'rank sum') to a
% known null, or we can model H0 directly using a permutation test. The
% simplest way to do this is just to shuffle the labels.

nObs = length(allData);
nPerm = 10000;
allW = zeros(nPerm,1);

for k = 1:nPerm
    shuffledLabels = allInd(randperm(nObs));
    allW(k) = sum(rnk(shuffledLabels == NA));
end

figure
histogram(allW);
hold on
ax = axis;
hl = line([wNA,wNA],[ax(3),ax(4)],'LineWidth',2,'Color','r');
xlabel('Rank sum under H0');
ylabel('#');

%% Compute a 2-tailed p-value

if wNA >= median(allW)
    pValPermTest = 2 * (sum(allW >= wNA) / nPerm);
else
    pValPermTest = 2 * (sum(allW < wNA) / nPerm);
end
display(pValPermTest);

%% Compare with the Wilcoxon Rank Sum Test

pValWRST = ranksum(natAm,cauc);
display(pValWRST);

% Moral #1: We get the same answer using a permuation test as we do with
% the "classical" Wilcoxon Rank Sum Test

%% Same question, but with the linear model on ranks
% This is only exact for N > 11.

% As for all 2-sample tests with the LM, it is the beta coefficient for
% the slope of the regression line that is relevant for the comparison.
[b,dev,stats] = glmfit(allInd,rnk);
pLM = stats.p(2);
display(pLM);

%% Plot results of LM

figure
plot(allInd,rnk,'bo');
hold on

% fix axes
ax = axis;
axis([-1,2,ax(3),ax(4)]);

% make it pretty:
xVals = get(gca,'XTick');
xLabels = get(gca,'XTickLabel');

for k = 1:length(xVals)
    if xVals(k) == NA
        xLabels{k} = 'Native Am.';
    elseif xVals(k) == C
        xLabels{k} = 'Caucasian';
    else
        xLabels{k} = ' ';
    end
end
set(gca,'XTickLabel',xLabels);
ylabel([xStr, ' Rank']);
title(['Comparing ',xStr, ' measurements']);

% Now plot regression line
X = [0,1];
Y = b(1) + (b(2) .* X);

% Not clear how one should interpret a standard error for ranks.
yEBhi = (b(1)+stats.se(1)) + ((b(2)+stats.se(2)).*X);
yEBlo = (b(1)-stats.se(1)) + ((b(2)-stats.se(2)).*X);
errorbar(X,Y,yEBlo,yEBhi,'r.');
plot(X,Y,'r-');

% Moral #2: We get very nearly the same p-value as the WRST using the LM on
% ranked data.