% etALS2Sample_ng_withAnswers.m
%
% RTB wrote it, working through Ch. 15 of E&T, Dec. 2017
% RTB converted to class exercise, 16 Sept. 2017; gray, damp morning;
% listening to Casals recording of Bach Cello Suites
%
% E&T Wisdom: "When there is something to permute, it is a good idea to do
% so . . . ."
%
% Source (E&T): Efron B, Tibshirani RJ (1993) "An introduction to the
% bootstrap" New York: Chapman and Hall.
%
% formerly named: etMouse2SamplePermutationTest_ng_withAnswers.m

% Concepts covered:
% 1. Using 'boxplot' to visualize grouped data
% 2. Jittering x-values to allow visualization of raw data
% 3. Permutation test for hypothesis testing
% 4. Comparison with 2-sample t-test
% 5. Bootstrap test for hypothesis testing
% 6. Permutation vs. bootstrap

% What to do: Login to learning catalytics and join the session for the
% module entitled "etALS2sample". You will answer a series of questions
% based on the guided programming below. Each section begins with a '%%'.
% Read through the comments and follow the instructions provided. In some
% cases you will be asked to answer a question, clearly indicated by
% 'QUESTION'. In other cases, you be asked to supply missing code,
% indicated by 'TODO'. The corresponding question in learning catalytics
% will be indicated in parentheses (e.g. Q1). If there is no 'Q#'
% accompanying a 'QUESTION' just type your answer into this script and
% discuss it with your team. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).

%% Read in and plot the data
myAlpha = 0.05;
nPerm = 10000;

% A pilot study is performed in which sixteen SOD1 mice (a model of
% Amyotrophic Lateral Sclerosis, also known as Lou Gehrig's disease) were
% randomly assigned to a treatment group (riluzole) or a control group.
% Values are their survival times, in days, from the beginning of
% treatment. Does the treatment improve survival time?

rxGrp = [94,197,16,38,99,141,23];           % Riluzole
ctrlGrp = [52,104,146,10,51,30,40,27,46];   % Placebo

nRx = length(rxGrp);
nCtrl = length(ctrlGrp);

% TODO: Calculate the mean difference in survival time between treatment
% and control.
muDiffHat = mean(rxGrp) - mean(ctrlGrp);
% QUESTION (Q1): What is the mean difference?

%% Plot the data
grpLabels = {'Rx';'Rx';'Rx';'Rx';'Rx';'Rx';'Rx'; ...
    'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl';'Ctrl'};
boxplot([rxGrp';ctrlGrp'],grpLabels);
hold on
xlabel('Treatment Group');
ylabel('Survival Time (days)');
tStr = sprintf('Mean Diff = %.2f days',muDiffHat);
text(1.5,180,tStr);

% NOTE: The way that the text string 'tStr' was created above is the way a
% 'C' programmer would do it, not a particularly MATLAB-y way. 
% QUESTION(Q2): How could you create the same text string without using 'sprintf'?
tStr = ['Mean Diff = ' num2str(muDiffHat) ' days'];

%% Superimpose raw data

% To do this, we need to know the numeric values associated with the two
% groups. This can be obtained using the 'gca' command
xVals = get(gca,'XTick');

% To prevent nearby values from superimposing, we can jitter the data along
% the x-axis. If we don't jitter too much, it will still be clear to which
% group the values belong.
jFactor = (xVals(2) - xVals(1)) / 25;
plot((ones(1,nRx) .* xVals(1)) + (randn(1,nRx).*jFactor),rxGrp,'ko');
plot((ones(1,nCtrl) .* xVals(2)) + (randn(1,nCtrl).*jFactor),ctrlGrp,'ko');

%% Use classical 2-sample t-test and a non-parametric test

% TODO: We are only interested in the alternative hypothesis that our drug
% prolongs survival. Perform the appropriate one-tailed test using both a
% parametric 2-sample t-test and the corresponding non-parametric test.
[hT,pValTtest,ci,statsT] = ttest2(rxGrp,ctrlGrp,'Tail','right');
% QUESTION (Q3): What is the p-value for the 2-sample t-test?

[pValWRST,hRS,statsRS] = ranksum(rxGrp,ctrlGrp,'Tail','right');
% QUESTION (Q4): What is the p-value for the non-parametric analog of the t-test?

%% Permutation test for this difference

% Key step: Our null hypothesis is that the two data sets came from the
% SAME underlying probability distribution. Put another way, "if H0 is true,
% then any of our observations (survival times) for any of the mice could
% have come equally well from either of the treatment groups." (E&T p.
% 205). So we just pool the data, shuffle it, then take the first nRx
% values to be our Rx group (under H0) and the remaining to be the ctrl
% group. Note that this works for any statistic. For example, compare
% looking at the difference with 'mean' vs. 'trimmean' (the trimmed mean).

% This is the key conceptual step: pool the data
H0data = [rxGrp, ctrlGrp];

nTotal = length(H0data);
muDiffPerm = zeros(nPerm,1);
rng default
for k = 1:nPerm
    % TODO: 'shuffle' the data (HINT: same as sampling without replacement)
    shuffledData = H0data(randperm(nTotal));
    %muDiffPerm(k) = trimmean(shuffledData(1:nRx),25) - trimmean(shuffledData(nRx+1:end),25);
    muDiffPerm(k) = mean(shuffledData(1:nRx)) - mean(shuffledData(nRx+1:end));
end
figure, hist(muDiffPerm,30);
hold on;
xlabel('Permuted differences: Rx - Ctrl'); ylabel('#');

% Draw a line to show our experimentally determined survival difference
ax = axis;
line([muDiffHat,muDiffHat],[ax(3),ax(4)],'Color','y');

% QUESTION (Q5): Stated in terms of 'muDiffPerm', what is our null value?

%% Calculate a one-tailed p-value

% This is just the # of permuted values that are greater-than-or-equal-to
% the one we observed experimentally divided by the total. E&T refer to
% this value as the "achieved significance level" or ASL.
pValPerm = sum(muDiffPerm >= muDiffHat) / nPerm;

% QUESTION (Q6): What is your p-value?

% add this to the figure title
tStr = sprintf('p-value = %0.3f',pValPerm);
title(tStr);

%% We can also use the bootstrap to do the hypothesis test

% Here the strategy is to resample WITH REPLACEMENT from each group
% separately, then look at how the difference behaves. In essence, we are
% creating the sampling distribution of the mean difference.

nBoot = 10000;
muDiffBoot = zeros(nBoot,1);

rng default
for k = 1:nBoot
    % TODO: Calculate a bootstrap replication of the mean difference
    muDiffBoot(k) = mean(rxGrp(unidrnd(nRx,nRx,1))) - ...
        mean(ctrlGrp(unidrnd(nCtrl,nCtrl,1)));
end

% TODO: Calculate the standard error of our mean difference:
seDiffBoot = std(muDiffBoot);
% QUESTION (Q7): What is the value of the standard error of the mean
% difference?

% Plot the sampling distribution:
figure, hist(muDiffBoot,30);
hold on;
xlabel('Bootstrap difference: Rx - Ctrl'); ylabel('#');
ax = axis;
line([muDiffHat,muDiffHat],[ax(3),ax(4)],'Color','y');
line([0,0],[ax(3),ax(4)],'Color','r');

%% Calculate a one-tailed p-value based on this distribution

pValBS = sum(muDiffBoot < 0) / nBoot;
tStr = sprintf('p-value = %0.3f',pValBS);
title(tStr);

% QUESTION (Q8): What is the p-value?

%% Final thoughts

% There is a beautiful and important section on the relationship between
% confidence intervals and hypothesis tests in section 15.4 of E&T (p.
% 214). 
%
% Look at your figure 3. Key question: What value of alpha will make the
% lower end of the bootstrap confidence interval equal to 0? It is just the
% probability mass of the bootstrap distribution to the left of 0:
% sum(muDiffBoot < 0) / nBoot
%
% Compare figures 2 and 3 and consider this statement from E&T (p. 216):
% "In this sense, the permutation p-value measures how far the observed
% estimate, theta_hat, is from 0, while the bootstrap p-value measures how
% far 0 is from theta_hat." (Remember that 0 is our null value for the
% difference in survival times.)
%
% Which test is preferred? In general, the permutation test is the most
% rigorous way to test H0, since we are directly instantiating it by
% "breaking" the relationship we seek to test. In this case, we are
% randomly assigning the observed data to either the treatment or the
% control group. As E&T point out (p. 216): "The permutation p-value is
% exact, while the bootstrap p-value is approximate." They conclude with
% what should become one of your statistical mantras (p. 218): "When there
% is something to permute, it is a good idea to do so . . . ."

%% Bonus: Using the linear model to make the same comparisons.

% The trick is to use an indicator variable for the two groups. We'll use
% the ctrl groups as the reference (iv = 0) and the treatment group gets a
% value of one (iv = 1).
iCtrl = zeros(size(ctrlGrp));
iRx = ones(size(rxGrp));

% For the 2-sample t-test, our 'Y' is the values our 'X' is the indicator
% variables:
Y = [ctrlGrp,rxGrp];
X = [iCtrl,iRx];

% Now just fit the linear model and look at the beta1 coeff for a p-value:
[~,~,stats] = glmfit(X,Y);
pLM1 = stats.p(2);
display(pLM1);

% This is intrinsically a 2-sided test, so it is appropriate that it is
% exactly double the value for our 1-sided t-test. If we run a 2-sided
% t-test . . .
[hT,pValTtest,ci,statsT] = ttest2(rxGrp,ctrlGrp,'Tail','both');
display(pValTtest);
% . . . we get exactly the same answer. Amazing!

% Now let's do the same thing, but on ranked data, and compare this to the
% 2-sided Wilcoxon Rank Sum Test. Note that the LM approximation will only
% be exact for n >= 11 in each group.
[rnkY,tiedJ] = tiedrank(Y);

[~,~,stats] = glmfit(X,rnkY);
pLM2 = stats.p(2);
display(pLM2);

[pValWRST,hRS,statsRS] = ranksum(rxGrp,ctrlGrp,'Tail','both');
display(pValWRST);