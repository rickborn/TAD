% etASAdemo_ng.m: 
%
% Bootstrap example for Stroke data from pp. 3-5 of Efron & Tibshirani
% 
% RTB wrote it 29 October 2016 (derived from BS_ex1.m)
% RTB modified it 30 Jan 2017: combined etASAhypoth.m and etASAstats.m into
% one file that is modular for my stats class.
% RTB modified it to emphasize stroke data (13 September 2018)

% The scenario: A study was done to see if low-dose aspirin would prevent
% heart attacks in healthy middle-aged men. The study design was optimal:
% controlled, randomized, double-blind. Subjects were randomly assigned to
% receive aspirin (ASA) or placebo. The summary statistics:
%
% aspirin group (n=11037): 104 heart attacks (MI); 10933 no MI
% placebo group (n=11034): 189 heart attacks; 10845 no MI
%
% Scientific question #1: Does aspirin help to prevent heart attacks?
%
% In the same study, the researchers also recorded the number of strokes in
% each group:
%
% aspirin group (n=11037): 119 strokes; 10918 without stroke
% placebo group (n=11034): 98 strokes; 10936 without stroke
%
% Scientific question #2: Does aspirin increase the risk of having a
% stroke?
%
% We will start by addressing the 2nd question regarding strokes. The code
% you generate here will then allow you to rapidly analyze the heart attack
% data.

% What to do: Login to learning catalytics and join the session for the
% module entitled "ASA Bootstrap". You will answer a series of
% questions based on the guided programming below. Each section begins with
% a '%%'. Read through the comments and follow the instructions provided.
% In some cases you will be asked to answer a question, clearly indicated
% by 'QUESTION'. In other cases, you be asked to supply missing code,
% indicated by 'TODO'. The corresponding question in learning catalytics
% will be indicated in parentheses (e.g. Q1). If there is no 'Q#'
% accompanying a 'QUESTION' just type your answer into this script and
% discuss it with your team. Once you have supplied the required code, you
% can execute that section by mouse-clicking in that section (The block
% will turn yellow.) and then simultaneously hitting the 'ctrl' and 'enter'
% keys (PC) or 'command' and 'enter' keys (Mac).
%
% NOTE: If you are using version R2018a (or later) of MATLAB, you won't be
% able to use the ctrl+enter feature, because it now checks the entire
% script for errors, rather than just the cell you are trying to execute.
% This is stupid, but we're stuck with it. What you can do instead is use
% the mouse to highlight the code you want to run, then hit the F9 key (PC)
% or Shift+F7 (Mac). You can always find what the "Evaluate" shortcut is
% for your system by right-clicking the highlighted code in your script
% window. If this is too confusing, you can just copy the section in your
% editor and then paste it into the command line.

%% Concepts covered:
% 1. Test for proportions: odds ratio
% 2. Comparing resampling tests with Fisher's exact test
% 3. Std. error and confidence intervals through bootstrapping
% 4. Relationship between CI and hypothesis test
% 5. Permutation test for strong test of H0.
% 6. One-tailed vs. two-tailed tests.
% 7. CI with 'bootci' and an anonymous function
% 8. Making data tables with 'table' command

%% Constants: these would normally be passed as arguments to a function

nBoot = 10000;
myAlpha = 0.05;

% useful numbers for stroke data
nRx = 11037;        % total number of patients in the treatment group (ASA)
nStrokeRx = 119;    % number of Strokes in the treatment group
nCtrl = 11034;      % total number of patients in the control group (placebo)
nStrokeCtrl = 98;   % number of Strokes in the control group
nTotal = nRx + nCtrl;

%% Calculate the actual ratio of rates of disease: an odds ratio

% This is defined as the ratio of 2 ratios: the numerator is the ratio of
% the number of subjects in the treatment group who had a stroke divided by
% the number who did not have a stroke. The denominator is the same, but
% for the control group.

% TODO: calculate the odds ratio for this study
orHat = ;

% QUESTION (Q1): What is your odds ratio to 2 decimal places?

%% Create a population from which to resample:

% The general approach in bootstrapping is to resample from our *original
% sample*. But you've been given only proportions, so you have to
% *re-create* the raw data based on the proportions. It's alarmingly
% simple, but you might have to think about it for a bit.

% TODO: Re-create the raw data.
% HINT: When you're done, the 'size' command for 'rxGrp' should return
% 11037,1. And 'size' for 'ctrlGrp' should return 11034,1
% HINT-HINT: You should be able to re-calculate your original odds ratio
% with this formula:
% orHat2 = (sum(rxGrp) / sum(~rxGrp)) / (sum(ctrlGrp) / sum(~ctrlGrp));
rxGrp = ;    % aspirin group for strokes
ctrlGrp = ;  % placebo group for strokes

%% Generate bootstrap replicates of the odds ratio

orStar = zeros(nBoot,1);    % holds each bootstrap calc. of the odds ratio
rng default
for k = 1:nBoot
    % TODO: Re-sample from each group WITH REPLACEMENT to create two new
    % samples: rxStar and ctrlStar. Then use these two bootstrap samples to
    % calculate an odds ratio and store it in orStar.
    rxStar = ;
    ctrlStar = ;
    orStar(k) = ;
end

%% Make a histogram of our bootstrap replicates of OR

figure
subplot(2,1,1);
hist(orStar,20);
hold on;
xlabel('OR^*'); ylabel('#');
title('Distribution of bootstrapped odds ratios');
% Draw a vertical line to see where our actual value lies within the
% distribution of bootstrap replications. Does it make sense?
ax = axis;
h1=line([orHat,orHat],[ax(3),ax(4)],'Color','g','LineWidth',2);

%% Calculate the standard error and the confidence intervals

% QUESTION (Q2): What is bootstrap estimate of the standard error of the
% odds ratio?
semBoot = ;

% TODO: Use the percentile method to determine the 95% confidence interval.
!!! More code here
confInterval = ;

% QUESTION (Q3): What is the 95% CI based on your bootstrap distribution?

% QUESTION (Q4): What is the null value of the odds ratio?

% QUESTION (Q5): How can we use the known null value of the odds ratio to
% perform a hypothesis test?

% QUESTION (Q6): Can we reject H0 at an alpha of 0.05?

%% Plot CIs on histogram

h2=line([confInterval(1),confInterval(1)],[ax(3),ax(4)],'Color','r');
line([confInterval(2),confInterval(2)],[ax(3),ax(4)],'Color','r');

legend([h1,h2],{'orHat','95% CI'},'Location','NorthEast');

%% Perform an explicit hypothesis test by modeling our OR under H0

% In this case, we will use a permutation test, where we resample WITHOUT
% replacement. The logic is that we are essentially randomly assigning each
% patient to the treatment or control group, then recalculating our odds
% ratio. Here, we are testing the most extreme version of H0, which is that
% the two distributions are the SAME.

% TODO: Perform resampling as though the patients all belonged to the
% same group (called H0data), shuffle this data, then arbitraily assign
% each patient to the treatment or control group and compute the odd ratio. 
% Store each bootstrapped odds ratio in orPerm

% Pool all the data:
H0data = [rxGrp;ctrlGrp];
% Place to store our results:
orPerm = zeros(nBoot,1);

rng default
for k = 1:nBoot
    % Shuffle and assign to rxStar and ctrlStar (ALL values are used!)
    !!! Your code here
    orPerm(k) = ;
end

%% Plot the distrubtion of permuted ORs

subplot(2,1,2);
hist(orPerm,20);
hold on
xlabel('Permuted ORs'); ylabel('#');
title('Distribution of ORs under H0');

% draw a line for our actual odds ratio
axis(ax);
h3=line([orHat,orHat],[ax(3),ax(4)],'Color','g','LineWidth',2);
legend(h3,'orHat','Location','NorthEast');

%% Calcualte a 1-tailed p-value from the orPerm values under H0

% TODO: Calculate a one-tailed p-value based on your permuted samples (i.e.
% orPerm) and store it in a variable called 'pVal1t'
pVal1t = ;

% QUESTION (Q7): What is our one-tailed p-value for the odds ratio?

%% Calculate a 2-tailed p-value

% TODO: Calculate a two-tailed p-value based on your permuted samples (i.e.
% orPerm) and store it in a variable called 'pVal2t'
pVal2t = ;

% QUESTION (Q8): What is our two-tailed p-value for the odds ratio?

% Philosophical interlude: The difference in implementation between
% 2-tailed and 1-tailed p-value is pretty clear. The philosophical
% difference, somewhat less. If you accidentally coded a 2-tailed test and
% get a p value of, say,0.06, and then remember "Oh! A 1-tailed test was
% actually more appropriate!" (and it really is in that instance, not for a
% "p-hacky" reason) and obtain p ~ 0.03, there's a sudden shift in
% perspective on the data. But it's the same data, and you're performing
% more or less the same analysis. Does this seem even remotely reasonable?
% This very subtle distinction would have a pretty heavy impact on a
% statistics-na�ve researcher. It can be helpful to think about edge cases
% like this, where our arbitrary thresholding statistical procedure leads
% to binarization of the same data into two categories which are
% interpreted in very different ways, and how we should consider data of
% this variety. Is it helpful to construct a new categorization, e.g.
% "statistically significant (p small)," "unlikely to produce statistical
% significance (p biggish)" and "of uncertain relationship (p kinda
% small?)" or does that just move the problem?

%% Compare with the Fisher Exact Test for Stroke data

% load the data into an object of type 'table'
strokeData = table([nStrokeRx;nStrokeCtrl],[nRx-nStrokeRx;nCtrl-nStrokeCtrl],...
    'VariableNames',{'Stroke','NoStroke'},'RowNames',{'ASA','NoASA'});

% TODO: Calculate a 2-tailed p-value and 95% confidence interval using
% Fisher's Exact Test

% QUESTION (Q9): What p-value does Fisher's Exact Test give?

% QUESTION (Q10): What is the lower bound of the 95% CI from Fisher's Exact
% Test?
%
% Be sure to compare this with your values from the bootstrap!

% QUESTION (Q11): Save your final figure for the stroke data as a jpeg and
% upload it to the LC site.

%% Repeat calculations for the heart attack data:

% useful numbers for heart attack data
nRx = 11037;    % number of patients in the treatment group (ASA)
nMIrx = 104;    % number of heart attacks (= MIs) in the treatment group
nCtrl = 11034;  % number of patients in the control group (placebo)
nMIctrl = 189;  % number of MIs in the control group

% Generate the raw data
rxGrp = ;       % aspirin group for MI
ctrlGrp = ;     % non-aspirin group for MI

% Odds ratio for MI data:
orHat = ;

% QUESTION (Q12): What is the odds ratio, orHat, for the heart attack data?

% TODO: Repeat the above analysis for the MI data. If you wrote
% everything in terms of 'rxGrp' and 'ctrlGrp', then once you've generated
% the corresponding raw values for the MI data, you should be able to just
% re-run everything without further edits to your code.

% QUESTION (Q13): What is the bootstrap estimate of the 95% CI of the odds
% ratio for the heart attack data?
%
% QUESTION (Q14): Can we reject H0 at an alpha of 0.05 for the heart attack data?
%
% QUESTION (Q15): Based on the permutation test, what is your 1-sided
% p-value for the heart attack data?
% 
% Be sure to compare the confidence intervals that you obtain
% through bootstrapping to those obtained with Fisher's Exact Test.