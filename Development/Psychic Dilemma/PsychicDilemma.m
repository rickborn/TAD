% PsychicDilemma.m
%
% RTB wrote it, 18 August 2018. New exercise for QMBC.
% RTB added Bayesian approach sections, 10 November 2021, at BWH for
% treatment of atrial fibrillation
%
% Tip o'the Pin to John K. Kruschke's, "Doing Bayesian Data Analysis"
% see Chapter 11, p. 297 on "Null Hypothesis Significance Testing"
%
% see also: pNflips2getZheads.m

%% An experiment in telekinesis

% You are a psychologist studying a subtle form of telekinesis that
% involves using thought to bias the outcomes of the flips of a fair coin.
% A woman claiming to have psychic powers tells you that by concentrating
% intensely on "negative thoughts about the outcome 'heads'" she can create
% runs of outcomes with fewer heads than one would predict by chance.
% Before you begin the experiment, the subject tells you that the psychic
% effort she must exert is extreme, so the most she can possibly do is 25
% trials. While you watch, she obtains the following sequence:
% 
%        TTHHTTHTTTHTTTTTTHTTHHTTH
% 
% That is, 8 heads on 25 flips of the coin. We have obtained the coin from
% the National Institute of Standards and Technology, and tested it
% extensively under non-psychic conditions, so we are abolutely certain
% that it is fair (i.e. pHeads = 0.5).
% 
% What is the probability that, out of 25 flips, she obtained 8 or fewer
% heads?

nHeads = 8;
nTrials = 25;
pNull = 0.5;
% With fixed 'N', the appropriate sample space is determined by the
% binomial distribution:
pFixedN = binocdf(nHeads,nTrials,pNull);
% p = 0.0539; Alas, she is not psychic.

%% Simulate for fixed N

nSims = 100000;
pHead = 0.5;
nReal = 25; % actual # of trials it took us
zHeads = 8;

allHeads = sum(round(rand(nReal,nSims)));
pSim = sum(allHeads <= zHeads) / nSims;

%% View the distribution

figure
pStr = sprintf('p = %.3f',pSim);
tStr = sprintf('Number of heads from %d tosses if p(Heads) = %.2f',nReal,pHead);
subplot(2,1,1)
histogram(allHeads);
xlabel('# of Heads');
ylabel('#');
title(tStr);
ax = axis;
line([zHeads,zHeads],[ax(3),ax(4)],'Color','r');
text(ax(2)*0.8,ax(4)*0.9,pStr);

%% A last-minute confession:

% After the experiment, however, your subject tells you that she was not
% completely truthful about her 'psychic effort' story. In fact, the only
% way she can exert her powers effectively is if one of her lucky
% numbers is involved. Since 8 is her luckiest number, she had determined
% beforehand that she would stop when she got the 8th head. This just
% happened to occur on the 25th flip.
% 
% Does this knowledge change your analysis or your conclusions? Ans.: Yes!
% This is a very different data-generating process and will give us a
% different p-value. We can calculate it by simulation (see below), or we
% can use the negative binomial distribution, y = nbincdf(x,R,p), which
% gives the probability of x-or-fewer failures *before* we get R successes
% given a single-trial probability of success of p. In our case, we want to
% know the probability of 17 or *more* failures before 8 successes. Since
% the function returns x-or-fewer, we need to use 1 - p(x-or-fewer - 1).
% In MATLAB, this would be:
pFixedZ = nbincdf(nTrials-nHeads-1,nHeads,pNull,'upper');
% p = 0.0320; She IS psychic! You begin to write your Nature paper.

% Teaching note: Most of the students will not know about the negative
% binomial distribution, but they can easily do the calculation by
% simulation.

%% Simulation for fixed Z

nSims = 100000;
pHead = 0.5;
nReal = 25; % actual # of trials it took us
zHeads = 8;

% Save the data in two different ways:
% 1. The # of coin flips it took us to get to zHeads on each sim trial:
allN = zeros(nSims,1);
% 2. The proportion of heads (zHeads/nFlips) we got on each sim trial:
allProps = zeros(nSims,1);

for k = 1:nSims
    nFlips = 0;
    nHeads = 0;
    while nHeads < zHeads
        nFlips = nFlips + 1;
        nHeads = nHeads + (rand < pHead);
    end
    allProps(k) = nHeads / nFlips;
    allN(k) = nFlips;
end
pVal = sum(allProps <= (zHeads ./ nReal)) ./ nSims;

%% Plot the results of our simulation:

% figure
% subplot(2,1,1)
% histogram(allProps);
% xlabel('Sample proportion (z/N)');
% ylabel('#');
% tStr = sprintf('Number of flips to get %d heads if p(Heads) = %.2f',...
%     zHeads,pHead);
% title(tStr);
% ax = axis;
% line([zHeads/nReal,zHeads/nReal],[ax(3),ax(4)],'Color','r');
% pStr = sprintf('p = %.3f',pVal);
% text(ax(2)*0.8,ax(4)*0.9,pStr);

tStr = sprintf('Number of flips to get %d heads if p(Heads) = %.2f',...
     zHeads,pHead);
subplot(2,1,2)
%figure
histogram(allN);
xlabel('# of Flips');
ylabel('#');
title(tStr);
ax = axis;
line([nReal,nReal],[ax(3),ax(4)],'Color','r');
pStr = sprintf('p = %.3f',pVal);
text(ax(2)*0.8,ax(4)*0.9,pStr);

%% Other teaching notes:
% 1. This problem, namely that the p-value we calculate depends on our
% stopping intention, has been touted as a fatal flaw in frequentist Null
% Hypothesis Significance Testing (NHST), and used to promote a Bayesian
% approach to statistics. In the latter, the statistician only considers
% the likelihood of the data actually obtained. For coin tosses, this is
% simply (pHeads^z)*[(1-pHeads)^(N-z)], for different values of pHeads.
%
% 2. The frequentist approach does not incorporate any prior information
% about the likelihood that psychic powers actually exist. That is, we are
% treating this experiment with this one subject as if it were the first
% and only such experiment. In reality, we might have a strong prior belief
% that psychic powers don't exist. The Bayesian approach gives us an
% elegant way to combine our prior beliefs with data from an experiment.
%
% 3. The arbitrary nature of using p < 0.05 as a cut-off for Truth. The
% two p-values we obtained are not that different--they just happen to fall
% on opposite sides of the sacred 0.05.

%% The Bayesian approach, dumbed down

% Let's imagine that we have only two specific hypotheses:
%   H0: She is NOT psychic and the behavior of the coin is pH = pT = 0.5
%   HA: She IS psychic, and can influence pH = 0.25 and pT = 0.75
%
% In this case, we first calculate the likelihood of our data under the two
% different hypotheses:
pH0 = 0.5;
pHA = 0.25;
nHeads = 8;
nTrials = 25;

pDataGivenH0 = (pH0.^nHeads).*((1-pH0).^(nTrials-nHeads));
pDataGivenHA = (pHA.^nHeads).*((1-pHA).^(nTrials-nHeads));

% A common way to compare these is with a ratio, known as a 'Bayes Factor'.
% In this case, we'll put the probability of the alternative hypothesis in
% the numerator, so that we effectively calculating the likelihood ratio in
% favor of our alternative hypothesis. This is typically notated as 'BF10'
% to indicate that the null is in the denominator. If we did it the
% opposite way (i.e. null in numerator), we'd write it as 'BF01'
BF10 = pDataGivenHA / pDataGivenH0;

% So we see that our data are roughly 4x more likely under the alternative
% hpothesis than under H0. Does this mean we stop here and declare our
% subject a psychic?

% NO! From an epistemological perspective, the Bayes Factor is an
% indication of how much our experiment should change our belief. In this
% case, whatever our beliefs about psychic powers were before our
% experiment, we should now be 4x more likely to believe that psychic
% powers exist. Bayes showed us how to make the calculation. But the trick
% is that we need to put a number on our prior belief in psychic powers.
% What should this be? Some people might say that this is zero, but this
% kind of begs the question. Let's start by being relatively lenient and
% saying that we think the prior probability that our subject is psychic is
% one in one thousand. Now we just plug values into Bayes Rule:
priorHA = 0.001;
pHAGivenData = (pDataGivenHA * priorHA) / ...
    ((pDataGivenHA * priorHA) + (pDataGivenH0 * (1 - priorHA)));

% Of course, we could have just multiplied our prior by the Bayes Factor.
% So our new belief is that the probability that our suject is psychic is
% 0.0038. Greater than when we started, but still a pretty remote
% possibility.
% [NOTE: To be precise, we should multiply our prior ODDS RATIO by the
% Bayes Factor. But for very small p, the odds = p. See my 'p2odds.m'

%% Maximum Bayes Factors

% Where did our HA that pHeads = 0.25 come from? We just pulled it out of
% thin air. But, in reality, we might consider the entire range of possible
% values of pHeads, from 0 to 1. In this way, we could find the value of
% pHeads that maximizes the likelihood of our data and choose to test the
% hypothesis that makes our Bayes Factor as large as possible:
pH = 0:0.01:1;
pDataGivenH = (pH.^nHeads).*((1-pH).^(nTrials-nHeads));
pHmax = pH(pDataGivenH == max(pDataGivenH));

% In this case, our most charitable estimate of the psychic's ability is
% that she can make pHead = 0.32 and pTail = 0.68

% So now we can re-run our calculation using this value:
priorHA = 0.001;
pH0 = 0.5;
pHA = pHmax;
pDataGivenH0 = (pH0.^nHeads).*((1-pH0).^(nTrials-nHeads));
pDataGivenHA = (pHA.^nHeads).*((1-pHA).^(nTrials-nHeads));
pHAGivenData = (pDataGivenHA * priorHA) / ...
    ((pDataGivenHA * priorHA) + (pDataGivenH0 * (1 - priorHA)));

% We now have a posterior belief of 0.0052, or, in other words, we found
% the largest possible Bayes Factor in favor of HA is 5.2
maxBF = pDataGivenHA / pDataGivenH0;

% Note that in the literature, it is more common to use the version of the
% Bayes Factor with pDataGivenH0 in the numerator, in which case we write
% it as BF01. In this case, we are looking for the "minimum Bayes Factor."

% There is an extensive literature on relating p-values to Bayes Factors
% and using them to adjust beliefs. For starters see:

% Goodman SN. Of P-values and Bayes: a modest proposal. Epidemiology. 2001
% May;12(3):295-7.
%
% Stern HS. A Test by Any Other Name: P Values, Bayes Factors, and
% Statistical Inference. Multivariate Behav Res. 2016;51(1):23-9.
%
% Johnson VE. Revised standards for statistical evidence. Proc Natl Acad
% Sci U S A. 2013 Nov 26;110(48):19313-7
%
% Nuzzo R. Scientific method: statistical errors. Nature. 2014 Feb
% 13;506(7487):150-2.

%% The full Bayesian approach: all H's considered


