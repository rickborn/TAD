% 4-choice probability: Answers to LC questions 1 to 5

% Q1: The figures below show two sets of 25 histograms. One set was
% generated with random draws from the uniform discrete distribution (i.e.
% 1,2,3 or 4 with equal probability); for the other set, some of the
% histograms were "tweaked" a little. Which set of histograms are truly
% random?
%
% A1: Panel 'A' is the truly random set, generated using 'unidrnd'. For set
% 'B', I also used 'unidrnd', but then I checked the values of the maximum
% and minimum bins, and, if the difference exceeded a certain threshold, I
% nudged them to be closer to each other. This had the effect of making
% this set look flatter, or more uniform, which might conform to our
% intuition about what a uniform random process should generate. The bottom
% line is that truly random processes often generate what look to our
% brains like interesting patterns. We need to use simulations and proper
% statistical analyses to check these (often faulty) intuitions.

% Q2: On August 7, 2017, 47 people picked a random number from 1 to 4
% (inclusive). Here is the distribution of responses: 4 people picked 1, 14
% picked 2, 20 picked 3, 9 picked 4.
%
% Do students prefer the number 3 over other numbers? Given 47 draws from a
% uniform discrete random distribution over [1,4], what is the probability
% of getting the number 3, 20 or more times?

nFolks = 47; R = 4; nSims = 100000; k = 20;

allSims = unidrnd(R,nFolks,nSims);
allBinned = hist(allSims,[1:R]);
allSuccesses = allBinned(3,:) >= k;     % only '3'
pVal = sum(allSuccesses) / nSims;
% p = 0.0064, so choice 'E' is correct

% Q3: What is your conclusion?
%
% A3: We might tentatively conclude that our class is biased towards the
% number '3', but this would be dangerous given that we looked at our data
% before we generated the hypothesis we tested. But we've been doing this
% demonstration for 14 years, and the class *always* is heavily biased
% towards the number '3'. But the correct thing to do would be to out and
% do a new experiment, this time stating our hypothesis ("People are biased
% towards the number '3'.") at the outset, before we collected more data.

% Q4: What if we are not interested in the number 3 in particular, but in
% the question of whether any number gets picked more often than the the
% criterion? What is your p-value then?

% A4: All we need to do is change our definition of 'success'. In this
% case, I used a neat logical trick: the 'any' function:
allSuccesses = any(allBinned >= k);    % 1, 2, 3 or 4
pVal = sum(allSuccesses) / nSims;
% p = 0.0255

% Q4b: What is the relationship between the p-value you obtained for this
% question and the one for question #2? Would you have needed to run this
% simulation to get the answer, knowing the answer from question #2?
%
% A4b: By the symmetry of the problem (i.e. We should get the same
% probability whether we look for outliers in any of the 4 bins.), our new
% p-value should be 4x our answer for Q2.

% Q5: For 65 people "randomly" picking a number from 1 to 4 (inclusive)
% what is the smallest number of picks, k, of any given number that would
% cause you to reject the null hypothesis (randomness) at p < 0.05? Note:
% In this case, we want to reject H0 if EITHER 1 or 2 or 3 or 4 is picked
% more than k times.
nFolks = 65; R = 4; nSims = 100000;

allSims = unidrnd(R,nFolks,nSims);
allBinned = hist(allSims,[1:R]);

% Note that at this point, we've done all of the hard work, and we
% effectively have the entire probability distribution. Let's look at our
% histogram:
figure, histogram(allBinned);

% Now all we really need to do is ask where we would draw a vertical line
% on our histogram such that the sum of the bins to the right (including
% our criterion bin) was only 5% of the total sum (= nSims)
pVal = 1;
k = -1;
while(pVal > 0.05)
    k = k+1;
    allSuccesses = any(allBinned >= k);    % 1, 2, 3 or 4
    pVal = sum(allSuccesses) / nSims;
end
disp(k);    % 25
disp(pVal); % 0.0445


