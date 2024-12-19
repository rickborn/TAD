% 4-choice probability

% For 65 people "randomly" picking a number from 1 to 4 (inclusive) what is 
% the smallest number of picks of any given number which would cause you to 
% reject the null hypothesis (randomness) at p < 0.05?
%
% For either 1,2,3 or 4, Ans. = 25.
% For 3 only, Ans. = 23.

nFolks = 65; R = 4; nSims = 1000000; k = 23;

allSims = unidrnd(R,nFolks,nSims);
allBinned = hist(allSims,[1:R]);
%allSuccesses = any(allBinned >= k);    % 1, 2, 3 or 4
allSuccesses = allBinned(3,:) >= k;     % only '3'
pSim = sum(allSuccesses) / nSims

% Note: You can solve the 'only 3' case using the binomial:
%pBino = binocdf(k-1,nFolks,1/R,'upper')
