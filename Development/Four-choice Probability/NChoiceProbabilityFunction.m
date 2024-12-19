function [pVal] = NChoiceProbabilityFunction(nFolks,R,k,laxFlag,nSims)

% NChoiceProbabilityFunction: probability of k or more choices of a given number
% when nFolks pick randomly from 1:R (inclusive)
%
% [pVal] = NChoiceProbabilityFunction(nFolks,R,k,laxFlag,nSims)
%
% Inputs:
% - nFolks: size of group doing the picking
% - R: range over which they chose (1:R, inclusive)
% - k: number who picked a given number (e.g. 3)
% - laxFlag: 1 (default) means no particular number is special
% - nSims: # of simulations to perform (default = 10,000)
%
% Outputs:
% - pVal: probability of k or more
%
% RTB wrote it, Friday 12 August 2016

% Set some defaults:
if nargin < 5, nSims = 10000; end
if nargin < 4, laxFlag = 1; end
if nargin < 3, k = 30; end
if nargin < 2, R = 4; end
if nargin < 1, nFolks = 65; end

% generate all of the random picks for all of the simulations
allSims = unidrnd(R,nFolks,nSims);

% use the histogram function to tally up the choices
% NOTE: The 2nd argument to 'hist' tells it to use bins of 1,2,3,4
allBinned = hist(allSims,[1:R]);

% determine for how many simulations we exceeded k choices
if laxFlag
    allSuccesses = any(allBinned >= k);    % 1, 2, 3 or 4
else
    % By symmetry, any bin will do
    allSuccesses = allBinned(3,:) >= k;     % only '3'
end

pVal = sum(allSuccesses) / nSims;
