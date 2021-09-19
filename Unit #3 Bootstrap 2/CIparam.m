function [ci] = CIparam(nHits,nAtBats,myAlpha,nBoot)
%
% use the bootstrap to calculate confidence intervals for binomial data

if nargin < 4, nBoot = 10000; end
if nargin < 3, myAlpha = 0.05; end

% Generate the estimate for pHit:
pHit = nHits / nAtBats;

% Run the bootstrap
allBA = zeros(nBoot,1);
for k = 1:nBoot
    thisSim = rand(nAtBats,1) <= pHit;
    allBA(k) = sum(thisSim) / nAtBats;
end

% Sort the data
allBAsorted = sort(allBA);

% Calculate the CI
idxLo = floor((myAlpha/2) * nBoot);    % index corresponding to lower bound
idxHi = ceil((1-(myAlpha/2)) * nBoot); % index corresponding to upper bound

ci = [allBAsorted(idxLo),allBAsorted(idxHi)];
end