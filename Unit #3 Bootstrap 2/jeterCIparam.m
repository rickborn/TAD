% jeterCIparam.m
%
% use parametric bootstrap to calculate confidence intervals for binomial data

% BA = 0.250 (12 hits in 48 at bats)
nAtBats = 48;
nHits = 12;
pHit = nHits / nAtBats;

% Run the bootstrap
nBoot = 10000;
allBA = zeros(nBoot,1);
for k = 1:nBoot
    thisSim = rand(nAtBats,1) <= pHit;
    allBA(k) = sum(thisSim) / nAtBats;
end

% Sort the data
allBAsorted = sort(allBA);

% Calculate the CI
myAlpha = 0.05;
idxLo = floor((myAlpha/2) * nBoot);    % index corresponding to lower bound
idxHi = ceil((1-(myAlpha/2)) * nBoot); % index corresponding to upper bound

ciLo = allBAsorted(idxLo);
ciHi = allBAsorted(idxHi);

% make a histogram:
figure, histogram(allBA);
hold on
ax = axis;
line([ciLo,ciLo],[ax(3),ax(4)],'Color','r');
line([ciHi,ciHi],[ax(3),ax(4)],'Color','r');
xlabel('Batting average');
ylabel('#');