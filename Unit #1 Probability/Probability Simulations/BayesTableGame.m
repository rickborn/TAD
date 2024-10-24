function [fairOdds] = BayesTableGame(nSims)
%
% Simulates the "table game" used by Thomas Bayes to discover his famous
% theorem. Inspired by a nice article:
%
% Eddy, Sean R. 2004. “What Is Bayesian Statistics?” Nature Biotechnology
% 22(9):1177–78. https://doi.org/10.1038/nbt0904-1177.
%
% RTB wrote it, 24 October 2024, gorgeous autumn day, breaking in Limmers

rng shuffle

% parameters: 1 million simulations takes about 15 seconds
if nargin < 1, nSims = 1000000; end

% counters:
nAliceWins = 0;
nBobWins = 0;

for k = 1:nSims
% Start game here; simulate rolling a ball down the center of a pool table
% to divide it into "halves":
    p = rand;
    
    % Simulate 8 rounds of the game where Alice wins with probability, p,
    % and Bob wins with probability, 1 - p. We're only interested in cases
    % where Alice has a 5 to 3 lead after 8 rounds. We'll code Alice
    % winning as 1 and Bob winning as 0
    if binornd(8,p) == 5
        if binornd(3,p) == 0
            nBobWins = nBobWins + 1;
        else
            nAliceWins = nAliceWins + 1;
        end
    end
end

pAliceWins = nAliceWins / (nAliceWins + nBobWins);
pBobWins = nBobWins / (nAliceWins + nBobWins);

fairOdds = pAliceWins / pBobWins;
%display(fairOdds);

end
