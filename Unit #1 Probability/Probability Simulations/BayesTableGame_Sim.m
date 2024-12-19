% BayesTableGame_Sim
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
nSims = 1000000;

% counters:
nAliceWins = 0;
nBobWins = 0;
allP = ones(nSims,1) .* NaN;

for k = 1:nSims
    % Start game here; simulate rolling a ball down the center of a pool table
    % to divide it into "halves":
    p = rand;
    
    % Simulate 8 rounds of the game where Alice wins with probability, p,
    % and Bob wins with probability, 1 - p. We're only interested in cases
    % where Alice has a 5 to 3 lead after 8 rounds. We'll code Alice
    % winning as 1 and Bob winning as 0
    if binornd(8,p) == 5
        allP(k) = p;
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
display(fairOdds);

% Plot a histogram of all p's that met our condition:
histogram(allP);
xlabel('p');
ylabel('#');

% If we've done things correctly, we should have modeled:
[pHat, pCI] = binofit(5,8);

% draw lines on histogram
ax = axis;
% expectation
h1 = line([pHat, pHat], [ax(3), ax(4)]);
set(h1,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',2);
% lower bound of 95% CI
h2 = line([pCI(1), pCI(1)], [ax(3), ax(4)]);
set(h2,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',1,'LineStyle','--');
% upper bound of 95% CI
h3 = line([pCI(2), pCI(2)], [ax(3), ax(4)]);
set(h3,'Color',[0.6350, 0.0780, 0.1840],'LineWidth',1,'LineStyle','--');

