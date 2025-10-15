% prob2pair.m
%
% Use a simulation to calculate the probability of drawing 2 pairs

rng shuffle
nSims = 100000;
nSuccesses = 0;
allHands = [];

% create a deck of cards
% for this calculation, only the ranks matter (not the suits)
myDeck = [1:13,1:13,1:13,1:13];
% number of cards to deal:
n = 5;

for k = 1:nSims
    % deal a hand of n cards
    myHand = datasample(myDeck,n,'Replace',false);

    % sort
    a = sort(myHand);
    b = diff(a);

    if sum(b == 0) == 2     % might be 2 pair, could be 3 of a kind
        c = diff(b);
        if sum(c == 0) == 0
            nSuccesses = nSuccesses + 1;
            allHands = [allHands;myHand];
        end
    end
end

prob = nSuccesses / nSims;
disp(prob);