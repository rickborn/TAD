% gptcointosssim.m
%
% On 14 June 2023, I asked Chat-GPT4 the following:
% Write a MATLAB script that uses a simulation to calculate the probability
% of getting "heads" on two successive tosses of a fair coin.
%
% It replied:
% Sure, you can simulate this scenario in MATLAB. Let's say you want to do
% this simulation 1 million times to get a reliable estimation of the
% probability. Here's an example of how you could do that:

% Set number of simulations
numSims = 1e6;

% Initialize a counter for successful outcomes
successes = 0;

% Iterate over the number of simulations
for i = 1:numSims
    
    % Simulate two coin tosses (1 is heads, 0 is tails)
    toss1 = rand > 0.5;
    toss2 = rand > 0.5;
    
    % Check if both tosses are heads
    if toss1 && toss2
        % If so, increment the success counter
        successes = successes + 1;
    end
end

% Calculate probability by dividing successful outcomes by total simulations
prob = successes / numSims;

% Display the result
disp(['The estimated probability of getting heads on two successive tosses is: ', num2str(prob)]);

% And then Chat-GPT4 provided the following explanation:
% This script first initializes a counter for successful outcomes (i.e.,
% two heads in a row). It then loops over the specified number of
% simulations, randomly simulating two coin tosses in each iteration. If
% both tosses result in heads, the script increments the counter for
% successful outcomes. Finally, it calculates the probability of getting
% two heads in a row as the ratio of successful outcomes to total
% simulations, and then it prints this probability.
%
% Note that we are using the rand function to simulate the coin toss. This
% function returns a random number between 0 and 1. If the number is
% greater than 0.5, we interpret that as a head; if it's less than or equal
% to 0.5, we interpret it as a tail. Since rand generates numbers uniformly
% distributed over the interval (0,1), this gives us a fair coin toss.