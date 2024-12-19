% GaltonBoard.m
% Francis Galton (1822-1911) developed a mechanical device that would
% generate a binomial distribution. Watch this short video to see how it
% works: https://www.youtube.com/watch?v=6YDHBFVIvIs
% 
% Use MATLAB to simulate 1000 balls going through a Galton board with 16
% rows of pins. Post a figure showing their final distribution.
%
% RTB wrote it, 01 August 2017, hot summer day listening to Tales from
% Topographic Oceans

% We can think of each row of pins as enforcing a binary decsision with
% equal probability of the ball turning to the right or the left. So, in
% essence, this is like flipping a fair coin 16 times and counting the
% number of heads (or tails). So getting 16 heads would be like the ball
% making 16 'rightwards' turns.

nRowsOfPins = 16;
nBalls = 100000;

% It can all be done with a single line of code. It's a good test of your
% MATLAB knowledge to try to parse such lines, but it's a terrible thing to
% do in practice. Why?
rng default
numRightTurns = sum(round(rand(nRowsOfPins,nBalls)));

% Plot a histogram
xBins = 0:nRowsOfPins;       % all possible #s of right-hand turns
binCounts = hist(numRightTurns,xBins);
figure
bar(xBins,binCounts);
xlabel('Number of right-hand turns');
ylabel('# of Balls');

% Questions #2-4: You probably noticed that your histogram for question
% #1 looked pretty close to a bell-shaped curve (a.k.a. Gaussian or Normal
% distribution).
% 
% Q2: Use 'binopdf' to calculate the probability of getting a ball in the
% middle bin (i.e. it makes 8 right and 8 left turns, in any order).
x=8; n=16; p=0.5;
pBino = binopdf(x,n,p); % 0.1964

% Q3: Use 'normpdf' to make the same calculation
pNorm = normpdf(x,n*p,sqrt(n*p*(1-p))); % 0.1995

% Q4: Which number is more accurate?
% Ans: The binomial probability is exact; the normal is an approximation
% that gets better as n gets large. Even with an n of 16, the normal
% approximation is only off by about 1.5%.

% Q5: What is your estimate of the same probability based on your
% simulation?
pSim = sum(numRightTurns == x) / nBalls;


