function [D,pHgH,pHgM] = GVTsim(nShots,kStreak,pHit,nSims,inverseSampleFlag,pFlag)

% Simulation of bias in counting streaks after Miller & Sanjurjo 2016
% "Surprised by the Gambler's and Hot Hand Fallacies? A Truth in the Law of
% Small Numbers"
% 
% [D,pHgH,pHgM] = GVTsim(100,3,0.5,10000,0,1);
% 
% Inputs:
% - nShots: total # of tries in each simulation (default = 100)
% - kStreak: # of hits that constitute a "streak" (default = 3)
% - pHit: probability of a hit on each try (default = 0.5)
% - nSims: # of simulations to run (default = 10,000)
% - inverseSampleFlag: use inverse sample to eliminate bias (default = 0)
% - pFlag: 1 = plot histogram (default)
% 
% Outputs:
% - D: Prob(hit|k hits) - Prob(hit | k misses)
% - pHgH: Prob(hit | k hits)
% - pHgM: Prob(hit | k misses)
% 
% Based on counting method and calculations used by:
% Gilovich, T., R. Vallone, and A. Tversky (1985) "The Hot Hand in
% Basketball: On the Misperception of Random Sequences," Cognitive
% Psychology, 17: 295-314.
% 
% This is a really nice simulation of a subtle counting bias that went
% undetected for >30 years. Perhaps the moral of the story is that any time
% you condition on 'strange things' (i.e. streaks), you can distort the
% sample space and introduce bias.
% 
% RTB wrote it, 09 November 2018, airplane trip from San Diego to Boston

if nargin < 6, pFlag = 1; end
if nargin < 5, inverseSampleFlag = 0; end
if nargin < 4, nSims = 10000; end
if nargin < 3, pHit = 0.5; end
if nargin < 2, kStreak = 3; end
if nargin < 1, nShots = 100; end

% shuffle the random number generator:
rng shuffle

% init variables to hold results of simulations:
D = zeros(nSims,1);
pHgH = D;
pHgM = D;
% filter for use with 'conv' to detect streaks of length kStreak
u = ones(1,kStreak);

% It occasionally happens that there will be no streaks (either hits or
% misses) of length 'k', in which case, we'll just re-generate a new random
% sequence. This is why I'm using a 'while' loop instead of a 'for' loop.
k = 1;
while k <= nSims
    % generate a random sequence of hits/misses:
    x = rand(1,nShots) <= pHit;
    
    % calculate p(hit|k hits)
    % We use the 'conv' trick to find the end of runs of >= kStreak.
    % NOTE: Using 'conv' without a 'shape' argument will pad 'x' with zeros
    % before doing the convolution and thus return a 'w' that is larger 
    % than 'x'. This works out so that the positions in 't' correspond to
    % the ends of a given streak. However, if you specify the 'shape'
    % parameter as 'same', then the indices in 't' will correspond to the
    % next-to-last position in the streak. You just have to be careful
    % about the indexing so that you will correctly identify the value of
    % the shot that follows a given streak.
    w = conv(x,u);
    
    % Different ways to do sampling: biased vs. "inverse sampling"
    % For the latter we don't count nested runs as multiple samples. That
    % is, we only sample the first instance of each run of kStreak.
    if inverseSampleFlag
        L = w >= kStreak;
        t = find(diff(L) == 1) + 1;
    else
        t = find(w >= kStreak);
    end
    
    if isempty(t)
        %warning('No hit streaks of length k')
        continue;
    end
        
    % Need to make sure we don't run off the end. That is, if 'x' ends in a
    % streak, 'conv' will find the last element as the end of a streak. But
    % we can't know the result of the next shot, because there isn't one.
    if max(t) == length(x)
        t = t(1:end-1);
    end
    
    pHitGivenKhits = sum(x(t+1)) / length(t);
    pHgH(k) = pHitGivenKhits;
    
    % calculate p(hit|k misses)
    % same but first convert misses to hits to find miss streaks
    w = conv(~x,u);
    if inverseSampleFlag
        L = w >= kStreak;
        t = find(diff(L) == 1) + 1;
    else
        t = find(w >= kStreak);
    end
    
    if isempty(t)
        %warning('No miss streaks of length k')
        continue;
    end
        
    if max(t) == length(x)
        t = t(1:end-1);
    end
    pHitGivenKmisses = sum(x(t+1)) / length(t);
    pHgM(k) = pHitGivenKmisses;
    
    % Gilovich et al. calculated the difference, which should be 0 under
    % the null hypothesis:
    D(k) = pHitGivenKhits - pHitGivenKmisses;
    k = k + 1;
end

if pFlag
    histogram(D);
    ylabel('#');
    xlabel(['P(Hit|',num2str(kStreak),' hits) - P(Hit|',num2str(kStreak),' misses)']);
    set(gca,'Fontsize',14);
    % draw vertical lines for the mean and the median:
    Dmean = mean(D);
    Dmedian = median(D);
    ax = axis;
    h1 = line([Dmean,Dmean],[ax(3),ax(4)],'Color',[0.7,0.2,0]);
    h2 = line([Dmedian,Dmedian],[ax(3),ax(4)],'Color',[0,0.7,0.2]);
    legend([h1,h2],{'mean','median'});
    
end
    