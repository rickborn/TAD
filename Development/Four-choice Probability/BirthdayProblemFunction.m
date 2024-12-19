function p = BirthdayProblemFunction(N, K, nSims)
%
% Given a group of N people, what is the probability that any two of them
% have the same birthday?
%
% See also: ProbHist.m
% See also: http://en.wikipedia.org/wiki/Birthday_problem
%
% RTB wrote it, 13 August 2014

if nargin < 3, nSims = 10000; end
if nargin < 2, K = 2; end
if nargin < 1, N = 20; end

% # of days in one year
R = 365;
p = sum(any(hist(randi(R,N,nSims),1:R) >= K)) / nSims;