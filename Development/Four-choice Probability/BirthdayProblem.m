% BirthdayProblem.m
%
% Given a group of N people, what is the probability that any two of them
% have the same birthday?
%
% See also: ProbHist.m
% See also: http://en.wikipedia.org/wiki/Birthday_problem
%
% RTB wrote it, 09 August 2013

R = 365; K = 2; nSims = 10000;

allP = zeros(1,100);
allN = 1:100;

for N = allN
    allP(N) = sum(any(hist(randi(R,N,nSims),1:R) >= K)) / nSims;
end

hp=plot(allN,allP,'r-');
set(hp,'LineWidth',3);
grid on
xlabel('# of people in group');
yStr = sprintf('probability that %d or more share a birthday',K);
ylabel(yStr);
title('The Birthday Problem');