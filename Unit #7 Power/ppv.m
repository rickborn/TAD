function [y] = ppv(pNull,alpha,power)

% ppv.m: positive predictive value
%
% PPV = ppv(pNull,alpha,power);
%
% Inputs:
% - pNull: prior probability of a null result (default = 0.9)
% - alpha: type I error (default = 0.05)
% - power: 1 - type II error (default = 0.8)
%
% Button KS, Ioannidis JP, Mokrysz C, Nosek BA, Flint J, Robinson ES,
% Munafò MR. Power failure: why small sample size undermines the
% reliability of neuroscience. Nat Rev Neurosci. 2013 May;14(5):365-76.
% doi: 10.1038/nrn3475. Epub 2013 Apr 10. Review. Erratum in: Nat Rev
% Neurosci. 2013 Jun;14(6):451. PubMed PMID: 23571845.
%
% Note that the default value for the prior probability of a null result is
% based on:
%
% Johnson VE, Payne RD, Wang T, Asher A, Mandal S (2017) On the
% reproducibility of psychological science. Journal of the American
% Statistical Association. 112(517):1-10.
%
% RTB wrote it 17 April 2018 for a slide in his Reproducibility lecture to postdocs

if nargin < 3, power = 0.8; end
if nargin < 2, alpha = 0.05; end
if nargin < 1, pNull = 0.9; end

R = p2odds(1-pNull);
y = (power.*R) ./ (power.*R + alpha);

end

