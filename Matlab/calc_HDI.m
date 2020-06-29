function hdi = calc_HDI(samples, q)
% function hdi = calc_HDI(samples, q)
%
% Calculates highest-density interval on chain of MCMC samples
% See Kruschke (2015) Doing Bayesian Data Analysis
%
% INPUTS
%
% samples - vector of MCMC samples
% q - credible mass, scalar between 0 and 1, default is 0.95
% for 95% HDI
%
% Based on R code by Kruschke (2015)
% http://www.indiana.edu/~kruschke/BEST/
% and Matlab translation by Nils Winter
% https://github.com/NilsWinter/matlab-bayesian-estimation/blob/master/mbe_hdi.m
%
% Steve Fleming 2015 stephen.fleming@ucl.ac.uk

if ~exist('q','var') %|| isempty(fninv)
    q = 0.95;
end

sortedVec = sort(samples);
ciIdx = ceil(q * length(sortedVec));
nCIs = length(sortedVec) - ciIdx;  % number of vector elements that make HDI

% Determine middle of HDI to get upper and lower bound
ciWidth = zeros(nCIs,1);
for ind = 1:nCIs
    ciWidth(ind) = sortedVec(ind + ciIdx) - sortedVec(ind);
end
[~,idxMin] = min(ciWidth);
HDImin = sortedVec(idxMin);
HDImax = sortedVec(idxMin + ciIdx);
hdi = [HDImin, HDImax];

%
% if ~exist('q','var') || isempty(fninv)
%
%     q = [0.025 0.975];
%
% end
%
% % hdi = quantile(samples, q);

