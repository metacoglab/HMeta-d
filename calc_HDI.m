function hdi = calc_HDI(samples, q)
% function hdi = calc_HDI(samples, q)
%
% Calculates highest-density interval on chain of MCMC samples
% See Kruschke (2010) DBDA
%
% INPUTS
%
% samples - vector of MCMC samples
% q - bounds of cumulative probabilities for HDI, default is [0.025 0.975]
% for 95% HDI
%
% Steve Fleming 2015 sf102@nyu.edu

if ~exist('q','var') || isempty(fninv)
    
    q = [0.025 0.975];
    
end

hdi = quantile(samples, q);

