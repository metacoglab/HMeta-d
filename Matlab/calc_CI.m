function hdi = calc_CI(samples, q)
% function hdi = calc_CI(samples, q)
%
% Calculates symmetric 95% confidence interval for MCMC samples
%
% INPUTS
%
% samples - vector of MCMC samples
% q - bounds on symmetric CI, default [0.025 0.975]
%
% Steve Fleming 2015 stephen.fleming@ucl.ac.uk

if ~exist('q','var') || isempty(fninv)

    q = [0.025 0.975];

end

hdi = quantile(samples, q);

