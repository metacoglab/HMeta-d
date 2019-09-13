% Demonstration of hierarchical model fits
%
% SF 2014

clear all
close all

Ntrials = 1000;
Nsub = 10;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;
sigma = 0.5;
noise = [1/3 2/3];

for i = 1:Nsub
    
    % Generate dprime
    d(i) = normrnd(group_d, sigma);
    
    % Generate data
    sim = type2_SDT_sim(d(i), noise, c, c1, c2, Ntrials);
    
    nR_S1{i} = sim.nR_S1;
    nR_S2{i} = sim.nR_S2;
    
end

% Get default parameters
mcmc_params = fit_meta_d_params;
% Change defaults to make response-conditional
mcmc_params.response_conditional = 1;

% Fit group data all at once
fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params);

% Plot output
plotSamples(exp(fit.mcmc.samples.mu_logMratio_rS1))
plotSamples(exp(fit.mcmc.samples.mu_logMratio_rS2))