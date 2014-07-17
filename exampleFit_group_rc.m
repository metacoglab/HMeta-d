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

% Fit the data using response-conditional model
% MCMC Parameters
mcmc_params.response_conditional = 1;
mcmc_params.nchains = 3; % How Many Chains?
mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
mcmc_params.nsamples = 10000;  %How Many Recorded Samples?
mcmc_params.nthin = 1; % How Often is a Sample Recorded?
mcmc_params.doparallel = 0; % Parallel Option
mcmc_params.dic = 1;
% Initialize Unobserved Variables
for i=1:mcmc_params.nchains
    S.mu_Mratio_rS1 = 1;
    S.lambda_Mratio_rS1 = 0.5;
    S.mu_Mratio_rS2 = 1;
    S.lambda_Mratio_rS2 = 0.5;
    S.cS1_raw = linspace(-1,0.2,length(c1));
    S.cS2_raw = linspace(0.2,1,length(c2));
    mcmc_params.init0(i) = S;
end

% Fit group data all at once
fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params);
