% Demonstration of hierarchical model fits
%
% SF 2014

clear all

Ntrials = 200;
Nsub = 30;
c = -0.3;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;
group_mratio = 0.8;
sigma = 0.5;

mcmc_params = fit_meta_d_params;
mcmc_params.estimate_dprime = 0;

for i = 1:Nsub
        
        % Generate dprime
        d(i) = normrnd(group_d, sigma);
        metad = group_mratio.*d(i);
        
        % Generate data
        sim = metad_sim(d(i), metad, c, c1, c2, Ntrials);
        
        nR_S1{i} = sim.nR_S1;
        nR_S2{i} = sim.nR_S2;
        
end

% Fit group data all at once
fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params);

% Call plotSamples to plot posterior of group Mratio
plotSamples(exp(fit.mcmc.samples.mu_logMratio))
hdi = calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
fprintf(['\n HDI on meta-d/d: ', num2str(hdi) '\n\n'])

metad_group_visualise
