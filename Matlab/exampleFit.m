% Example Bayesian meta-d fit (single subject)
clear all
close all

Ntrials = 1000;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];
d = 1.5;
meta_d = 1;

mcmc_params = fit_meta_d_params;
mcmc_params.estimate_dprime = 0;

% Generate data
sim = metad_sim(d, meta_d, c, c1, c2, Ntrials);

% Fit data
fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2, mcmc_params);
hdi = calc_HDI(fit.mcmc.samples.meta_d(:));
fprintf(['\n HDI: ', num2str(hdi) '\n\n'])

% Visualise single-subject fits
metad_visualise