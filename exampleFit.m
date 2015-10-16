% Example Bayesian meta-d fit (single subject)
clear all
close all

Ntrials = 1000;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];
d = 2;
noise = 0.5;

% Generate data
sim = type2_SDT_sim(d, noise, c, c1, c2, Ntrials);

% Fit data
fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2)

% Visualise fits
metad_visualise