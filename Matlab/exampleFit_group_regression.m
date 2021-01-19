% Example script for estimating a regression model
%
% Steve Fleming 2018

clear all

Ntrials = 100;
Nsub = 50;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;
group_baseline_mratio = 0.8;
sigma = 0.5;
sigma_beta = 0.2;
gen_beta = 0.5;
cov = rand(Nsub,1)';

mcmc_params = fit_meta_d_params;
mcmc_params.estimate_dprime = 0;

for i = 1:Nsub
        
        % Generate dprime
        d(i) = normrnd(group_d, sigma);
        beta(i) = normrnd(gen_beta, sigma_beta);
        metad(i) = (group_baseline_mratio + beta(i).*cov(i)).*d(i);    % note this is assuming no unmodelled variance in beta; previous line would do this
        
        % Generate data
        sim = metad_sim(d(i), metad(i), c, c1, c2, Ntrials);
        
        nR_S1{i} = sim.nR_S1;
        nR_S2{i} = sim.nR_S2;
        
end

%% Regression fit
% fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params);
fit = fit_meta_d_mcmc_regression(nR_S1, nR_S2, cov, mcmc_params);

% % Call plotSamples to plot posterior of group Mratio
plotSamples(exp(fit.mcmc.samples.mu_logMratio))
hdi = calc_HDI(exp(fit.mcmc.samples.mu_logMratio(:)));
fprintf(['\n HDI on meta-d/d: ', num2str(hdi) '\n\n'])

plotSamples(fit.mcmc.samples.mu_beta1)
hdi = calc_HDI(fit.mcmc.samples.mu_beta1(:));
fprintf(['\n HDI on beta1: ', num2str(hdi) '\n\n'])
