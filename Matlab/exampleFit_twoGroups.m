% Demonstration of comparison between two independent groups
%
% SF 2014

clear all
close all

Ntrials = 500;
Nsub = 30;

c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;    % same dprime across groups
sigma = 0.5;

group_mratio(1) = 1;
group_mratio(2) = 0.6; % group 2 has worse metacognition than group 1

j=1;
for group = 1:2
    for i = 1:Nsub
        
        % Generate dprime
        d(j) = normrnd(group_d, sigma);
        
        % Generate data
        metad = group_mratio(group).*d;
        sim = metad_sim(d(j), metad, c, c1, c2, Ntrials);
        
        DATA(group).nR_S1{i} = sim.nR_S1;
        DATA(group).nR_S2{i} = sim.nR_S2;
        
        j = j+1;
    end
    
end

% Fit group 1
fit1 = fit_meta_d_mcmc_group(DATA(1).nR_S1, DATA(1).nR_S2);

% Fit group 2
fit2 = fit_meta_d_mcmc_group(DATA(2).nR_S1, DATA(2).nR_S2);

% Compute HDI of difference
sampleDiff = fit1.mcmc.samples.mu_logMratio - fit2.mcmc.samples.mu_logMratio;
hdi = calc_HDI(sampleDiff(:));
fprintf(['\n HDI on difference in log(meta-d''/d''): ', num2str(hdi) '\n\n'])

% Plot group Mratio and the difference
plotSamples(exp(fit1.mcmc.samples.mu_logMratio))
plotSamples(exp(fit2.mcmc.samples.mu_logMratio))
plotSamples(sampleDiff)
