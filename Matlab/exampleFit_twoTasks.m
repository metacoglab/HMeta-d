% Demonstration of estimating covariance between paired measures of meta-d'
% (e.g. two tasks from the same individuals)
%
% SF 2019

clear all

Ntrials = 400;
Nsub = 100;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;
group_mratio = 0.8;
type1_sigma = 0.2;
rho = 0.6;
type2_sigma = 0.5;

for i = 1:Nsub
        
        % Generate Mratios for this subject
        bigSigma = [type2_sigma^2 rho.*type2_sigma^2; rho.*type2_sigma^2 type2_sigma^2];
        mratios(i,:) = mvnrnd([group_mratio group_mratio], bigSigma); 
        
        %% Task 1
        % Generate dprime
        d = normrnd(group_d, type1_sigma);
        metad = mratios(i,1).*d;
        
        % Generate data
        sim = metad_sim(d, metad, c, c1, c2, Ntrials);
        
        nR_S1(1).counts{i} = sim.nR_S1;
        nR_S2(1).counts{i} = sim.nR_S2;
        
        %% Task 2
        d = normrnd(group_d, type1_sigma);
        metad = mratios(i,2).*d;
        
        % Generate data
        sim = metad_sim(d, metad, c, c1, c2, Ntrials);
        
        nR_S1(2).counts{i} = sim.nR_S1;
        nR_S2(2).counts{i} = sim.nR_S2;        
        
end

% Fit group data all at once
mcmc_params = fit_meta_d_params;
fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params);
plotSamples(fit.mcmc.samples.mu_logMratio(:,:,1))
plotSamples(fit.mcmc.samples.mu_logMratio(:,:,2))

% Compute HDI of difference between tasks 
sampleDiff = fit.mcmc.samples.mu_logMratio(:,:,1) - fit.mcmc.samples.mu_logMratio(:,:,2);
hdi = calc_HDI(sampleDiff(:));
fprintf(['\n HDI on difference in log(meta-d''/d''): ', num2str(hdi) '\n\n'])

% Plot difference in meta-d/d ratio between two tasks 
plotSamples(sampleDiff)

% Plot estimate of correlation
h1 = figure;
set(gcf, 'Position', [200 200 400 300])
h= histogram(fit.mcmc.samples.rho(:), 'Normalization', 'probability');
xlabel('\rho');
ylabel('Posterior density');
line([rho rho],[0 max(h.Values)+0.015], 'LineWidth', 2, 'Color', 'k', 'LineStyle', '--')
ci = calc_CI(fit.mcmc.samples.rho(:));
line([ci(1) ci(2)],[0.002 0.002], 'LineWidth', 3, 'Color', [1 1 1])
box off
set(gca, 'FontSize', 14, 'XLim', [-1 1])