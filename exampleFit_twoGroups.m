% Demonstration of hierarchical model fits
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

groupIndex = [];
j=1;
for group = 1:2
    for i = 1:Nsub
        
        % Generate dprime
        d(j) = normrnd(group_d, sigma);
        
        % Generate data
        metad = group_mratio(group).*d;
        sim = metad_sim(d(j), metad, c, c1, c2, Ntrials);
        
        nR_S1{j} = sim.nR_S1;
        nR_S2{j} = sim.nR_S2;
        
        DATA(group).nR_S1{i} = sim.nR_S1;
        DATA(group).nR_S2{i} = sim.nR_S2;
        
        groupIndex = [groupIndex group-1];
        j = j+1;
    end
    
end

% Fit group 1
fit1 = fit_meta_d_mcmc_group(DATA(1).nR_S1, DATA(1).nR_S2);

% Fit group 2
fit2 = fit_meta_d_mcmc_group(DATA(2).nR_S1, DATA(2).nR_S2);

% Fit group data all at once to get difference in Mratio between tasks
fit = fit_meta_d_mcmc_group_twoGroups(nR_S1, nR_S2, groupIndex);

% Compute HDIs for two methods
sampleDiff = fit1.mcmc.samples.mu_Mratio(:) - fit2.mcmc.samples.mu_Mratio(:);
hdi_separate = calc_HDI(sampleDiff);
hdi_combined = calc_HDI(fit.mcmc.samples.mu_MratioG(:));

% Make some trace plots
figure;
subplot(1,2,1);
plot(fit.mcmc.samples.mu_Mratio');
xlabel('Sample');
ylabel('meta-d/d''');
box off

subplot(1,2,2);
plot(fit.mcmc.samples.mu_MratioG');
xlabel('Sample');
ylabel('Difference');
box off

% Calculate the HDI on the difference, test whether it overlaps zero
hdi = calc_HDI(fit.mcmc.samples.mu_MratioG(:));
fprintf(['\n HDI on difference = ', num2str(hdi) '\n\n'])
