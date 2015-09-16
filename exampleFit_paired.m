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

group_d = 2;    % same dprime across tasks
sigma = 0.5;
noise(1) = 0.2;
noise(2) = 0.5; % task 2 has worse metacognition than task 1 
% (in this simulation the noise terms are independent but the model assumes correlated noise, i.e. within-subject design)

for i = 1:Nsub
    
    for task = 1:2
        % Generate dprime
        d(i) = normrnd(group_d, sigma);
        
        % Generate data
        sim = type2_SDT_sim(d(i), noise(task), c, c1, c2, Ntrials);
        
        nR_S1{i}(task,:) = sim.nR_S1;
        nR_S2{i}(task,:) = sim.nR_S2;
        
    end
        
end

% Fit group data all at once to get difference in Mratio between tasks
fit = fit_meta_d_mcmc_group_paired(nR_S1, nR_S2);

% Make some trace plots
figure;
subplot(1,2,1);
plot(fit.mcmc.samples.mu_Mratio');
xlabel('Sample');
ylabel('meta-d/d');
box off

subplot(1,2,2);
plot(fit.mcmc.samples.mu_diff');
xlabel('Sample');
ylabel('Difference');
box off

% Calculate the HDI on the difference, test whether it overlaps zero
hdi = calc_HDI(fit.mcmc.samples.mu_diff(:));
fprintf(['\n HDI on difference = ', num2str(hdi) '\n\n'])
