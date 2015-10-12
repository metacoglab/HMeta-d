% Demonstration of hierarchical model fits
%
% SF 2014

clear all
close all

Ntrials = 100;
Nsub = 20;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

group_d = 2;    % same mean dprime across tasks
sigma = 0.5;
group_noise = 0.5;
r = 0.8;        % correlation coefficient

for i = 1:Nsub
    
    % Generate dprime
    d(i) = normrnd(group_d, sigma);
    
    % Generate correlated noise terms (metacognitive sensitivity)
    noise(i,:) = mvnrnd([group_noise group_noise], [sigma^2 r*sigma^2; r*sigma^2 sigma^2]);
    
    for task = 1:2
        % Generate data
        sim = type2_SDT_sim(d(i), noise(i,task), c, c1, c2, Ntrials);
        
        nR_S1{i}(task,:) = sim.nR_S1;
        nR_S2{i}(task,:) = sim.nR_S2;
        
    end
end

% Fit group data all at once to get difference in Mratio between tasks
fit = fit_meta_d_mcmc_group_correl(nR_S1, nR_S2);

% Make some trace plots
figure;
plot(fit.mcmc.samples.r');
xlabel('Sample');
ylabel('r');
box off

% Calculate the HDI on correlation coefficient
hdi = calc_HDI(fit.mcmc.samples.r(:));
fprintf(['\n HDI on correlation coefficient = ', num2str(hdi) '\n\n'])
