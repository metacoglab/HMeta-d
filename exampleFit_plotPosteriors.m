% Example Bayesian meta-d fits with plots of posteriors
%
% SF 2014

Ntrials = 10000;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];
d = 2;
noise = 0;

% Generate data
sim = type2_SDT_sim(d, noise, c, c1, c2, Ntrials);
fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2);

% Plot posteriors and true values for criteria and meta-d'
h1 = figure(1);
maxSamp = [];
subplot(1,2,1);
for i = 1:length(fit.t2ca_rS1)
    [n x] = hist(fit.mcmc.samples.cS1(1,:,i));
    bar(x, n, 'edgecolor','b','facecolor',[1 1 1]);
    hold on
    maxSamp = [maxSamp max(n)];
end
for i = 1:length(fit.t2ca_rS2)
    [n x] = hist(fit.mcmc.samples.cS2(1,:,i));
    bar(x, n, 'edgecolor','b','facecolor',[1 1 1]);
    hold on
    maxSamp = [maxSamp max(n)];
    
end
for i = 1:length(c1)
    line([c1(i) c1(i)], [0 max(maxSamp)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
end
for i = 1:length(c1)
    line([c2(i) c2(i)], [0 max(maxSamp)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
end
line([c c], [0 max(maxSamp)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);

subplot(1,2,2);
[n x] = hist(fit.mcmc.samples.meta_d(1,:));
bar(x, n, 'edgecolor','b','facecolor',[1 1 1]);
hold on
line([d d], [0 max(n)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);

Ntrials = 10000;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];
d = 2;
noise = 0.5;

% Generate data
sim = type2_SDT_sim(d, noise, c, c1, c2, Ntrials);
fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2);

% Plot posteriors and true values for criteria and meta-d'
h1 = figure(1);
maxSamp = [];
subplot(1,2,1);
for i = 1:length(fit.t2ca_rS1)
    [n x] = hist(fit.mcmc.samples.cS1(1,:,i));
    bar(x, n, 'edgecolor','r','facecolor',[1 1 1]);
    hold on
    maxSamp = [maxSamp max(n)];
end
for i = 1:length(fit.t2ca_rS2)
    [n x] = hist(fit.mcmc.samples.cS2(1,:,i));
    bar(x, n, 'edgecolor','r','facecolor',[1 1 1]);
    hold on
    maxSamp = [maxSamp max(n)];
    
end
for i = 1:length(c1)
    line([c1(i) c1(i)], [0 max(maxSamp)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
end
for i = 1:length(c1)
    line([c2(i) c2(i)], [0 max(maxSamp)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);
end
line([c c], [0 max(maxSamp)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);

subplot(1,2,2);
[n x] = hist(fit.mcmc.samples.meta_d(1,:));
bar(x, n, 'edgecolor','r','facecolor',[1 1 1]);
hold on
line([d d], [0 max(n)], 'color', 'k', 'linestyle', '--', 'linewidth', 1.5);