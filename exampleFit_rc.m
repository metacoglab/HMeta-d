% Comparison of vanilla and response-conditional meta-d' model fit to
% simulated data with response-conditional noise

clear all

Ntrials = 10000;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];

d = 2;
noise = [1/4 3/4];

% Generate data (passing in a 2-vector of noise terms generates
% differential RC meta-d')
sim = type2_SDT_sim(d, noise, c, c1, c2, Ntrials);

% Fit the data using vanilla model
mcmc_params = fit_meta_d_params;
vanilla_fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2, mcmc_params);

% Fit the data using response-conditional model
mcmc_params = fit_meta_d_params;
mcmc_params.response_conditional = 1;
rc_fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2, mcmc_params);

% Visualise fits of both models
h1 = figure;
set(gcf, 'Position', [500 500 1000 500]);

subplot(1,2,1);
plot(vanilla_fit.obs_FAR2_rS1, vanilla_fit.obs_HR2_rS1, 'ko-','linewidth',1.5,'markersize',12);
hold on
plot(vanilla_fit.est_FAR2_rS1, vanilla_fit.est_HR2_rS1, '+-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
plot(rc_fit.est_FAR2_rS1, rc_fit.est_HR2_rS1, '+--','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
ylabel('HR2, "S1"');
xlabel('FAR2, "S1"');
line([0 1],[0 1],'linestyle','--','color','k');
axis square
box off

subplot(1,2,2);
plot(vanilla_fit.obs_FAR2_rS2, vanilla_fit.obs_HR2_rS2, 'ko-','linewidth',1.5,'markersize',12);
hold on
plot(vanilla_fit.est_FAR2_rS2, vanilla_fit.est_HR2_rS2, '+-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
plot(rc_fit.est_FAR2_rS2, rc_fit.est_HR2_rS2, '+--','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
ylabel('HR2, "S2"');
xlabel('FAR2, "S2"');
line([0 1],[0 1],'linestyle','--','color','k');
axis square
box off
legend('Data','Vanilla fit','RC fit','Location','SouthEast');