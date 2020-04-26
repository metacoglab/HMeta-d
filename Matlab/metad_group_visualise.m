%% Generate mean and 95% CI for group-level ROC

Nsub = length(fit.d1);
ts = tinv([0.05/2,  1-0.05/2],Nsub-1);

if any(isnan(fit.obs_FAR2_rS1(:))) || any(isnan(fit.obs_HR2_rS1(:))) || any(isnan(fit.obs_FAR2_rS2(:))) || any(isnan(fit.obs_HR2_rS2(:)))
    warning('One or more subjects have NaN entries for observed confidence rating counts; these will be omitted from the plot')
end

mean_obs_FAR2_rS1 = nanmean(fit.obs_FAR2_rS1);
mean_obs_HR2_rS1 = nanmean(fit.obs_HR2_rS1);
mean_obs_FAR2_rS2 = nanmean(fit.obs_FAR2_rS2);
mean_obs_HR2_rS2 = nanmean(fit.obs_HR2_rS2);

CI_obs_FAR2_rS1(1,:) = ts(1).*(nanstd(fit.obs_FAR2_rS1)./sqrt(Nsub));
CI_obs_FAR2_rS1(2,:) = ts(2).*(nanstd(fit.obs_FAR2_rS1)./sqrt(Nsub));
CI_obs_HR2_rS1(1,:) = ts(1).*(nanstd(fit.obs_HR2_rS1)./sqrt(Nsub));
CI_obs_HR2_rS1(2,:) = ts(2).*(nanstd(fit.obs_HR2_rS1)./sqrt(Nsub));
CI_obs_FAR2_rS2(1,:) = ts(1).*(nanstd(fit.obs_FAR2_rS2)./sqrt(Nsub));
CI_obs_FAR2_rS2(2,:) = ts(2).*(nanstd(fit.obs_FAR2_rS2)./sqrt(Nsub));
CI_obs_HR2_rS2(1,:) = ts(1).*(nanstd(fit.obs_HR2_rS2)./sqrt(Nsub));
CI_obs_HR2_rS2(2,:) = ts(2).*(nanstd(fit.obs_HR2_rS2)./sqrt(Nsub));

mean_est_FAR2_rS1 = nanmean(fit.est_FAR2_rS1);
mean_est_HR2_rS1 = nanmean(fit.est_HR2_rS1);
mean_est_FAR2_rS2 = nanmean(fit.est_FAR2_rS2);
mean_est_HR2_rS2 = nanmean(fit.est_HR2_rS2);

CI_est_FAR2_rS1(1,:) = ts(1).*(nanstd(fit.est_FAR2_rS1)./sqrt(Nsub));
CI_est_FAR2_rS1(2,:) = ts(2).*(nanstd(fit.est_FAR2_rS1)./sqrt(Nsub));
CI_est_HR2_rS1(1,:) = ts(1).*(nanstd(fit.est_HR2_rS1)./sqrt(Nsub));
CI_est_HR2_rS1(2,:) = ts(2).*(nanstd(fit.est_HR2_rS1)./sqrt(Nsub));
CI_est_FAR2_rS2(1,:) = ts(1).*(nanstd(fit.est_FAR2_rS2)./sqrt(Nsub));
CI_est_FAR2_rS2(2,:) = ts(2).*(nanstd(fit.est_FAR2_rS2)./sqrt(Nsub));
CI_est_HR2_rS2(1,:) = ts(1).*(nanstd(fit.est_HR2_rS2)./sqrt(Nsub));
CI_est_HR2_rS2(2,:) = ts(2).*(nanstd(fit.est_HR2_rS2)./sqrt(Nsub));

%% Observed and expected type 2 ROCs for S1 and S2 responses
h1 = figure(1);
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.2 0.2 0.5 0.5]);

subplot(1,2,1);
errorbar([1 mean_obs_FAR2_rS1 0], [1 mean_obs_HR2_rS1 0], [0 CI_obs_HR2_rS1(1,:) 0], [0 CI_obs_HR2_rS1(2,:) 0], ...
    [0 CI_obs_FAR2_rS1(1,:) 0], [0 CI_obs_FAR2_rS1(2,:) 0], 'ko-','linewidth',1.5,'markersize', 12);
hold on
errorbar([1 mean_est_FAR2_rS1 0], [1 mean_est_HR2_rS1 0], [0 CI_est_HR2_rS1(1,:) 0], [0 CI_est_HR2_rS1(2,:) 0], ...
    [0 CI_est_FAR2_rS1(1,:) 0], [0 CI_est_FAR2_rS1(2,:) 0], 'd-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
ylabel('HR2');
xlabel('FAR2');
line([0 1],[0 1],'linestyle','--','color','k');
axis square
box off
legend('Data', 'Model', 'Location', 'SouthEast')
legend boxoff
title('Response = S1')

subplot(1,2,2);
errorbar([1 mean_obs_FAR2_rS2 0], [1 mean_obs_HR2_rS2 0], [0 CI_obs_HR2_rS2(1,:) 0], [0 CI_obs_HR2_rS2(2,:) 0], ...
    [0 CI_obs_FAR2_rS2(1,:) 0], [0 CI_obs_FAR2_rS2(2,:) 0], 'ko-','linewidth',1.5,'markersize', 12);
hold on
errorbar([1 mean_est_FAR2_rS2 0], [1 mean_est_HR2_rS2 0], [0 CI_est_HR2_rS2(1,:) 0], [0 CI_est_HR2_rS2(2,:) 0], ...
    [0 CI_est_FAR2_rS2(1,:) 0], [0 CI_est_FAR2_rS2(2,:) 0], 'd-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
ylabel('HR2');
xlabel('FAR2');
line([0 1],[0 1],'linestyle','--','color','k');
axis square
box off
legend('Data', 'Model', 'Location', 'SouthEast')
legend boxoff
title('Response = S2')
