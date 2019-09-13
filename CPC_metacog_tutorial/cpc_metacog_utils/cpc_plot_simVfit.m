%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CPC METACOGNITION TUTORIAL 2019: PLOT TYPE2ROC %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to plot type2 ROC curve from observed ± estimated data fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cpc_plot_simVfit(sim, fit, parameter, fitType)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear est

% Pull out parameter of interest
if strcmp(fitType, 'single') == 1
    for a = 1:length(sim.params.meta_d)
        for b = 1:sim.params.Nsims
            sim.meta_d(b,a) = sim.values{a}.meta_d(b);
            sim.Mratio(b,a) = sim.values{a}.meta_d(b) / sim.values{a}.d(b);
            est.meta_d(b,a) = fit{b,a}.meta_d;
            est.Mratio(b,a) = fit{b,a}.M_ratio;
        end
    end
elseif strcmp(fitType, 'group') == 1
    for a = 1:length(sim.params.meta_d)
        for b = 1:sim.params.Nsims
            sim.meta_d(b,a) = sim.values{a}.meta_d(b);
            sim.Mratio(b,a) = sim.values{a}.meta_d(b) / sim.values{a}.d(b);
            est.meta_d(b,a) = fit{a}.meta_d(b);
            est.Mratio(b,a) = fit{a}.Mratio(b);
        end
    end
end

figure
hold on
if strcmp(parameter, 'meta-d') == 1
    for a = 1:sim.params.Nsims
        scatter(sim.meta_d(a,:), est.meta_d(a,:), 'k');
    end
    xlabel('Simulated meta-d');
    ylabel('Recovered meta-d');
    title('META-D PARAMETER RECOVERY');
    rf1 = refline(1, 0);
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 1;
    rf2 = refline(0, 0);
    rf2.LineStyle = ':';
    rf2.Color = 'k';
    rf2.LineWidth = 1.5;
elseif strcmp(parameter, 'Mratio') == 1
    for a = 1:sim.params.Nsims
        scatter(sim.Mratio(a,:), est.Mratio(a,:), 'k');
    end
    xlabel('Simulated Mratio');
    ylabel('Recovered Mratio');
    title('MRATIO PARAMETER RECOVERY');
    rf1 = refline(1, 0);
    rf1.LineStyle = '--';
    rf1.Color = 'k';
    rf1.LineWidth = 1;
    rf2 = refline(0, 0);
    rf2.LineStyle = ':';
    rf2.Color = 'k';
    rf2.LineWidth = 1.5;
elseif strcmp(fitType, 'regression') == 1 && strcmp(parameter, 'log(Mratio)') == 1
    S1 = scatter(sim.values.covZscored, log(sim.fit.bayesGroupMean.Mratio));
    refline
    S2 = scatter(sim.values.covZscored, log(sim.fit.bayesGroupRegression.Mratio));
    refline
    xlabel('Covariate');
    ylabel('Recovered log(Mratio)');
    title('REGRESSION COMPARISON');
    legend([S1 S2], 'HMeta-d','RHMeta-d','Location','SouthEast');
end

set(gca, 'FontSize', 12)

end