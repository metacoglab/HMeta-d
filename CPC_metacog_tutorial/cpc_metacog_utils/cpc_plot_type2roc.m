%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% CPC METACOGNITION TUTORIAL 2019: PLOT TYPE2ROC %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to plot type2 ROC curve from observed ± estimated data fit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cpc_plot_type2roc(data, titleText, toPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(toPlot, 'obs') == 1
    nRatings = length(data.responses.nR_S1) / 2;
    % Find incorrect observed ratings...
    I_nR_rS2 = data.responses.nR_S1(nRatings+1:end);
    I_nR_rS1 = data.responses.nR_S2(nRatings:-1:1);
    I_nR = I_nR_rS2 + I_nR_rS1;
    % Find correct observed ratings...
    C_nR_rS2 = data.responses.nR_S2(nRatings+1:end);
    C_nR_rS1 = data.responses.nR_S1(nRatings:-1:1);
    C_nR = C_nR_rS2 + C_nR_rS1;
    % Calculate type 2 hits and false alarms
    for i = 1:nRatings
        obs_FAR2_rS2(i) = sum( I_nR_rS2(i:end) ) / sum(I_nR_rS2);
        obs_HR2_rS2(i)  = sum( C_nR_rS2(i:end) ) / sum(C_nR_rS2);
        obs_FAR2_rS1(i) = sum( I_nR_rS1(i:end) ) / sum(I_nR_rS1);
        obs_HR2_rS1(i)  = sum( C_nR_rS1(i:end) ) / sum(C_nR_rS1);
        obs_FAR2(i) = sum( I_nR(i:end) ) / sum(I_nR);
        obs_HR2(i)  = sum( C_nR(i:end) ) / sum(C_nR);
    end
    obs_FAR2(nRatings+1) = 0;
    obs_HR2(nRatings+1)  = 0;
    plot(obs_FAR2, obs_HR2, 'ko-', 'linewidth', 1.5, 'markersize', 12);
    hold on
    legend('Observed','Location','SouthEast');
    title(['TYPE2 ROC CURVE: ', titleText]);
    ylabel('TYPE2 P(CORRECT)');
    xlabel('TYPE2 P(INCORRECT)');
elseif strcmp(toPlot, 'est') == 1
    figure
    set(gcf, 'Units', 'normalized');
    set(gcf, 'Position', [0.2 0.2 0.5 0.43]);
    subplot(1,2,1);
    plot(data.fit.obs_FAR2_rS1, data.fit.obs_HR2_rS1, 'ko-','linewidth',1.5,'markersize',12);
    hold on
    plot(data.fit.est_FAR2_rS1, data.fit.est_HR2_rS1, '+-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
    legend('Observed','Estimated','Location','SouthEast');
    text(0.5, 1.15, [titleText, ': ', sprintf('fit meta-d = %.2f\n', data.fit.meta_d)], 'FontSize', 25, 'FontWeight', 'bold');
    title('STIMULUS 1');
    ylabel('TYPE2 HR');
    xlabel('TYPE2 FAR');
end   
set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
line([0 1],[0 1],'linestyle','--','color','k','HandleVisibility','off');
axis square

if strcmp(toPlot, 'est') == 1
    subplot(1,2,2);
    if strcmp(toPlot, 'obs') == 1
        plot(obs_FAR2_rS2, obs_HR2_rS2, 'ko-', 'linewidth', 1.5, 'markersize', 12);
        hold on
        legend('Observed','Location','SouthEast');
    elseif strcmp(toPlot, 'est') == 1
        plot(data.fit.obs_FAR2_rS2, data.fit.obs_HR2_rS2, 'ko-','linewidth',1.5,'markersize',12);
        hold on
        plot(data.fit.est_FAR2_rS2, data.fit.est_HR2_rS2, '+-','color',[0.5 0.5 0.5], 'linewidth',1.5,'markersize',10);
        legend('Observed','Estimated','Location','SouthEast');
    end   
    set(gca, 'XLim', [0 1], 'YLim', [0 1], 'FontSize', 16);
    ylabel('TYPE2 HR');
    xlabel('TYPE2 FAR');
    line([0 1],[0 1],'linestyle','--','color','k','HandleVisibility','off');
    title('STIMULUS 2');
    axis square

    hold off
end

end