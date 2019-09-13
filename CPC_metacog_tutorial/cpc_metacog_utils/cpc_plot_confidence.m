%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% CPC METACOGNITION TUTORIAL 2019: PLOT CONFIDENCE %%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to plot confidence scores from simulated and fit data

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function cpc_plot_confidence(data, titleText, toPlot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Show model fit in terms of conditional probability of confidence ratings, collapse
% over responses
if strcmp(toPlot, 'est') == 1
    if isfield(data.fit, 'meta_d')
        meta_d = data.fit.meta_d;
    elseif isfield(data.fit, 'meta_da')
        meta_d = data.fit.meta_da;
    end
    mu1 = meta_d./2;
    mean_c2 = (data.fit.t2ca_rS2 + abs(data.fit.t2ca_rS1(end:-1:1)))./2;
    I_area = 1-normcdf(0,-mu1,1);
    C_area = 1-normcdf(0,mu1,1);
    allC = [0 mean_c2 Inf];
    for i = 1:length(allC)-1
        I_prop_model(i) = (normcdf(allC(i+1), -mu1, 1) - normcdf(allC(i), -mu1, 1))./I_area;
        C_prop_model(i) = (normcdf(allC(i+1), mu1, 1) - normcdf(allC(i), mu1, 1))./C_area;
    end
end

% Get the number of ratings
Nrating = length(data.responses.nR_S1) / 2;

% Collapse across two stimuli for correct and incorrect responses
obsCount = (data.responses.nR_S1 + data.responses.nR_S2(end:-1:1)); % this gives flipped corrects followed by incorrects
C_prop_data = fliplr(obsCount(1:Nrating))./sum(obsCount(1:Nrating));
I_prop_data = obsCount(Nrating+1:2*Nrating)./sum(obsCount(Nrating+1:2*Nrating));

% Plot responses
if strcmp(toPlot, 'obs') == 1
    bar([1:Nrating]-0.2, I_prop_data, 0.3, 'r', 'LineWidth', 2);
    hold on
    bar([1:Nrating]+0.2, C_prop_data, 0.3, 'g', 'LineWidth', 2);
    legend('Obs Incorrect','Obs Correct','Location','NorthEast');
    title(titleText);
elseif strcmp(toPlot, 'est') == 1
    bar([1:Nrating]-0.2, I_prop_data, 0.3, 'r', 'LineWidth', 2);
    hold on
    plot([1:Nrating]-0.2, I_prop_model, 'ro ', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor','r', 'MarkerFaceColor', [1 1 1]);
    bar([1:Nrating]+0.2, C_prop_data, 0.3, 'g', 'LineWidth', 2);
    plot([1:Nrating]+0.2, C_prop_model, 'go ', 'MarkerSize', 8, 'LineWidth', 2, 'MarkerEdgeColor','g', 'MarkerFaceColor', [1 1 1]);
    legend('Obs Incorrect','Est Incorrect','Correct','Est Correct','Location','NorthEast');
    title([titleText, ': ', sprintf('meta-d = %.2f\n', meta_d)]);
end
ylabel('P(conf = y) | outcome)');
xlabel('Confidence rating');
set(gca, 'FontSize', 14, 'XTick', [1:Nrating], 'YLim', [0 0.7])
box off

end
