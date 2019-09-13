% Shows underlying SDT model and relationship to multinomial probabilities
% for fitted type 2 data
%
% Requires fit, nR_S1 and nR_S2 objects from either subject-level HMM or MLE meta-d fit to be in the
% workspace
%
% SF 2014

figure;
set(gcf, 'Units', 'normalized');
set(gcf, 'Position', [0.2 0.2 0.5 0.3]);
base = linspace(-4,4,500);

subplot(1,2,1);
mu1 = fit.meta_d./2;
S = normpdf(base, mu1, 1);
N = normpdf(base, -mu1, 1);
plot(base, S, 'k', 'LineWidth', 2);
hold on
plot(base, N, 'k--', 'LineWidth', 2);
set(gca, 'YLim', [0 0.5], 'FontSize',12);
line([fit.c1 fit.c1],[0 0.5],'Color','k','LineWidth',1);
for i = 1:length(fit.t2ca_rS1)
    line([fit.t2ca_rS1(i) fit.t2ca_rS1(i)],[0 0.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
    line([fit.t2ca_rS2(i) fit.t2ca_rS2(i)],[0 0.5], 'Color', 'k', 'LineStyle', '--', 'LineWidth', 1);
end
xlabel('X','FontSize', 14)
ylabel('f(x|S)', 'FontSize', 14);

% Show model fit in terms of proportion total confidence ratings, collapse
% over response
mean_c2 = (fit.t2ca_rS2 + abs(fit.t2ca_rS1(end:-1:1)))./2;
I_area = 1-normcdf(0,-mu1,1);
C_area = 1-normcdf(0,mu1,1);
allC = [0 mean_c2 Inf];
for i = 1:length(allC)-1
    I_prop(i) = (normcdf(allC(i+1), -mu1, 1) - normcdf(allC(i), -mu1, 1))./I_area;
    C_prop(i) = (normcdf(allC(i+1), mu1, 1) - normcdf(allC(i), mu1, 1))./C_area;
end

subplot(1,2,2);
modelProp = [C_prop(end:-1:1) I_prop];
obsCount = (nR_S1 + nR_S2(end:-1:1));
Nrating = length(nR_S1)./2;
obsProp(1:Nrating) = obsCount(1:Nrating)./sum(obsCount(1:Nrating));
obsProp(Nrating+1:2*Nrating) = obsCount(Nrating+1:2*Nrating)./sum(obsCount(Nrating+1:2*Nrating));
bar(obsProp);
hold on
plot(1:length(modelProp), modelProp, 'ro ', 'MarkerSize', 10, 'LineWidth', 2);
set(gca, 'YLim', [0 1], 'XTick', [1:8], 'XTickLabel', {'4','3','2','1','1','2','3','4'},'FontSize',12);
line([4.5 4.5], [0 1], 'Color', 'k', 'LineStyle', '--');
ylabel('P(conf = y) | outcome)','FontSize',14);
xlabel('Confidence rating','FontSize',14);
text(8, 0.7, 'Correct','FontSize',14)
text(1, 0.7, 'Error','FontSize',14);