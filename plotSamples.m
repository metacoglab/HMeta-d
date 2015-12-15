function plotSamples(samples)
% function out = plotSamples(samples)
%
% Plots chains and histograms of samples to check mixing
%
% INPUTS
%
% samples - vector of MCMC samples
%
% Steve Fleming 2015 stephen.fleming@ucl.ac.uk

figure;
set(gcf, 'Position', [200 200 800 300])

subplot(1,2,1)
plot(samples');
xlabel('Sample');
ylabel('Parameter');
box off
set(gca, 'FontSize', 14)

subplot(1,2,2)
histogram(samples(:))
xlabel('Parameter');
ylabel('Sample count');
box off
set(gca, 'FontSize', 14)

