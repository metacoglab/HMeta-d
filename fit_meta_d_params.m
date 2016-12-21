function mcmc_params = fit_meta_d_params()
% mcmc_params = fit_meta_d_params()
%
% Returns default mcmc_params for fit_meta_d_mcmc_group and fit_meta_d_mcmc
%
% SF 2015

mcmc_params.response_conditional = 0;
mcmc_params.estimate_dprime = 0;
mcmc_params.nchains = 3; % How Many Chains?
mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
mcmc_params.nsamples = 10000;  %How Many Recorded Samples?
mcmc_params.nthin = 1; % How Often is a Sample Recorded?
mcmc_params.doparallel = 0; % Parallel Option
mcmc_params.dic = 1;
%% These lines are now part of main functions to avoid issues when changing chain number
% for i=1:mcmc_params.nchains
%     mcmc_params.init0(i) = struct;
% end