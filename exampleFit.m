% Example Bayesian meta-d fit
%
% SF 2014

Ntrials = 10000;
c = 0;
c1 = [-1.5 -1 -0.5];
c2 = [0.5 1 1.5];
d = 2;
noise = 0.5;

% Generate data
sim = type2_SDT_sim(d, noise, c, c1, c2, Ntrials);

% Fit the data
mcmc_params.response_conditional = 0;
mcmc_params.nchains = 3; % How Many Chains?
mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
mcmc_params.nsamples = 10000;  %How Many Recorded Samples?
mcmc_params.nthin = 1; % How Often is a Sample Recorded?
mcmc_params.doparallel = 0; % Parallel Option
mcmc_params.dic = 1;
init_d = [1.5 2 2.5];
init_cs1 = [linspace(-1,-0.2,length(c1)); linspace(-2,-0.5,length(c1)); linspace(-1.5,-0.3,length(c1))];
init_cs2 = [linspace(0.2,1,length(c1)); linspace(0.5,2,length(c1)); linspace(0.3,1.5,length(c1))];
for i=1:mcmc_params.nchains
    S.meta_d = init_d(i);
    S.cS1_raw = init_cs1(i,:);
    S.cS2_raw = init_cs2(i,:);
    mcmc_params.init0(i) = S;
end

fit = fit_meta_d_mcmc(sim.nR_S1, sim.nR_S2, mcmc_params)

% Visualise fits
metad_visualise