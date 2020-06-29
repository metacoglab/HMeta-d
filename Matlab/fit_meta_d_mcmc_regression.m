function fit = fit_meta_d_mcmc_regression(nR_S1, nR_S2, cov, mcmc_params, fncdf, fninv, name)
% HMeta-d for between-subjects regression on meta-d'/d'
% cov is a n x s matrix of covariates, where s=number of subjects, n=number
% of covariates 
% See fit_meta_d_mcmc_group for full details
%
% Steve Fleming 2017
%
% 16.09.2019
% Ofaull added options to specify models for multiple covariates (up to 5)

fprintf('\n')
disp('----------------------------------------')
disp('Hierarchical meta-d'' regression model')
disp('https://github.com/smfleming/HMeta-d')
disp('----------------------------------------')
fprintf('\n')

cwd = pwd;

% Select model file and parameters to monitor
if size(cov, 1) == 1
    model_file = 'Bayes_metad_group_regress_nodp.txt';
    monitorparams = {'d1', 'c1', 'mu_logMratio', 'sigma_logMratio', 'mu_c2', 'sigma_c2', 'mu_beta1', 'Mratio', 'cS1', 'cS2'};
elseif size(cov, 1) == 2
    model_file = 'Bayes_metad_group_regress_nodp_2cov.txt';
    monitorparams = {'d1', 'c1', 'mu_logMratio', 'sigma_logMratio', 'mu_c2', 'sigma_c2', 'mu_beta1', 'mu_beta2', 'Mratio', 'cS1', 'cS2'};
elseif size(cov, 1) == 3
    model_file = 'Bayes_metad_group_regress_nodp_3cov.txt';
    monitorparams = {'d1', 'c1', 'mu_logMratio', 'sigma_logMratio', 'mu_c2', 'sigma_c2', 'mu_beta1', 'mu_beta2', 'mu_beta3', 'Mratio', 'cS1', 'cS2'};
elseif size(cov, 1) == 4
    model_file = 'Bayes_metad_group_regress_nodp_4cov.txt';
    monitorparams = {'d1', 'c1', 'mu_logMratio', 'sigma_logMratio', 'mu_c2', 'sigma_c2', 'mu_beta1', 'mu_beta2', 'mu_beta3', 'mu_beta4', 'Mratio', 'cS1', 'cS2'};
elseif size(cov, 1) == 5
    model_file = 'Bayes_metad_group_regress_nodp_5cov.txt';
    monitorparams = {'d1', 'c1', 'mu_logMratio', 'sigma_logMratio', 'mu_c2', 'sigma_c2', 'mu_beta1', 'mu_beta2', 'mu_beta3', 'mu_beta4', 'mu_beta5', 'Mratio', 'cS1', 'cS2'};
else
    error('Too many covariates specified: Max = 5')
end

findpath = which(model_file);
if isempty(findpath)
    error('Please add HMetaD directory to the path')
else
    hmmPath = fileparts(findpath);
    cd(hmmPath)
end

if ~exist('fncdf','var') || isempty(fncdf)
    fncdf = @normcdf;
end

if ~exist('fninv','var') || isempty(fninv)
    fninv = @norminv;
end

if ~exist('name','var') || isempty(name)
    tmpfolder = 'tmpjags';
else
    tmpfolder = name;
    mkdir(tmpfolder);
end

Nsubj = length(nR_S1);
nRatings = length(nR_S1{1})/2;

for n = 1:Nsubj
    if length(nR_S1{n}) ~= nRatings*2 || length(nR_S2{n}) ~= nRatings*2
        error('Subjects do not have equal numbers of response categories');
    end
    % Get type 1 SDT parameter values
    counts(n,:) = [nR_S1{n} nR_S2{n}];
    nTot(n) = sum(counts(n,:));
    % Adjust to ensure non-zero counts for type 1 d' point estimate (not
    % necessary if estimating d' inside JAGS)
    adj_f = 1/length(nR_S1{n});
    nR_S1_adj = nR_S1{n} + adj_f;
    nR_S2_adj = nR_S2{n} + adj_f;
    
    ratingHR  = [];
    ratingFAR = [];
    for c = 2:nRatings*2
        ratingHR(end+1) = sum(nR_S2_adj(c:end)) / sum(nR_S2_adj);
        ratingFAR(end+1) = sum(nR_S1_adj(c:end)) / sum(nR_S1_adj);
    end
    
    t1_index = nRatings;
    
    d1(n) = fninv(ratingHR(t1_index)) - fninv(ratingFAR(t1_index));
    c1(n) = -0.5 .* (fninv(ratingHR(t1_index)) + fninv(ratingFAR(t1_index)));
end

%% Sampling
if ~exist('mcmc_params','var') || isempty(mcmc_params)
    % MCMC Parameters
    mcmc_params.response_conditional = 0;   % response-conditional meta-d?
    mcmc_params.estimate_dprime = 0;    % also estimate dprime in same model?
    mcmc_params.nchains = 3; % How Many Chains?
    mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
    mcmc_params.nsamples = 10000;  %How Many Recorded Samples?
    mcmc_params.nthin = 1; % How Often is a Sample Recorded?
    mcmc_params.doparallel = 0; % Parallel Option
    mcmc_params.dic = 1;
end
% Ensure init0 is correct size
if ~isfield(mcmc_params, 'init0')
    for i=1:mcmc_params.nchains
        mcmc_params.init0(i) = struct;
    end
end

datastruct = struct('d1', d1, 'c1', c1,'nsubj',Nsubj,'counts', counts,'cov', cov, 'nratings', nRatings, 'nTot', nTot, 'Tol', 1e-05);

% Use JAGS to Sample
try
    tic
    fprintf( 'Running JAGS ...\n' );
    [samples, stats] = matjags( ...
        datastruct, ...
        fullfile(pwd, model_file), ...
        mcmc_params.init0, ...
        'doparallel' , mcmc_params.doparallel, ...
        'nchains', mcmc_params.nchains,...
        'nburnin', mcmc_params.nburnin,...
        'nsamples', mcmc_params.nsamples, ...
        'thin', mcmc_params.nthin, ...
        'dic', mcmc_params.dic,...
        'monitorparams', monitorparams, ...
        'savejagsoutput' , 0 , ...
        'verbosity' , 1 , ...
        'cleanup' , 1 , ...
        'workingdir' , tmpfolder );
    toc
catch ME
    % Remove temporary directory if specified
    if exist('name','var')
        if exist(['../', tmpfolder],'dir')
            rmdir(['../', tmpfolder], 's');
        end
    end
    % Print the error message
    rethrow(ME);
end

% Remove temporary directory if specified
if exist('name','var')
    if exist(tmpfolder,'dir')
        rmdir(tmpfolder, 's');
    end
end

% Package group-level output
if isrow(stats.mean.cS1)
    stats.mean.cS1 = stats.mean.cS1';
    stats.mean.cS2 = stats.mean.cS2';
end
fit.t2ca_rS1  = stats.mean.cS1;
fit.t2ca_rS2  = stats.mean.cS2;
fit.d1 = stats.mean.d1;
fit.c1 = stats.mean.c1;

fit.mu_logMratio = stats.mean.mu_logMratio;
fit.sigma_logMratio = stats.mean.sigma_logMratio;
fit.mu_beta1 = stats.mean.mu_beta1;
if size(cov, 1) > 1
    fit.mu_beta2 = stats.mean.mu_beta2;
end
if size(cov, 1) > 2
    fit.mu_beta3 = stats.mean.mu_beta3;
end
if size(cov, 1) > 3
    fit.mu_beta4 = stats.mean.mu_beta4;
end
if size(cov, 1) > 4
    fit.mu_beta5 = stats.mean.mu_beta5;
end

fit.Mratio = stats.mean.Mratio;
fit.meta_d   = fit.Mratio.*stats.mean.d1;

fit.mcmc.dic = stats.dic;
fit.mcmc.Rhat = stats.Rhat;
fit.mcmc.samples = samples;
fit.mcmc.params = mcmc_params;

for n = 1:Nsubj
    
    
    %% Data is fit, now package output
    I_nR_rS2 = nR_S1{n}(nRatings+1:end);
    I_nR_rS1 = nR_S2{n}(nRatings:-1:1);
    
    C_nR_rS2 = nR_S2{n}(nRatings+1:end);
    C_nR_rS1 = nR_S1{n}(nRatings:-1:1);
    
    for i = 2:nRatings
        obs_FAR2_rS2(i-1) = sum( I_nR_rS2(i:end) ) / sum(I_nR_rS2);
        obs_HR2_rS2(i-1)  = sum( C_nR_rS2(i:end) ) / sum(C_nR_rS2);
        
        obs_FAR2_rS1(i-1) = sum( I_nR_rS1(i:end) ) / sum(I_nR_rS1);
        obs_HR2_rS1(i-1)  = sum( C_nR_rS1(i:end) ) / sum(C_nR_rS1);
    end
    
    
    % Calculate fits based on either vanilla or response-conditional model
    s = 1;
    %% find estimated t2FAR and t2HR
    meta_d = fit.meta_d(n);
    S1mu = -meta_d/2; S1sd = 1;
    S2mu =  meta_d/2; S2sd = S1sd/s;
    
    C_area_rS2 = 1-fncdf(fit.c1(n),S2mu,S2sd);
    I_area_rS2 = 1-fncdf(fit.c1(n),S1mu,S1sd);
    
    C_area_rS1 = fncdf(fit.c1(n),S1mu,S1sd);
    I_area_rS1 = fncdf(fit.c1(n),S2mu,S2sd);
    
    t2c1 = [fit.t2ca_rS1(n,:) fit.t2ca_rS2(n,:)];
    
    for i=1:nRatings-1
        
        t2c1_lower = t2c1(nRatings-i);
        t2c1_upper = t2c1(nRatings-1+i);
        
        I_FAR_area_rS2 = 1-fncdf(t2c1_upper,S1mu,S1sd);
        C_HR_area_rS2  = 1-fncdf(t2c1_upper,S2mu,S2sd);
        
        I_FAR_area_rS1 = fncdf(t2c1_lower,S2mu,S2sd);
        C_HR_area_rS1  = fncdf(t2c1_lower,S1mu,S1sd);
        
        
        est_FAR2_rS2(i) = I_FAR_area_rS2 / I_area_rS2;
        est_HR2_rS2(i)  = C_HR_area_rS2 / C_area_rS2;
        
        est_FAR2_rS1(i) = I_FAR_area_rS1 / I_area_rS1;
        est_HR2_rS1(i)  = C_HR_area_rS1 / C_area_rS1;
        
    end
    
    fit.est_HR2_rS1(n,:)  = est_HR2_rS1;
    fit.obs_HR2_rS1(n,:)  = obs_HR2_rS1;
    
    fit.est_FAR2_rS1(n,:) = est_FAR2_rS1;
    fit.obs_FAR2_rS1(n,:) = obs_FAR2_rS1;
    
    fit.est_HR2_rS2(n,:)  = est_HR2_rS2;
    fit.obs_HR2_rS2(n,:)  = obs_HR2_rS2;
    
    fit.est_FAR2_rS2(n,:) = est_FAR2_rS2;
    fit.obs_FAR2_rS2(n,:) = obs_FAR2_rS2;
end
cd(cwd);
