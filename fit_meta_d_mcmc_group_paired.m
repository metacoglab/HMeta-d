function fit = fit_meta_d_mcmc_group_paired(nR_S1, nR_S2, mcmc_params, fncdf, fninv)
% fit = fit_meta_d_mcmc_paired(nR_S1, nR_S2, mcmc_params, s, fncdf, fninv)
%
% Estimate difference between two conditions in within-subject design
%
% nR_S1 and nR_S2 must be cell arrays of 2xNratings, where each row
% corresponds to a different condition and each cell to a different subject
%
% See fit_meta_d_mcmc_group for further details about input and output

if ~exist('fncdf','var') || isempty(fncdf)
    fncdf = @normcdf;
end

if ~exist('fninv','var') || isempty(fninv)
    fninv = @norminv;
end

Nsubj = length(nR_S1);
nRatings = length(nR_S1{1})/2;

for n = 1:Nsubj
    
    if length(nR_S1{n}) ~= nRatings*2 || length(nR_S2{n}) ~= nRatings*2
        error('Subjects do not have equal numbers of response categories');
    end
    
    for task = 1:size(nR_S1{n},1)
        % Organise data
        counts(n,:,task) = [nR_S1{n}(task,:) nR_S2{n}(task,:)];
        nTot(n,task) = sum(counts(n,:,task));
        
        % Adjust to ensure non-zero counts for type 1 d' point estimate (not
        % necessary if estimating d' inside JAGS)
        adj_f = 1/length(nR_S1{n}(task,:));
        nR_S1_adj = nR_S1{n}(task,:) + adj_f;
        nR_S2_adj = nR_S2{n}(task,:) + adj_f;
        
        ratingHR  = [];
        ratingFAR = [];
        for c = 2:nRatings*2
            ratingHR(end+1) = sum(nR_S2_adj(c:end)) / sum(nR_S2_adj);
            ratingFAR(end+1) = sum(nR_S1_adj(c:end)) / sum(nR_S1_adj);
        end

        t1_index = nRatings;

        d1(n,task) = fninv(ratingHR(t1_index)) - fninv(ratingFAR(t1_index));
        c1(n,task) = fninv(ratingHR(t1_index)) + fninv(ratingFAR(t1_index));
    end
    
end

%% Sampling
if ~exist('mcmc_params','var') || isempty(mcmc_params)
    % MCMC Parameters
    mcmc_params.response_conditional = 0;
    mcmc_params.estimate_dprime = 0;    % also estimate dprime in same model?
    mcmc_params.nchains = 3; % How Many Chains?
    mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
    mcmc_params.nsamples = 10000;  %How Many Recorded Samples?
    mcmc_params.nthin = 1; % How Often is a Sample Recorded?
    mcmc_params.doparallel = 0; % Parallel Option
    mcmc_params.dic = 1;
    % Initialize Unobserved Variables
    for i=1:mcmc_params.nchains
        S.mu_Mratio = 1;
        mcmc_params.init0(i) = S;
    end
end
% Assign variables to the observed nodes
switch mcmc_params.estimate_dprime
    case 1
        datastruct = struct('nsubj',Nsubj,'counts', counts, 'nratings', nRatings, 'nTot', nTot, 'Tol', 1e-05);
    case 0
        datastruct = struct('d1', d1, 'c1', c1, 'nsubj',Nsubj,'counts', counts, 'nratings', nRatings, 'nTot', nTot, 'Tol', 1e-05);
end

% Select model file and parameters to monitor

model_file = 'Bayes_metad_group_paired.txt';
monitorparams = {'d1', 'c', 'mu_Mratio','lambda_Mratio','mu_diff','lambda_diff','Mratio','diff','cS1','cS2'};

% Use JAGS to Sample
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
    'workingdir' , 'tmpjags' );
toc

% Package group-level output
fit.d1 = stats.mean.d1;
fit.c1 = stats.mean.c;
fit.mu_diff = stats.mean.mu_diff;
fit.mu_Mratio = stats.mean.mu_Mratio;
fit.diff = stats.mean.mu_diff;
fit.lambda_Mratio = stats.mean.lambda_Mratio;
fit.lambda_diff = stats.mean.lambda_diff;
fit.Mratio = stats.mean.Mratio;
fit.diff = stats.mean.diff;
fit.meta_d(:,1) = fit.Mratio(1).*fit.d1(:,1);
fit.meta_d(:,2) = fit.Mratio(2).*fit.d1(:,2);

if isrow(stats.mean.cS1)
    stats.mean.cS1 = stats.mean.cS1';
    stats.mean.cS2 = stats.mean.cS2';
end

fit.t2ca_rS1  = stats.mean.cS1;
fit.t2ca_rS2  = stats.mean.cS2;

fit.mcmc.dic = stats.dic;
fit.mcmc.Rhat = stats.Rhat;
fit.mcmc.samples = samples;
fit.mcmc.params = mcmc_params;

