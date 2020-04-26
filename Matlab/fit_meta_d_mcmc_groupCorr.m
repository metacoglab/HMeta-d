function fit = fit_meta_d_mcmc_groupCorr(nR_S1, nR_S2, mcmc_params, fncdf, fninv, name)
% fit = fit_meta_d_mcmc_groupCorr(nR_S1_1, nR_S2_1, nR_S1_2, nR_S2_2, mcmc_params, fncdf, fninv)
%
% Estimates correlation coefficient between metacognitive effiency
% estimates between two domains. See fit_meta_d_mcmc_group for full details
% of the basic model
%
% To enter data from the two tasks use the following format:
% 
% nR_S1(1).counts, nR_S2(1).counts, nR_S1(2).counts, nR_S2(2).counts
%
% 5/8/2016 Steve Fleming www.metacoglab.org
% Parts of this code are adapted from Brian Maniscalco's meta-d' toolbox
% which can be found at http://www.columbia.edu/~bsm2105/type2sdt/

fprintf('\n')
disp('----------------------------------------')
disp('Hierarchical meta-d'' model')
disp('https://github.com/smfleming/HMeta-d')
disp('----------------------------------------')
fprintf('\n')

cwd = pwd;
findpath = which('Bayes_metad_group.txt');
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

Nsubj = length(nR_S1(1).counts);
nRatings = length(nR_S1(1).counts{1})/2;

if length(nR_S1(1).counts) ~= length(nR_S1(2).counts) || length(nR_S2(1).counts) ~= length(nR_S2(2).counts)
    error('There are different numbers of subjects across the two conditions')
end

for n = 1:Nsubj
    for task = 1:2
        
        if length(nR_S1(task).counts{n}) ~= nRatings*2 || length(nR_S2(task).counts{n}) ~= nRatings*2
            error('Subjects do not have equal numbers of response categories');
        end
        % Get type 1 SDT parameter values
        counts{task}(n,:) = [nR_S1(task).counts{n} nR_S2(task).counts{n}];
        nTot(n,task) = sum(counts{task}(n,:));
        % Adjust to ensure non-zero counts for type 1 d' point estimate (not
        % necessary if estimating d' inside JAGS)
        adj_f = 1/length(nR_S1(task).counts{n});
        nR_S1_adj = nR_S1(task).counts{n} + adj_f;
        nR_S2_adj = nR_S2(task).counts{n} + adj_f;
        
        ratingHR  = [];
        ratingFAR = [];
        for c = 2:nRatings*2
            ratingHR(end+1) = sum(nR_S2_adj(c:end)) / sum(nR_S2_adj);
            ratingFAR(end+1) = sum(nR_S1_adj(c:end)) / sum(nR_S1_adj);
        end
        
        t1_index = nRatings;
        
        d1(n,task) = fninv(ratingHR(t1_index)) - fninv(ratingFAR(t1_index));
        c1(n,task) = -0.5 .* (fninv(ratingHR(t1_index)) + fninv(ratingFAR(t1_index)));
    end
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
% Assign variables to the observed nodes

datastruct = struct('d1', d1, 'c1', c1, 'nsubj',Nsubj,'counts1', counts{1}, 'counts2', counts{2}, 'nratings', nRatings, 'Tol', 1e-05);

% Select model file and parameters to monitor
model_file = 'Bayes_metad_group_corr.txt';
monitorparams = {'d1', 'c1', 'mu_logMratio', 'sigma_logMratio', 'rho', 'Mratio'};

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
        'workingdir' , 'tmpjags' );
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

fit.mu_logMratio = stats.mean.mu_logMratio;
fit.sigma_logMratio = stats.mean.sigma_logMratio;
fit.rho = stats.mean.rho;
fit.Mratio = stats.mean.Mratio;
fit.d1 = stats.mean.d1;
fit.c1 = stats.mean.c1;

fit.mcmc.dic = stats.dic;
fit.mcmc.Rhat = stats.Rhat;
fit.mcmc.samples = samples;
fit.mcmc.params = mcmc_params;

cd(cwd);