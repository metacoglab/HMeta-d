function fit = fit_meta_d_mcmc_group_paired(nR_S1, nR_S2, mcmc_params, s, fncdf, fninv)
% fit = fit_meta_d_mcmc_group(nR_S1, nR_S2, mcmc_params, s, fncdf, fninv)


if ~exist('s','var') || isempty(s)
    s = 1;
end

if ~exist('fncdf','var') || isempty(fncdf)
    fncdf = @normcdf;
end

if ~exist('fninv','var') || isempty(fninv)
    fninv = @norminv;
end

Nsubj = length(nR_S1);
nRatings = length(nR_S1{1})/2;
c1_index = nRatings;
padFactor = 1/(2*nRatings);

for n = 1:Nsubj
    
    if length(nR_S1{n}) ~= nRatings*2 || length(nR_S2{n}) ~= nRatings*2
        error('Subjects do not have equal numbers of response categories');
    end
    
    for task = 1:size(nR_S1{n},1)

        % Get type 1 SDT parameter values
        counts(n,:,task) = [nR_S1{n}(task,:) nR_S2{n}(task,:)];
        nTot(n,task) = sum(counts(n,:,task));

        pad_nR_S1 = nR_S1{n}(task,:) + padFactor;
        pad_nR_S2 = nR_S2{n}(task,:) + padFactor;

        j=1;
        for c = 2:nRatings*2
            ratingHR(j) = sum(pad_nR_S2(c:(nRatings*2)))/sum(pad_nR_S2);
            ratingFAR(j) = sum(pad_nR_S1(c:(nRatings*2)))/sum(pad_nR_S1);
            j=j+1;
        end

        % Get type 1 estimate (from pair of middle HR/FAR ratings)
        d1(n,task) = norminv(ratingHR(c1_index))-norminv(ratingFAR(c1_index));
        c1(n,task) = -0.5 * (norminv(ratingHR(c1_index)) + norminv(ratingFAR(c1_index)));
    end
    
end

%% Sampling
if ~exist('mcmc_params','var') || isempty(mcmc_params)
    % MCMC Parameters
    mcmc_params.response_conditional = 0;
    mcmc_params.nchains = 3; % How Many Chains?
    mcmc_params.nburnin = 1000; % How Many Burn-in Samples?
    mcmc_params.nsamples = 10000;  %How Many Recorded Samples?
    mcmc_params.nthin = 1; % How Often is a Sample Recorded?
    mcmc_params.doparallel = 0; % Parallel Option
    mcmc_params.dic = 1;
    % Initialize Unobserved Variables
    for i=1:mcmc_params.nchains
        if mcmc_params.response_conditional
            S.mu_Mratio_rS1 = 1;
            S.lambda_Mratio_rS1 = 0.5;
            S.mu_Mratio_rS2 = 1;
            S.lambda_Mratio_rS2 = 0.5;
        else % Not response conditional
            S.mu_Mratio(1) = 1;
            S.lambda_Mratio(1) = .5;
            S.mu_Mratio(2) = 1;
            S.lambda_Mratio(2) = .5;
        end
        S.cS1_raw = linspace(-1,0.2,nRatings);
        S.cS2_raw = linspace(0.2,1,nRatings);
        mcmc_params.init0(i) = S;
    end
end
% Assign variables to the observed nodes
%                      # subjects   responses    type1dprime type1criterion     # levels      # trials total      tolerance                   
datastruct = struct('nsubj',Nsubj,'counts', counts, 'd1', d1, 'c', c1, 'nratings', nRatings, 'nTot', nTot, 'Tol', 1e-05);

% Select model file and parameters to monitor

switch mcmc_params.response_conditional
    case 0
        model_file = 'Bayes_metad_group_correl.txt';
        monitorparams = {'mu_Mratio','lambda_Mratio','Mratio','r','cS1','cS2'};
    case 1
        model_file = 'Bayes_metad_rc_group.txt';
        monitorparams = {'mu_Mratio_rS1','mu_Mratio_rS2','lambda_Mratio_rS1','lambda_Mratio_rS2','Mratio_rS1','Mratio_rS2','cS1','cS2'};
end

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

if ~mcmc_params.response_conditional
    
    fit.mu_Mratio(1) = stats.mean.mu_Mratio(1);
    fit.mu_Mratio(2) = stats.mean.mu_Mratio(2);
    fit.lambda_Mratio(1) = stats.mean.lambda_Mratio(1);
    fit.lambda_Mratio(2) = stats.mean.lambda_Mratio(2);
    fit.Mratio = stats.mean.Mratio;
    fit.r = stats.mean.r;
    fit.meta_d(:,1) = fit.Mratio(1).*d1(:,1);
    fit.meta_d(:,2) = fit.Mratio(2).*d1(:,2);
    fit.meta_da(:,1) = sqrt(2/(1+s^2)) * s * fit.meta_d(:,1);
    fit.meta_da(:,2) = sqrt(2/(1+s^2)) * s * fit.meta_d(:,2);

else
    
    fit.mu_Mratio_rS1 = stats.mean.mu_Mratio_rS1;
    fit.mu_Mratio_rS2 = stats.mean.mu_Mratio_rS2;
    fit.lambda_Mratio_rS1 = stats.mean.lambda_Mratio_rS1;
    fit.lambda_Mratio_rS2 = stats.mean.lambda_Mratio_rS2;
    fit.Mratio_rS1 = stats.mean.Mratio_rS1;
    fit.Mratio_rS2 = stats.mean.Mratio_rS2;
    fit.meta_d_rS1   = fit.Mratio_rS1.*d1;
    fit.meta_d_rS2   = fit.Mratio_rS2.*d1;
    fit.meta_da_rS1 = sqrt(2/(1+s^2)) * s * fit.meta_d_rS1;
    fit.meta_da_rS2 = sqrt(2/(1+s^2)) * s * fit.meta_d_rS2;

end
if isrow(stats.mean.cS1)
    stats.mean.cS1 = stats.mean.cS1';
    stats.mean.cS2 = stats.mean.cS2';
end

fit.da        = sqrt(2/(1+s^2)) .* s .* d1;
fit.s         = s;
fit.meta_ca   = ( sqrt(2).*s ./ sqrt(1+s.^2) ) .* c1;
fit.t2ca_rS1  = ( sqrt(2).*s ./ sqrt(1+s.^2) ) .* stats.mean.cS1;
fit.t2ca_rS2  = ( sqrt(2).*s ./ sqrt(1+s.^2) ) .* stats.mean.cS2;



fit.mcmc.dic = stats.dic;
fit.mcmc.Rhat = stats.Rhat;
fit.mcmc.samples = samples;
fit.mcmc.params = mcmc_params;

for n = 1:Nsubj
    
    
    %% Data is fit, now package output
    I_nR_rS2 = nR_S1{n}(:,nRatings+1:end);
    I_nR_rS1 = nR_S2{n}(:,nRatings:-1:1);
    
    C_nR_rS2 = nR_S2{n}(:,nRatings+1:end);
    C_nR_rS1 = nR_S1{n}(:,nRatings:-1:1);
    
    for i = 2:nRatings
        obs_FAR2_rS2(i-1,:) = sum( I_nR_rS2(:,i:end) ) ./ sum(I_nR_rS2,2)';
        obs_HR2_rS2(i-1,:)  = sum( C_nR_rS2(:,i:end) ) ./ sum(C_nR_rS2,2)';
        
        obs_FAR2_rS1(i-1,:) = sum( I_nR_rS1(:,i:end) ) ./ sum(I_nR_rS1,2)';
        obs_HR2_rS1(i-1,:)  = sum( C_nR_rS1(:,i:end) ) ./ sum(C_nR_rS1,2)';
    end
    
    for task = 1:size(nR_S1{n},1)
        % Calculate fits based on either vanilla or response-conditional model
        switch mcmc_params.response_conditional

            case 0
            

                %% find estimated t2FAR and t2HR
                meta_d = fit.meta_d(n,task);
                S1mu = -meta_d/2; S1sd = 1;
                S2mu =  meta_d/2; S2sd = S1sd/s;

                C_area_rS2 = 1-fncdf(c1(n,task),S2mu,S2sd);
                I_area_rS2 = 1-fncdf(c1(n,task),S1mu,S1sd);

                C_area_rS1 = fncdf(c1(n,task),S1mu,S1sd);
                I_area_rS1 = fncdf(c1(n,task),S2mu,S2sd);

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
            
            
            case 1

                %% find estimated t2FAR and t2HR
                S1mu_rS1 = -fit.meta_d_rS1(n)/2; S1sd = 1;
                S2mu_rS1 =  fit.meta_d_rS1(n)/2; S2sd = S1sd/s;
                S1mu_rS2 = -fit.meta_d_rS2(n)/2/2;
                S2mu_rS2 =  fit.meta_d_rS2(n)/2;

                C_area_rS2 = 1-fncdf(c1(n),S2mu_rS2,S2sd);
                I_area_rS2 = 1-fncdf(c1(n),S1mu_rS2,S1sd);

                C_area_rS1 = fncdf(c1(n),S1mu_rS1,S1sd);
                I_area_rS1 = fncdf(c1(n),S2mu_rS1,S2sd);

                t2c1 = [fit.t2ca_rS1(n,:) fit.t2ca_rS2(n,:)];

                for i=1:nRatings-1

                    t2c1_lower = t2c1(nRatings-i);
                    t2c1_upper = t2c1(nRatings-1+i);

                    I_FAR_area_rS2 = 1-fncdf(t2c1_upper,S1mu_rS2,S1sd);
                    C_HR_area_rS2  = 1-fncdf(t2c1_upper,S2mu_rS2,S2sd);

                    I_FAR_area_rS1 = fncdf(t2c1_lower,S2mu_rS1,S2sd);
                    C_HR_area_rS1  = fncdf(t2c1_lower,S1mu_rS1,S1sd);


                    est_FAR2_rS2(i) = I_FAR_area_rS2 / I_area_rS2;
                    est_HR2_rS2(i)  = C_HR_area_rS2 / C_area_rS2;

                    est_FAR2_rS1(i) = I_FAR_area_rS1 / I_area_rS1;
                    est_HR2_rS1(i)  = C_HR_area_rS1 / C_area_rS1;

                end
            
        end
    fit.est_HR2_rS1(n,:,task)  = est_HR2_rS1;
    fit.obs_HR2_rS1(n,:,task)  = obs_HR2_rS1;
    
    fit.est_FAR2_rS1(n,:,task) = est_FAR2_rS1;
    fit.obs_FAR2_rS1(n,:,task) = obs_FAR2_rS1;
    
    fit.est_HR2_rS2(n,:,task)  = est_HR2_rS2;
    fit.obs_HR2_rS2(n,:,task)  = obs_HR2_rS2;
    
    fit.est_FAR2_rS2(n,:,task) = est_FAR2_rS2;
    fit.obs_FAR2_rS2(n,:,task) = obs_FAR2_rS2;
    end
end
