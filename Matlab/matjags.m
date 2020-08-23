function [samples, stats, structArray] = matjags(dataStruct, jagsFilenm, initStructs , varargin)
% MATJAGS, a Matlab interface for JAGS
% Version 1.3.1. Tested on JAGS 3.3.0, Windows 64-bit version
%
% This code has been adapted from MATBUGS that was written by Kevin Murphy and Maryam Mahdaviani
%
% [samples, stats] = matjags(dataStruct,  bugsModelFileName,  initStructs, ...)
%
% INPUT:
% dataStruct contains values of observed variables.
% jagsFilenm is the name of the model file
% initStructs contains initial values for the latent variables (unlike
% matbugs, this is a required variable)
%
% Note: variables with names 'a.b' in the bugs model file
% should be called 'a_b' in the matlab data structure.
%
% Optional arguments passed as 'string', value pairs [default in brackets, case
% insensitive]:
% 'monitorParams' - cell array of field names (use 'a_b' instead of 'a.b')
%                   [defaults to *, which currently does nothing...]
% 'nChains'  - number of chains [3]
% 'nBurnin'  - num samples for burn-in per chain [1000]
% 'nSamples' - num samples to keep after burn-in [5000]
% 'thin'     - keep every n'th step [1]
% 'dic'      - [1] read out DIC values
% 'workingDir' - directory to store temporary data/init/coda files [pwd/tmp]. Note that total number of iterations  = nBurnin + nSamples * thin.
% 'savejagsoutput' - 0/1 = do/do not produce text files with output from JAGS
% 'verbosity'  -
%    0 = no messages during runtime;
%    1 = minimum number of messages (e.g. which chain is being executed);
%    2 = maximum number of messages
% 'cleanup' - 0/1 -- do we want to remove PREVIOUS temporary files?
% 'dotranspose' - Set to 0 (default) if you want to insure compatibility with matbugs/WinBUGS
%
% OUTPUT
% S contains the samples; each field may have a different shape:
%  S.theta(c, s)       is the value of theta in sample s, chain c
%                      (scalar variable)
%  S.theta(c, s, i)    is the value of theta(i) in sample s, chain c
%                      (vector variable)
%  S.theta(c, s, i, j) is the value of theta(i,j) in sample s, chain c
%                      (matrix variable)
%
% stats contains various statistics, currently:
%    stats.mean, stats.std and stats.Rhat, stats.DIC.
% Each field may have a different shape:
%    stats.mean.theta
%    stats.mean.theta(i)
%    stats.mean.theta(i,j)
%
% Rhat is the "estimated potential scale reduction" statistic due to
%     Gelman and Rubin.
% Rhat values less than 1.1 mean the chain has probably converged for
%     this variable.
%
% Example
%
%     [samples, stats ] = matjags( ...
%         datastruct, ...
% 		  fullfile(pwd, 'Gaussian_3.txt'), ...
% 		  init0, ...
%         'doparallel' , 0, ...
%         'nchains', nchains,...
% 		  'nburnin', nburnin,...
%         'nsamples', nsamples, ...
% 		  'thin', 1, ...
%         'dic' , 1,...
% 		  'monitorparams', {'mu','sigma'}, ...
% 		  'savejagsoutput' , 1 , ...
%         'verbosity' , 1 , ...
%         'cleanup' , 0 , ...
%         'showwarnings' , 1 , ...
%         'workingdir' , 'tmpjags' );
%
% For Windows users:
% The JAGS executable should be placed in the windows path
% In Windows 7, go to Control Panel, System and Security, System
% and click on "Advanced System Settings" followed by "Environment Variables"
% Under System variables, click on Path, and add the jags path to the string
% This could look something like "C:\Program Files\JAGS\JAGS-3.3.0\x64\bin"
%
% For MAC users:
% Mike Kalish provided some suggestions to make matjags work on a Mac. These were implemented
% in the current version but this code is not thouroughly tested yet. Please send me any changes in
% the code needed to make matjags run on a Mac/Linux/Unix system.
%
% Written by Mark Steyvers (mark.steyvers@uci.edu) based on the code
% MATBUGS that was written by Maryam Mahdaviani (maryam@cs.ubc.ca)
% and Kevin Murphy (murphyk@cs.ubc.ca)

% Changes in version 1.3:
% * Bug fix. Added lines marked by "% GP" as suggested by Ganesh Padmanabhan

% Changes in version 1.1:
% * The warnings produces by JAGS are now suppressed by default. This removes
% any message about adaptation being incomplete due to a small number of
% burnin iterations
% * Changed the default working directory for JAGS to make it platform
% independent

if ispc
    defaultworkingDir = 'c:\temp\jagstmp'; % change this to fit your particular needs
elseif ismac | isunix
    defaultworkingDir = [ '..' filesep 'jagstmp' ]; % change this to fit your particular needs
end

% Get the parameters
[ nChains, workingDir, nBurnin, nSamples, ...
    monitorParams, thin, dodic, doParallel, savejagsoutput, verbosity, cleanup, showwarnings,dotranspose,doboot ] =  ...
    process_options(...
    varargin, ...
    'nChains', 1, ...
    'workingDir', defaultworkingDir, ...
    'nBurnin', 1000, ...
    'nSamples', 5000, ...
    'monitorParams', {}, ...
    'thin', 1, ...
    'dic' , 1, ...
    'doParallel' , 0, ...
    'savejagsoutput' , 1, ...
    'verbosity' , 0,...
    'cleanup' , 0, ...
    'showwarnings' , 0 , ...
    'dotranspose' , 0 );

if length( initStructs ) ~= nChains
    error( 'Number of structures with initial values should match number of chains' );
end

[ whdir , jagsModelBase , modelextension ] = fileparts( jagsFilenm );
jagsModel = [ jagsModelBase modelextension ];

% get the current directory
curdir = pwd;
cd( whdir );
%[status,result] = dos( 'dir' )

% Does the temporary directory exist? If not, create it
if ~exist( workingDir , 'dir' )
    [SUCCESS,MESSAGE,MESSAGEID] = mkdir(workingDir);
    if SUCCESS == 0
        error( MESSAGE );
    end
end

cd( workingDir );

% Do we want to cleanup files before we start?
if cleanup==1
    delete( 'CODA*' );
    delete( 'jag*' );
end

% Create the data file
jagsData = 'jagsdata.R';
dataGenjags(dataStruct, jagsData , '', dotranspose );


nmonitor = length( monitorParams );
if nmonitor == 0
    error( 'Please specify at least one node name to monitor' );
end

% Develop a separate JAGS script for each chain
for whchain=1:nChains
    jagsScript   = sprintf( 'jagscript%d.cmd' , whchain );
    codastem     = sprintf( 'CODA%d' , whchain );
    InitData     = sprintf( 'jagsinit%d.R' , whchain );
    
    % Create the jags script for this chain
    [ fid , message ] = fopen( jagsScript , 'wt' );
    if fid == -1
        error( message );
    end
    
    %fprintf( fid , 'model in "..%s%s"\n' , filesep , jagsModel' );
    if dodic
        fprintf( fid , 'load dic\n' );
    end
        
    if ~isempty(whdir) && (strcmp(whdir(1),filesep) || (length(whdir) > 2 && whdir(2) == ':'))
        % Case when a full path string is specified for the jagsModel
        fprintf( fid , 'model in "%s"\n' , fullfile(whdir , jagsModel));
    else
        % Case when a relative path string is specified for the jagsModel
        fprintf( fid , 'model in "%s"\n' , fullfile(curdir , whdir, jagsModel));
    end
    
    fprintf( fid , 'data in %s\n' , jagsData' );
    fprintf( fid , 'compile, nchains(1)\n' );
    fprintf( fid , 'parameters in %s\n' , InitData );
    fprintf( fid , 'initialize\n' );
    fprintf( fid , 'update %d\n' , nBurnin );
    for j=1:nmonitor
        fprintf( fid , 'monitor set %s, thin(%d)\n' , monitorParams{ j } , thin );
    end
    if dodic
        %fprintf( fid , 'monitor deviance, thin(%d)\n' , thin );
        fprintf( fid , 'monitor deviance\n' );
        %fprintf( fid , 'monitor pD\n' );
    end
    fprintf( fid , 'update %d\n' , nSamples * thin );
    fprintf( fid , 'coda *, stem(''%s'')\n' , codastem );
    fclose( fid );
    
    % Create the init file
    addlines = { '".RNG.name" <- "base::Mersenne-Twister"' , ...
        sprintf( '".RNG.seed" <- %d' , whchain ) };
    dataGenjags( initStructs, InitData , addlines, dotranspose );
end

% Do we use the Matlab parallel computing toolbox?
if doParallel==1
    numworkers = parpool('size');
    if (numworkers == 0)
        error( 'Matlab pool of workers not initialized. Use command "matlabpool open 7" for example to open up a pool of 7 workers' );
    end
    
    status = cell( 1,nChains );
    result = cell( 1,nChains );
    if ispc
        parfor whchain=1:nChains
            jagsScript   = sprintf( 'jagscript%d.cmd' , whchain );
            cmd = sprintf( 'jags %s' , jagsScript );
            if verbosity > 0
                fprintf( 'Running chain %d (parallel execution)\n' , whchain  );
            end
            [status{ whchain },result{whchain}] = dos( cmd );
        end
    elseif ismac | isunix
        parfor whchain=1:nChains
            jagsScript   = sprintf( 'jagscript%d.cmd' , whchain );
            jagsPrefix = sprintf('/usr/local/bin/');
            cmd = sprintf( '%sjags %s' ,jagsPrefix, jagsScript );
            if verbosity > 0
                fprintf( 'Running chain %d (parallel execution)\n' , whchain  );
            end
            [status{ whchain },result{whchain}] = dos( cmd );
        end
    end
    
else % Run each chain serially
    status = cell( 1,nChains );
    result = cell( 1,nChains );
    for whchain=1:nChains
        jagsScript   = sprintf( 'jagscript%d.cmd' , whchain );
        if ispc
            cmd = sprintf( 'jags %s' , jagsScript );
        elseif ismac | isunix
            jagsPrefix = sprintf('/usr/local/bin/');
            cmd = sprintf( '%sjags %s' ,jagsPrefix, jagsScript );
        end
        if verbosity > 0
            fprintf( 'Running chain %d (serial execution)\n' , whchain );
        end
        [status{ whchain },result{whchain}] = dos( cmd );
    end
end

% Save the output from JAGS to a text file?
if savejagsoutput==1
    for whchain=1:nChains
        filenm = sprintf( 'jagoutput%d.txt' , whchain );
        [ fid , message ] = fopen( filenm , 'wt' );
        if fid == -1
            error( message );
        end
        resultnow = result{whchain};
        fprintf( fid , '%s' , resultnow );
        fclose( fid );
    end
end

%% Do some error checking.
% For each chain, check if the output contains some error or warning message.
for whchain=1:nChains
    resultnow = result{whchain};
    statusnow = status{ whchain };
    if status{whchain} > 0
        cd( curdir );
        error( [ 'Error from dos environment: ' resultnow ] );
    end
    
    % Do we get an error message anywhere from JAGS --> produce an error
    pattern = [ 'can''t|RUNTIME ERROR|syntax error|failed' ];
    errstr = regexpi( resultnow , pattern , 'match' );
    if ~isempty( errstr )
        cd( curdir );
        fprintf( 'Error encountered in jags (chain %d). Check output from JAGS below:\n' , whchain  );
        fprintf( 'JAGS output for chain %d\n%s\n' , whchain , resultnow );
        error( 'Stopping execution because of jags error' );
    end
    
    % Do we get a warning message anywhere from JAGS --> produce a matlab warning
    if showwarnings ~= 0
        pattern = [ 'WARNING' ];
        errstr = regexpi( resultnow , pattern , 'match' );
        if ~isempty( errstr )
            warning( 'JAGS produced a warning message. Check the output below produced by the JAGS run' );
            fprintf( 'JAGS output for chain %d\n%s\n' , whchain , resultnow );
        end
    end
    
    if verbosity == 2
        fprintf( 'JAGS output for chain %d\n%s\n' , whchain , resultnow );
    end
    
    
    % NOTE: if the error is "jags is not recognized as an internal or external
    % command, then the jags bin folder is not on the windows path"
end

%% Extract information from the output files so we can pass it back to Matlab
codaIndex = 'CODA1index.txt'; % the index files are identical across chains, just pick first one
for i=1:nChains
    codaF = [ 'CODA' , num2str(i) , 'chain1.txt' ];
    
    S = bugs2mat(codaIndex, codaF);
    structArray(i) = S;
end
samples = structsToArrays(structArray);
stats = computeStats(samples,doboot);


cd( curdir );

%% DIC calculation
if dodic
    dbar = mean( samples.deviance(:));
    dhat = min( samples.deviance(:));
    pd = dbar - dhat;
    stats.dic = pd + dbar;
end

end

%% ----- nested functions -----------------
function dataGenjags(dataStruct, fileName, addlines, dotranspose )
% This is a helper function to generate data or init files for JAGS
% Inputs:
%   fileName: name of the text file containing initial values. for each
%             chain we'll fileName_i where 'i' is the chain number,
%   dataStruct: is a Struct with name of params(consistant in the same
%               order with paramList) are fields and intial values are functions

if nargin<3
    error(['This function needs three arguments']);
end

fieldNames = fieldnames(dataStruct);
Nparam = size(fieldNames, 1);

%fileName = [fileName, '.txt'];
fid = fopen(fileName, 'w');
if fid == -1
    error(['Cannot open ', fileName ]);
end

%fprintf(fid,'list(');
for i=1:Nparam
    fn = fieldNames(i);
    fval = fn{1};
    val = getfield(dataStruct, fval);
    [sfield1, sfield2]= size(val);
    
    msfield = max(sfield1, sfield2);
    newfval = strrep(fval, '_', '.');
    newfval = [ '"' newfval '"' ];
    
    if ((sfield1 == 1) && (sfield2 == 1))  % if the field is a singleton
        fprintf(fid, '%s <-\n%G',newfval, val);
        
        %
        % One-D array:
        %   beta = c(6, 6, ...)
        %
        % 2-D or more:
        %   Y=structure(
        %     .Data = c(1, 2, ...), .Dim = c(30,5))
        %
    elseif ((length(size(val)) == 2) && ((sfield1 == 1) || (sfield2 == 1)))
        fprintf(fid, '%s <-\nc(',newfval);
        for j=1:msfield
            if (isnan(val(j)))
                fprintf(fid,'NA');
            else
                % format for winbugs
                fprintf(fid,wb_strval(val(j)));
            end
            if (j<msfield)
                fprintf(fid, ', ');
            else
                fprintf(fid, ')');
            end
        end
    else
        % non-trivial 2-D or more array
        valsize    = size(val);
        alldatalen = prod(valsize);
        %alldata = reshape(val', [1, alldatalen]);
        %alldata = alldata(:)';
        
        %Truccolo-Filho, Wilson <Wilson_Truccolo@brown.edu>
        if length(valsize)<3
            if dotranspose==0
                alldata = reshape(val, [1, alldatalen]);
            else
                alldata = reshape(val', [1, alldatalen]);
                valsize = size( val' );
            end
        elseif length(valsize)==3
            clear valTransp
            if dotranspose==1
                for j=1:valsize(3)
                    valTransp(j,:,:)=val(:,:,j)';%need a new variable, since val might be rectangular
                end
                alldata=valTransp(:)';
            else % GP 
                alldata = reshape(val, [1, alldatalen]); % GP
            end
        else
            ['Error: 4D and higher dimensional arrays not accepted']
            return
        end
        
        %fprintf(fid, '%s <-\nstructure(.Data=c(', newfval);
        fprintf(fid, '%s <-\nstructure(c(', newfval);
        for j=1:alldatalen
            if (isnan(alldata(j)))
                fprintf(fid,'NA');
            else
                % format for winbugs
                fprintf(fid,wb_strval(alldata(j)));
            end
            if (j < alldatalen)
                fprintf(fid,',');
            else
                fprintf(fid,'), .Dim=c(', alldata(j));
            end
        end
        
        for j=1:length(valsize)
            if (j < length(valsize))
                fprintf(fid, '%G,', valsize(j));
            else
                fprintf(fid, '%G))', valsize(j));
            end
        end
    end
    if (i<Nparam)
        %fprintf(fid, ', ');
        fprintf(fid, '\n');
    else
        %fprintf(fid, ')\n');
        fprintf(fid, '\n');
    end
end

if length( addlines ) > 0
    nextra = length( addlines );
    for j=1:nextra
        fprintf( fid , '%s\n' , addlines{ j } );
    end
end

fclose(fid);

%%%%%%%%
end

function s = wb_strval(v)
% Converts numeric value to a string that is acceptable by winbugs.
% This is most problematic for exponent values which must have at least 1
% decimal and two digits for the exponent. So 1.0E+01 rather than 1E+001
% Note that only Matlab on PC does 3 digits for exponent.
s = sprintf('%G', v);
if strfind(s, 'E')
    if length(strfind(s, '.')) == 0
        s = strrep(s, 'E', '.0E');
    end
    s = strrep(s, 'E+0', 'E+');
    s = strrep(s, 'E-0', 'E-');
end

%%%%%%%%
end
function f = fullfileKPM(varargin)
% fullfileKPM Concatenate strings with file separator, then convert it to a/b/c
% function f = fullfileKPM(varargin)

f = fullfile(varargin{:});
f = strrep(f, '\', '/');

%%%%%%%%
end
function A = structsToArrays(S)
% Suppose S is this struct array
%
% S(c).X1(s)
% S(c).X2(s,i)
% S(c).X3(s,i,j)
%
% where s=1:N in all cases
%
% Then we return
% A.X1(c,s)
% A.X2(c,s,i)
% A.X3(c,s,i,j)

C = length(S);
fld = fieldnames(S);
A = [];
for fi=1:length(fld)
    fname = fld{fi};
    tmp = getfield(S(1), fname);
    sz = size(tmp);
    psz = prod(sz);
    data = zeros(C, psz);
    for c=1:C
        tmp = getfield(S(c), fname);
        %data = cat(1, data, tmp);
        data(c,:) = tmp(:)';
    end
    if sz(2) > 1 % vector or matrix variable
        data = reshape(data, [C sz]);
    end
    A = setfield(A, fname, data);
end

end
%%%%%%%%%%%%

function [Rhat, m, s] = EPSR(samples)
%
% function [R, m, s] = EPSR(samples)
% "estimated potential scale reduction" statistics due to Gelman and Rubin.
% samples(i,j) for sample i, chain j
%
% R = measure of scale reduction - value below 1.1 means converged:
%                                  see Gelman p297
% m = mean(samples)
% s = std(samples)

% This is the same as the netlab function convcalc(samples')

[n m] = size(samples);
meanPerChain = mean(samples,1); % each column of samples is a chain
meanOverall = mean(meanPerChain);

% Rhat only works if more than one chain is specified.
if m > 1
    % between sequence variace
    B = (n/(m-1))*sum( (meanPerChain-meanOverall).^2);
    
    % within sequence variance
    varPerChain = var(samples);
    W = (1/m)*sum(varPerChain);
    
    vhat = ((n-1)/n)*W + (1/n)*B;
    Rhat = sqrt(vhat/(W+eps));
else
    Rhat = nan;
end

m = meanOverall;
s = std(samples(:));

%%%%%%%%%
end
function stats = computeStats(A,doboot)

fld = fieldnames(A);
N = length(fld);
stats = struct('Rhat',[], 'mean', [], 'std', [], 'ci_low' , [] , 'ci_high' , [] );
for fi=1:length(fld)
    fname = fld{fi};
    samples = getfield(A, fname);
    sz = size(samples);
    clear R m s
    % samples(c, s, i,j,k)
    Nchains = sz(1);
    Nsamples = sz(2);
    
    st_mean_per_chain = mean(samples, 2);
    st_mean_overall   = mean(st_mean_per_chain, 1);
    
    
    % "estimated potential scale reduction" statistics due to Gelman and
    % Rubin.
    if Nchains > 1
        B = (Nsamples/Nchains-1) * ...
            sum((st_mean_per_chain - repmat(st_mean_overall, [Nchains,1])).^2);
        varPerChain = var(samples, 0, 2);
        W = (1/Nchains) * sum(varPerChain);
        vhat = ((Nsamples-1)/Nsamples) * W + (1/Nsamples) * B;
        Rhat = sqrt(vhat./(W+eps));
    else
        Rhat = nan;
    end
    
    % reshape and take standard deviation over all samples, all chains
    samp_shape = size(squeeze(st_mean_overall));
    % padarray is here http://www.mathworks.com/access/helpdesk/help/toolbox/images/padarray.html
    %reshape_target = padarray(samp_shape, [0 1], Nchains * Nsamples, 'pre');
    reshape_target = [Nchains * Nsamples, samp_shape]; % fix from Andrew Jackson  a.jackson@tcd.ie
    reshaped_samples = reshape(samples, reshape_target);
    st_std_overall = std(reshaped_samples);
    
    % get the 95% interval of the samples (not the mean)
    ci_samples_overall = prctile( reshaped_samples , [ 2.5 97.5 ] , 1 );
    ci_samples_overall_low = ci_samples_overall( 1,: );
    ci_samples_overall_high = ci_samples_overall( 2,: );
    
    if ~isnan(Rhat)
        stats.Rhat = setfield(stats.Rhat, fname, squeeze(Rhat));
    end
    
    % special case - if mean is a 1-d array, make sure it's long
    squ_mean_overall = squeeze(st_mean_overall);
    st_mean_size = size(squ_mean_overall);
    if (length(st_mean_size) == 2) && (st_mean_size(2) == 1)
        stats.mean = setfield(stats.mean, fname, squ_mean_overall');
    else
        stats.mean = setfield(stats.mean, fname, squ_mean_overall);
    end
    
    stats.std = setfield(stats.std, fname, squeeze(st_std_overall));
    
    stats.ci_low = setfield(stats.ci_low, fname, squeeze(ci_samples_overall_low));
    stats.ci_high = setfield(stats.ci_high, fname, squeeze(ci_samples_overall_high));
end
end


%%%%%%%%%%%%
% Andrew Jackson - function getDICstats added a.jackson@tcd.ie
% Used to retrieve the DIC statistics from the log file
% if matbugs input DICstatus = 1
% This code is probably a bit clunky but it does the job.
% Note that the values read from the log file have 3 decimal places but
% matlab decides to give them 4 - with a zero tagged on the end. The
% precision is obviously only to 3 decimal places. sscanf.m does not seem
% to recognise the field-width and precision codes that sprintf.m does.
function DICstats = getDICstats(workingDir)
DICstats = [];
FIDlog = fopen([workingDir '\log.txt'],'r');
ct = 0;
test = 0;
endloop = 0;
while 1
    
    tline = fgets(FIDlog);
    
    if tline == -1; break; end
    if endloop; break; end
    
    if strfind(tline,'dic.set cannot be executed');
        DICstats.error = 'DIC monitor could not be set by WinBUGS';
    end
    
    if size(tline,2)>6
        % The string 'total' in the log file denotes the end of the DIC
        % stats so the loop can be ended in the next iteration.
        if strcmp(tline(1:5),'total'); endloop = 1; end;
    end
    
    if size(tline,2)>2
        % locate the DIC string identifier in the log file
        if strcmp(tline(1:3),'DIC'); test = 1; end
    end
    
    if test
        ct=ct+1;
        % DIC results are located 3 lines after the DIC string identifier
        % in the log file.
        if ct >= 4
            A = sscanf(tline,'%*s %f %f %f %f');
            S = sscanf(tline, '%s %*f %*f %*f %*f');
            % Cheng-Ta, Yang suggested this change
            %DICstats = setfield(DICstats,S,'Dbar',A(1));
            %%DICstats = setfield(DICstats,S,'Dhat',A(2));
            %DICstats = setfield(DICstats,S,'pD',A(3));
            %DICstats = setfield(DICstats,S,'DIC',A(4));
            DICstats.S.Dbar = A(1);
            DICstats.S.Dhat = A(2);
            DICstats.S.pD = A(3);
            DICstats.S.DIC = A(4);
            DICstats.S.Dhat = A(2);
        end
    end
end

fclose(FIDlog)
end

%%%%%%%%%%%%

function S=bugs2mat(file_ind,file_out,dir)
%BUGS2MAT  Read (Win)BUGS CODA output to matlab structure
%
% S=bugs2mat(file_ind,file_out,dir)
%  file_ind - index file (in ascii format)
%  file_out - output file (in ascii format)
%  dir      - directory where the files are found (optional)
%  S        - matlab structure, with CODA variables as fields
%
% The samples are stored in added 1'st dimension,
% so that 2 x 3 variable R with 1000 samples would be
% returned as S.R(1000,2,3)
%
% Note1: the data is returned in a structure that makes extraction
% of individual sample sequencies easy: the sequencies are
% directly Nx1 double vectors, as for example S.R(:,1,2).
% The computed statistics must, however, be squeezed,
% as mean(S.R,1) is a 1x2x2 matrix.
%
% Note2: in variable names "." is replaced with "_"

% To change the output structure, edit the 'eval' line in the m-file.
% For example, to return all samples as a cell, wich possibly varying
% number of samples for elements of a multidimensional variable,
% cange the 'eval' line to
%    eval(['S.' varname '={samples};']);
% Then the samples of R(2,1) would be returned as cell S.R(2,1)

% (c) Jouko.Lampinen@hut.fi, 2000
% 2003-01-14 Aki.Vehtari@hut.fi - Replace "." with "_" in variable names
% slightly modified by Maryam Mahdaviani, August 2005 (to suppress redundant output)

if nargin>2,
    file_ind=[dir '/' file_ind];
    file_out=[dir '/' file_out];
end

ind=readfile(file_ind);

data=load(file_out);

Nvars=size(ind,1);
S=[];
for k=1:Nvars
    [varname,indexstr]=strtok(ind(k,:));
    varname=strrep(varname,'.','_');
    indices=str2num(indexstr);
    if size(indices)~=[1 2]
        error(['Cannot read line: [' ind(k,:) ']']);
    end
    sdata = size(data);
    %indices
    samples=data(indices(1):indices(2),2);
    varname(varname=='[')='(';
    varname(varname==']')=')';
    leftparen=find(varname=='(');
    outstruct=varname;
    if ~isempty(leftparen)
        outstruct=sprintf('%s(:,%s',varname(1:leftparen-1),varname(leftparen+1:end));
    end
    eval(['S.' outstruct '=samples;']);
end
end


function T=readfile(filename)
f=fopen(filename,'r');
if f==-1, fclose(f); error(filename); end
i=1;
while 1
    clear line;
    line=fgetl(f);
    if ~isstr(line), break, end
    n=length(line);
    T(i,1:n)=line(1:n);
    i=i+1;
end
fclose(f);
end

% PROCESS_OPTIONS - Processes options passed to a Matlab function.
%                   This function provides a simple means of
%                   parsing attribute-value options.  Each option is
%                   named by a unique string and is given a default
%                   value.
%
% Usage:  [var1, var2, ..., varn[, unused]] = ...
%           process_optons(args, ...
%                           str1, def1, str2, def2, ..., strn, defn)
%
% Arguments:
%            args            - a cell array of input arguments, such
%                              as that provided by VARARGIN.  Its contents
%                              should alternate between strings and
%                              values.
%            str1, ..., strn - Strings that are associated with a
%                              particular variable
%            def1, ..., defn - Default values returned if no option
%                              is supplied
%
% Returns:
%            var1, ..., varn - values to be assigned to variables
%            unused          - an optional cell array of those
%                              string-value pairs that were unused;
%                              if this is not supplied, then a
%                              warning will be issued for each
%                              option in args that lacked a match.
%
% Examples:
%
% Suppose we wish to define a Matlab function 'func' that has
% required parameters x and y, and optional arguments 'u' and 'v'.
% With the definition
%
%   function y = func(x, y, varargin)
%
%     [u, v] = process_options(varargin, 'u', 0, 'v', 1);
%
% calling func(0, 1, 'v', 2) will assign 0 to x, 1 to y, 0 to u, and 2
% to v.  The parameter names are insensitive to case; calling
% func(0, 1, 'V', 2) has the same effect.  The function call
%
%   func(0, 1, 'u', 5, 'z', 2);
%
% will result in u having the value 5 and v having value 1, but
% will issue a warning that the 'z' option has not been used.  On
% the other hand, if func is defined as
%
%   function y = func(x, y, varargin)
%
%     [u, v, unused_args] = process_options(varargin, 'u', 0, 'v', 1);
%
% then the call func(0, 1, 'u', 5, 'z', 2) will yield no warning,
% and unused_args will have the value {'z', 2}.  This behaviour is
% useful for functions with options that invoke other functions
% with options; all options can be passed to the outer function and
% its unprocessed arguments can be passed to the inner function.

% Copyright (C) 2002 Mark A. Paskin
% GNU GPL
function [varargout] = process_options(args, varargin)

% Check the number of input arguments
n = length(varargin);
if (mod(n, 2))
    error('Each option must be a string/value pair.');
end

% Check the number of supplied output arguments
if (nargout < (n / 2))
    error('Insufficient number of output arguments given');
elseif (nargout == (n / 2))
    warn = 1;
    nout = n / 2;
else
    warn = 0;
    nout = n / 2 + 1;
end

% Set outputs to be defaults
varargout = cell(1, nout);
for i=2:2:n
    varargout{i/2} = varargin{i};
end

% Now process all arguments
nunused = 0;
for i=1:2:length(args)
    found = 0;
    for j=1:2:n
        if strcmpi(args{i}, varargin{j})
            varargout{(j + 1)/2} = args{i + 1};
            found = 1;
            break;
        end
    end
    if (~found)
        if (warn)
            warning(sprintf('Option ''%s'' not used.', args{i}));
            args{i}
        else
            nunused = nunused + 1;
            unused{2 * nunused - 1} = args{i};
            unused{2 * nunused} = args{i + 1};
        end
    end
end

% Assign the unused arguments
if (~warn)
    if (nunused)
        varargout{nout} = unused;
    else
        varargout{nout} = cell(0);
    end
end
end
