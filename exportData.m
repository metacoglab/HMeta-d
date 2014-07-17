% Loops over analysis.m for different subjects
%
% SF 2012

clear all
close all

cwd = pwd;

root_direc = '~/Dropbox/Published/FlemingScience2010/Code & analysis/Behaviour/Data/';

subjects = { ...
    'hb', 1,1,6;
    'an',2,1,6;
    'am',3,1,6;
    'pw',4,1,6;
    'cm',5,1,6;
    'wa',8,1,6;
    'sr',9,1,6;
    'js',10,1,6;
    'rm',11,1,6;
    'ec',12,1,6;
    'ob',14,1,6;
    'tg',15,1,6;
    'so',16,1,6;
    'rg',17,1,6;
    'cb',18,1,6;
    'ak',19,1,6;
    'md',20,1,6;
    'kc',20,1,6;
    % 'cb',22,1,6;
    'hb',21,1,6;
    'kj', 24,1,6;
    'ms',23,1,6;
    'gp',25,1,6;
    'sb',26,1,6;
    'im',27,1,6;
    % 'km',28,1,6;
    'sm',29,1,6;
    'rr',30,1,6;
    'js',31,1,6;
    'cy',32,1,6;
    'js',33,1,6;
    'mk',34,1,6;
    'em',35,1,6;
    'lj',36,1,6;
    };

all_dataframe = [];

for i = 1:length(subjects)
    
    subj_direc=[root_direc subjects{i,1} '_' num2str(subjects{i,2})  '_' num2str(subjects{i,3})];
    cd (subj_direc)
    
    datafile = ['metacog_sub_' subjects{i,1} '_' num2str(subjects{i,2}) '_' num2str(subjects{i,3}) '_' num2str(subjects{i,4})];
    load(datafile)
    
    ntrials = DATA.stim.ntrials;
    nblock = 5;
    start = 2;
    finish = 6;
    
    typeIcode = reshape(DATA.responses.typeI.code(start:finish,:)',nblock*ntrials,1);
    typeIIcon = reshape(DATA.responses.typeII.con(start:finish,:)',nblock*ntrials,1);
    typeIRT = reshape(DATA.responses.typeI.RT(start:finish,:)',nblock*ntrials,1);
    typeIIRT = reshape(DATA.responses.typeII.RT(start:finish,:)',nblock*ntrials,1);
    
    err = typeIIcon > 6;
    
    stimID = (typeIcode(~err) == 1 | typeIcode(~err) == 3); % S1 presented
    response = typeIcode(~err) < 3;   % S1 response
    rating = typeIIcon(~err);
    nRatings = 6;

    % get tallies of "S1" rating responses for S1 and S2 stim
    for r = 1:nRatings
        nR_S1(r) = sum(stimID==0 & response==0 & rating==nRatings+1-r);
        nR_S2(r) = sum(stimID==1 & response==0 & rating==nRatings+1-r);
    end
    
    % get tallies of "S2" rating responses for S1 and S2 stim
    for r = 1:nRatings
        nR_S1(r+nRatings) = sum(stimID==0 & response==1 & rating==r);
        nR_S2(r+nRatings) = sum(stimID==1 & response==1 & rating==r);
    end
    
    nR_S1_all{i} = nR_S1;
    nR_S2_all{i} = nR_S2;
end

cd(cwd);

% Fit data using hierarchical model
fit = fit_meta_d_mcmc_group(nR_S1_all, nR_S2_all);


