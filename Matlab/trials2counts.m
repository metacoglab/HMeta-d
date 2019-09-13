function [nR_S1 nR_S2] = trials2counts(stimID,response,rating,nRatings,cellpad)

% [nR_S1 nR_S2] = trials2counts(stimID,response,rating,nRatings,cellpad)
%
% Convert trial by trial experimental information for N trials into response counts.
%
% INPUTS
% stimID:   1xN vector. stimID(i) = 0 --> stimulus on i'th trial was S1.
%                       stimID(i) = 1 --> stimulus on i'th trial was S2.
%
% response: 1xN vector. response(i) = 0 --> response on i'th trial was "S1".
%                       response(i) = 1 --> response on i'th trial was "S2".
%
% rating:   1xN vector. rating(i) = X --> rating on i'th trial was X.
%                       X must be in the range 1 <= X <= nRatings.
%
% nRatings: total # of available subjective ratings available for the
%           subject. e.g. if subject can rate confidence on a scale of 1-4,
%           then nRatings = 4
%
% cellpad: if set to 1, each response count in the output has the value of
%          1/(2*nRatings) added to it. This is desirable if trial counts of
%          0 interfere with model fitting.
%          if set to 0, trial counts are not manipulated and 0s may be
%          present.
%          default value is 0.
%
% OUTPUTS
% nR_S1, nR_S2
% these are vectors containing the total number of responses in
% each response category, conditional on presentation of S1 and S2.
%
% e.g. if nR_S1 = [100 50 20 10 5 1], then when stimulus S1 was
% presented, the subject had the following response counts:
% responded S1, rating=3 : 100 times
% responded S1, rating=2 : 50 times
% responded S1, rating=1 : 20 times
% responded S2, rating=1 : 10 times
% responded S2, rating=2 : 5 times
% responded S2, rating=3 : 1 time

nR_S1 = [];
nR_S2 = [];

% S1 responses
for r = nRatings : -1 : 1
    nR_S1(end+1) = sum(stimID==0 & response==0 & rating==r);
    nR_S2(end+1) = sum(stimID==1 & response==0 & rating==r);
end

% S2 responses
for r = 1 : nRatings
    nR_S1(end+1) = sum(stimID==0 & response==1 & rating==r);
    nR_S2(end+1) = sum(stimID==1 & response==1 & rating==r);
end

% cell pad
if ~exist('cellpad','var') || isempty(cellpad), cellpad = 0; end

if cellpad
    
    padFactor = 1/(2*nRatings);
    
    nR_S1 = nR_S1 + padFactor;
    nR_S2 = nR_S2 + padFactor;

end