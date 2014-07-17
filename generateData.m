% Generate model recovery data for d' model without noise
%
% SF 2014

% Parameters
Ntrials = 1000;
s = 1;
meta_d = 2;
c1 = 0;
t2ca_rS1 = [-2 -1 -0.5];
t2ca_rS2 = [0.5 1 2];
nRatings = length(t2ca_rS1)+1;

S1mu = -meta_d/2; S1sd = 1;
S2mu =  meta_d/2; S2sd = S1sd/s;

C_area_rS2 = 1-normcdf(c1,S2mu,S2sd);
I_area_rS2 = 1-normcdf(c1,S1mu,S1sd);

C_area_rS1 = normcdf(c1,S1mu,S1sd);
I_area_rS1 = normcdf(c1,S2mu,S2sd);

t2c1 = [t2ca_rS1 t2ca_rS2];

t2c1x = [-Inf t2c1(1:nRatings-1) c1 t2c1(nRatings:end) Inf];

% No division by type 1 response as we want actual probabilities of each confidence bin for this dprime not
% conditional on previous type 1 responses
for i = 1:nRatings
    prC_rS1(i) = ( normcdf(t2c1x(i+1),S1mu,S1sd) - normcdf(t2c1x(i),S1mu,S1sd) );
    prI_rS1(i) = ( normcdf(t2c1x(i+1),S2mu,S2sd) - normcdf(t2c1x(i),S2mu,S2sd) );
    
    prC_rS2(i) = ( (1-normcdf(t2c1x(nRatings+i),S2mu,S2sd)) - (1-normcdf(t2c1x(nRatings+i+1),S2mu,S2sd)) );
    prI_rS2(i) = ( (1-normcdf(t2c1x(nRatings+i),S1mu,S1sd)) - (1-normcdf(t2c1x(nRatings+i+1),S1mu,S1sd)) );
end

nC_rS1 = round(prC_rS1.*Ntrials);
nI_rS1 = round(prI_rS1.*Ntrials);
nC_rS2 = round(prC_rS2.*Ntrials);
nI_rS2 = round(prI_rS2.*Ntrials);
sim.nR_S1 = [nC_rS1 nI_rS2];
sim.nR_S2 = [nI_rS1 nC_rS2];