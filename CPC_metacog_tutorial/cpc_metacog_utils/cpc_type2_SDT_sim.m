function sim = cpc_type2_SDT_sim(d, noise, c, Nratings, Ntrials)
% Type 2 SDT simulation with variable noise
% sim = type2_SDT_sim(d, noise, c, c1, c2, Ntrials)
%
% INPUTS
% d - type 1 dprime
% noise - standard deviation of noise to be added to type 1 internal
% response for type 2 judgment. If noise is a 1 x 2 vector then this will
% simulate response-conditional type 2 data where noise = [sigma_rS1
% sigma_rS2]
%
% c - type 1 criterion
% c1 - type 2 criteria for S1 response
% c2 - type 2 criteria for S2 response
% Ntrials - number of trials to simulate
%
% OUTPUT
%
% sim - structure containing nR_S1 and nR_S2 response counts
%
% SF 2014

% Specify the confidence criterions based on the number of ratings
c1 = c + linspace(-1.5, -0.5, (Nratings - 1));
c2 = c + linspace(0.5, 1.5, (Nratings - 1));

if length(noise) > 1
    rc = 1;
    sigma1 = noise(1);
    sigma2 = noise(2);
else
    rc = 0;
    sigma = noise;
end

S1mu = -d/2;
S2mu = d/2;

% Initialise response arrays
nC_rS1 = zeros(1, length(c1)+1);
nI_rS1 = zeros(1, length(c1)+1);
nC_rS2 = zeros(1, length(c2)+1);
nI_rS2 = zeros(1, length(c2)+1);

for t = 1:Ntrials
    s = round(rand);
    
    % Type 1 SDT model
    if s == 1
        x = normrnd(S2mu, 1);
    else
        x = normrnd(S1mu, 1);
    end
    
    % Add type 2 noise to signal
    if rc % add response-conditional noise
        if x < c
            if sigma1 > 0
                x2 = normrnd(x, sigma1);
            else
                x2 = x;
            end
        else
            if sigma2 > 0
                x2 = normrnd(x, sigma2);
            else
                x2 = x;
            end
        end
    else
        if sigma > 0
            x2 = normrnd(x,sigma);
        else
            x2 = x;
        end
    end
    
    % Generate confidence ratings
    if s == 0 && x < c      % stimulus S1 and response S1
        pos = (x2 <= [c1 c]);
        [y ind] = find(pos);
        i = min(ind);
        nC_rS1(i) = nC_rS1(i) + 1;
        
    elseif s == 0 && x >= c   % stimulus S1 and response S2
        pos = (x2 >= [c c2]);
        [y ind] = find(pos);
        i = max(ind);
        nI_rS2(i) = nI_rS2(i) + 1;
        
    elseif s == 1 && x < c  % stimulus S2 and response S1
        pos = (x2 <= [c1 c]);
        [y ind] = find(pos);
        i = min(ind);
        nI_rS1(i) = nI_rS1(i) + 1;
        
    elseif s == 1 && x >= c % stimulus S2 and response S2
        pos = (x2 >= [c c2]);
        [y ind] = find(pos);
        i = max(ind);
        nC_rS2(i) = nC_rS2(i) + 1;
    end
    
end

sim.nR_S1 = [nC_rS1 nI_rS2];
sim.nR_S2 = [nI_rS1 nC_rS2];
