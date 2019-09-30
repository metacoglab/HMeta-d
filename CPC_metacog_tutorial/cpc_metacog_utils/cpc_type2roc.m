%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CPC METACOGNITION TUTORIAL 2019: CALCULATE TYPE2ROC %%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Function to calculate the area under the type2 ROC curve

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function auroc2 = cpc_type2roc(nR_S1, nR_S2, Nratings)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
flipped_nR_S1 = fliplr(nR_S1);
flipped_nR_S2 = fliplr(nR_S2);

for c = 1:Nratings
    S1_H2(c) = nR_S1(c) + 0.5;
    S2_H2(c) = flipped_nR_S2(c) + 0.5;
    S1_FA2(c) = flipped_nR_S1(c) + 0.5;
    S2_FA2(c) = nR_S2(c) + 0.5;
end

H2 = S1_H2 + S2_H2;
FA2 = S1_FA2 + S2_FA2;
 
H2 = H2./sum(H2);
FA2 = FA2./sum(FA2);
cum_H2 = [0 cumsum(H2)];
cum_FA2 = [0 cumsum(FA2)];
 
i=1;
for c = 1:Nratings
        k(i) = (cum_H2(c+1) - cum_FA2(c))^2 - (cum_H2(c) - cum_FA2(c+1))^2;
        i = i+1;
end
auroc2 = 0.5 + 0.25*sum(k);


end