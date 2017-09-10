% A short test for FDR procedures 
P = rand(1,1000); % generate a vector of p-values, of which a 100 are taken from the alternative hypothesis
P(1:500) = P(1:500) ./ 33; % some p-values are smaller
R = 150; % number of desired rejected hypothesis
q_BH = FDR_R_to_q(P, R, 'bh')  % Compute the FDR of the lowest 150 p-values using standard BH procedure
q_IBH = FDR_R_to_q(P, R, 'ibh')  % Compute the FDR of the lowest 150 p-values using our imptoved BH procedure

q = 0.05;
R_BH = FDR_q_to_R(P, q, 'bh')  % Compute number of rejected hypothesis at q=0.05 using standard BH procedure
R_IBH = FDR_q_to_R(P, q, 'ibh_up')  % Compute number of rejected hypothesis at q=0.05 using our imptoved BH procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


