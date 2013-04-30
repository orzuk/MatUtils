% A little example showing that the FDR is not monotonic 
eps = 0.05; %  a small constant 
n = 10000; % number of simulations

P = rand(3,n); P(2:3,:) = min(P(2:3,:), eps); 

% q_BH = FDR_R_to_q(P, 1, 'bh')  % Compute the FDR of the lowest 1 p-values using standard BH procedure

[R1 V1] = FDR_q_to_R(P', 0.34, 'min_k', [], 1); % reject the minimal p-value ( a more conservative procedure) 
FDR1 = mean(V1 ./ max(R1,1))

[R2 V2] = FDR_q_to_R(P', 0.67, 'min_k', [], 1); % reject the two minimal p-values (a more liberal procedure)
FDR2 = mean(V2 ./ max(R2,1))

FDR_diff = FDR2 - FDR1 %  liberal  - conservative. The 'surprise' here is that this value is negative 

% Now check if step-up vs. step-down preserves monotinicity
q = eps/2; % P = rand(3,n);
[R_down V_down] = FDR_q_to_R(P(1:3,:)', q, 'bh_down', [], 1); % standard BH step-down ( a more conservative procedure) 
FDR_down = mean(V_down ./ max(R_down,1))

[R_up V_up] = FDR_q_to_R(P(1:3,:)', q, 'bh95', [], 1); % standard BH step-down ( a more conservative procedure) 
FDR_up = mean(V_up ./ max(R_up,1))

FDR_diff = FDR_up - FDR_down %  liberal  - conservative. The 'surprise' here is that this value is negative 

