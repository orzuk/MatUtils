% function test_rare_variants()

gamma = 1; % assume ALL variants are functional 
f_input=0.01;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ahituv
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1=0.1,t2=0.95,n1=378,n2=379,r1=26,r2=46
beta_ahituv_from_eliana = 5.40401375 
f_ahituv_from_eliana = 0.06422145

[beta_ahituv f_ahituv V_ahituv min_gamma_ahituv] = ...
    ratio_QTL_to_var_explained(n1, n2, r1, r2, t1, t2, gamma) % Use shamil's formula
[beta_ahituv0_1, ~, V_ahituv_0_1, min_gamma_ahituv] = ...
    ratio_QTL_to_var_explained(n1, n2, r1, r2, t1, t2, gamma, f_input) % Use shamil's formula

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Romeo:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t2 = 0.75, t1 = 0.25, n2 = 878, n1 = 897, r2 = 8, r1 = 36

beta_romeo_from_eliana = -0.32 
beta_romeo_from_eliana = 0.017
0.33 % variance explained

[beta_romeo f_romeo V_romeo min_gamma_romeo] = ...
    ratio_QTL_to_var_explained(n1, n2, r1, r2, t1, t2, gamma) % Use shamil's formula
[beta_romeo0_1, ~, V_romeo, min_gamma_romeo] = ...
    ratio_QTL_to_var_explained(n1, n2, r1, r2, t1, t2, gamma, f_input) % Use shamil's formula


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Simulate data: test that reconstruction gives the same result
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
beta = 0.01; % effect size for rare allele carrier 
f = 0.1; % combined frequency of rare alleles 
V = beta^2*f*(1-f);
gamma=1; % fraction of functional alleles
t1 = 0.1; % 
t2 = 0.9; % 
n1 = 100; % individuals in left tail
n2 = 200; % individuals in right tail

[n1 r1 n2 r2] = simulate_rare_variants_tails(beta, f, gamma, t1, t2, n1, n2)
[beta_hat f_hat V_hat min_gamma_hat] = ...
    ratio_QTL_to_var_explained(n1, n2, r1, r2, t1, t2, gamma) % Use shamil's formula
beta_error = beta - beta_hat
f_error = f-f_hat
V_error = V-V_hat




