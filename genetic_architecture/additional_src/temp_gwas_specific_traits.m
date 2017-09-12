% Temp script for analyzing height and T2D
% function temp_gwas_specific_traits()

% Deal with height: input effect size is in beta
H = [0.28	0.019
0.48	-0.049
0.59	0.031
0.26	-0.036
0.67	0.02
0.75	-0.025
0.12	0.045
0.61	0.02
0.37	-0.027
0.24	-0.048
0.58	-0.048
0.73	-0.045
0.43	-0.016
0.64	-0.042
0.53	-0.02
0.71	-0.019
0.47	-0.023
0.77	0.033
0.54	0.021
0.28	-0.026
0.91	-0.061
0.27	0.02
0.35	0.036
0.23	-0.059
0.67	-0.029
0.8	-0.033
0.08	0.028
0.79	0.044
0.9	0.028
0.55	0.023
0.19	-0.051
0.24	-0.018
0.9	-0.031
0.94	-0.048
0.55	-0.028
0.46	-0.032
0.88	0.028
0.22	-0.033
0.21	-0.035
0.25	-0.017
0.56	-0.075
0.31	0.03
0.39	0.026
0.36	0.026
0.85	0.081
0.2	0.028
0.47	0.038
0.68	-0.05
0.49	0.028
0.16	-0.08
0.24	-0.019
0.4	0.032
0.7	-0.026
0.56	0.04
0.07	-0.037
0.47	-0.019
0.4	-0.032
0.73	0.029
0.2	-0.038
0.25	-0.024
0.36	-0.032
0.22	-0.028
0.39	-0.01
0.54	-0.035
0.75	-0.037
0.7	0.051
0.39	-0.041
0.45	0.051
0.51	-0.037
0.92	-0.072
0.02	-0.068
0.22	-0.019
0.89	-0.033
0.52	0.02
0.68	-0.051
0.58	0.016
0.42	-0.026
0.5	0.028
0.76	-0.047
0.29	-0.055
0.4	-0.037
0.06	-0.045
0.3	-0.042
0.18	0.023
0.22	-0.04
0.3	0.038
0.32	-0.025
0.31	0.062
0.74	-0.017
0.69	-0.028
0.25	-0.019
0.87	0.064
0.72	-0.023
0.2	-0.056
0.6	0.024
0.32	0.017
0.11	-0.04
0.77	-0.021
0.53	0.028
0.24	0.037
0.44	0.033
0.04	0.07
0.92	0.05
0.23	0.025
0.25	-0.026
0.72	0.021
0.64	0.024
0.33	0.024
0.44	-0.02
0.49	-0.027
0.38	-0.036
0.11	0.05
0.55	-0.019
0.35	0.019
0.69	-0.02
0.34	0.02
0.18	0.026
0.06	-0.057
0.35	-0.023
0.14	0.035
0.41	0.02
0.62	0.028
0.64	-0.029
0.33	0.036
0.68	0.026
0.93	-0.058
0.49	0.073
0.35	0.042
0.22	0.05
0.46	-0.034
0.78	-0.035
0.62	-0.02
0.29	-0.063
0.4	0.019
0.29	0.038
0.58	-0.029
0.2	-0.026
0.36	-0.037
0.36	-0.016
0.05	-0.049
0.47	-0.017
0.91	0.034
0.97	-0.051
0.54	-0.031
0.48	-0.047
0.03	-0.124
0.88	0.062
0.68	-0.015
0.74	-0.039
0.46	0.04
0.34	0.033
0.79	0.015
0.61	0.017
0.33	-0.021
0.39	-0.032
0.15	0.017
0.45	-0.024
0.35	-0.037
0.3	0.013
0.34	0.018
0.65	-0.034
0.33	0.04
0.73	-0.052
0.34	0.026
0.79	0.056
0.58	-0.039
0.76	-0.035
0.6	-0.027
0.24	0.049
0.74	-0.035
0.46	-0.034
0.74	0.018
0.36	0.037
0.63	-0.016
0.65	-0.04
0.23	-0.042
0.58	-0.061
0.21	0.053
0.65	0.024
0.84	0.027]; 

beta = H(:,2); p = H(:,1);
sign_vec = (sign(beta) + 1) / 2;
p = p .* sign_vec + (1-p) .* (1-sign_vec);
h_liab = 2 .* beta.^2 .* p .* (1-p);
total_h_liab_height = sum(h_liab) % shoulf be ~8.8%
figure; plot(p, abs(beta), '.'); title('RAF vs. \beta for Height');
xlabel('Risk-Allele-Frequency'); ylabel('\beta');


GRR_vec_height = heritability_to_genetic_relative_risk(h_liab, 'liability', p, 0.1);
[lambda_s_vec lambda_s_add lambda_mz_add max_h_add max_V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat ] = ...
    genetic_relative_risk_to_heritability( ...
    mat_into_vec(repmat(p, 1, 2)'), mat_into_vec(repmat(GRR_vec_height, 1, 2)'), 0.1);


figure; plot(p, GRR_vec_height, '.'); title('RAF vs. GRR for Height (\mu = 0.1)');
xlabel('Risk-Allele-Frequency'); ylabel('GRR');

lambda_s_vec2 = lambda_s_vec(1:2:end).^2;
figure; plot(cumsum(h_liab), cumprod(lambda_s_vec2), '.')
title('Variance explained for Height'); 
xlabel('Variance explained '); ylabel('Sibling Relative Risk Explained (Tall)');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Deal with T2D: input effect size is GRR 
T = [0.46 	1.08 
0.26	1.07 
0.55	1.06 
0.48	1.05 
0.93	1.08 
0.52	1.07 
0.88	1.14 
0.1	1.08 
0.85	1.05 
0.6	1.05 
0.22	1.06 
0.21	1.32 
0.64	1.10 
0.28	1.08 ];



f_vec = mat_into_vec(repmat(T(:,1), 1, 2)');
GRR_vec = mat_into_vec(repmat(T(:,2), 1, 2)');
[lambda_s_vec lambda_s_add lambda_mz_add h_add V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat h_liab] = ...
    genetic_relative_risk_to_heritability(f_vec, GRR_vec, 0.1);
total_lambda_s = prod(lambda_s_vec)

% Try 32 loci: 
read_voight_diabetes_loci;
voights_diabetes_loci_outfile = '../../common_disease_model/data/gwas_database_catalog/Voight_T2D_32_loci.txt';

% X = rmfield(X, 'pos_vec'); % why??
diabetes_prevalence = 0.08;
WriteDataFile(X, voights_diabetes_loci_outfile); 

MODY_freq = 0.008; MODY_GRR = 6;
X.GRR_vec = [X.GRR_vec' MODY_GRR]'; 
X.RAF_min = [X.RAF_min' MODY_freq]';
X.RAF_max = [X.RAF_max' MODY_freq]';

[min_lambda_s_vec lambda_s_add lambda_mz_add min_h_add min_V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat min_h_liab] = ...
    genetic_relative_risk_to_heritability( ...
    mat_into_vec(repmat(X.RAF_min, 1, 2)'), mat_into_vec(repmat(X.GRR_vec, 1, 2)'), diabetes_prevalence);
total_min_lambda_s = prod(min_lambda_s_vec)

[max_lambda_s_vec lambda_s_add lambda_mz_add max_h_add max_V_add ...
    lambda_s_mult lambda_mz_mult h_mult V_mult sib_freq_mat max_h_liab] = ...
    genetic_relative_risk_to_heritability( ...
    mat_into_vec(repmat(X.RAF_max, 1, 2)'), mat_into_vec(repmat(X.GRR_vec, 1, 2)'), diabetes_prevalence);
total_max_lambda_s = prod(max_lambda_s_vec)

figure; plot(max_V_add ./ V_mult, max_lambda_s_vec, '.')
figure; plot(cumsum(max_V_add .* max_h_liab ./ max_h_add) ./ V_mult, cumprod(max_lambda_s_vec), '.')
title('Variance explained for Type 2 Diabetes'); 
xlabel('Variance explained (liability)'); ylabel('Sibling Relative Risk Explained');

prevalence = 0.08;
sqrt(heritability_to_mz_twin_risk( heritability_scale_change(0.26, 'binary', prevalence), prevalence))
sqrt(heritability_to_mz_twin_risk( heritability_scale_change(0.997, 'binary', prevalence), prevalence))

%% Try to replicate heritability estimates from twin pairs. Should get 26% !!! 
MZ_DZ_table = [XX YY; ZZ TT]; % counts of number of affected twins 
% First row: MZ, second row: DZ. First column: healthy, second column: sick
% Total of column: 5.9% * 303 = 18
% Total of table: 303 

[h h_add thresh mean_liability_given_affected] = ...
    family_segregation_to_heritability(family_tree, family_trait_values, scale, varargin)

h = mz_twin_risk_to_heritability(); % compute heritabiity using ...
h_x = heritability_scale_change(h, 'liability', prevalence) % get heritability on liability scale 




T2D_NHGRI = {'rs4430796','rs1470579','rs13081389','rs5945326','rs8042680','rs11634397','rs7957197','rs1552224', ...
    'rs231362','rs896854','rs972283','rs11642841','rs4760790','rs7903146','rs5015480','rs10965250','rs1387153', ...
    'rs7578326','rs1531343','rs13292136','rs243021','rs1801214','rs3802177','rs4457053','rs10440833','rs849134'};


T2D_V32 = {'rs10923931','rs11899863','rs13081389','rs6795735','rs1470579','rs1801214', ...
'rs10440833','rs849134','rs3802177','rs10965250','rs12779790','rs5015480','rs7903146','rs163184', ... 
'rs5215','rs4760790','rs11642841','rs4430796','rs1387153','rs7578326','rs243021','rs4457053', ...
'rs972283','rs896854','rs13292136','rs231362','rs1552224','rs1531343','rs7957197','rs11634397','rs8042680','rs5945326'};


