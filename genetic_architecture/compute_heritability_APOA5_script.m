% Compute heritability for APOA5 Data for early Myocardial Infraction from Ron Do
% function compute_heritability_APOA5_script()

APOA5_variants_data_file = '../../common_disease_model/data/RonDo_APOA5/APOA5_summary_Tables_forOR.txt'; 
APOA5_samples_data_file = '../../common_disease_model/data/RonDo_APOA5/APOA5_Samples_summary_Table_forOR.txt';

R = ReadDataFile(APOA5_variants_data_file, [], [], 2, sprintf('\t')); % First load data. Use tab-delimiter 
R_samples = ReadDataFile(APOA5_samples_data_file, [], [], [], sprintf('\t')); % Load samples data. Use tab-delimiter 


num_cases = R_samples.N_cases(end-1)
num_controls = R_samples.N_controls(end-1)
num_carriers_cases = R.casetotal(end)
num_carriers_controls = R.controltotal(end)
alpha_vec=[0.1 0.3 0.5 0.7 1]; % assume ~50% of alleles are null
prevalence_vec=[0.001 0.01 0.03 0.05 0.1 0.2 0.5]; % assume 10% prevalence

f_cases = num_carriers_cases / (2*num_cases)
f_controls = num_carriers_controls / (2*num_controls)

for i=1:length(alpha_vec)
    for j=1:length(prevalence_vec)
        [beta(i,j) f(i,j) V(i,j) min_alpha] = ratio_QTL_to_var_explained(num_cases, num_controls, num_carriers_cases, num_carriers_controls, ...
            prevalence_vec(j), prevalence_vec(j), alpha_vec(i));
    end
end



% New: add PFKM single variant
PFKM_OR = 0.61; PFKM_Control_AF = 0.0192; PFKM_Case_AF = 0.0120; MI_prevalence = 0.03; 
PFKM_MAF = PFKM_Control_AF * (1-MI_prevalence) + PFKM_Case_AF * MI_prevalence; 
PFKM_GRR = odds_ratio_to_genetic_relative_risk(PFKM_OR, PFKM_Control_AF, MI_prevalence)
PFKM_beta = genetic_relative_risk_to_beta(PFKM_MAF, PFKM_GRR, MI_prevalence)
PFKM_Var_Explained = genetic_relative_risk_to_variance_explained(PFKM_MAF, PFKM_GRR, MI_prevalence, 'diploid', 'liability')


% New data for LDLR gene:
APOA5_mat = [93 42 6721 6711]; % [carriers-cases carriers-controls cases controls
LDLR_mat = [88 27 2743 2465];  % [carriers-cases carriers-controls cases controls  enriched-alleles
LDLR_all_mat = [192 103 2743 2465]; %  [carriers-cases carriers-controls cases controls all-nonsynonymous alleles

LDLR_new_mat = [285 208 4703 5090];  % NEW DATA FOR LDLR: 11/2013

for i=1:length(alpha_vec)
    for j=1:length(prevalence_vec)
        [APOA5_beta(i,j) APOA5_f(i,j) APOA5_V(i,j) APOA5_min_alpha] = ...
            ratio_QTL_to_var_explained(APOA5_mat(3), APOA5_mat(4), APOA5_mat(1), APOA5_mat(2), ...
            prevalence_vec(j), prevalence_vec(j), alpha_vec(i));
        [LDLR_beta(i,j) LDLR_f(i,j) LDLR_V(i,j) LDLR_min_alpha] = ...
            ratio_QTL_to_var_explained(LDLR_mat(3), LDLR_mat(4), LDLR_mat(1), LDLR_mat(2), ...
            prevalence_vec(j), prevalence_vec(j), alpha_vec(i));
        [LDLR_all_beta(i,j) LDLR_all_f(i,j) LDLR_all_V(i,j) LDLR_all_min_alpha] = ...
            ratio_QTL_to_var_explained(LDLR_all_mat(3), LDLR_all_mat(4), LDLR_all_mat(1), LDLR_all_mat(2), ...
            prevalence_vec(j), prevalence_vec(j), alpha_vec(i));
        [LDLR_new_beta(i,j) LDLR_new_f(i,j) LDLR_new_V(i,j) LDLR_new_min_alpha] = ...
            ratio_QTL_to_var_explained(LDLR_new_mat(3), LDLR_new_mat(4), LDLR_new_mat(1), LDLR_new_mat(2), ...
            prevalence_vec(j), prevalence_vec(j), alpha_vec(i));
    end
    
end





