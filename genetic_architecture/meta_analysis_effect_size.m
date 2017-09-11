% Determine combine effect size from effect sizes of each seperate study
% 
% Input: 
% num_cases_vec 
% num_controls_vec
% f_vec 
% grr_vec 
% mu 
% trait_type
%
% Output: 
% combined_f_vec - allele frequency in population
% combined_grr_vec - effect size in population 
%
function [combined_f_vec combined_grr_vec] = ...
    meta_analysis_effect_size(num_cases_vec, num_controls_vec, f_vec, grr_vec, mu, trait_type)


num_studies = length(num_cases_vec);
num_snps = size(f_vec, 1)

combined_f_vec = sum(f_vec .* repmat(num_cases_vec, num_snps, 1),2) ./ sum(num_cases_vec); % take average for QTL 

switch trait_type    
    case 'QTL' % for QTL just compute a weighted average of effect sizes 
        combined_grr_vec = sum(grr_vec .* num_cases_vec) ./ sum(num_cases_vec);
        
    case 'Binary' % More complicated: compute for each study the 2x2 table 
        p_mat = zeros(num_snps, 4); 
        
        for i=1:num_studies
            study_p_mat = genetic_relative_risk_to_p_z_x_marginal(f_vec(:,i), grr_vec(:,i), mu);
            study_p_mat = pop_prob_to_case_control_prob(study_p_mat, num_cases_vec(i), num_controls_vec(i)); 
            p_mat = p_mat + study_p_mat .* (num_cases_vec(i) + num_controls_vec(i));             
        end
        p_mat = p_mat ./ (sum(num_cases_vec) + sum(num_controls_vec));
        p_mat = case_control_prob_to_pop_prob(p_mat, mu); % transfer back to population prob.        
        [~, combined_grr_vec] = p_z_x_marginal_to_genetic_relative_risk(p_mat); % convert back to grr
end        
   

        
        