% Parse data from Eliana for different effect size under different
% envirounments. Compute the underlying heritabilities and their gap 
%function parse_GxE_data()

GxE_file_name = '../../common_disease_model/data/GxE/glucosevraffinose.txt';
R = loadcellfile(GxE_file_name);

title_inds_vec = [7 24]; num_vals = 9; % Find lambda_s and prevalence

for j=1:length(title_inds_vec)
    title_ind = title_inds_vec(j);
    mu = cell2vec(R(title_ind+1:title_ind+num_vals,3));
    lambda_s_top = cell2vec(R(title_ind+1:title_ind+num_vals,4));
    lambda_s_bottom = cell2vec(R(title_ind+1:title_ind+num_vals,5));
    if(j==1)
        R{title_ind,6} = [R{title_ind,6} ' (%)'];
        R{title_ind,7} = [R{title_ind,7} ' (%)'];
        R{title_ind,8} = 'h^2 phantom (%)';
    end
    for i=1:num_vals
        run_h_i = i
        h_top(i) = familial_risk_to_heritability(lambda_s_top(i),  'liability', mu(i), 0.5);
        h_bottom(i) = familial_risk_to_heritability(lambda_s_bottom(i),  'liability', mu(i), 0.5);
        R{title_ind+i,6} = 100*h_top(i);
        R{title_ind+i,7} = 100*h_bottom(i);
        R{title_ind+i,8} = 100*(h_top(i)-h_bottom(i))/h_top(i);
    end
end

R = num2str_cell(R, 3); 
savecellfile(R, [remove_suffix_from_file_name(GxE_file_name) '_with_h.txt'], [], 1);
