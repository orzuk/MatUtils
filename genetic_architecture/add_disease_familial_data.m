% Add information from familial studies to disease data struct
function data = add_disease_familial_data(data, disease_data_from_literature_file)

prevalence_base = 0.01; % we give this to diseases where we don't know the prevalence
lambda_s_familial_base = -1; % this is when we don't know lambda_s
num_filtered =length(data.Trait);
disease_literature_from_file = 1;
if(disease_literature_from_file)
    disease_data_from_literature_vec = ...
        ReadDataFile(disease_data_from_literature_file, [], -1); % read data
    disease_data_from_literature_vec.Lifetime_prevalence = ...
        empty_cell_to_numeric_val(disease_data_from_literature_vec.Lifetime_prevalence, '-1');
    disease_data_from_literature_vec.Point_prevalence = ...
        empty_cell_to_numeric_val(disease_data_from_literature_vec.Point_prevalence,  num2str(prevalence_base*100));
    disease_data_from_literature_vec.Genetic_heritability__twin_study_  = ...
        empty_cell_to_numeric_val(disease_data_from_literature_vec.Genetic_heritability__twin_study_,  num2str(-1*100));
    disease_data_from_literature_vec.Sibling_relative_risk___  = ...
        empty_cell_to_numeric_val(disease_data_from_literature_vec.Sibling_relative_risk___,  num2str(-1));
else % just take what's here in the .mat source
    disease_data_from_literature_vec = ... % Data doesn't contain prevalence - we have to collect it ourselves from different sources
        {'Type 1 diabetes', 0.004, 15;
        'Type 2 diabetes', 0.04, [1.8 3];
        'Alzheimer''s disease', 0.01, 1.31; %  [0.58 0.79];
        'Parkinson''s disease', 0.01, 2.3;
        'Crohn''s disease', 0.001, [17 35];
        'Breast cancer', 0.04, 2.3;
        'Colorectal cancer', 0.02, [1.72 1.78];
        'Lung cancer', 0.02, 1.85;
        'Prostate cancer', 0.06, 2.87;
        'Schizophrenia', 0.01, 6.99;
        'Asthma', 0.082, 2.6;
        'Psoriasis', 0.02, 8.8;
        'Autism', prevalence_base, 22;
        };
end

[intersect_literature I J] = intersect_all(lower(data.Trait), ...
    lower(disease_data_from_literature_vec.Disease));
data.Prevalence = zeros(num_filtered,1) + prevalence_base; % Here all are pre-set to be of prevalence base
%data.Prevalence(I) = max( ... % input is in percents
tmp_prevalence = ... %%%
    cell2mat(str2nums_cell(disease_data_from_literature_vec.Point_prevalence(J),1)) ./ 100;
good_prevalence_inds = find(tmp_prevalence ~= prevalence_base); % we don't have enough lifetime prevalence ...
tmp_prevalence2 = ...
    cell2mat(str2nums_cell(disease_data_from_literature_vec.Lifetime_prevalence(J),1)) ./ 100;  % we don't have enough lifetime prevalence ..
good_prevalence_inds2 = find(tmp_prevalence2 ~= prevalence_base);
intersect_good_inds = intersect(good_prevalence_inds, good_prevalence_inds2);
diff_good_inds = setdiff(good_prevalence_inds2, good_prevalence_inds);
tmp_prevalence(intersect_good_inds) = max(tmp_prevalence(intersect_good_inds), ...
    tmp_prevalence2(intersect_good_inds));
tmp_prevalence(diff_good_inds) = tmp_prevalence2(diff_good_inds);
tmp_prevalence(tmp_prevalence < 0) = prevalence_base;
data.Prevalence(I) = tmp_prevalence;
data.base_prevalence_inds = ones(length(data.Prevalence),1);
data.base_prevalence_inds(union(intersect_good_inds, diff_good_inds)) = 0; 

data.h_familial = zeros(num_filtered,1) - 1;
data.h_familial(I) = ...
    cell2mat(str2nums_cell(disease_data_from_literature_vec.Genetic_heritability__twin_study_(J),1)) ./ 100;  % we don't have enough lifetime prevalence ...
data.lambda_s_familial = zeros(num_filtered,1) - 1;
data.lambda_s_familial(I) = ...
    cell2mat(str2nums_cell(disease_data_from_literature_vec.Sibling_relative_risk___(J),1));  % we don't have enough lifetime prevalence ...
