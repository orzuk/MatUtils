% Get a set of constraint on parameters for various types of diseases
% 
% Input: 
% disease_type_str - string representing type of disease model
% 
% Output: 
% intervals_struct - allowable ranges for various disease parameters 
% 
function intervals_struct = disease_type_to_constraint(disease_type_str)

switch disease_type_str % determine constraint based on disease type
    case 'rare'
        intervals_struct.ratio_interval = [0.001 0.3]; % highlight any ratios below this threshold
        intervals_struct.h_interval = [0.5 0.8]; % allowable genetic heretability content
        intervals_struct.freq_interval = [0.01 0.10]; % allowable frequency of disease in the population
        intervals_struct.penetrance_interval = [1.2 100.8]; % allowable lambda_mz (relative risk for a twin)
        intervals_struct.h_add_interval = [0.1 0.92]; % allowable fraction of variance explained by additive effects
        intervals_struct.lods_interval = [1 1.6]; % allowable lods-ratio for disease between 0 and 1 markers
    case 'common'
        intervals_struct.ratio_interval = [0.01 0.25]; % highlight any ratios below this threshold
        intervals_struct.h_interval = [0.25 0.8]; % allowable genetic heretability content
        intervals_struct.freq_interval = [0.15 0.30]; % allowable frequency of disease in the population
        intervals_struct.penetrance_interval = [1.2 8]; % allowable lambda_mz (risk for a twin)
        intervals_struct.h_add_interval = [0.00002 0.92]; % allowable fraction of variance explained by additive effects
        intervals_struct.lods_interval = [1 1.5]; % allowable lods-ratio for disease between 0 and 1 markers
    case 'high-heret'
        intervals_struct.ratio_interval = [0.01 0.25]; % highlight any ratios below this threshold
        intervals_struct.h_interval = [0.5 1]; % allowable genetic heretability content
        intervals_struct.freq_interval = [0.01 0.25]; % allowable frequency of disease in the population
        intervals_struct.penetrance_interval = [5 100]; % allowable lambda_mz (risk for a twin)
        intervals_struct.h_add_interval = [0.00002 0.92]; % allowable fraction of variance explained by additive effects
        intervals_struct.lods_interval = [1 1.5]; % allowable lods-ratio for disease between 0 and 1 markers
    case 'high-ratio'
        intervals_struct.ratio_interval = [0 0.10]; % highlight any ratios below this threshold
        intervals_struct.h_interval = [0.4 0.8]; % allowable genetic heretability content
        intervals_struct.freq_interval = [0.01 0.10]; % allowable frequency of disease in the population
        intervals_struct.penetrance_interval = [2 100]; % allowable penetrance (risk for a twin)
        intervals_struct.h_add_interval = [0.00002 0.1]; % allowable fraction of variance explained by additive effects
        intervals_struct.lods_interval = [1 1.5]; % allowable lods-ratio for disease between 0 and 1 markers
end
