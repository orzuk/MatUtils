%function run_load_and_norm_data()
function run_load_and_norm_data()

sample_names = {'HD78_9_d'; 'HD78_9_n'; 'TEL74_5_d'; 'TEL74_5_n'};
array_names = {'D:\Program Files\Affymetrix\GeneChip\Affy_Data\Data\78B_Xba.CEL'; ...
    'D:\Program Files\Affymetrix\GeneChip\Affy_Data\Data\79B_Xba.CEL'; ...
    'D:\Program Files\Affymetrix\GeneChip\Affy_Data\Data\74B_Xba.CEL'; ...
    'D:\Program Files\Affymetrix\GeneChip\Affy_Data\Data\75B_Xba.CEL'};

working_path = 'E:\Libi\tools\SNP_tool\data\Leukemia_our_norm\';
chip = 'xba';
diploid_ind = get_diploid_ind(sample_names);
load_and_norm_data(sample_names, array_names, working_path, diploid_ind, chip);