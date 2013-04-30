%function create_and_save_norm_summary_file(sample_names, median_intensity_mat, user_dir, chip, array_outlier_vec)
function create_and_save_norm_summary_file(sample_names, median_intensity_mat, user_dir, chip, array_outlier_vec)

num_samples = length(sample_names);

if(size(array_outlier_vec,1) ==1) array_outlier_vec = array_outlier_vec'; end
norm_summary_cell = cell(num_samples+1,3); 
norm_summary_cell{1,1} = 'Sample'; norm_summary_cell{1,2} = 'Median Intensity (unnormalized)';
norm_summary_cell{1,3} = 'Normalized Median Intensity';
norm_summary_cell{1,4} = '% Array outlier';
% fill the cell
if(size(sample_names,1)==1) sample_names = sample_names'; end

norm_summary_cell(2:end,1) = sample_names;
norm_summary_cell(2:end,2) = num2cell(median_intensity_mat(:,1));
norm_summary_cell(2:end,3) = num2cell(median_intensity_mat(:,2));
norm_summary_cell(2:end,4) = num2cell(array_outlier_vec);

saveCellFile(norm_summary_cell, fullfile(user_dir, ['Normalization_and_HMM_summary_' lower(chip) '.txt']));


