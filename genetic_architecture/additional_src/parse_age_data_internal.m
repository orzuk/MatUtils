% Convert allelic age data to .mat format 
function parse_age_data_internal(age_file_name)

run_file_name = age_file_name
R = loadcellfile(age_file_name);
num_iters = size(R,1);
num_alleles_vec = cell2mat(R(:,1));
allele_freq_vec = cell(num_iters,1);
allele_age_vec = cell(num_iters,1);

for j=1:num_iters % loop on different iterations
    allele_freq_vec{j} = cell2mat(R(j,2:2:end));
    allele_age_vec{j} = cell2mat(R(j,3:2:end));
end
save(fullfile(dir_from_file_name(age_file_name), 'mat', ...
    [remove_dir_from_file_name(age_file_name) '.mat']), ...
    'num_iters', 'allele_freq_vec', 'allele_age_vec', 'num_alleles_vec');

