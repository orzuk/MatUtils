%function [average_copy_num_mat, chr_num_snps_vec, data_snp_ids] =
%load_hmm_genotype_mat(in_sample_names, user_dir, chip_type)
function [average_copy_num_mat, chr_num_snps_vec, data_snp_ids] = load_hmm_avg_copy_mat(in_sample_names, user_dir, chip_type)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Load hmm output for all samples
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f_name = fullfile(user_dir, 'display', ['all_samples_hmm_out_' chip_type '.mat']);
load(f_name, 'chr_num_snps_vec', 'average_copy_num_mat', 'sample_names');

% find wanted samples
[C, IA, IB] = intersect_order_by_first_gr(in_sample_names, sample_names);

average_copy_num_mat = single(average_copy_num_mat(:, IB));

if(nargout>2)
    load(f_name, 'data_snp_ids');
end