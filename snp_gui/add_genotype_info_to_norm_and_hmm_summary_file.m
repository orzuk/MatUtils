% Add genotype inforation to summary file 
%function add_genotype_info_to_norm_and_hmm_summary_file(user_dir, chip_type, genotype_mat, sample_names)
function add_genotype_info_to_norm_and_hmm_summary_file(user_dir, chip_type, genotype_mat, sample_names)

% first load the file created by the normalization function
norm_summary_f = fullfile(user_dir, ['Normalization_and_HMM_summary_' lower(chip_type) '.txt']);
if(exist(norm_summary_f,'file'))
    table = loadCellFile_str(norm_summary_f);
    % extract the sample_names from the file
    f_sample_names = table(2:end,1);
    num_samples = length(f_sample_names);
    add_table = cell(num_samples+1,4);
    add_table{1,1} = '%AA'; add_table{1,2} = '%AB'; add_table{1,3} = '%BB';
    add_table{1,4} = '%Identity';
    [sample_pairs_cell, pairs_samples_ind] = samples_into_pairs(f_sample_names);
    [C, IA, IB] = intersect_order_by_first_gr(f_sample_names, sample_names);
    genotype_mat = genotype_mat(:, IB);
    [AA, AB, BB, NoCall, call_cell] = genotype_call_into_num();
    genotype_mat = hmm_geno_into_old(genotype_mat);
    num_snps = size(genotype_mat, 1);
    for i = 1:num_samples
        add_table{i+1,1} = ceil((length(find(genotype_mat(:,i)==AA))/num_snps)*100);
        add_table{i+1,2} = ceil((length(find(genotype_mat(:,i)==AB))/num_snps)*100);
        add_table{i+1,3} = ceil((length(find(genotype_mat(:,i)==BB))/num_snps)*100);
    end
    num_samples_pairs = size(sample_pairs_cell, 1);
    if(num_samples_pairs)
        for i = 1:num_samples_pairs
            disease_ind = pairs_samples_ind(i,2);
            normal_ind = pairs_samples_ind(i,1);
            identity_frac = (length(find(genotype_mat(:,disease_ind)==genotype_mat(:,normal_ind)))/num_snps)*100;
            add_table{disease_ind+1,4} = num2str(identity_frac,4);
        end
    end
    table = concat_cells(table, add_table, 1);
    saveCellFile(table, norm_summary_f);
end