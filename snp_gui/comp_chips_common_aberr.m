function comp_chips_common_aberr(params_str, chips)

stretches_flag = 1;
genome_assembly = get_genome_assembly();
load(fullfile('..','database',['refgenes_' genome_assembly '.mat']), 'gene_symbols');

num_genes_all = length(gene_symbols);


if(nargin<2)
    chips = {'Xba';'Hind'};
end

if(nargin<1)
    params_str = '';
end


%app_folder = 'E:\Libi\tools\SNP_tool\data\Leukemia_all\output\Leukemia_all_ucsc\';
app_folder = 'E:\Libi\tools\SNP_tool\data\Leukemia_all\output\';

num_chips = length(chips);

hmm_del_cell_d = cell(num_chips, 1);
hmm_amp_cell_d = cell(num_chips, 1);

num_pairs = num_chips*(num_chips-1)/2;
pair_first_line = {'Chip1' 'Chip2' 'Del1' 'Del2' 'Del intersect' 'Amp1' 'Amp2' 'Amp intersect'};

for i = 1:num_chips
    hmm_del_f_d = ['DEL_genes_' chips{i} '_d' params_str '.txt'];
    hmm_amp_f_d = ['AMP_genes_' chips{i} '_d' params_str '.txt'];
    hmm_del_cell_d{i} = loadCellFile_str([app_folder hmm_del_f_d]);
    hmm_amp_cell_d{i} = loadCellFile_str([app_folder hmm_amp_f_d]);
end

for i = 1:num_chips
    hmm_del_f_n = ['DEL_genes_' chips{i} '_n' params_str '.txt'];
    hmm_amp_f_n = ['AMP_genes_' chips{i} '_n' params_str '.txt'];
    hmm_del_cell_n{i} = loadCellFile_str([app_folder hmm_del_f_n]);
    hmm_amp_cell_n{i} = loadCellFile_str([app_folder hmm_amp_f_n]);
end

fdr_vec = [0.01:0.01:0.2:0.1:1];
for i = 1:num_chips
    for j = i+1:num_chips
        chip1 = char(chips{i});
        chip2 = char(chips{j});
        comp_chips_common_aberr_graph(chip1, chip2, app_folder, '_d', 'DEL', params_str, fdr_vec);
        comp_chips_common_aberr_graph(chip1, chip2, app_folder, '_d', 'AMP', params_str, fdr_vec);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % First compare disease results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ret_cell = cell(0,length(pair_first_line));

        hmm_del1_d = hmm_del_cell_d{i};
        hmm_del2_d = hmm_del_cell_d{j};
        [C_d, IA_d, IB_d] = intersect(hmm_del2_d(2:end,1), hmm_del1_d(2:end,1));
        del_len1_d = size(hmm_del1_d, 1)-1;
        del_len2_d = size(hmm_del2_d, 1)-1;
        del_intersect_d = length(C_d);
        del_intersect_p_d = give_one_cluster_TFs_pvalues(del_len1_d, del_intersect_d, num_genes_all, del_len2_d, 1);

        hmm_amp1_d = hmm_amp_cell_d{i};
        hmm_amp2_d = hmm_amp_cell_d{j};
        [C_d, IA_d, IB_d] = intersect(hmm_amp2_d(2:end,1), hmm_amp1_d(2:end,1));
        amp_len1_d = size(hmm_amp1_d, 1)-1;
        amp_len2_d = size(hmm_amp2_d, 1)-1;
        amp_intersect_d = length(C_d);
        amp_intersect_p_d = give_one_cluster_TFs_pvalues(amp_len1_d, amp_intersect_d, num_genes_all, amp_len2_d, 1);

        pair_cell = cell(2,length(pair_first_line));
        pair_cell(1,:) = pair_first_line;
        pair_cell{2,1} = chip1; pair_cell{2,2} = chip2;
        pair_cell{2,3} = del_len1_d; pair_cell{2,4} = del_len2_d;
        pair_cell{2,5} = [num2str(del_intersect_d) ' ( ' num2str(del_intersect_p_d) ' )'];
        pair_cell{2,6} = amp_len1_d; pair_cell{2,7} = amp_len2_d;
        pair_cell{2,8} = [num2str(amp_intersect_d) ' ( ' num2str(amp_intersect_p_d) ' )'];
        ret_cell_d = concat_cells(ret_cell, pair_cell, 2);
        comp_chips_common_aberr_graph(chip1, chip2, app_folder, '_n', 'DEL', params_str, fdr_vec);
        comp_chips_common_aberr_graph(chip1, chip2, app_folder, '_n', 'AMP', params_str, fdr_vec);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Compare normals results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        ret_cell = cell(0,length(pair_first_line));

        hmm_del1_n = hmm_del_cell_n{i};
        hmm_del2_n = hmm_del_cell_n{j};
        [C_n, IA_n, IB_n] = intersect(hmm_del2_n(2:end,1), hmm_del1_n(2:end,1));
        del_len1_n = size(hmm_del1_n, 1)-1;
        del_len2_n = size(hmm_del2_n, 1)-1;
        del_intersect_n = length(C_n);
        del_intersect_p_n = give_one_cluster_TFs_pvalues(del_len1_n, del_intersect_n, num_genes_all, del_len2_n, 1);

        hmm_amp1_n = hmm_amp_cell_n{i};
        hmm_amp2_n = hmm_amp_cell_n{j};
        [C_n, IA_n, IB_n] = intersect(hmm_amp2_n(2:end,1), hmm_amp1_n(2:end,1));
        amp_len1_n = size(hmm_amp1_n, 1)-1;
        amp_len2_n = size(hmm_amp2_n, 1)-1;
        amp_intersect_n = length(C_n);
        amp_intersect_p_n = give_one_cluster_TFs_pvalues(amp_len1_n, amp_intersect_n, num_genes_all, amp_len2_n, 1);

        pair_cell = cell(2,length(pair_first_line));
        pair_cell(1,:) = pair_first_line;
        pair_cell{2,1} = chip1; pair_cell{2,2} = chip2;
        pair_cell{2,3} = del_len1_n; pair_cell{2,4} = del_len2_n;
        pair_cell{2,5} = [num2str(del_intersect_n) ' ( ' num2str(del_intersect_p_n) ' )'];
        pair_cell{2,6} = amp_len1_n; pair_cell{2,7} = amp_len2_n;
        pair_cell{2,8} = [num2str(amp_intersect_n) ' ( ' num2str(amp_intersect_p_n) ' )'];
        ret_cell_n = concat_cells(ret_cell, pair_cell, 2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % compare disease to normals
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [C_d_n, IA_d_n, IB_d_n] = intersect(hmm_del1_d(2:end,1), hmm_del1_n(2:end,1));
        del_intersect_d_n1 = length(C_d_n);
        del_intersect_p_d_n1 = give_one_cluster_TFs_pvalues(del_len1_d, del_intersect_d_n1, num_genes_all, del_len1_n, 1);
        [C_d_n, IA_d_n, IB_d_n] = intersect(hmm_amp1_d(2:end,1), hmm_amp1_n(2:end,1));
        amp_intersect_d_n1 = length(C_d_n);
        amp_intersect_p_d_n1 = give_one_cluster_TFs_pvalues(amp_len1_d, amp_intersect_d_n1, num_genes_all, amp_len1_n, 1);

        [C_d_n, IA_d_n, IB_d_n] = intersect(hmm_del2_d(2:end,1), hmm_del2_n(2:end,1));
        del_intersect_d_n2 = length(C_d_n);
        del_intersect_p_d_n2 = give_one_cluster_TFs_pvalues(del_len2_d, del_intersect_d_n2, num_genes_all, del_len2_n, 1);
        [C_d_n, IA_d_n, IB_d_n] = intersect(hmm_amp2_d(2:end,1), hmm_amp2_n(2:end,1));
        amp_intersect_d_n2 = length(C_d_n);
        amp_intersect_p_d_n2 = give_one_cluster_TFs_pvalues(amp_len2_d, amp_intersect_d_n2, num_genes_all, amp_len2_n, 1);

        d_n_cell = cell(1, size(ret_cell_n,2)+1); d_n_cell(:) = {' '};
        d_n_cell{1,1} = 'Intersect Disease-Normal';
        d_n_cell{1,4} = [num2str(del_intersect_d_n1) ' ( ' num2str(del_intersect_p_d_n1) ' )']; % chip1 del
        d_n_cell{1,5} = [num2str(del_intersect_d_n2) ' ( ' num2str(del_intersect_p_d_n2) ' )']; % chip2 del
        d_n_cell{1,7} = [num2str(amp_intersect_d_n1) ' ( ' num2str(amp_intersect_p_d_n1) ' )']; % chip1 amp        
        d_n_cell{1,8} = [num2str(amp_intersect_d_n2) ' ( ' num2str(amp_intersect_p_d_n2) ' )']; % chip2 amp
    end
end


ret_cell = concat_cells(ret_cell_d, ret_cell_n, 2);
num_rows = size(ret_cell, 1);
first_column = cell(num_rows, 1);
first_column(1:num_rows/2) = {'Disease'};
first_column(num_rows/2+1:end) = {'Normals'};
ret_cell = concat_cells(first_column, ret_cell, 1);
ret_cell = concat_cells(ret_cell, d_n_cell, 2);

if(stretches_flag)
    s_ret_cell = comp_d_n_stretches(params_str, chips, app_folder);
    if(size(s_ret_cell,2) < size(ret_cell,2))
        empty_cell = cell(size(s_ret_cell, 1), size(ret_cell,2)-size(s_ret_cell,2));
        empty_cell(:) = {' '};
        s_ret_cell = concat_cells(s_ret_cell, empty_cell, 1);
    end
    ret_cell = concat_cells(ret_cell, s_ret_cell, 2);
end

saveCellFile(ret_cell, [app_folder 'Comp_chips_common_aberrations' params_str '.txt']);

t=1;