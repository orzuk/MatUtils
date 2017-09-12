%function comp_proj_common_aberr(folder1, folder2, proj1, proj2, params_str, chips, n_d)
function comp_proj_common_aberr(folder1, folder2, proj1, proj2, params_str, chips, n_d)

stretches_flag = 1;
genome_assembly = get_genome_assembly();
load(fullfile('..','database',['refgenes_' genome_assembly '.mat']), 'gene_symbols');

num_genes_all = length(gene_symbols);

if(nargin<5)
    params_str = '';
end
if(nargin<6)
    chips = {'Xba';'Hind'};
end
if(nargin<7)
    n_d = '_d';
end
if(strcmp(n_d,'_d')
    n_d_str = '_d';
else
    n_d_str = '_n';
end

num_chips = length(chips);

hmm_del_cell_n_d = cell(num_chips, 2);
hmm_amp_cell_n_d = cell(num_chips, 2);

pair_first_line = {'proj1' 'proj2' 'Del1' 'Del2' 'Del intersect' 'Amp1' 'Amp2' 'Amp intersect'};

for i = 1:num_chips
    hmm_del_f_n_d1 = ['DEL_genes_' chips{i} n_d_str params_str '.txt'];
    hmm_amp_f_n_d1 = ['AMP_genes_' chips{i} n_d_str params_str '.txt'];
    hmm_del_cell_n_d{i, 1} = loadCellFile_str([folder1 hmm_del_f_n_d1]);
    hmm_amp_cell_n_d{i, 1} = loadCellFile_str([folder1 hmm_amp_f_n_d1]);
    hmm_del_f_n_d2 = ['DEL_genes_' chips{i} n_d_str params_str '.txt'];
    hmm_amp_f_n_d2 = ['AMP_genes_' chips{i} n_d_str params_str '.txt'];
    hmm_del_cell_n_d{i, 2} = loadCellFile_str([folder2 hmm_del_f_n_d2]);
    hmm_amp_cell_n_d{i, 2} = loadCellFile_str([folder2 hmm_amp_f_n_d2]);
end

for i = 1:num_chips
    chip = char(chips{i});
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First compare disease results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ret_cell = cell(0,length(pair_first_line));

    hmm_del1_n_d = hmm_del_cell_n_d{i, 1};
    hmm_del2_n_d = hmm_del_cell_n_d{i, 2};
    [C_n_d, IA_n_d, IB_n_d] = intersect(hmm_del2_n_d(2:end,1), hmm_del1_n_d(2:end,1));
    del_len1_n_d = size(hmm_del1_n_d, 1)-1;
    del_len2_n_d = size(hmm_del2_n_d, 1)-1;
    del_intersect_n_d = length(C_n_d);
    del_intersect_p_n_d = give_one_cluster_TFs_pvalues(del_len1_n_d, del_intersect_n_d, num_genes_all, del_len2_n_d, 1);

    hmm_amp1_n_d = hmm_amp_cell_n_d{i, 1};
    hmm_amp2_n_d = hmm_amp_cell_n_d{i, 2};
    [C_n_d, IA_n_d, IB_n_d] = intersect(hmm_amp2_n_d(2:end,1), hmm_amp1_n_d(2:end,1));
    amp_len1_n_d = size(hmm_amp1_n_d, 1)-1;
    amp_len2_n_d = size(hmm_amp2_n_d, 1)-1;
    amp_intersect_n_d = length(C_n_d);
    amp_intersect_p_n_d = give_one_cluster_TFs_pvalues(amp_len1_n_d, amp_intersect_n_d, num_genes_all, amp_len2_n_d, 1);

    pair_cell = cell(2,length(pair_first_line));
    pair_cell(1,:) = pair_first_line;
    pair_cell{2,1} = erase_backslash(proj1); pair_cell{2,2} = erase_backslash(proj2);
    pair_cell{2,3} = del_len1_n_d; pair_cell{2,4} = del_len2_n_d;
    pair_cell{2,5} = [num2str(del_intersect_n_d) ' ( ' num2str(del_intersect_p_n_d) ' )'];
    pair_cell{2,6} = amp_len1_n_d; pair_cell{2,7} = amp_len2_n_d;
    pair_cell{2,8} = [num2str(amp_intersect_n_d) ' ( ' num2str(amp_intersect_p_n_d) ' )'];
    ret_cell_n_d = concat_cells(ret_cell, pair_cell, 2);

end

num_rows = size(ret_cell_n_d, 1);

% if(stretches_flag)
%     s_ret_cell = comp_d_n_stretches(params_str, chips, app_folder);
%     if(size(s_ret_cell,2) < size(ret_cell,2))
%         empty_cell = cell(size(s_ret_cell, 1), size(ret_cell,2)-size(s_ret_cell,2));
%         empty_cell(:) = {' '};
%         s_ret_cell = concat_cells(s_ret_cell, empty_cell, 1);
%     end
%     ret_cell = concat_cells(ret_cell, s_ret_cell, 2);
% end

saveCellFile(ret_cell_n_d, [folder1  'Comp_' proj1 '_' proj2 '_common_aberrations' params_str '.txt']);

