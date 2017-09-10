% save fitted gene-specific statistics to file. Tmp: seperate file 
% 
% Input: 
% GeneStruct - structure with original and fitted gene parameters
% output_data_dir - where to save output
% 
% Output: 
% None - output is saved to files
% 
function internal_save_gene_stats(GeneStruct, ExonsGeneStruct, output_data_dir, gene_struct_input_file, exome_struct)
 
save(fullfile(output_data_dir, ...
    [remove_suffix_from_file_name(remove_dir_from_file_name(gene_struct_input_file)) '_fitted_stats.mat']), '-struct', 'GeneStruct'); % save .mat file 

R_header = ['Gene' exome_struct.pop_str(1:end-1) exome_struct.pop_str(1:end-1)]; % TMP! remove last one 

% R = cell(GeneStruct.num_genes, 10); 
R = [GeneStruct.gene_names GeneStruct.s_MLE_vec GeneStruct.alpha_MLE_vec];
R = [R_header' R']'; 

%R = [GeneStruct.gene_names GeneStruct.s_MLE_vec
%for i_pop=1:num_poulations
%    R = [R num2str(GeneStruct.s_MLE_vec(:,i_pop))]; 
%end

savecellfile(R, fullfile(output_data_dir, ...
    [remove_suffix_from_file_name(remove_dir_from_file_name(gene_struct_input_file)) '_fitted_stats.txt'])); % save .txt file
