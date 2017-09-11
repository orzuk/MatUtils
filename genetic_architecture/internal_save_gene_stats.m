% save fitted gene-specific statistics to file. Tmp: seperate file 
% 
% Input: 
% GeneStruct - structure with original and fitted gene parameters
% output_data_dir - where to save output
% exome_struct - structure with annotation information
% gene_fit_I - indices of fitted genes (save only these!)
%
% Output: 
% None - output is saved to files
% 
function internal_save_gene_stats(GeneStruct, ExonsGeneStruct, output_data_dir, gene_struct_input_file, exome_struct, fit_genes_I)

GeneStruct.fit_genes_I = fit_genes_I; 
save(fullfile(output_data_dir, ...
    [remove_suffix_from_file_name(remove_dir_from_file_name(gene_struct_input_file)) '_fitted_stats.mat']), '-struct', 'GeneStruct'); % save .mat file 

R_header = ['Gene' strcat('s_MLE.', exome_struct.pop_str(1:2)) strcat('alpha_MLE.', exome_struct.pop_str(1:2))]; % TMP! remove last one 
R = [GeneStruct.gene_names(fit_genes_I) num2str_cell(num2cell(GeneStruct.s_MLE_vec(fit_genes_I,:))) num2str_cell(num2cell(GeneStruct.alpha_MLE_vec(fit_genes_I,:)))];
R = [R_header' R']'; 
savecellfile(R, fullfile(output_data_dir, ...
    [remove_suffix_from_file_name(remove_dir_from_file_name(gene_struct_input_file)) '_fitted_stats.txt']), [], 1); % save .txt file
