% Save expression data matrix in .txt format
%
% Input: 
% labels - sample names
% affy_id - affy ids of genes
% gene_symbol - gene symbols of genes
% data - data matrix
% output_file - where to save the data
%
function save_expression_mat_as_text(labels, affy_id, gene_symbol, data, output_file)

n = length(gene_symbol);
m = length(labels);
R = cell(n+1,m+2); R{1,1} = 'affy ids'; R{1,2} = 'gene symbols';
R(1,3:end) = labels; % first line is labels
R(2:end,1) = affy_id; % first column is affy ids
R(2:end,2) = gene_symbol; % second column is gene symbol
for i=1:n
    for j=1:m
        R{i+1,j+2} = data(i,j);
    end
end
savecellfile(R, output_file);




