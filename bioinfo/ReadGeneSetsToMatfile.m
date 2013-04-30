% Convert GSEA gmt format gene sets to .mat files 
%
% Input: 
% genesets_file - a gmt file with gene sets name
% output_file - where to save the output
% 
function ReadGeneSetsToMatfile(genesets_file, output_file, varargin)

addutils;

if(~exist('output_file', 'var'))
    output_file = [genesets_file(1:end-4) '.mat'];
end

%R = textread(genesets_file, '%s', 'delimiter', ' ');  % read the genesets 
R = loadcellfile(genesets_file); % textread doesn't work well here (overflow)

geneset_name = R(:,1);
geneset_description = R(:,2);
num_sets = size(R,1);
gene_symbol = cell(num_sets,1); 
for i=1:num_sets
    last_full_ind = find(~isempty_cell(R(i,:)), 1, 'last');
    gene_symbol{i} = R(i,3:last_full_ind);
end

% Now getting the type (trickiest part)

type_vec = strsplit(remove_dir_from_file_name(genesets_file), '.'); % NEW! use matlab strsplit (first string, then delimiter)

switch type_vec{1} % determine how was the gene set constructed
    case 'c1'
        geneset_type = 'positional';
    case 'c2'
        geneset_type = 'curated';
    case 'c3'
        geneset_type = 'motif';
    case 'c4'
        geneset_type = 'expression';
    case 'c5'
        geneset_type = 'GO';
    otherwise 
        geneset_type = type_vec{1}; % unknown new type
end

geneset_type = [geneset_type '_' type_vec{2}]; % add a more specific sub-type 
    

save(output_file, 'geneset_type', 'geneset_name', 'geneset_description', 'gene_symbol', 'num_sets'); 



