% This function reads a generic parameters file in a specified standard form.
%
% Input: 
% params_struct - a structure containing all the parameters appearing in the input file
% params_file - file name to save parameters  
% 
% The first field gives the name of the variable. Then the rest of the
% fields, until the first '%' (indicating a comment) are read into this variable.
% Here is an example for a 'default' input parameters file which we read
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% data_dir:               /seq/orzuk/24mammals/data/ncRNA     % directory with (which kind of?) data 
% regions_dir:            /seq/orzuk/24mammals/data/ncRNA     % directory with regions file
% regions_file:           ESMEFMLFNPCNovelK4.unique_mm9.mat   % file containing the regions
% background_genome_dir:  /seq/orzuk/24mammals/data/mm9/  % directory of beckground genome
% genome_version:         mm9            % genome version (irrelevant?) 
% method:                 KNOWN_PWMS     % method of finding motifs 
% background_method       ENRICHMENT     % method of comparing to background 
% k_vec:                  8              % length of kmers to use 
% alpha_fdr:              0.05           % allowed fdr level for found motifs         
% top_k:                  1000           % maximum number of motifs to keep
% outfile:                motifs.mat     % output file with found motifs  
% do_counting:            1              % if to count kmers statistics 
% compute_pval_flag:      1              % if to compute p-values for kmers statistics
% conservaion:            0.05           % cut-off above which we consider things as conserved
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function WriteParametersFile(params_struct, params_file) 

addutils; % add several utilities (specially for manipulating strings) 
Assign24MammalsGlobalConstants;

% text_vec = textread(params_file, '%s', 'delimiter', '\n'); % read text into strings vector 
str_vecs = fieldnames(params_struct); 
num_fields = length(str_vecs); 
fid = fopen(params_file, 'w'); 
for i=1:num_fields
    eval_str = ['str_vecs{i} = [str_vecs{i} '':'' tab num2str(params_struct.' str_vecs{i} ')];'];
    eval(eval_str); 
    fprintf(fid, '%s\n', str_vecs{i}); 
end
fclose(fid); 


