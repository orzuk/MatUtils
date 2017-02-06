% Read a generic parameters file in a specified standard form.
%
% Input: 
% params_file - input file in a specific format
% 
% Output: 
% params_struct - a structure containing all the parameters appearing in the input file
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

function params_struct = ReadParametersFile(params_file) % read parameters

addutils; % add several utilities (specially for manipulating strings) 
Assign24MammalsGlobalConstants;

params_struct = []; % just put an empty variable at the beginning in order not to upset matlab
text_vec = textread(params_file, '%s', 'delimiter', '\n'); % read text into strings vector 
num_fields = size(text_vec, 1); % number of fields
str_vecs = cell(num_fields,1); 
for i=1:num_fields
    str_vecs{i} = strread(text_vec{i}, '%s'); % seperate to words in each line 
end

% Here try to make the function generic: The first field is the variable
% name, the second field is the variable value (could be many) 

for i=1:num_fields
    comment_start = strmatch('%', str_vecs{i}); % parse the strings in str_vecs
    if(isempty(comment_start))
        in_str = str_vecs{i}(2:end);
    else
        in_str = str_vecs{i}(2:comment_start(1)-1);
    end    
    if(~isempty(in_str))
        in_num = str2num_cell(in_str); % convert whenever possible strings to nums
        numeric_flag = 1;
        for j=1:length(in_str)
            if(~isempty(in_num{j}))
                in_str{j} = in_num{j};
            else
                numeric_flag = 0; % at least one non-numeric variable
            end
        end
        if(numeric_flag) % if all are numeric (each has one number!)
            in_str = cell2mat(in_str); 
        end
        if(iscell(in_str) && (length(in_str) == 1))
            eval_str = ['params_struct.' strdiff(str_vecs{i}{1}, ':') ' = in_str{1};'];
        else
            eval_str = ['params_struct.' strdiff(str_vecs{i}{1}, ':') ' = in_str;'];
        end
        eval(eval_str);
    end
end

