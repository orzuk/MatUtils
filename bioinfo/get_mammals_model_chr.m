% Get all chromosome-specific model files
% 
% Input: 
% genome_version - what genome is used
% organism_flag - decides if we take local model or mammals model
% encode_flag - decides if to use models with encode namings (slightly different species names)
%
% Output: 
% mammals_model_chr_file - model file for a specific chromosome
% 
function mammals_model_chr_file = get_mammals_model_chr(genome_version, ...
    organism_flag, encode_flag, varargin)

Assign24MammalsGlobalConstants;
mammals_model_chr_file = cell(organism_chr+1, 1);

if(exist('organism_flag', 'var') && organism_flag)
    model_dir_used = local_model_dir;
else
    model_dir_used = mammals_model_dir;
end
if(exist('encode_flag', 'var') && encode_flag)
    model_dir_used = fullfile(model_dir_used, 'encode'); 
end
for i_chr=1:organism_chr % all autosomal and X
    mammals_model_chr_file{i_chr} = fullfile(model_dir_used, ['chr' chr_num2str(i_chr, genome_version) '_4D.mod']);
end
mammals_model_chr_file{organism_chr+1} = fullfile(model_dir_used, 'autosomal.mod'); % Y chromosome (don't have yet so use autosomal)

