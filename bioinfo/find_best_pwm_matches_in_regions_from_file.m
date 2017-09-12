% Look for motif hits given coordinates, with input from a file
% and then just call the function find_best_pwm_matches_in_sequences_pack
% and store the output in an output file. Supported are two modes of
% operation: 
% 1. find the best hit in each region 
% 2. find the set of best hits over all 
%
% Input: 
% regions_file - file with regions coordinates
% pwms_file - file with pwms 
% params_file - file with the following parameters: 
%       genome_version: what version to use for coordinates   
%       do_log: if to do log-transform to the pwm
%       conservation_flag: if to consider only conserved binding sites
%       alpha: fraction of conservation cutoff
%       bs_output_file: where to write the binding sites  (optional, can be .mat, .txt or none)
% thresholds - cutoff above which a motif score is considered a 'hit'
% bs_output_file - where to save the binding sites 
% pwm_ind - we can run only a certain pwm 
%
% Output: 
% BS_regions - regions indices
% BS_positions - position within region 
% BS_scores - score for matching the binding sites
% BS_seqs - the sequence of the binding sites themselves
% BS_strand - strand in which site was found
% BS_pwms - indices of the pwms for each site
% BS_loglike - the loglikelihood conservation score for each instance found
%
function [BS_regions BS_positions BS_scores BS_seqs BS_strand BS_pwms BS_loglike] = ...
    find_best_pwm_matches_in_regions_from_file(regions_file, pwms_file, params_file, ...
    thresholds, bs_output_file, pwms_ind, varargin)

if(ischar(params_file))
    params_struct = ReadParametersFile(params_file); % read parameters. We need: genome_version, do_log, strand, conservation_flag, alpha, bs_output_file, alignment_flag
else
    params_struct = params_file;
end
params_struct_is = params_struct
genome_version = params_struct.genome_version; 

Assign24MammalsGlobalConstants; % assign constants

switch regions_file(end-2:end) % load regions
    case 'mat' % Read regions file to get the regions to work on
        regions = load(regions_file); % assume format is: chr_vec, pos_start_vec, pos_end_vec
    case '.fa' % fasta file containing already the regions themselves
        [Header SEQS.seqs] = fastaread(regions_file);
    otherwise % .txt file with the regions. Assume a .txt file with format: chr_vec, pos_start_vec, pos_end_vec
        regions = load(regions_file, '-ASCII');
end
if(ischar(pwms_file)) % enable input from a variable
    if(strmatch(pwms_file(end-3:end), '.mat', 'exact')) % read pwms
        load(pwms_file); % must have a variable named 'pwms'
    else % read a .txt file
        pwms = load_pwms_from_txt_file(pwms_file);
    end
else
    pwms = pwms_file;
end
num_pwms = size(pwms,1);
% regions = load(regions_file); 
% load(pwms_file);  % load pwms

if(exist('pwms_ind', 'var')) % read command-line input (over-ride parameters file)
    params_struct.pwms_ind = pwms_ind;
end
if(isfield(params_struct, 'pwms_ind')) % take only a subset of the pwms
    if(~isempty(params_struct.pwms_ind))
        if(~isnumeric(params_struct.pwms_ind))
            params_struct.pwms_ind = strmatch(params_struct.pwms_ind, pwms(:,1));
        end
        pwms = pwms(params_struct.pwms_ind,:); 
    end
end
seq_flag = 1; % this is needed
positions_flag = 0; tree_flag = 0; branch_length_flag = 0; % all these are not needed
seqs_dir = fullfile('../data/', params_struct.genome_version);
if(~exist('SEQS', 'var'))
    tic;
%    chr_is = regions.chr_vec
%    start_is = regions.pos_start_vec 
%    end_is = regions.pos_end_vec
%    dir_is = seqs_dir
%    ver_is = params_struct.genome_version
    SEQS = ExtractSeqsByPositions(regions.chr_vec, regions.pos_start_vec, regions.pos_end_vec, ...
        seqs_dir, params_struct.genome_version, ...
        seq_flag, params_struct.conservation_flag, positions_flag, tree_flag, branch_length_flag); % This could be a problem if we take a whole genome ???
    if(isfield(regions, 'score_vec')) % also copy the scores
        SEQS.score_vec = regions.score_vec(SEQS.sort_perm); 
    end
    toc % see how much time does sequence extraction take
end
num_seqs = length(SEQS.seqs);
if(~isfield(SEQS, 'packed_seqs'))
    [SEQS.packed_seqs, SEQS.seqs_lens] = pack_seqs(SEQS.seqs); % pack sequences (time consuming)
end
if(iscell(SEQS.seqs_lens))
    SEQS.seqs_lens = cell2mat(SEQS.seqs_lens);
end

% Call the function to do the binding site search
%%smart_thresholds = 1; % allow choosing thresholds automatically
diff_lengths_flag = 1; % assume always different lengths
if(~isfield(SEQS, 'loglike')) % check if we have conservation information 
    SEQS.loglike = [];
end
if(~isfield(params_struct, 'alignment_flag'))
    params_struct.alignment_flag = 0;
end
for i=1:size(pwms,1) % perform derichlet correction in case some of the pwms have zeros within them 
    if(min(pwms{i,2}(:)) == 0)
        pwms{i,2} = derich_correct(pwms{i,2}, derich_alpha);
    end
end
if(isfield(params_struct, 'background')) % here correct according to a background model
    for i=1:size(pwms,1)
        if(params_struct.do_log) % correct multiplicatively
            pwms{i,2} = pwms{i,2} ./ repmat(vec2column(params_struct.background), 1, size(pwms{i,2},2));
        else % correct additively
            pwms{i,2} = pwms{i,2} - repmat(vec2column(params_struct.background), 1, size(pwms{i,2},2));
        end            
    end
end
tic;
if(~isfield(params_struct, 'bs_overlap'))
	params_struct.bs_overlap = [];
end 
if(~isfield(params_struct, 'save_format'))
	params_struct.save_format = {'mat'};
end 
if(~isfield(params_struct, 'background_model'))
	params_struct.background_model = [];
end 

call_pwms_packed_pwms_is = pwms
[BS_regions BS_positions BS_scores BS_seqs BS_strand BS_pwms BS_loglike] = ...
    find_best_pwm_matches_in_sequences_pack(SEQS.packed_seqs, SEQS.seqs_lens, pwms, params_struct.top_sites, ...
    params_struct.threshold_type, params_struct.do_log, params_struct.background_model, params_struct.strand, diff_lengths_flag, ... 
    params_struct.bs_overlap, SEQS.loglike, params_struct.conservation_flag, params_struct.alignment_flag);
num_BS = size(BS_regions) % see if we really got enough sites 
toc % see how much time does motif finding take 

if(~isempty(strmatch('mat', lower(params_struct.save_format))))  % save as mat file
    my_mkdir(dir_from_file_name(bs_output_file));
    save(file_name_to_mat(bs_output_file), 'SEQS', 'BS_regions', 'BS_positions', 'BS_scores', 'BS_seqs', 'BS_strand', ...
        'BS_pwms', 'BS_loglike', 'pwms', 'genome_version'); % save the results into a .mat file
    if(isfield(params_struct, 'background')) 
        background = params_struct.background;
        save(bs_output_file, 'background', '-append');
    end
%    Dummy = AddRegionTypeToBSFile(bs_output_file, params_struct.genome_version); % quickly just add the regions type 
end
if(~isempty(strmatch('txt',  lower(params_struct.save_format))))  % save also txt file
	if(isfield(SEQS, 'score_vec'))
    score_vec = SEQS.score_vec; SEQS = rmfield(SEQS, 'score_vec'); 
else 
score_vec = zeros(length(BS_regions{1}),1);
	end
%    BS_binary = double( <= length(BS_regions{1})/2); % Temp Patch!!! 
    Dummy = save_bs_mat_as_text(SEQS, BS_regions, BS_positions, BS_scores, BS_seqs, ...
        BS_strand, BS_loglike, pwms, file_name_to_txt(bs_output_file), ...
        score_vec(BS_regions{1}), {'BS-region'}); % save results into also .txt file
end

