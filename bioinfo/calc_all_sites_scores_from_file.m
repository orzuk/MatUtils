% Scan a set of sequences given as an input file,
% and compute the simple pwm match score for each one of them
%
% Input:
% pwms_file - file containing the pwms
% regions_file - file containing the regions
% seqs_dir - directory where genomic sequence is present (optional)
% genome_version - what genome version to use (optional)
% scores_outfile - where to save the scores output file 
% do_log - if to do log to the pwms (optional)
% background_model - enable a background model subtraction (default: NONE!)
% strand - which strand to run on (default is only positive strand)
%
% Output: 
% calc_time - the time it took to scan the motifs on the regions (without all pre/post-processing)
%
function calc_time = calc_all_sites_scores_from_file(pwms_file, regions_file, seqs_dir, ...
    genome_version, scores_outfile, do_log, background_model, strand, varargin)

Assign24MammalsGlobalConstants();

if((~exist('seqs_dir', 'var')) || isempty(seqs_dir)) % set default seqs dir (it needs to contain matlab files)
    seqs_dir = fullfile('/seq/orzuk/24mammals/data', genome_version);
end
if(~exist('background_model', 'var'))
    background_model = [];
end
if(~exist('strand', 'var'))
    strand = 0;
end
if(strmatch(pwms_file(end-3:end), '.mat', 'exact')) % read pwms
    load(pwms_file); % must have a variable named 'pwms'
else % read a .txt file
    pwms = load_pwms_from_txt_file(pwms_file);
end
num_pwms = size(pwms,1);

switch regions_file(end-2:end)
    case 'mat' % Read regions file to get the regions to work on
        load(regions_file); % assume format is: chr_vec, pos_start_vec, pos_end_vec
    case '.fa' % fasta file containing already the regions themselves
        [Header SEQS.seqs] = fastaread(regions_file);
    otherwise % .txt file with the regions. Assume a .txt file with format: chr_vec, pos_start_vec, pos_end_vec
        load(regions_file, '-ASCII');
end

if(~exist('SEQS', 'var'))
    seq_flag = 1; loglike_flag = 0; positions_flag = 0; tree_flag = 0; branch_length_flag = 0; % get only the sequence
    SEQS = ExtractSeqsByPositions(chr_vec, pos_start_vec, pos_end_vec, seqs_dir, genome_version, ...
        seq_flag, loglike_flag, positions_flag, tree_flag, branch_length_flag);  % Reading input sequence
end

calc_time = cputime % count time including sequence packing 
[SEQS.packed_seqs, SEQS.seqs_len] = pack_seqs(SEQS.seqs);
if(iscell(SEQS.seqs_len))
    SEQS.seqs_len = cell2mat(SEQS.seqs_len);
end

% Calling the function to run the scanning
% do both strands
[sites_scores rev_comp_sites_scores] = ...
    calc_all_sites_scores_matlab(pwms(:,2), SEQS.packed_seqs, SEQS.seqs_len, do_log, background_model, strand); % run searching for sites

calc_time = cputime - calc_time % count time including sequence packing 

% if(~exist('chr_vec', 'var'))
%     chr_vec = Header;
%     pos_start_vec = Header;
%     pos_end_vec = Header;
% end

if(strmatch(scores_outfile(end-3:end), '.mat', 'exact')) % save scores to output .mat file
    save(scores_outfile, 'sites_scores', 'rev_comp_sites_scores');
    if(exist('chr_vec', 'var'))
        save(scores_outfile, 'chr_vec', 'pos_start_vec', 'pos_end_vec', '-append');
    end
    if(exist('Header', 'var'))
        save(scores_outfile, 'Header', '-append');
    end
else % save as txt file: one .txt per region
    if(exist('chr_vec', 'var'))
        len = length(chr_vec);
        for i=1:len % loop on different regions
            cur_file_name{i} = [scores_outfile(1:end-4) '.chr' chr_num2str(chr_vec(i), genome_version) ':' ...
                num2str(pos_start_vec(i)) '-' num2str(pos_end_vec(i))];
        end
    else
        len = 1;
        cur_file_name{1} = [scores_outfile(1:end-4) '.head.' Header];
    end
    for i=1:len
        for j=1:length(sites_scores) % loop on different pwms
            if(iscell(sites_scores{j}))
                temp_scores = [ vec2column(double(sites_scores{j}{i})) vec2column(double(rev_comp_sites_scores{j}{i})) ];
            else
                temp_scores = [ vec2column(double(sites_scores{j})) vec2column(double(rev_comp_sites_scores{j})) ];
            end
%             temp_scores = num2str(temp_scores); 
            save([cur_file_name{i} '.' pwms{j,1} '.txt'], 'temp_scores', '-ASCII');
        end
    end
end









