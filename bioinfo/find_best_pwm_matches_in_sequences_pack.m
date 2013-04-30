% Find matches to position weight matrices in a given set of sequences.
% seqs is an array of the promoters sequences, represented as digits from 1
% to 4 or in a packed format. pwms is a cell array of the matrices we use in order to look for
% binding sites. thresholds is a vector of threshold determining when to
% accept a 'hit' as a binding site. (For each pwm we put a different
% threshold). smart_threshold is a flag, if it is one the program
% calculates the best threshold automatically and ignores the threshold
% vetor given by the user. do_log is a flag saying if to take the logaritm
% of the pwms. If they are for example probabilities, we need the do_log
% flag to be one.
%
% Input:
% seqs - a vector of sequences (in packed form)
% seqs_lens - the length of the sequences
% pwms - a set of pwms
% thresholds - which threshold to use to declare a hit
% smart_thresholds - flag saying if to determine the threshold alone
%   The options for smart threshold are the following: (see function get_scores_above_threshold_Matlab)
%       -1 - pick the single best instance of loc_scores (in each member of the cell array)
%       0 - pick the input threshold given from outside : 'threshold'
%       1 - choose threshold in a smart way. Treat 'threshold' as the top fraction above which we keep sites
%       2 - choose threshold in a smart way, such that exactly k scores pass, where k is the 'threshold' parameter
% do_log - flag saying if to do log to the pwms
% background_model - enable a background model subtraction (default: NONE!)
% strand - which strand to look for (could be both strand)
% diff_lengths_flag - flag saying if each sequence is of different length or not
% loglike_mat - used if we also use conservation cutoffs when calling binding sites
% alpha - the cutoff for declaring a site as 'conserved'
% alignments_flag - also include alignments of the bs
%
% Output:
% BS_regions - indices of the regions
% BS_positions - position within each region. The sign represents the strand
% BS_scores - scores of matching the pwm
% BS_seqs - the sequence at the BS
% BS_strand - the strand of the found binding sites
% BS_pwms - which pwms index correspond to each binding site
% BS_loglike - log-likelihood of omega conservation (optional)
% BS_alignments - multiple local alignment of the binding sites (optional)
%
function [BS_regions BS_positions BS_scores BS_seqs BS_strand BS_pwms BS_loglike BS_alignments] = ...
    find_best_pwm_matches_in_sequences_pack(seqs, seqs_lens, pwms, thresholds, ...
    smart_thresholds, do_log, background_model, strand, diff_lengths_flag, min_overlap, ...
    loglike_mat, alpha, alignment_flag, varargin)

Assign24MammalsGlobalConstants; % assign constants
tomy = cputime
if(iscell(pwms)) % adjust pwms as a column vector
    if(length(pwms) > 4)
        pwms = vec2column(pwms);
    else
        pwms = vec2row(pwms);
    end
    ends
pwms_is = pwms
pwms_size_is = size(pwms)
TFs = size(pwms, 1); % The number of pwms

overlap_flag = 0;
if(exist('min_overlap', 'var') && (~isempty(min_overlap)))
    if(isnumeric(min_overlap))
        overlap_flag = 1;
    end
end
if(exist('loglike_mat', 'var') && isempty(loglike_mat)) % check if to include conservation
    clear loglike_mat;
end
if(~exist('alpha', 'var') || isempty(alpha)) % default: don't take consrevation into account
    alpha = 1;
end
if(~exist('alignment_flag', 'var') || isempty(alignment_flag))
    alignment_flag = 0;
end


BS_pwms = []; % for the case where there is only one TF
BS_loglike = []; % for the case where no conservation information is given

% transfer cell to one long array
do_it_all = 1;  % do all regions at once and report the best score in each region
if((~do_it_all) && iscell(seqs))
    [seqs cumsum_seqs_lens_words] = cell2vec(seqs);
end

TOL = 0.00000000001;
% up_down = 'up';   % 'up' - an upperbound , 'down' - a lower bound on the score
% MAX_BS = 10;   % Maximal putative binding sites for one gene
% my_pval_thresh = pval; %%0.1;   % A threshold below it to take ........ need to be taken care of more carefully

N = max(seqs_lens); % N is the maximal length of the sequences. Here we give it from outside since the sequence is packed
num_genes = length(seqs_lens); % The number of genes prom other sequences


% Work chunks in order to avoid memory problems
MAX_SEQS_TOGETHER = min(5000,num_genes);  % Never take a block larger than the number of genes !!!
num_blocks = ceil(num_genes / MAX_SEQS_TOGETHER); % The number of blocks to work with

local_pwms = pwms % Copy to a local array
if(do_log)
    for t = 1:TFs
        local_pwms{t,2} = log(local_pwms{t,2}); % Perform log transform on the matrices if neccessary
    end
end
if(exist('background_model', 'var') && ~isempty(background_model))
    for t = 1:TFs
        local_pwms{t,2} = local_pwms{t,2} - repmat(log(background_model), 1, size(local_pwms{t,2},2));
    end
end


counter = 1;

% We need to be clever here, since seqs_lens is in nucleotides, and the
% seqs are packed so we need them in words
cum_seqs_lens = zeros(1, length(seqs_lens));
cum_seqs_lens(1) = floor((seqs_lens(1)-1)/16)+1;
for i=2:length(seqs_lens)
    cum_seqs_lens(i) = cum_seqs_lens(i-1) +  floor((seqs_lens(i)-1)/16)+1;
end

% Find best site in each promomter (?) Make generic place-holders for outputs
n = length(seqs_lens);
if(TFs > 1)
    BS_scores = cell(TFs,1); BS_regions = cell(TFs,1); BS_positions = cell(TFs,1);
    BS_rev_scores = cell(TFs,1); BS_rev_positions = cell(TFs,1); BS_rev_regions = cell(TFs,1);
    BS_strand = cell(TFs,1); BS_rev_strand = cell(TFs,1);
    BS_seqs = cell(TFs,1);
    BS_pwms = cell(TFs,1);
    if(~isempty(thresholds))
        if(length(thresholds) == 1)
            thresholds = repmat(thresholds, TFs, 1);
        end
    end
end
if(exist('loglike_mat', 'var')) % save also the conservation of each binding site
    BS_loglike = cell(TFs,1); % zeros(n,TFs);
end

num_TFs_to_search = TFs 
for t=1:TFs % loop over all the possible pwms, one by one
    sprintf('Doing TF %ld out of %ld', t, TFs)
    L = size(local_pwms{t,2}, 2);  % The length of the matrix binding site
    if(~overlap_flag) % set default overlap as L
        min_overlap = L;
    end
    loc_scores_lengths = max(seqs_lens-L+1, 0); % Avoid negative scores length. If sequence is shorter than pwm, put length=0
    
    do_it_all_flag = do_it_all
    if(do_it_all)  % do all regions at once
        [loc_scores rev_comp_loc_scores] = ...
            calc_all_sites_scores_matlab(local_pwms{t,2}, seqs, seqs_lens, 0, [], strand); % assume they're all equal and take the first
        if(exist('loglike_mat', 'var'))
            smoothed_loglike_mat = cell(length(loc_scores),1);
            for i=1:length(loc_scores)
                smoothed_loglike_mat{i} = my_smooth(loglike_mat{i}, single(L));
                if(isempty(smoothed_loglike_mat{i}))
                    smoothed_loglike_mat{i} = single(-L);
                end
            end
            if(alpha < 1) % threhsold conservation scores
                quantile_smoothed_loglike = quantile(cell2vec(smoothed_loglike_mat), 1-alpha);
                for i=1:length(loc_scores)
                    loc_scores{i}((smoothed_loglike_mat{i} < quantile_smoothed_loglike)) = single(-99999999999);
                    rev_comp_loc_scores{i}((smoothed_loglike_mat{i} < quantile_smoothed_loglike)) = single(-99999999999);
                end
            end
        end
        if(~iscell(loc_scores))
            temp = loc_scores; loc_scores = {}; loc_scores{1} = temp;
            temp = rev_comp_loc_scores; rev_comp_loc_scores = {}; rev_comp_loc_scores{1} = temp;
        end
        for i=1:length(loc_scores)
            if(isempty(loc_scores{i}))
                loc_scores{i} = -999999999;
                seqs_lens(i) = 1; % temporary patch
            end
            if(strand == 2) % do reverse strand
                if(isempty(rev_comp_loc_scores{i}))
                    rev_comp_loc_scores{i} = -999999999;
                    seqs_lens(i) = 1; % temporary patch
                end
            end
        end
        if(isempty(thresholds))
            local_threshold = 0.05;
        else
            local_threshold = thresholds(t);
        end
        [BS_regions{t} BS_positions{t} BS_scores{t} thresh_out] = ...
            get_scores_above_threshold_Matlab( loc_scores, loc_scores_lengths, local_threshold,...
            min_overlap, smart_thresholds, diff_lengths_flag); % This new code should replace the lines below
        num_scores = length(BS_regions{t}); % how many scores did we get
        BS_strand{t} = repmat(single(POS_STRAND), num_scores,1); % set default strand (positive)
        if(exist('loglike_mat', 'var')) % save also the conservation of each binding site
            BS_loglike{t} = zeros(num_scores,1, 'single');
            for i=1:num_scores
                BS_loglike{t}(i) = smoothed_loglike_mat{BS_regions{t}(i)}(BS_positions{t}(i));
            end
        end
        if(strand == 2) % deal also with the reverse strand
            [BS_rev_regions{t} BS_rev_positions{t} BS_rev_scores{t} thresh_out] = ...
                get_scores_above_threshold_Matlab( rev_comp_loc_scores, loc_scores_lengths, local_threshold,...
                min_overlap, smart_thresholds, diff_lengths_flag); % This new code should replace the lines below
            BS_rev_strand{t} = repmat(single(REV_STRAND), length(BS_rev_regions{t}),1); % set default strand (positive)
            
            BS_scores{t} = vec2column([vec2row(BS_scores{t}) vec2row(BS_rev_scores{t})]);
            BS_positions{t} = vec2column([vec2row(BS_positions{t}) vec2row(BS_rev_positions{t})]);
            BS_strand{t} = vec2column([vec2row(BS_strand{t}) vec2row(BS_rev_strand{t})]);
            BS_regions{t} = vec2column([vec2row(BS_regions{t}) vec2row(BS_rev_regions{t})]);
            
            %            if(smart_thresholds ~= 0) % in this case we just keep all the scores above a given threshold
            switch smart_thresholds
                case 0 % do nothing
                    
                case -1 % here we need to duplicate smartly (within each region)
                    [BS_scores{t} BS_strand{t}] = ...
                        max_with_inds(BS_scores{t}(num_scores+1:end), BS_scores{t}(1:num_scores));
                    BS_strand{t} = double(BS_strand{t});
                    BS_regions{t} = BS_regions{t}(1:num_scores);
                    BS_positions{t} = BS_positions{t}(1:num_scores) .* BS_strand{t} + ...
                        BS_positions{t}(num_scores+1:end) .* (1-BS_strand{t});
                otherwise
                    [sorted_scores sort_perm] = sort(BS_scores{t}, 'descend');
                    sort_perm = sort_perm(1:num_scores); % take only the top scores
                    BS_scores{t} = BS_scores{t}(sort_perm);
                    BS_positions{t} = BS_positions{t}(sort_perm);
                    BS_strand{t} = BS_strand{t}(sort_perm);
                    BS_regions{t} = BS_regions{t}(sort_perm);
            end
        end
        
        L = size(local_pwms{t,2},2);         % Find at which WORD each sequence starts
        BS_seqs{t} = positions_to_seqs(seqs, seqs_lens, ...
            BS_positions{t}, BS_regions{t}, L, diff_lengths_flag, BS_strand{t}, N);
        if(TFs > 1)
            BS_pwms{t} = t + zeros(length(BS_positions{t}),1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% This Part Doesn't Work %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else % go region by region and ... (doesn't work now (???))
        do_it_all_flag = do_it_all
        
        % We must here to use the mex functions
        % Pack the sequences ? Assume seqs are already packed!
        local_threshold = thresholds(t);
        local_smart_thresholds = smart_thresholds;
        
        for cur_block = 1:num_blocks % divide the work into blocks
            cur_block_size = min(MAX_SEQS_TOGETHER, num_genes-(cur_block-1)*MAX_SEQS_TOGETHER );
            
            % New ! (4.7.05) - We need to deal with different lengths !!!
            if(diff_lengths_flag == 0)
                is_double = isa(local_pwms{t,2}, 'double');
                loc_scores = calc_all_sites_scores(local_pwms{t,2}, ...
                    seqs((cur_block-1)*MAX_SEQS_TOGETHER+1:(cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size,:), seqs_lens(1), -10, is_double); % assume they're all equal and take the first
            else
                if(cur_block == 1)
                    start_index = 1;
                    end_index = cum_seqs_lens(cur_block_size);
                else
                    start_index = cum_seqs_lens((cur_block-1)*MAX_SEQS_TOGETHER)+1;
                    end_index = cum_seqs_lens((cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size);
                end
                
                % Note : Now the scores are returned as one big chain !!!
                if(size(seqs_lens,2) == 1)
                    loc_scores = calc_all_sites_scores_diff_lengths(local_pwms{t,2}, ...
                        seqs(start_index:end_index), ...
                        seqs_lens((cur_block-1)*MAX_SEQS_TOGETHER+1:(cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size)', -10); % New ! we send only one sequence, which contains
                else
                    loc_scores = calc_all_sites_scores_diff_lengths(local_pwms{t,2}, ...
                        seqs(start_index:end_index), ...
                        seqs_lens((cur_block-1)*MAX_SEQS_TOGETHER+1:(cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size), -10); % New ! we send only one sequence, which contains
                end
            end
            scanning_time = cputime - tomy;
            sprintf('Finished scanning block %ld with time %f on TF %d, now do smart thresholding ....\n', ...
                cur_block, scanning_time, t)
            
            [cur_BS_regions cur_BS_positions cur_BS_scores thresh_out] = get_scores_above_threshold_Matlab( loc_scores, ...
                loc_scores_lengths((cur_block-1)*MAX_SEQS_TOGETHER+1:(cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size), ...
                local_threshold, min_overlap, local_smart_thresholds, diff_lengths_flag );
            local_threshold = thresh_out;
            local_smart_thresholds = 0; % For the next blocks keep the same threshold !!!!
            len = length(cur_BS_positions); % see how many binding sites were found
            
            cur_BS_regions = cur_BS_regions + (cur_block-1)*MAX_SEQS_TOGETHER;   % Adjust genes to the current block
            % Avoid empty files. Put one that is like the threshold
            if(len == 0)
                len = 1;
                cur_BS_positions = 1;
                cur_BS_regions = 1;
                cur_BS_scores = thresh_out-TOL; %thresholds(t);
            end
            cur_BS_seqs = positions_to_seqs(seqs, seqs_lens, ...
                cur_BS_positions, cur_BS_regions, L, diff_lengths_flag, [], N);  % get bs sequences
            
            switch strand % Move location to the middle of the binding sites
                case 1 % positive strand
                    cur_BS_positions = cur_BS_positions + floor(L/2);
                    
                case 0 % reverse strand
                    cur_BS_positions = seqs_lens(cur_BS_regions) - cur_BS_positions - floor(L/2) + 1;
            end
            if(TFs > 1)
                cur_BS_pwms = t + zeros(length(cur_BS_scores),1);
            end
            if(cur_block == 1)
                BS_positions = cur_BS_positions;
                BS_regions = cur_BS_regions;
                BS_scores = cur_BS_scores;
                BS_seqs = cur_BS_seqs;
                if(TFs > 1)
                    BS_pwms = t + zeros(length(cur_BS_scores),1);
                end
            else
                BS_positions = [BS_positions cur_BS_positions];
                BS_regions = [BS_regions cur_BS_regions];
                BS_scores = [BS_scores cur_BS_scores];
                BS_seqs = [BS_seqs cur_BS_seqs];
                %            BS_strand = [BS_strand cur_BS_strand];
                if(TFs > 1)
                    BS_pwms = [BS_pwms cur_BS_pwms];
                end
            end
            counter = counter + len;             % update the counter
        end % loop on blocks
        threshold_used_was = local_threshold
        
        % Remove all those dummies
        [BS_unique u] = unique([BS_positions BS_regions BS_scores BS_seqs], 'rows');
        BS_positions = BS_positions(u);
        BS_regions = BS_regions(u);
        BS_scores = BS_scores(u);
        BS_seqs = BS_seqs(u,:);
        
        size_BS = size(BS_scores, 1)         % Remove the threshold if we have more than one !!
        if(size_BS > 1)
            dummy_ind = find(BS_scores == thresh_out-TOL);
            if(~isempty(dummy_ind))
                ind_vec = [1:dummy_ind-1,dummy_ind+1:size_united];
                BS_positions = BS_positions(ind_vec);
                BS_regions = BS_regions(ind_vec);
                BS_scores = BS_scores(ind_vec);
                BS_seqs = BS_seqs(ind_vec,:);
            end
        end
    end  % else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Finished Part Doesn't Work %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    sort_by_scores = 1; % Sort positions by scores
    if(sort_by_scores)
        if(iscell(BS_scores))
            [BS_scores{t} sort_perm] = sort(BS_scores{t}, 'descend');
            BS_seqs{t} = BS_seqs{t}(sort_perm,:);
            BS_regions{t} = BS_regions{t}(sort_perm);
            BS_positions{t} = BS_positions{t}(sort_perm);
            if(exist('BS_strand', 'var'))
                BS_strand{t} = BS_strand{t}(sort_perm);
            end
            if(exist('BS_loglike', 'var'))
                if(~isempty(BS_loglike))
                    BS_loglike{t} = BS_loglike{t}(sort_perm);
                end
            end
        else
            [BS_scores(:,t) sort_perm] = sort(BS_scores(:,t), 'descend');
            BS_seqs{t} = BS_seqs{t}(sort_perm,:);
            BS_regions(:,t) = BS_regions(sort_perm,t);
            BS_positions(:,t) = BS_positions(sort_perm,t);
            if(exist('BS_strand', 'var'))
                BS_strand(:,t) = BS_strand(sort_perm,t);
            end
            if(exist('BS_loglike', 'var'))
                if(~isempty(BS_loglike))
                    BS_loglike(:,t) = BS_loglike(sort_perm,t);
                end
            end
        end
    end
    BS_alignments = [];
    if(alignment_flag)
        seq_flag=0; loglike_flag=0; positions_flag=0; tree_flag=0; branch_length_flag=0;
        BS_alignments  = ExtractAligmentByPositions(chr_vec, pos_start_vec, pos_end_vec, mammals_alignment_dir, ...
            seq_flag, loglike_flag, positions_flag, tree_flag, branch_length_flag)% Extract alignments
    end
end      % loop on TFs

if(TFs == 1) % get rid of cell array
    BS_seqs = BS_seqs{1};
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A 'helper' function that is used to pull out the sequences out of a set of packed sequence
% Input:
% seqs - sequences
% seqs_lens - lengths of sequences
% BS_positions - indices for the binding sites positions
% BS_regions - indices of regions with binding sites
% L - binding site length
% diff_lengths_flag
% BS_strand - strand of binding sites (currently not used)
%
% Output:
% BS_seqs - the sequences comprising the binding sites
function BS_seqs = positions_to_seqs(seqs, seqs_lens, BS_positions, BS_regions,...
    L, diff_lengths_flag, BS_strand, N, varargin)

% Find at which WORD each sequence starts
cum_seqs_lens_in_words = zeros(1, length(seqs_lens)+1);
ind=1;
for i=1:length(seqs_lens)
    cum_seqs_lens_in_words(i) = ind;
    ind = ind + floor((seqs_lens(i)-1)/16+1);
end
cum_seqs_lens_in_words(length(seqs_lens)+1) = ind;
BS_words = ceil(L/16);
BS_seqs = zeros(length(BS_positions), BS_words);
base_vec = 4 .^ [0:16]; % powers of four
% Now do it in numerical arrays
prev_region_ind = -9.9; % index for the previous region
for k=1:length(BS_positions)
    kkk = BS_positions(k);
    if(BS_regions(k) ~= prev_region_ind) % extract and expand the relevant sequence if needed
        if(diff_lengths_flag == 0) % here seqs is a two-dim array
            cur_seq = unpack_seqs(seqs(BS_regions(k),:),N);
        else  % here seqs can a cell-array
            if(iscell(seqs))
                cur_seq = unpack_seqs(seqs{BS_regions(k)}, seqs_lens(BS_regions(k)));
            else
                cur_seq = unpack_seqs(seqs(cum_seqs_lens_in_words(BS_regions(k)):cum_seqs_lens_in_words(BS_regions(k)+1)-1),...
                    seqs_lens(BS_regions(k)));
            end
        end
    end
    prev_region_ind = BS_regions(k);
    cur_bs_seq = cur_seq(kkk:min(kkk+L-1, length(cur_seq)));
    if(length(cur_bs_seq) < L)
        cur_bs_seq(length(cur_bs_seq)+1:L) = 1; % add dummy 'A's
    end
    if(L <= 16)
        BS_seqs(k) = sum((cur_bs_seq-1) .* base_vec(1:L));   % Do packing to less than 32 'online'
    else
        % Deal with all words except last
        for j=0:BS_words-2
            BS_seqs(k,j+1) = sum((cur_bs_seq(j*16+1:(j+1)*16)-1) .* base_vec(1:16));   % Do packing to less than 32 'online'
        end
        Remainder =  mod_max(L, 16);   % Deal with last word
        BS_seqs(k,BS_words) = sum((cur_bs_seq((BS_words-1)*16+1:end)-1) .* base_vec(1:Remainder));   % Do packing to less than 32 'online'
    end
end
% if(exist('BS_strand', 'var')) % reverse the ones with opposite strands ? Currently we've decided do not reverse them!
%     for k=1:length(BS_positions)
%         if(BS_strand(k) == 0)
%             BS_seqs(k,:) = seqrcomplement(BS_seqs(k,:));
%         end
%     end
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

