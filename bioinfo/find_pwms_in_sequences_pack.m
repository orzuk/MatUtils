% Find the matches to position weight matrices in a given set of sequences.
%
% Note: this function is not so useful anymore and barely used. Please use
% the function 'find_best_pwm_matches_in_sequences_pack' which offers
% similar functionality and much more user friendly interface
%
% seqs is an array of the promoters sequences, represented as digits from 1
% to 4. pwms is a cell array of the matrices we use in order to look for
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
% do_log - flag saying if to do log to the pwms
% strand - which strand to look for (???)
%
% diff_lengths_flag - flag saying if each sequence is of different length or not 
% 
% Output: 
% BS_united - a structure containig all the binding sites which passed the thresholds
%
function BS_united = find_pwms_in_sequences_pack(seqs, seqs_lens, pwms, thresholds, ...
    smart_thresholds, do_log, strand, diff_lengths_flag)
% Add a path to libi's functions
%path(path,'E:\Research\sequences\libis_software');
%path(path,'libis_software');

tomy = cputime;

% transfer cell to a long array
if(iscell(seqs))
    seqs2 = [];
    for i=1:length(seqs)
        seqs2 = [seqs2 seqs{i}];
    end
end
seqs = seqs2;


TOL = 0.00000000001;
% Here DO NOT calculate the nucleotide distribution. It is given from
% outside !!!

up_down = 'up';   % 'up' - an upperbound , 'down' - a lower bound on the score
% my_pval_thresh = pval; %%0.1;   % A threshold below it to take ........ need to be taken care of more carefully

% N is the maximal length !!
N = max(seqs_lens); % size(seqs,2);% The length of the sequences. Here we give it from outside since the sequence is packed
%packed_N = size(seqs,2);   % The packed length can be obtained from seqs
num_genes = length(seqs_lens); % The number of genes prom oter sequences
MAX_BS = 10;   % Maximal putative binding sites for one gene

% Do chunks in order to avoid memory problems
MAX_SEQS_TOGETHER = min(5000,num_genes);  % Never take a block larger than the number of genes !!!
num_blocks = ceil(num_genes / MAX_SEQS_TOGETHER); % The number of blocks to work with
TFs = size(pwms, 1); % The number of pwms

% Copy to a local array
local_pwms = pwms;

% Perform log transform on the matrices if neccessary
if(do_log)
    for t = 1:TFs
        local_pwms{t,2} = log(local_pwms{t,2});
    end
end

timtim = cputime;

counter = 1;

% We need to be clever here, since seqs_lens is in nucleotides, and the
% seqs are packed so we need them in words
cum_seqs_lens = zeros(1, length(seqs_lens));
cum_seqs_lens(1) = floor((seqs_lens(1)-1)/16)+1;
for i=2:length(seqs_lens)
    cum_seqs_lens(i) = cum_seqs_lens(i-1) +  floor((seqs_lens(i)-1)/16)+1;
end

% Go over all the possible pwms, one by one
for t=1:TFs %%%%% Do it for now just for one TF !!!!!
    L = size(local_pwms{t,2}, 2);  % The length of the matrix binding site

    % We must here to use the mex functions
    % Pack the sequences ? Assume seqs are already packed !
    local_threshold = thresholds(t);
    local_smart_thresholds = smart_thresholds;    
    loc_scores_lengths = seqs_lens-L+1;
    loc_scores_lengths = max(loc_scores_lengths, 0); % Avoid negative scores length. If sequence is shorter than pwm, put length=0

    for cur_block = 1:num_blocks % divide the work into blocks
        cur_block_size = min(MAX_SEQS_TOGETHER, num_genes-(cur_block-1)*MAX_SEQS_TOGETHER );
        
        % New ! (4.7.05) - We need to deal with different lengths !!!
        if(diff_lengths_flag == 0)
            is_double = isa(local_pwms{t,2}, 'double');
            loc_scores = calc_all_sites_scores(local_pwms{t,2}, ...
                seqs((cur_block-1)*MAX_SEQS_TOGETHER+1:(cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size,:), seqs_lens(1), -10); % assume they're all equal and take the first
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

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Here try to remove the overlaps AFTER we got the indexes ..

        %     loc_scores_size = size(loc_scores)
        %    smart_thresholds = 0; % Temp ! must remove later !!!!

        
        [BS_genes BS_places BS_scores thresh_out] = get_scores_above_threshold_Matlab( loc_scores, ...
            loc_scores_lengths((cur_block-1)*MAX_SEQS_TOGETHER+1:(cur_block-1)*MAX_SEQS_TOGETHER+cur_block_size), ...
            local_threshold, L, local_smart_thresholds, diff_lengths_flag );
        local_threshold = thresh_out;
        local_smart_thresholds = 0; % For the next blocks keep the threshold !!!! 
        len = length(BS_places); % see how many binding sites were found

        % Adjust genes to the current block
        BS_genes = BS_genes + (cur_block-1)*MAX_SEQS_TOGETHER;
        
        % Avoid empty files. Put one that is like the threshold
        if(len == 0)
            len = 1;
            BS_places = 1;
            BS_genes = 1;
            BS_scores = thresh_out-TOL; %thresholds(t);
        end

        % Find at which WORD each sequence starts
        cum_seqs_lens_in_words = zeros(1, length(seqs_lens)+1);
        ind=1;
        for i=1:length(seqs_lens)
            cum_seqs_lens_in_words(i) = ind;
            ind = ind + floor((seqs_lens(i)-1)/16+1);
        end
        cum_seqs_lens_in_words(length(seqs_lens)+1) = ind;
        BS_words = ceil(L/16);
        BS_seqs = zeros(len, BS_words);
        base_vec = 4 .^ [0:16]; % powers of four
        % Now do it in numerical arrays
        for k=1:len
            kkk = BS_places(k);
            if(diff_lengths_flag == 0)
                temp = unpack_seqs(seqs(BS_genes(k),:),N);
            else
                temp = unpack_seqs(seqs(cum_seqs_lens_in_words(BS_genes(k)):cum_seqs_lens_in_words(BS_genes(k)+1)-1),seqs_lens(BS_genes(k)));
            end
            temp = temp(kkk:kkk+L-1);
            if(L <= 16)
                BS_seqs(k) = sum((temp-1) .* base_vec(1:L));   % Do packing to less than 32 'online'
            else
                % Deal with all words except last
                for j=0:BS_words-2
                    BS_seqs(k,j+1) = sum((temp(j*16+1:(j+1)*16)-1) .* base_vec(1:16));   % Do packing to less than 32 'online'
                end
                % Deal with last word
                Remainder =  mod(L, 16);
                if(Remainder == 0) then
                    Remainder = 16;
                end
                BS_seqs(k,BS_words) = sum((temp((BS_words-1)*16+1:end)-1) .* base_vec(1:Remainder));   % Do packing to less than 32 'online'
            end
        end

        % Move location to the middle of the binding sites
        if(strand==1) % coding strand
            BS_places = BS_places + floor(L/2);
        end
        if(strand==0) % opposite strand
            BS_places = seqs_lens(BS_genes) - BS_places(k) - floor(L/2) + 1;
%            for k=1:len
%                BS_places(k) = seqs_lens(BS_genes(k))-floor(L/2)+1-BS_places(k);
%            end
        end
        BS_places = 2^40*strand + BS_places; % Add the strand ...
        cur_BS_united = [(2^20)*BS_genes + BS_places  BS_scores  BS_seqs];     % Note that these are already packed !!! We do here packing 'online' !!!
        % Note !!! NEW PACKING !!! STRAND IS WITH THE FIRST WORD
        
        if(cur_block == 1)
            BS_united = cur_BS_united;
        else
            BS_united = [BS_united' cur_BS_united']';
        end

        % update the counter
        counter = counter + len;
    end % loop on blocks
    threshold_used_was = local_threshold
    
    % Remove all those dummies
    BS_united = unique(BS_united, 'rows');
    
    % Remove the threshold if we have more than one !!
    size_united = size(BS_united, 1);
    if(size_united > 1)
       dummy_ind = find(BS_united(:,2) == thresh_out-TOL);
       if(isempty(dummy_ind) == 0)
          ind_vec = [1:dummy_ind-1,dummy_ind+1:size_united];
          BS_united = BS_united(ind_vec,:);
       end        
    end
end      % loop on TFs


