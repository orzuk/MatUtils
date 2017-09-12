% A matlab interface to the c function calc_all_sites_scores.
% In here we enable to look for multiple regions and/or multiple pwms
% Input:
% pwms - the pwms to look for
% packed_seqs - the sequences to look for
% seqs_lens - the lengths of the sequences
% do_log - flag saying if to perform log (natural base) on the pwms
% background_model - enable a background model subtraction (default: NONE!)
% strand - which strand to look at (we enable both)
% max_block_size - how big a block to send the c function each time 
%
% Output:
% sites_scores - an array/cell-array of all the scores
% rev_comp_sites_scores - we also include the scores of the opposite strand (in case both strands are computed)
%
function [sites_scores rev_comp_sites_scores] = ...
    calc_all_sites_scores_matlab(pwms, packed_seqs, seqs_lens, do_log, background_model, ...
    strand, max_block_size, varargin) % assume they're all equal and take the first

rev_comp_sites_scores = []; 
if(~exist('do_log', 'var')) % set default parameters
    do_log = 0; 
end
if(~exist('background_model', 'var'))
    background_model = [];
end
if(~exist('strand', 'var'))
    strand = 1;
end
if(strand == 2) % work on both strands. New! answer is returned in single format!!! 
    sites_scores = calc_all_sites_scores_matlab(pwms, packed_seqs, seqs_lens, do_log, background_model, 1); % positive strand
    rev_comp_sites_scores = calc_all_sites_scores_matlab(pwms, packed_seqs, seqs_lens, do_log, background_model, 0); % negative strand
    return;
end
if(iscell(pwms))  % enumerate all pwms
    sites_scores = cell(1, length(pwms));
    for i=1:length(pwms)
        sites_scores{i} = calc_all_sites_scores_matlab(pwms{i}, packed_seqs, seqs_lens, do_log, background_model, strand);
    end
    return;
end
if(do_log) % log-transform
    pwms = log(pwms);    
end
if(~isempty(background_model)) % correct for background
    pwms = pwms - repmat(log(vec2column(background_model)), 1, size(pwms,2));
end
if(strand == 0) % revesre strand
    pwms = pwmrcomplement(pwms);
end
pwms = single(pwms); is_double = 0; % MAKE EVERYTHING SINGLE TO SAVE TIME&SPACE !!! 
if(~exist('max_block_size', 'var'))
    max_block_size = 160000; % set default block size. Must be a multiplication of 16 (due to packed seqs) 
end
L = size(pwms, 2);
overlap = ceil((L+1)/16)*16; % take possibly an extra block of overlap to be sure
if(iscell(packed_seqs))
    sites_scores = cell(1, length(packed_seqs));
    for i=1:length(packed_seqs)
        if(seqs_lens(i) >= L)
            [block_starts block_ends block_lengths block_num] = divide_region_to_blocks(1, seqs_lens(i), max_block_size);
            block_starts(2:end) = block_starts(2:end) - overlap; % take two-word overlap
            block_lengths = block_ends-block_starts+1;
            for j=1:block_num
                sites_scores{i}(max(block_starts(j),1) :block_ends(j)-L+1) = ...
                    calc_all_sites_scores(pwms, ...
                    packed_seqs{i}((block_starts(j)-1)/16+1: ceil(block_ends(j)/16)), block_lengths(j), -10, is_double); % here we do one by one (single)!
                % Note: There used to be 6 vars here (a zero before is_double) - why????
            end
            %%%%%            sites_scores{i} = single(calc_all_sites_scores(pwms, packed_seqs{i}, seqs_lens(i), -10, is_double)); % New: singles to save space
        else
            sites_scores{i} = []; % leave an empty vector
        end
    end
else % here we got one sequence
    if(size(packed_seqs, 2) == 1)
        packed_seqs = packed_seqs';
    end
    if(seqs_lens >= L)
        [block_starts block_ends block_lengths block_num] = divide_region_to_blocks(1, seqs_lens, max_block_size);
        for j=1:block_num
            sites_scores(block_starts(j):block_ends(j)-L+1) = ...
                calc_all_sites_scores(pwms, packed_seqs((block_starts(j)-1)/16+1:block_ends(j)/16), block_lengths(j), -10, 0, is_double); % here we do one by one !
        end
    else
        sites_scores = [];
    end
end


