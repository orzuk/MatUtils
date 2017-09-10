% Divide a region into blocks (where the last block is a 'remainder' and typically smaller)
% New: enable multiple blocks and save block inds.
%
% Input: 
% region_start - range start index
% region_end - range end index
%
% Output: 
% block_starts - where each block starts
% block_ends - where each block ends
% block_lengths - lengths of blocks (except last block)
% block_num - number of blocks
% block_ind - index of block (what region did it come from)
%
function [block_starts block_ends block_lengths block_num block_ind] = ...
    divide_region_to_blocks(region_start, region_end, block_size)

ctr=0;
for i=1:length(region_start) % determine total blocks num
    ctr = ctr + ceil((region_end(i) - region_start(i)+1) / block_size);
end

block_starts = zeros(1,ctr); block_ends = zeros(1,ctr); block_ind = zeros(1,ctr); ctr=1;
for i=1:length(region_start) % new: enable multiple regions
    block_num = ceil((region_end(i) - region_start(i)+1) / block_size);
    
    block_starts(ctr:ctr+block_num-1) = region_start(i) + (block_size .* [0:block_num-1]);
    block_ends(ctr:ctr+block_num-1) = block_starts(ctr:ctr+block_num-1) + block_size - 1; 
    block_ends(ctr+block_num-1) = region_end(i);
    
    block_ind(ctr:ctr+block_num-1) = i; 
    ctr=ctr+block_num; 
end
block_lengths = block_ends - block_starts + 1;
block_num = length(block_lengths);
