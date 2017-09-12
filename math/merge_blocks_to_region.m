% Divide a region into blocks (where the last block is a 'remainder' and typically smaller)
%
% Input: 
% block_starts - where each block starts
% block_ends - where each block ends
%
% Output: 
% regions_start - range start index
% regions_end - range end index
% block_lengths - lengths of blocks (except last block)
% block_num - number of blocks
%
function [region_start region_end block_size block_lengths block_num] = ...
    merge_blocks_to_region(block_starts, block_ends)

block_num = length(block_starts); 
region_start = min(block_starts); 
region_end = max(block_ends);
block_size = block_ends(1) - block_starts(1) + 1;
block_lengths = block_ends - block_starts + 1;

