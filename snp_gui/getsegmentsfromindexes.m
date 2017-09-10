% The function finds the start and end of each segment of ones.
% It stores a binary vectoy in a 'run-length' code.
%
% Input:
% Inds - a vector of zeros and ones
%
% Output: 
% SegsStarts - Indexes of one's segments starts
% SegsEnds - Indexes of one's segments ends
%
% For example, if the input is : 1 1 1 1 0 1 0 1 0 0 1 0 0 
% Then the output will be [1 6 8 11], [4 6 8 11]
function [SegsStarts SegsEnds] = GetSegmentsFromIndexes(Inds)

Inds = [0 Inds 0]; % add artificially zeros on both sides
Segs = Inds(1:end-1) - Inds(2:end); % detect start and end points
SegsStarts = find(Segs == -1); SegsEnds = find(Segs == 1)-1;
