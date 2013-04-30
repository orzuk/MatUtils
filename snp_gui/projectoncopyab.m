% Written by Or Zuk 5/2007
%
% This function'projects' the probability matrix with 36 states
% on the different joint values for (A,B) copy. There are 13 possible such
% values.
%
% Input:
% GammaProbs - A table of 36*seq_len of states probs
%
% Output :
% ProjectedABGammaProbs - A matrix containing the projected probabilities
% AB_copy - A vector of the joint AB copy number values
function [ProjectedABGammaProbs AB_copy] = ProjectOnCopyAB(GammaProbs)

x_dim = 3; x_dim_eq = 2*x_dim-1;

% Tables for indexing the A and B copy
A_copy_tab = [ 0, 0, 1, 0, 2, 0, 0, 0, 1, ...
    0, 2, 0, 1, 1, 2, 1, 3, 1, ...
    0, 0, 1, 0, 2, 0, 2, 2, 3, ...
    2, 4, 2, 0, 0, 1, 0, 2, 0];
B_copy_tab = [ 0, 0, 0, 1, 0, 2, 0, 0, 0, ...
    1, 0, 2, 0, 0, 0, 1, 0, 2, ...
    1, 1, 1, 2, 1, 3, 0, 0, 0, ...
    1, 0, 2, 2, 2, 2, 3, 2, 4];
AB_copy_tab = A_copy_tab + x_dim_eq*B_copy_tab; % here put A in the lsb (base 5) and B in the  msb (another 5)

for j=1:36
    IndsTabAB{j} = find(AB_copy_tab == j-1);
end

NonEmptyInds = unique(AB_copy_tab)+1;
ProjectedABGammaProbs = zeros(length(NonEmptyInds),  size(GammaProbs,2));
k=1;
for j=1:36
    if(~isempty(IndsTabAB{j}))
        ProjectedABGammaProbs(k,:) = sum(GammaProbs(IndsTabAB{j},:),1);
        k=k+1;
    end
end
AB_copy = unique(AB_copy_tab);


A_copy = mod(AB_copy,5); B_copy = (AB_copy - A_copy)/5; total_copy = A_copy+B_copy; figure; 
imagesc([A_copy' B_copy' total_copy']); title('A copy                   B copy             total copy'); colorbar;
