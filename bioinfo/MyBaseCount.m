% Similar to Matlab Bioinformatics Toolbox's 'basecount', 
% except that it: 
% 1. Handles a matrix of sequences
% 2. Takes sequences in a packed (numeric) form 
function bc = MyBaseCount(packed_seqs, len)

seqs = unpack_seqs(packed_seqs, len); 
bc = zeros(size(seqs,1), 4); 

for nuc=1:4
    bc(:,nuc) = sum(seqs == nuc,2);
end
