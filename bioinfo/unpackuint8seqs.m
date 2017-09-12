% Pack the sequences such that every 8-bit uint contains 
% 3 basepairs. We work on base 5, where zero is unknown/missing nucleotide
% Note : It is assumed that all sequences here are of the same length and
% is divisible by three! 
function [Seqs] = UnPackUint8Seqs(PackedSeqs)

Seqs = zeros(size(PackedSeqs,1), 3*size(PackedSeqs,2), 'uint8');
Seqs(:,3:3:end) = mod(PackedSeqs,5);
Seqs(:,2:3:end) = mod((PackedSeqs - Seqs(:,3:3:end)) ./5, 5);
Seqs(:,1:3:end) = mod( (PackedSeqs - 5*Seqs(:,2:3:end) - Seqs(:,3:3:end)) ./25, 5);


