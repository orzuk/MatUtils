% Pack the sequences such that every 8-bit uint contains 
% 3 basepairs. We work on base 5, where zero is unknown/missing nucleotide
% Note : It is assumed that all sequences here are of the same length and
% is divisible by three! 
function PackedSeqs = PackUint8Seqs(Seqs)

PackedSeqs = mod(Seqs, 5); % make 5 as zero
PackedSeqs = PackedSeqs(:,1:3:end)*25 + ...
    PackedSeqs(:,2:3:end)*5 + ...
    PackedSeqs(:,3:3:end);

