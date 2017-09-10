% Compute a weighted average of aligned sequences + pseudocounts to form a pwm
%
% Input:
% seqs-aligned sequences
% counts-weight for each sequence
% pseudocount-scalar, weight of pseudocounts to be added to each letter in
% each position of the pwm
% Output:
% pwm-weighted average of all seqs for all letters and postions in the pwm;
% columns add to 1.
%
function pwm = my_seqlogov2(seqs,counts,pseudocounts)
if iscell(seqs)
    N = legnth(seqs);
    L = length(seqs{1});
    intSeq = zeros(N,L);
    for i = 1:N
        intSeq(i,:) = nt2int(seqs{i});
    end
else
    [N,L] = size(seqs);
    intSeq = nt2int(seqs);
end

pwm = zeros(4,L);
for i = 1:L
    for j = 1:4
        ii_j = intSeq(:,i) == j;
        pwm(j,i) = sum(counts(ii_j));
    end
    pwm(:,i) = (pwm(:,i) + pseudocounts)/sum(pwm(:,i) + pseudocounts);
end


 