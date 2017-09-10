% Convert a consensus sequence to a pwm (no derichlet correction)
% 
% Input: 
% cons_seq - nucleotide sequence e.g. 'AGGAGTCCACA'
% 
% Output: 
% pwm - a 4xL pwm representing the consensus sequence
% 
function pwm = consensus_to_pwm(cons_seq)

if(~isnumeric(cons_seq))
    cons_seq = nt2int(cons_seq);
end
L = size(cons_seq,2);  % seqs length
n = size(cons_seq,1); % number of seqs

if(n == 1) % just one pwm
    pwm = zeros(4,L);
    for i=1:4
        pwm(i,cons_seq==i) = 1;
    end
else
    pwm = zeros(4,L,n);
    for j=1:n % costly loop - should be replaced in the future
        for i=1:4
            pwm(i,cons_seq(j,:) == i,j) = 1;
        end
    end
end

