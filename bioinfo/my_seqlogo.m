% Compute seqlogo without displaying, fill in missing letters and normalize
function pwm = my_seqlogo(seqs)
pwm = seqlogo(seqs, 'displaylogo',0);

L = size(pwm{2},2); P = zeros(4,L);
for i=1:length(pwm{1})
    P(nt2int(pwm{1}(i)),:) = pwm{2}(i,:);
end
pwm = P ./ repmat(sum(P), size(P,1), 1);

 