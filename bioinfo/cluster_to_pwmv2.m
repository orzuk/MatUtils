% Generate a pwm from a set of sequences (cluster)
%
% Input:
% c - set of sequences from cluster
% orientation_vec - the best orientation for each kmer (the one matching the alignment)
%
% Output:
% pwm - counts after aligning all sequences
% num_seqs - total number of sequences in cluster
%
function [pwm SeqsAligned CountsAligned num_seqs] = cluster_to_pwmv2(c, orientation_vec, counts, pseudocounts)

num_seqs = length(c);
CountsAligned = counts;
for i=vec2row(find(orientation_vec))
    c{i} = seqrcomplement(c{i});
end
switch num_seqs
    case 1
        SeqsAligned = c{1}; pwm = my_seqlogov2(SeqsAligned,counts,pseudocounts); return;
    case 2
        [score SeqsAligned] = nwalign(c{1}, c{2}); SeqsAligned = SeqsAligned([1 3],:);
    otherwise
        SeqsAligned = multialign(c,'terminalGapAdjust',true,'GAPOPEN',@(sm,sx,len1,len2) 10*sm,'EXTENDGAP',@(sm,sx,len1,len2) sm/3);%prevent gaps
end
pwm = my_seqlogov2(SeqsAligned,counts,pseudocounts);

