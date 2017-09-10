% Get the consensus sequence of a pwm
function cons_seq = pwm_to_consensus(P)

[dummy seq] = max(P); 
cons_seq = int2nt(seq);
