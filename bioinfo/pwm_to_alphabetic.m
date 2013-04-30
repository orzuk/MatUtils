% Transfers the PWM into alphabetical first consensus lexicographic ..
function P_alphabetic = pwm_to_alphabetic(P)

cons_seq = pwm_to_consensus(P);
rev_comp_cons_seq = seqrcomplement(cons_seq);

if(strlexcmp(cons_seq, rev_comp_cons_seq) == 1)
    P_alphabetic = pwmrcomplement(P);
else
    P_alphabetic = P;
end