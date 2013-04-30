% Compute the reverse complement when sequences are packed
function rev_comp_seq = reverse_complement_for_packed(seq, seq_len)
% A = 1; 
% C = 2; 
% G = 3; 
% T = 4; 

% Do in a 'dumb' way 
temp = unpack_seqs(seq, seq_len);
rev_comp_seq = temp;

rev_comp_seq(find(temp == 1)) = 4; % transfer A-T C-G - Do the OLD style !
rev_comp_seq(find(temp == 2)) = 3;
rev_comp_seq(find(temp == 3)) = 2;
rev_comp_seq(find(temp == 4)) = 1;

rev_comp_seq = rev_comp_seq(:,end:-1:1); % reverse   


%%%%rev_comp_seq = 5-rev_comp_seq(end:-1:1);
[rev_comp_seq dumb] = pack_seqs(rev_comp_seq);
