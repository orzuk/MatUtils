% Compute all the possible maximal sequences
% for a whole range of parameters eps_vec, p_vec.
% What do we maximize? according to the flag :
% If flag = find_Y   (0) : find Y maximizing P(Y).
% If flag = finx_X   (1)  : find X maximizing P(X,Y) for a given Y
% If flag = find_XY (2):  find (X,Y) maximizing P(X,Y)
function [min_seqs, min_vals, max_seqs, max_vals] = FindAllMaximalYSeqs(eps_vec, p_vec, N, seqs_flag, Y, varargin)

min_seqs = zeros(length(p_vec), length(eps_vec)); 
min_vals = zeros(length(p_vec), length(eps_vec)); 
max_seqs = zeros(length(p_vec), length(eps_vec)); 
max_vals = zeros(length(p_vec), length(eps_vec)); 

% we loop over p and send to the P_Y function
switch seqs_flag 

    case 0 % maximize over all Y's

        for p=1:length(p_vec)
            p_is = p
            cur_p_vec = p_vec(p) + zeros(size(eps_vec, 1), size(eps_vec, 2));
            P_Y = HMP_ProbY(eps_vec, cur_p_vec, N);

            [max_vals(p,:) max_seqs(p,:)] = max(P_Y(:,1:2^(N-1)),[],2);     
            [min_vals(p,:) min_seqs(p,:)] = min(P_Y(:,1:2^(N-1)),[],2);     
        end


    case 1 % maximize over all X's for a given Y


    case 2 % maximize over all pairs (X,Y)

end

