% Compute all the probabilities for all Y's for complex epsilon parameters
% for a whole range of parameters eps_vec, p_vec.
% What do we maximize? according to the flag :
% If flag = find_Y   (0) : find Y maximizing P(Y).
% If flag = finx_X   (1)  : find X maximizing P(X,Y) for a given Y
% If flag = find_XY (2):  find (X,Y) maximizing P(X,Y)
function [min_real_seqs min_real_vals min_imag_seqs min_imag_vals  ...
    max_real_seqs max_real_vals max_imag_seqs max_imag_vals ] =  PlotComplexProbY(eps_real_vec, eps_imag_vec, p, N, seqs_flag, Y)

min_real_seqs = zeros(length(eps_real_vec), length(eps_imag_vec)); 
min_real_vals = zeros(length(eps_real_vec), length(eps_imag_vec)); 
min_imag_seqs = zeros(length(eps_real_vec), length(eps_imag_vec)); 
min_imag_vals = zeros(length(eps_real_vec), length(eps_imag_vec)); 

max_real_seqs = zeros(length(eps_real_vec), length(eps_imag_vec)); 
max_real_vals = zeros(length(eps_real_vec), length(eps_imag_vec)); 
max_imag_seqs = zeros(length(eps_real_vec), length(eps_imag_vec)); 
max_imag_vals = zeros(length(eps_real_vec), length(eps_imag_vec)); 


% we loop over p and send to the P_Y function
switch seqs_flag 

    case 0 % maximize over all Y's

        p_vec = p + zeros(size(eps_real_vec, 1), size(eps_real_vec, 2));

        for eps_ind=1:length(eps_imag_vec)
            cur_eps_vec = eps_real_vec + i.* eps_imag_vec(eps_ind); % Make a complex vector
            P_Y = HMP_ProbY(cur_eps_vec, p_vec, N);

            [max_real_vals(eps_ind,:) max_real_seqs(eps_ind,:)] = max(real(P_Y(:,1:2^(N-1))),[],2);
            [max_imag_vals(eps_ind,:) max_imag_seqs(eps_ind,:)] = max(imag(P_Y(:,1:2^(N-1))),[],2);
            [min_real_vals(eps_ind,:) min_real_seqs(eps_ind,:)] = min(real(P_Y(:,1:2^(N-1))),[],2);
            [min_imag_vals(eps_ind,:) min_imag_seqs(eps_ind,:)] = min(imag(P_Y(:,1:2^(N-1))),[],2);
        end
            

    case 1 % maximize over all X's for a given Y


    case 2 % maximize over all pairs (X,Y)

end


