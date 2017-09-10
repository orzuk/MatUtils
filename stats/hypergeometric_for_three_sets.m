% Give hypergeometric p-value for 3-way intersection, given pairwise intersections
%
% The formula is as follows : Sum over all n's  such that ...
%
%   sum_{n=ABC}^{max of ABC}   (AB n) * (B-AB BC-n) * (A-AB AC-n) * (N-A-B+AB C-AC-BC+n)
%
%
% Input: 
% N - universe size
% A - first set size
% B - second set size
% C - third set size
% AB - A&B intersection size
% AC - A&C intersection size
% BC - B&C intersection size
% ABC - three sets intersection size
%
% Output: 
% over_pval - prob. to get at random a larger intersection
% under_pval - prob. to get at random a smaller intersection
%
function [over_pval, under_pval] = hypergeometric_for_three_sets(N, A, B, C, AB, AC, BC, ABC)


max_sum = min( min(AB, AC), BC); % Determine up to what we need to sum
min_sum = max(max(max(AC+BC-C, 0), AB+BC-B), AB+AC-A); % Determine the minimal value we need to sum

if(ABC > max_sum)
    sprintf('Sorry, illegal input, ABC is too large ..\n')
    over_pval = -1; under_pval = -1;
    return
end
if(ABC < min_sum)
    sprintf('Sorry, illegal input, ABC is too small ..\n')
    over_pval = -1; under_pval = -1;
    return
end
if(A+B+C-AB-AC-BC+ABC > N)
    sprintf('Sorry, illegal input, N is too small ..\n')
    over_pval = -1; under_pval = -1;
    return
end

sum_len = max_sum-min_sum+1;

% Note : We keep all the vectors in a log form, in order to avoid
% underflows and overflows 

% Determine the vecotrs 
AB_over_n_vec = log_binom(repmat(AB, 1, sum_len), [min_sum:max_sum]);
B_minus_AB_over_BC_minus_n_vec = log_binom( repmat(B-AB, 1, sum_len), BC-[min_sum:max_sum]);
A_minus_AB_over_AC_minus_n_vec = log_binom( repmat(A-AB, 1, sum_len), AC-[min_sum:max_sum]);
N_minus_A_minus_B_plus_AB_over_C_minus_AC_minus_BC_plus_n_vec = log_binom(repmat(N-A-B+AB, 1, sum_len), C-AC-BC+[min_sum:max_sum]);


% Now we can calculate the p-value

% Note that the vectors are in log form. We want to find the maximal
% exponent and divide by it
sum_log_vec = AB_over_n_vec +  B_minus_AB_over_BC_minus_n_vec + A_minus_AB_over_AC_minus_n_vec + ...
    N_minus_A_minus_B_plus_AB_over_C_minus_AC_minus_BC_plus_n_vec;
max_sum_log = max(sum_log_vec);

% Needed to avoid overflow or underflow
sum_log_vec = sum_log_vec - max_sum_log;
exp_log_vec = exp(sum_log_vec);
over_pval = sum(exp_log_vec(ABC-min_sum+1:end)) / sum(exp_log_vec);
under_pval = sum(exp_log_vec(1:ABC-min_sum+1)) / sum(exp_log_vec);

