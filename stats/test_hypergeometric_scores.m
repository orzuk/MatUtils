% Test all hypergeometric functions
% % % N = 9000; 
% % % A = 4400; 
% % % B = 3030; 
% % % C = 3000; 
% % % AB = 1500;
% % % AC = 1450;
% % % BC = 1000;
% % % ABC = 521;

% Michal's example
N = 15427; 
A = 1964; 
B = 1557; 
C = 1858; 
AB = 367+317;
AC = 383+317;
BC = 317+165;
ABC = 317;

N = 10487;
A = 822;
B = 1997;
AB = 322;

[AB_over_pval AB_under_pval] = hypergeometric_for_two_sets(N, A, B, AB)
[AC_over_pval AC_under_pval] = hypergeometric_for_two_sets(N, A, C, AC)
[BC_over_pval BC_under_pval] = hypergeometric_for_two_sets(N, B, C, BC)
[ABC_over_pval ABC_under_pval] = hypergeometric_for_three_sets(N, A, B, C, AB, AC, BC, ABC)
[again_ABC_over_pval again_ABC_under_pval] = hypergeometric_for_three_sets(N, B, C, A, BC, AB, AC, ABC);  % permute A B and C to see if results remain unchanged


should_be_approximately_zero = [ABC_over_pval - again_ABC_over_pval, ABC_under_pval - again_ABC_under_pval]
% Give p-value for the hypergeometric score of intersection of three sets, 
% given the intersection of every couple

% The formula is as follows : Sum over all n's  such that ...
%
%   sum_{n=ABC}^{max of ABC}   (AB n) * (B-AB BC-n) * (A-AB AC-n) * (N C-AC-BC+n)
%
%

