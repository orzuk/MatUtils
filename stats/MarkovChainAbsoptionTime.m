% Compute time (# of generations) to spend in each state until we reach any state of a set A, 
% QU: starting from ?? 
% Input: 
% M - transition matrix of Markov chain
% A - set to be absorped at 
% 
% Output: 
% T_A - vector of times, with T_A(i) = E [ min_T s.t. X_T in A | X_0 = i]
% T_A_mat - matrix of times, with T_A_mat(i) = E [ \sum_{t=1}^T_A 1_{X_t = j} | X_0 = i]
% 
% Note: we also need: T(i,j) - matrix of times, with T(i,j) = E [ time spent at j before in A | X_0 = i]
function [T_A, T_A_mat] = MarkovChainAbsoptionTime(M, A)

% scalar equation: T(i) = 1 + sum_j M(i,j) T(j) not in A, T(i) = 0 in A 
% vector equation: T = ONES + M' * T -> T * (I - M) = ONES -> T = ONES \ (I-M)
n = size(M, 1); M2 = M; M2(:,A)=0; T_A =  (eye(n)-M2) \ ones(n,1);

% Matrix version: Now we want T(i,j) for j not in A 
% scalar equation: T(i,j) = \sum_k M(i,k) T(k,j) 
% Solve for each j separately 
% Consistency condition: We should have: \sum_{j \notin A} T(i,j) = T(i) 
T_A_mat = (eye(n)-M2) \ M; 

