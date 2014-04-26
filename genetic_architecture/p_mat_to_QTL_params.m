% Encoding of one SNP QTL in a matrix 
%
% Input: 
% p_mat - 4xN matrix with QTL parameters
% 
% Output: 
% f_vec - allele freq.
% beta - effect size 
% V - total trait's variance 
%
function [f_vec beta V] = p_mat_to_QTL_params( p_mat )

f_vec = p_mat(:,1); % encoding convention: first value is the frequency
beta = p_mat(:,2);
V = p_mat(:,4); % encoding convension: last value is traits variance


