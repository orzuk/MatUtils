% Encoding of one SNP QTL in a matrix
%
% Input: 
% f_vec - allele frequency 
% beta - effect size 
% V - toyal trait's variance
% 
% Output: 
% p_mat - matrix of width 4 encoding the parameters
% 
function p_mat = QTL_params_to_p_mat(f_vec, beta, V)
p_mat = zeros(size(f_vec,1), 4);
p_mat(:,1) = f_vec; % encoding convention: first value is the frequency
p_mat(:,2) = beta ;
p_mat(:,4) = V; % encoding convension: last value is traits variance


