% Exponential Ratio used in Fisher-Wright model
% 
% Input: 
% x - allele frequency
% S - effective selective coefficient (=4Ns)
% 
% Output: 
% r - ratio [1-e^(-S(1-x))] / [1-e^(-S)]
% 
function r = h_s(x,S)
r = (1-exp(-S*(1-x))) / (1-exp(-S));

