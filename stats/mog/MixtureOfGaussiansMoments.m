% Compute mean and variance for a Mixture of Gaussians distribution.
%
% Input:
% P - the probs. of each gaussian
% M - the means
% S - the st.ds 
%
% Output:
% mu - mean
% V - variance
% 
% Written by Or Zuk 10/2010
%
function [mu V] = MixtureOfGaussiansMoments(P, M, S) 

mu = sum(P .* M); % mean
V = sum(P .* (M.^2 + S.^2)) - mu^2; % variance 

