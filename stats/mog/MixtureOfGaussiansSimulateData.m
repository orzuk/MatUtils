% Simulate data from a Mixture of Gaussians distribution.
%
% Input:
% P - the probs. of each gaussian
% M - the means
% S - the st.ds 
% n - length of data to generate
%
% Output:
% R - the data vector
% G - (optional) the Gaussian chosen for each data point 
% 
% Written by Or Zuk 10/2006
%
function [R G] = MixtureOfGaussiansSimulateData(P, M, S, n) 

num_gaussians = length(M);
G = weighted_rand(P, n); % First determine which gaussian generates each data point
R = randn(1,n); % Draw from the Gaussian distributions 
for i=1:num_gaussians
    R(G == i) = R(G == i) .* S(i) + M(i);
end
