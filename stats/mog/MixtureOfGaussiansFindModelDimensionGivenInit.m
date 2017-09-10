% Fit a Mixture of Gaussians model to a given data and find model dimension given initial guess.
% It tries different dimensions (# of Gaussians), until it finds the best one.
% The score used is the MDL score, thus giving a likelihood
% term and a model complexity term.
% The method for maximum likelihood (for each dimension) is using EM algorithm. 
%
% Input:
% x - input data vec
% num_of_gaussians - # gaussians in the model
% num_of_iterations - # iterations to perform in the EM algorithm
% num_starting_points - # different (random) starting points for the EM algorithm
%
% Output: 
% P - the probs. of each gaussian
% M - the means
% S - the st.ds 
% LogLike - the loglikelihood score of the fit (higher is better)
%
% Written by Or Zuk 10/2006
%
function [P,M,S, Dim, LogLike]=MixtureOfGaussiansFindModelDimensionGivenInit(x,max_num_of_Gaussians,num_of_iterations, num_starting_points)

N=length(x); % # data points
max_MDL_score = -9999999999999999; max_D = -1;
for D=1:max_num_of_Gaussians
    [P{D},M{D},S{D}, LogLike{D}]=MixtureOfGaussiansEM(x,D,num_of_iterations, num_starting_points);
%     cur_S = S{D}
%     cur_M = M{D}
%     cur_P = P{D} 
%     cur_LogLike = LogLike{D}
    MMM{D} = (LogLike{D} - 3*D*log2(N)/2)
    if(  (LogLike{D} - 3*D*log2(N)/2) > max_MDL_score )
        max_D = D;
        max_MDL_score = (LogLike{D} - 3*D*log2(N)/2); % Here model length is in bits
    end
end

% return the optimal parameters
Dim = max_D;
S = S{Dim}; M = M{Dim}; P = P{Dim}; LogLike = LogLike{Dim};
    
