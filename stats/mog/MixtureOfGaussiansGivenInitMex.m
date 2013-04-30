% Fit a Mixture of Gaussian model to a given data - use a compuled c mex file.
% The method is using EM algorithm.
% The function just calls a mex file to do the job. Both single and double
% formats are supported, but note: we make ALL the inputs of the same type
% (either all 'single' or all 'double') and the type is set according to
% the data type (x)
%
% Input:
% x - input data vec
% num_of_gaussians - # gaussians in the model
% num_of_iterations - # iterations to perform in the EM algorithm
% num_starting_points - # different (random) starting points for the EM algorithm
%
% Output:
% S - the st.ds
% M - the means
% P - the probs. of each gaussian
% LogLike - the loglikelihood score of the fit (higher is better)
%
function [P,M,S, LogLike]=MixtureOfGaussiansGivenInitMex(x,num_of_Gaussians,num_of_iterations, INIT_P, INIT_M, INIT_S)

if(isa(x, 'single')) % here assume single
    [P,M,S, LogLike]=MixtureOfGaussiansGivenInitSingle(x,single(num_of_Gaussians),single(num_of_iterations), ...
        single(INIT_P), single(INIT_M), single(INIT_S));
else  % here assume double
    [P,M,S, LogLike]=MixtureOfGaussiansGivenInit(x,double(num_of_Gaussians),double(num_of_iterations), ...
        double(INIT_P), double(INIT_M), double(INIT_S));
end
