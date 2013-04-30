% Fit a Mixture of Gaussians model to a given data and find the number of Gaussians.
% It tries different dimensions (# of Gaussians), until it finds the best one.
% The score used is the MDL score, thus giving a likelihood
% term and a model complexity term.
% The method for maximum likelihood (for each dimension) is using EM algorithm.
% In EM there is always the question of what starting point to take, how
% many starting points etc. We are experimenting with this. Current
% approach: Go 'bottom-up' (from one gaussian to many), and at each stage
% choose as a starting point the gaussians of the previous stage and an
% additional gaussian which is between the gaussians, where it can be the
% most right, the next to most right etc. so if we had in the previous
% stage k gaussians, then there are currently k+1 places to put the next
% gaussian:
%  | - previous gaussians. x - possibilies for putting the current gaussian
%  x  |  x    |    x        |     x   | x
%
% Input:
% x - input data vec
% num_of_gaussians - # gaussians in the model
% num_of_iterations - # iterations to perform in the EM algorithm
% num_starting_points - # different (random) starting points for the EM algorithm
%
% Output:
% P - the probs. of each Gaussian
% M - the means
% S - the st.ds
% LogLike - the loglikelihood score of the fit (higher is better)
%
% Written by Or Zuk 10/2006
%
function [P,M,S, Dim, LogLike]=MixtureOfGaussiansFindModelDimension(x,max_num_of_Gaussians,num_of_iterations, num_starting_points)

N=length(x); % # data points
max_MDL_score = -9999999999999999; max_D = -1;

ttt = cputime;
for D=1:max_num_of_Gaussians
    if(D == 1)
        [P{D},M{D},S{D}, LogLike{D}]=MixtureOfGaussiansEM(x,D,num_of_iterations, num_starting_points);
    else
        LogLike{D} = -99999999999999;
        for j=1:D % loop over where to put the next gaussian
            INIT_P(setdiff(1:D, j)) = P{D-1} .* (D-1)/D; INIT_P(j) = 1/D;
            INIT_M(setdiff(1:D, j)) = M{D-1};
            switch j
                case 1
                    INIT_M(j) = mean(x(x <= max(M{D-1}(j), min(x))));
                case D
                    INIT_M(j) = mean(x(x >= min(M{D-1}(j-1),max(x))));
                otherwise
                    INIT_M(j) = mean(x(  (x >= M{D-1}(j-1)) & (x <= M{D-1}(j)) ));
            end


            INIT_S(setdiff(1:D, j)) = S{D-1}; INIT_S(j) = mean(S{D-1});
            %              ttt_mat = cputime;
            %              [tmp_P, tmp_M, tmp_S, tmp_LogLike]=MixtureOfGaussiansGivenInitMatlab(x,D,num_of_iterations, INIT_P, INIT_M, INIT_S);
            %              ttt_mat = cputime - ttt_mat
            ttt_c = cputime;
            [tmp_P, tmp_M, tmp_S, tmp_LogLike]=MixtureOfGaussiansGivenInit(double(x),D,num_of_iterations, INIT_P, double(INIT_M), double(INIT_S));
            ttt_c = cputime - ttt_c;

            if(tmp_LogLike >= LogLike{D})
                S{D} = tmp_S; M{D} = tmp_M; P{D} = tmp_P; LogLike{D} = tmp_LogLike;
            end
        end
    end


    %   [P{D},M{D},S{D}, LogLike{D}]=MixtureOfGaussiansEM(x,D,num_of_iterations, num_starting_points); Old version: no starting points
    MMM{D} = (LogLike{D} - 3*D*log2(N)/2);
    if(  (LogLike{D} - 3*D*log2(N)/2) > max_MDL_score )
        max_D = D;
        max_MDL_score = (LogLike{D} - 3*D*log2(N)/2); % Here model length is in bits
    end
end

cputime - ttt

% return the optimal parameters
Dim = max_D;
S = S{Dim}; M = M{Dim}; P = P{Dim}; LogLike = LogLike{Dim};

