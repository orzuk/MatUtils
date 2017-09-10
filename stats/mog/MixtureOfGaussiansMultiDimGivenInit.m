% Fit a Mixture of Gaussian model to a given multi-dimensional data and initial guess.
% The method is using EM algorithm.
%
% Input:
% x - input data matrix 
% num_of_gaussians - # gaussians in the model
% num_of_iterations - # iterations to perform in the EM algorithm
% num_starting_points - # different (random) starting points for the EM algorithm
% INIT_P - initial guess for the prob. of each Gaussian
% INIT_M - initial guess for the means of each Gaussian
% INIT_S - initial guess for the st.d. of each Gaussian
%
% Output:
% P - the probs. of each gaussian
% M - the means vector
% S - the covariance matrix
% LogLike - the loglikelihood score of the fit (higher is better)
%
% Written by Liat Ein-Dor and Or Zuk 10/2006
%
function [P,M,S, LogLike] = ...
    MixtureOfGaussiansMultiDimGivenInit(x,num_of_Gaussians,num_of_iterations, INIT_P, INIT_M, INIT_S)

dim=size(x,1); % Data dimension
N=size(x,2); % # data points
epsilon = 0.00000000000000000000000001; % used to avoid underflows

if(num_of_Gaussians == 1) % Deal first with the trivial case of one gaussian
    M = mean(x'); S = cov(x'); P = 1; S_inv = inv(S);
    LogLike = -0.5* (N*log((2*pi)^dim*det(S)) + sum(sum( ((x'-repmat(M,N,1)) * S_inv) .* (x'-repmat(M,N,1)))));
else

    p=ones(N,num_of_Gaussians)./num_of_Gaussians; % start with uniform probs.

    %two_pi_sqrt = sqrt(2*pi); % A constant to save time
    %OldKL = -9999999999999999;
    
    TOL=0.000000001; % tolerance for score improvement
    maxKL=-1000000000; % start with very very low log-likelihood

    sigma=INIT_S; mu = INIT_M; prior = INIT_P;
    for itt=1:num_of_iterations % Start EM iterations
        % E-step: update counts
        for m=1:num_of_Gaussians
            if(sum(sum(sigma{m}.^2)))
                p(:,m)=exp(- sum(((x'-repmat(mu(m,:),N,1)) * inv(sigma{m})  .* (x'-repmat(mu(m,:),N,1)   ))  ./ 2,2 ) ) ...
                    ./ (sqrt((2*pi)^dim.*det(sigma{m})));
            else
                p(:,m)=1;
            end
        end
        p = p + epsilon; % avoid zero probabilities by adding a small thing
        for m=1:num_of_Gaussians % Normalization
            z(:,m)=(p(:,m).*prior(m))./((p*prior'));
        end

        % M-step: update parameters
        sum_z=sum(z); sum_z(sum_z==0)=1;
        for m=1:num_of_Gaussians
            sigma{m}= (repmat(z(:,m),1,dim)' .*(x' - repmat(mu(m,:),N,1))') * (x' - repmat(mu(m,:),N,1)) ./ sum_z(m);      %x is a column vector
            
            if(det(sigma{m}) < TOL)                
                sigma{m} = sigma{m} .* (TOL/det(sigma{m})); % Regularization: Do not allow too small sigmas
            end
            mu(m,:)=sum(repmat(z(:,m),1,dim).*x')./sum_z(m);
        end
        prior=sum_z/N;

        % Calculate the score, and check if improvement is not negligible
        for m=1:num_of_Gaussians
            y(m,:)=prior(m)*exp(-sum ((x'-repmat(mu(m,:),N,1)) * inv(sigma{m}) .* (x'-repmat(mu(m,:),N,1))  ./2, 2)) ...
                ./ (sqrt((2*pi)^dim*det(sigma{m})));
        end

        KL=sum(log(sum(y+epsilon,1)))
        % Stop if improvement in log-likelihood to small
        %         if(KL - OldKL < TOL)
        %             iter = itt
        %             break;
        %         end
        OldKL = KL;
    end

    % Check for improvement in score (with respect to previous starting
    % points) and update parameters
    if(KL>maxKL)
        LogLike = KL % print current score to the screen
        maxKL=KL;
        P=prior; S=sigma; M=mu;
    end

    % Option for future: Sort the means in increasing order 
    % %     [M IndM] = sort(M);
    % %     P = P(IndM); S = S(IndM);
    % return also the best log-likelihood
    LogLike = maxKL;
end