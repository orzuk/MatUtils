% Fit a Mixture of Gaussians model to a given data.
% The method is using EM algorithm.
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
% S - the st.d.s
% LogLike - the loglikelihood score of the fit (higher is better)
%
% Written by Liat Ein-Dor and Or Zuk 10/2006
%
function [P,M,S, LogLike]=MixtureOfGaussiansEM(x,num_of_Gaussians,num_of_iterations, num_starting_points)

N=length(x); % # data points
if(num_of_Gaussians == 1) % Deal first with the trivial case of one gaussian
    M = mean(x); S = std(x); P = 1;
    LogLike = -N*log(sqrt(2*pi)*S) - sum((x-M).^2 ./ (2*S^2));
else

    min_x=min(x);
    max_x=max(x);
    gap_x=max_x-min_x;
    step_x=gap_x/(num_of_Gaussians-1);
    mu=min_x:step_x:max_x; % start with mus evenly spaced across the support of the data
    prior=ones(1,num_of_Gaussians)/num_of_Gaussians;

    p=ones(N,num_of_Gaussians)/num_of_Gaussians; % start with uniform probs.

    two_pi_sqrt = sqrt(2*pi); % A constant to save time


    TOL=0.000000001; % tolerance for score improvement

    maxLL=-1000000000; % start with very very low log-likelihood
    for replica=1:num_starting_points % Start picking random starting points
        OldLL = -9999999999999999;
        sigma=ones(1,num_of_Gaussians).*rand(1,num_of_Gaussians)*step_x; % Randomize only the st.ds. here
        mu = mu + (rand(1,num_of_Gaussians)-0.5) .* (max(mu)-min(mu));
        for itt=1:num_of_iterations % Start EM iterations
            % E-step: update counts
            for m=1:num_of_Gaussians
                if(sigma(m))
                    p(:,m)=exp(- (((x-mu(m))./(sigma(m))).^2) ./ 2 ) ./ (two_pi_sqrt.*sigma(m));
                else
                    p(:,m)=1;
                end
            end
            for m=1:num_of_Gaussians
                z(:,m)=(p(:,m).*prior(m))./((p*prior'));
            end
            
            % M-step: update parameters
            sum_z=sum(z); sum_z(sum_z==0)=1;
            sigma=sqrt(sum(z.*(repmat(x,1,num_of_Gaussians)-repmat(mu,N,1)).^2)./sum_z);  %x is a column vector
            sigma = max(sigma,TOL); % Regularization: Do not allow too small sigmas
            mu=(z'*x)'./sum_z;
            prior=sum_z/N;

            % Calculate the score, and check if improvement is not negligible
            for m=1:num_of_Gaussians
                y(m,:)=prior(m)*1/(two_pi_sqrt*sigma(m))*exp(-( (x-mu(m))./sigma(m)).^2./2);
            end
            LL=sum(log(sum(y,1)));
            if(LL - OldLL < TOL)
                iter = itt
                break;
            end
            OldLL = LL;
        end

        % Check for improvement in score (with respect to previous starting
        % points) and update parameters
        if(LL>maxLL)
            LogLike = LL % print current score to the screen
            maxLL=LL;
            P=prior; S=sigma; M=mu;
        end
    end
    % Sort the means in increasing order
    [M IndM] = sort(M);
    P = P(IndM); S = S(IndM);
    % return also the best log-likelihood
    LogLike = maxLL;
    
end