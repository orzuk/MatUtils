% Fit a Mixture of Gaussian model to a given data given an initial guess
% The method is using EM algorithm.
% New: 3/2009: Enables to 'freeze' some of the parameters (the one whose value the
% user thinks he knows and doesn't want to change them)
%
% Input:
% x - input data vec
% num_of_gaussians - # gaussians in the model
% num_of_iterations - # iterations to perform in the EM algorithm
% num_starting_points - # different (random) starting points for the EM algorithm
% INIT_P - initial guess for the prob. of each Gaussian
% INIT_M - initial guess for the means of each Gaussian
% INIT_S - initial guess for the st.d. of each Gaussian
% freeze_P - binary vector saying which probs. should be freezed (optional)
% freeze_M - binary vector saying which means should be freezed (optional)
% freeze_S - binary vector saying which st.d. should be freezed (optional)
%
% Output:
% P - the probs. of each gaussian
% M - the means
% S - the st.d.s
% LogLike - the loglikelihood score of the fit (higher is better)
%
% Written by Liat Ein-Dor and Or Zuk 10/2006
%
function [P,M,S, LogLike]=MixtureOfGaussiansGivenInit(x,num_of_Gaussians,num_of_iterations, ...
    INIT_P, INIT_M, INIT_S, ...
    freeze_P, freeze_M, freeze_S, varargin)

epsilon = 0.00000000000000000000000001; % used to avoid underflows
MAX_EXP = 80; % exonent of less then ~100 gives zero for single (~700 for double)
use_C_flag = 0;
if(use_C_flag) % just call a c function to do the work for us
    [S,M,P, LogLike]=MixtureOfGaussiansGivenInit_C(x,num_of_Gaussians,num_of_iterations, ...
        INIT_P, INIT_M, INIT_S); % Need to check that the C function takes it in the right order !!!
else
    N=length(x); % # data points
    if(num_of_Gaussians == 1) % Deal first with the trivial case of one gaussian
        M = mean(x); S = std(x); P = 1;
        LogLike = -N*log(sqrt(2*pi)*S) - sum((x-M).^2 ./ (2*S^2));
    else
        freeze_flag = 0;
        if(exist('freeze_S', 'var'))
            freeze_flag = 1;
            freeze_S_inds = find(freeze_S);
        end
        if(exist('freeze_M', 'var'))
            freeze_M_inds = find(freeze_M);
        end
        if(exist('freeze_P', 'var'))
            freeze_P_inds = find(freeze_P);
            learn_P_inds = setdiff(1:num_of_Gaussians, freeze_P_inds);
            sum_freeze_P = sum(INIT_P(freeze_P_inds));
        end

        p=ones(N,num_of_Gaussians)/num_of_Gaussians; % start with uniform probs.
        two_pi_sqrt = sqrt(2*pi); % A constant to save time
        TOL=0.000000001; % tolerance for score improvement
        maxKL=-1000000000; % start with very very low log-likelihood
        OldKL = -9999999999999999;
        sigma=INIT_S; mu = INIT_M; prior = INIT_P;
        y = zeros(num_of_Gaussians, N);     z=p; % just allocate memory to save some time ..
        for itt=1:num_of_iterations % Start EM iterations
            for m=1:num_of_Gaussians         % E-step: update counts
                if(sigma(m))
                    p(:,m)=exp(- (((x-mu(m))./(sigma(m))).^2) ./ 2 ) ./ (two_pi_sqrt.*sigma(m));
                else
                    p(:,m)=1;
                end
            end
            p = p + epsilon; % avoid zero probabilities by adding a small thing
            for m=1:num_of_Gaussians
                z(:,m)=(p(:,m).*prior(m))./((p*prior'));
            end

            sum_z=sum(z); sum_z(sum_z==0)=1;         % M-step: update parameters
            sigma=sqrt(sum(z.*(repmat(x,1,num_of_Gaussians)-repmat(mu,N,1)).^2)./sum_z);          %x is a column vector
            sigma = max(sigma,TOL); % Regularization: Do not allow too small sigmas
            mu=(z'*x)'./sum_z;
            prior=sum_z/N;
            if(freeze_flag) % keep the original parameters
                sigma(freeze_S_inds) = INIT_S(freeze_S_inds);
                mu(freeze_M_inds) = INIT_M(freeze_M_inds);
                prior(freeze_P_inds) = INIT_P(freeze_P_inds);
                prior(learn_P_inds) = prior(learn_P_inds) .* (1-sum_freeze_P) ./ sum(prior(learn_P_inds));
            end

            for m=1:num_of_Gaussians   % Calculate score, check if improvement is not negligible
                y(m,:)=max(epsilon, ...
                    prior(m)*1/(two_pi_sqrt*sigma(m))*exp(max(-( (x-mu(m))./sigma(m)).^2./2, -MAX_EXP)));
            end
            KL=sum(log(sum(y,1)));
            if(KL - OldKL < TOL)
                mog_iter = itt
                break;
            end
            OldKL = KL;
        end

        % Check for improvement in score (with respect to previous starting
        % points) and update parameters
        if(KL>maxKL)
            LogLike = KL % print current score to the screen
            maxKL=KL;
            P=prior; S=sigma; M=mu;
        end
        
        if(exist('M', 'var'))
            [M IndM] = sort(M);     % Sort the means in increasing order
            P = P(IndM); S = S(IndM);
            LogLike = maxKL;     % return also the best log-likelihood
        else
            sprintf('Problem! never got likelihood improvement. Probably skewed data!');
            error(99);
        end
    end
end % use C flag

