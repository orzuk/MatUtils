% Test Penalized Likelihood Ratio Test when Covariance matrix is estimated using Ridge regression

p = 6; % dimension (divide into two blocks)
n = 1; % samples
iters = 400;
epsilon = 0.001; % tolerance 
options.regularized = 1; 

kron_flag = 1; % 1 - run a Kronecker model !!!

P =  5*eye(p) + ones(p); 
P((p/2)+1:end, 1:(p/2)) = 2; 
P(1:(p/2), (p/2)+1:end) = 2;

P0 = P;
P0((p/2)+1:end, 1:(p/2)) = 0;
P0(1:(p/2), (p/2)+1:end) = 0;

if(kron_flag)
    q = 15; % set q
    Q = eye(q) + 2.*ones(q); % set second part of sigma
    Sigma = kron(P, Q); % take Kronecker product
    Sigma0 = kron(P0, Q); % Take Kronecker product for null model
    d = p*q; % set total dimension
else
    q=1;
    Sigma = P; % eye(p) + ones(p); % [1 0.5; 0.5 1];
    Sigma0 = P0;
    d = p;
end


%Sigma0 = diag(diag(Sigma));

df = (p/2)^2;

lambda_vec = logspace(-1, 3, 10);
lambda = 50; % regularization parameter !!!
alpha = 0.05; % set significance level
power_vec = zeros(size(lambda_vec));
power_vec_known_Q = zeros(size(lambda_vec));



for i_lambda =1:length(lambda_vec) % loop on lamda !!!
    lambda = lambda_vec(i_lambda); 
    options.lambda = lambda; 

    LRT = zeros(iters, 2);
    PLRT_Ridge = zeros(iters, 2);
    LRT_known_Q = zeros(iters, 2);
    for i=1:iters
        if(mod(i, 200) == 0)
            run_i = i
        end
        
        for H=0:1 % simulate under null and alternative hypothesis
            if(H==0)
                Z = mvnrnd(zeros(d, 1), Sigma0, n); % Sample according to H0
            else
                Z = mvnrnd(zeros(d, 1), Sigma, n); % Sample according to H1
            end
            
            S = (Z'*Z)./ n; % Compute sample covaraince matrix: MLE
            
            
            if(~kron_flag) % estimate Sigma
                S0 = S;     % Sample cvariance for null model - here need to estimate differently
                S0((p/2)+1:end, 1:(p/2)) = 0;
                S0(1:(p/2), (p/2)+1:end) = 0;
                LRT(i,H+1) = (n/2) * ( log(det(S)) - log(det(S0)) ); %  -(1/2) * ( trace(S*
                
                % Now add penalized version!!!
                
                S_Ridge0 = zeros(p);
                for j=1:3 % estimate two parts using Ridge
                    if(j < 3)
                        [U, D, V] = svd(Z(:, (1:(p/2)) + (p/2)*(j-1)));
                    else
                        [U, D, V] = svd(Z);
                    end
                    if(size(D, 1) > 1)
                        eig_vals = diag(D);
                    else
                        eig_vals = D(1,1); % only one eigenvalue
                    end
                    if(n < p / min(4-j,2))
                        eig_vals = [eig_vals' zeros(1,  p / min(4-j,2)-n)]';
                    end
                    theta = (eig_vals.^2 + sqrt( eig_vals.^4 + 16*n*lambda)) ./ (2*n);
                    if(j < 3)
                        S_Ridge0((1:(p/2)) + (p/2)*(j-1), (1:(p/2)) + (p/2)*(j-1)) = V*diag(theta)*V';
                    else
                        S_Ridge = V*diag(theta)*V';
                    end
                end
            else % Estimate under Kronecker product. Only works for n=1. Only for regularization !!!
%                Z_mat = vec2mat(Z, p);
                Z_mat = vec2mat(Z, q)'; % Try transpose !!
                
                P_Ridge0 = zeros(p);
                P_known_Q0 = zeros(p); 
                for j=1:3 % estimate two parts using Ridge
                    if(j < 3)
                        [U, D, V] = svd(Z_mat(:, (1:(p/2)) + (p/2)*(j-1))); r = rank(Z_mat(:, (1:(p/2)) + (p/2)*(j-1))); % min(p, q)
                        pp = p/2;
                    else
                        [U, D, V] = svd(Z_mat); r = rank(Z_mat); % min(p, q)
                        pp = p;
                        
                    end
                    eig_vals = diag(D);
                    % %                 if(n < p / min(4-j,2)) % require similar thing for 'padding' eigenvalues?
                    % %                     eig_vals = [eig_vals' zeros(1,  (p / min(4-j,2))-n)]';
                    % %                 end
                    
                    % Compute beta, theta with penalties
                    c1 = -4*lambda*pp^2;
                    c2 = 32*lambda^2*pp + eig_vals.^4 .* (q-pp);
                    c3 = 4*lambda* (eig_vals.^4 - 16*lambda^2);
                    
                    beta = repmat(2*sqrt(lambda/pp), q, 1);
                    beta(1:r) = sqrt( (-c2 - sqrt(c2.^2 - 4.*c1.*c3))  ./ (2.*c1) );
                    theta = repmat(2*sqrt(lambda/q), pp, 1);
                    theta(1:r) = eig_vals.^2 .* beta(1:r) ./ (pp .* beta(1:r).^2 - 4*lambda);
                    
                    if(j < 3)
                        P_Ridge0((1:(p/2)) + (p/2)*(j-1), (1:(p/2)) + (p/2)*(j-1)) = V*diag(theta)*V';
                        Q_Ridge0 = U*diag(beta)*U'; % Why set twice to same value ???

                        P_known_Q0((1:(p/2)) + (p/2)*(j-1), (1:(p/2)) + (p/2)*(j-1)) = Z_mat(:, (1:(p/2)) + (p/2)*(j-1))'*inv(Q)*Z_mat(:, (1:(p/2)) + (p/2)*(j-1)) ./ q;
                    
                    else
                        P_Ridge = V*diag(theta)*V';
                        Q_Ridge = U*diag(beta)*U';
                        P_known_Q = Z_mat'*inv(Q)*Z_mat ./ q; % Fit for known Q ('tree')
                   end
                    
                    if(j == 2)
                        S_Ridge0 = kron(P_Ridge0, Q_Ridge0);
%                        [P_Ridge_iter, Q_Ridge_iter] = ...
%                            KroneckerFlipFlop(Z_mat, P_Ridge0, Q_Ridge0, epsilon, options)
                        S_known_Q0 = kron(P_known_Q0, Q); 
                    end
                    if(j == 3)
                        S_Ridge = kron(P_Ridge, Q_Ridge);
                        S_known_Q = kron(P_known_Q, Q); 
                    end
                end % loop on which model to estimate
                
            end
%             PLRT_Ridge(i,H+1) = (n/2) * ( logdet(P_Ridge) - logdet(P_Ridge0) ) + ... % should be 0 or not??
%                 (n/2) * ( logdet(Q_Ridge) - logdet(Q_Ridge0) ) + ...
%                 (n/2) * ( trace(S/S_Ridge) - trace(S/S_Ridge0) ) + ... % new ! Add differences in penalties
%                 lambda * ( norm(inv(P_Ridge), 'fro') - norm(inv(P_Ridge0), 'fro') + norm(inv(Q_Ridge), 'fro') - norm(inv(Q_Ridge0), 'fro') ); 
            
            PLRT_Ridge(i,H+1) = KroneckerLogLike(Z, P_Ridge, Q_Ridge, lambda) - KroneckerLogLike(Z, P_Ridge0, Q_Ridge0, lambda);
%            LRT_known_Q(i,H+1) = (n/2) * ( sum(log(eig(P_known_Q))) - sum(log(eig(P_known_Q0))) ) + (n/2) * ( trace(S/S_known_Q) - trace(S/S_known_Q0) );  % Should always be positive ?? 
            LRT_known_Q(i,H+1) = KroneckerLogLike(Z, P_known_Q, Q, 0) - KroneckerLogLike(Z, P_known_Q0, Q, 0);
        end % loop on hypothesis
    end
    
    % Compute threshold and power 
    T_alpha = my_quantile(2*PLRT_Ridge(:,1), 1-alpha); 
    power_vec(i_lambda) = mean(2*PLRT_Ridge(:,2) > T_alpha);
    T_alpha = my_quantile(2*LRT_known_Q(:,1), 1-alpha); 
    power_vec_known_Q(i_lambda) = mean(2*LRT_known_Q(:,2) > T_alpha);
end % loop on lambda

figure; semilogx(lambda_vec, power_vec); hold on; 
semilogx(lambda_vec, power_vec_known_Q);
title(['Test Covariance, Power as function of \lambda n=' num2str(n) ', p=' num2str(p) ', q=' num2str(q)]); 
legend({'Unknown', 'Known-Q'}); legend('boxoff'); 
xlabel('\lambda'); ylabel('Power'); 


num_bins = min(1000, iters/20);

for cum_flag = 0:1
    figure; hold on; x_vec = 0:0.001:(5*df);
    
    if(cum_flag)
        plot(x_vec, chi2cdf(x_vec, df), 'r', 'linewidth', 2);
        plot(sort(2*PLRT_Ridge(:,1)), (1:iters) ./ iters, 'g');
        plot(sort(2*PLRT_Ridge(:,2)), (1:iters) ./ iters, 'c');
        plot(sort(2*LRT_known_Q(:,1)), (1:iters) ./ iters, 'g--');
        plot(sort(2*LRT_known_Q(:,2)), (1:iters) ./ iters, 'c--');

        ylabel('Cumulative');
    else
        plot(x_vec, chi2pdf(x_vec, df), 'r', 'linewidth', 2);
        hist_density(2*PLRT_Ridge(:,1), num_bins, 'g');
        hist_density(2*PLRT_Ridge(:,2), num_bins, 'c');
        hist_density(2*LRT_known_Q(:,1), num_bins, 'k');
        hist_density(2*LRT_known_Q(:,2), num_bins, 'y');
        ylabel('Density');
    end
    
    xlabel('x');
    
    legend_vec = {['$\chi^2(' num2str(df) ')$'], ...
        ['H0-Ridge($\lambda=' num2str(round(lambda,2)) '$)'], 'H1-LRT-Ridge'};
    if(~kron_flag)
        legend_vec = [legend_vec {'H0-LRT', 'H1-LRT'}];
        hist_density(2*LRT(:,1), num_bins);
        hist_density(2*LRT(:,2), num_bins, 'm');
    end
    legend(legend_vec, 'interpreter', 'latex'); legend('boxoff');
    xlim([0 max(5*df, 0.1+max(-2*PLRT_Ridge(:)))]);
    title(['Test Covariance, n=' num2str(n) ', p=' num2str(p) ', q=' num2str(q)]);
end

