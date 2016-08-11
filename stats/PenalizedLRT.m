% Test Penalized Likelihood Ratio Test when Covariance matrix is estimated using Ridge regression

p_vec = 2:2:12;
lambda_vec = 10.5345; % logspace(-1, 3, 10);
lambda = 50; % regularization parameter !!!
alpha = 0.05; % set significance level

equal_parts_flag = 1; % here p equals q


p = 6; % dimension (divide into two blocks)
n = 1; % samples
iters = 400;
epsilon = 0.001; % tolerance
options.regularized = 1;
kron_flag = 1; % 1 - run a Kronecker model !!!

power_vec = zeros(length(lambda_vec), length(p_vec));
power_vec_known_Q = zeros(length(lambda_vec), length(p_vec));
power_vec_known_Q_Ridge = zeros(length(lambda_vec), length(p_vec));


for i_p = 1:length(p_vec)
    p = p_vec(i_p)
    
    if(equal_parts_flag) % two parts have same size (e.g. two proteins)
        p1 = p/2;
        p2 = p/2;
    else % one has only one coordinate (e.g. protein vs. phenotype)
        p1 = p-1;
        p2 = 1;
    end
    
    P = eye(p) + 2.*ones(p); % set second part of sigma
    P = posdefrnd(p);
    %     P =  4*(p-1)*eye(p) + ones(p);
    %     P((p1+1):end, 1:p1) = 2;
    %     P(1:p1, (p1+1):end) = 2;
    
    P0 = P;
    P0((p1+1):end, 1:p1) = 0;
    P0(1:p1, (p1+1):end) = 0;
    
    if(kron_flag)
        q = p+2; % set q
        Q = eye(q) + 2.*ones(q); % set second part of sigma
        Q = posdefrnd(q);
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
    
    df = p1*p2; % (p/2)^2;
    
    
    
    All_Z_mat = zeros(q, iters);
    
    for i_lambda =1:length(lambda_vec) % loop on lamda !!!
        lambda = lambda_vec(i_lambda);
        options.lambda = lambda; options.p1 = p1;
        
        LRT = zeros(iters, 2);
        PLRT_Ridge = zeros(iters, 2);
        LRT_known_Q = zeros(iters, 2);
        PLRT_known_Q_Ridge = zeros(iters, 2);

        for i=1:iters
            if(mod(i, 10) == 0)
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
                    S0((p1+1):end, 1:p1) = 0;
                    S0(1:p1, (p1+1):end) = 0;
                    LRT(i,H+1) = (n/2) * ( log(det(S)) - log(det(S0)) ); %  -(1/2) * ( trace(S*
                    
                    % Now add penalized version!!!
                    
                    S_Ridge0 = zeros(p);
                    for j=1:3 % estimate two parts using Ridge
                        if(j < 3)
                            [U, D, V] = svd(Z(:, (1:(p1*(2-j)+p2*(j-1))) + p1*(j-1)));
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
                            if(j == 1)
                                S_Ridge0(1:p1, 1:p1) = V*diag(theta)*V';
                            else
                                S_Ridge0((p1+1):end, (p1+1):end) = V*diag(theta)*V';
                            end
                        else
                            S_Ridge = V*diag(theta)*V';
                        end
                    end
                else % Estimate under Kronecker product. Only works for n=1. Only for regularization !!!
                    %                Z_mat = vec2mat(Z, p);
                    Z_mat = vec2mat(Z, q)'; % Try transpose !!
                    All_Z_mat(:,i) = Z_mat(:,2); % save one column for testing later
                    
                    P_Ridge0 = zeros(p);
                    P_known_Q0 = zeros(p);
                    for j=3:-1:1 % estimate two parts using Ridge
                        if(j < 3)
                            [U, D, V] = svd(Z_mat(:, (1:(p1*(2-j)+p2*(j-1))) + p1*(j-1)));
                            r = rank(Z_mat(:, (1:(p1*(2-j)+p2*(j-1))) + p1*(j-1)));
                            %                         if(j == 1)
                            %                             pp = p1;
                            %                         else
                            %                             pp = p2;
                            %                         end
                            options.h = 'h0'; % fit seperately
                            
                            if(j==1)
                                P_known_Q0(1:p1, 1:p1) = Z_mat(:, 1:p1)'*inv(Q)*Z_mat(:, 1:p1) ./ q;
                            else
                                P_known_Q0((p1+1):end, (p1+1):end) = Z_mat(:, (p1+1):end)'*inv(Q)*Z_mat(:, (p1+1):end) ./ q;
                                [P_Ridge0, Q_Ridge0] = ...
                                    KroneckerFlipFlop(Z_mat, P_Ridge, Q_Ridge, epsilon, options);
                                P_known_Q_Ridge0 = KroneckerFlop(Z_mat, Q, options); % just fit
                            end
                            
                        else % j==3. Fit under H1
                            P_known_Q = Z_mat'*inv(Q)*Z_mat ./ q; % Fit for known Q ('tree')

                            [U, D, V] = svd(Z_mat);
                            r = rank(Z_mat); pp = p;
                            eig_vals = diag(D); % closed-form solution from Tibshirani's paper
                            
                            % Compute beta, theta with penalties
                            c1 = -4*lambda*pp^2;
                            c2 = 32*lambda^2*pp + eig_vals.^4 .* (q-pp);
                            c3 = 4*lambda* (eig_vals.^4 - 16*lambda^2);
                            
                            beta = repmat(2*sqrt(lambda/pp), q, 1);
                            beta(1:r) = real(sqrt( (-c2(1:r) - sqrt(c2(1:r).^2 - 4.*c1.*c3(1:r)))  ./ (2.*c1) ));
                            theta = repmat(2*sqrt(lambda/q), pp, 1);
                            theta(1:r) = eig_vals(1:r).^2 .* beta(1:r) ./ (pp .* beta(1:r).^2 - 4*lambda);
                            theta(pp .* beta(1:r).^2 - 4*lambda == 0) = 2*sqrt(lambda/q); % fix eigenvalues near zero
                            theta = max(theta(1:r), 2*sqrt(lambda/q)); 
                            P_Ridge = V*diag(theta)*V';
                            Q_Ridge = U*diag(beta)*U';
                        
                            %                        grad_P = q*P_Ridge - Z_mat' * inv(Q_Ridge) * Z_mat - 4*lambda * inv(P_Ridge) % New! COmpute Gradient !!
                            %                        grad_Q = p*Q_Ridge - Z_mat * inv(P_Ridge) * Z_mat' - 4*lambda * inv(Q_Ridge)
                            options.h = 'h1'; % fit jointly
                            P_known_Q_Ridge = KroneckerFlop(Z_mat, Q, options); % just fit 
                            
                        end
                        
                        % Compute large matrix
                        if(j == 2)
                            S_Ridge0 = kron(P_Ridge0, Q_Ridge0);
                            S_known_Q0 = kron(P_known_Q0, Q);
                            S_known_Q_Ridge0 = kron(P_known_Q_Ridge0, Q);
                        end
                        if(j == 3)
                            S_Ridge = kron(P_Ridge, Q_Ridge);
                            S_known_Q = kron(P_known_Q, Q);
                            S_known_Q_Ridge = kron(P_known_Q_Ridge, Q);
                        end
                    end % loop on which model to estimate
                    
                end
                %             PLRT_Ridge(i,H+1) = (n/2) * ( logdet(P_Ridge) - logdet(P_Ridge0) ) + ... % should be 0 or not??
                %                 (n/2) * ( logdet(Q_Ridge) - logdet(Q_Ridge0) ) + ...
                %                 (n/2) * ( trace(S/S_Ridge) - trace(S/S_Ridge0) ) + ... % new ! Add differences in penalties
                %                 lambda * ( norm(inv(P_Ridge), 'fro') - norm(inv(P_Ridge0), 'fro') + norm(inv(Q_Ridge), 'fro') - norm(inv(Q_Ridge0), 'fro') );
                %            LRT_known_Q(i,H+1) = (n/2) * ( sum(log(eig(P_known_Q))) - sum(log(eig(P_known_Q0))) ) + (n/2) * ( trace(S/S_known_Q) - trace(S/S_known_Q0) );  % Should always be positive ??
                
                PLRT_Ridge(i,H+1) = KroneckerLogLike(Z, P_Ridge, Q_Ridge, lambda) - ...
                    KroneckerLogLike(Z, P_Ridge0, Q_Ridge0, lambda);
                if(~isreal(PLRT_Ridge(i,H+1)))
                    xxxx = 124
                end
                if(PLRT_Ridge(i,H+1)<0)
                    xxxx = 1241234214
                end
                
                    
                LRT_known_Q(i,H+1) = KroneckerLogLike(Z, P_known_Q, Q, 0) - ...
                    KroneckerLogLike(Z, P_known_Q0, Q, 0);
                PLRT_known_Q_Ridge(i,H+1) = KroneckerLogLike(Z, P_known_Q_Ridge, Q, lambda) - ...
                    KroneckerLogLike(Z, P_known_Q_Ridge0, Q, lambda);
                
                %             x_vec = -0.01:0.0001:0.01;
                %             for jj = 1:length(x_vec)
                %                  KRON(jj) = KroneckerLogLike(Z, P_Ridge+x_vec(jj)*eye(p), Q_Ridge, lambda);
                %             end
                %             figure; plot(x_vec, KRON);
            end % loop on hypothesis
        end
        
        % Compute threshold and power
        T_alpha = my_quantile(2*PLRT_Ridge(:,1), 1-alpha);
        power_vec(i_lambda, i_p) = mean(2*PLRT_Ridge(:,2) > T_alpha);
        T_alpha = my_quantile(2*LRT_known_Q(:,1), 1-alpha);
        power_vec_known_Q(i_lambda, i_p) = mean(2*LRT_known_Q(:,2) > T_alpha);
        T_alpha = my_quantile(2*PLRT_known_Q_Ridge(:,1), 1-alpha);
        power_vec_known_Q_Ridge(i_lambda, i_p) = mean(2*PLRT_known_Q_Ridge(:,2) > T_alpha);
        
        num_bins = min(1000, iters/20);
        for cum_flag = 2:0 % 1:1 % 0
            figure; x_vec = 0:0.001:(5*df);
            
            if(cum_flag)
                semilogx(x_vec, chi2cdf(x_vec, df), 'r', 'linewidth', 2);  hold on; % plot on log-scale
                semilogx(sort(2*PLRT_Ridge(:,1)), (1:iters) ./ iters, 'g');
                semilogx(sort(2*PLRT_Ridge(:,2)), (1:iters) ./ iters, 'c');
                semilogx(sort(2*LRT_known_Q(:,1)), (1:iters) ./ iters, 'g--');
                semilogx(sort(2*LRT_known_Q(:,2)), (1:iters) ./ iters, 'c--');
                
                ylabel('Cumulative');
            else
                plot(x_vec, chi2pdf(x_vec, df), 'r', 'linewidth', 2); hold on;
                hist_density(2*PLRT_Ridge(:,1), num_bins, 'g');
                hist_density(2*PLRT_Ridge(:,2), num_bins, 'c');
                hist_density(2*LRT_known_Q(:,1), num_bins, 'k');
                hist_density(2*LRT_known_Q(:,2), num_bins, 'y');
                ylabel('Density');
            end
            
            xlabel('x');
            
            legend_vec = {['$\chi^2(' num2str(df) ')$'], ...
                ['$H_0$-LRT-Ridge($\lambda=' num2str(round(lambda,2)) '$)'], ...
                ['$H_1$-LRT-Ridge($\lambda=' num2str(round(lambda,2)) '$)'], ...
                '$H_0-LRT-known-Q$', '$H_1-LRT-known-Q$'};
            if(~kron_flag)
                legend_vec = [legend_vec {'H0-LRT', 'H1-LRT'}];
                hist_density(2*LRT(:,1), num_bins);
                hist_density(2*LRT(:,2), num_bins, 'm');
            end
            legend(legend_vec, 'interpreter', 'latex'); legend('boxoff');
            xlim([0 max( max(5*df, 0.1+max(2*PLRT_Ridge(:))), 0.1+max(2*LRT_known_Q(:))) ]);
            title(['Test Covariance, n=' num2str(n) ', p=' num2str(p) ', q=' num2str(q)]);
        end
        
        
        
    end % loop on lambda
    
end % loop on p

figure; semilogx(lambda_vec, power_vec); hold on;
semilogx(lambda_vec, power_vec_known_Q);
title(['Test Covariance, Power as function of \lambda n=' num2str(n) ', p=' num2str(p) ', q=' num2str(q)]);
legend({'Unknown', 'Known-Q'}); legend('boxoff');
xlabel('\lambda'); ylabel('Power');


figure; plot(p_vec, power_vec); hold on;
plot(p_vec, power_vec_known_Q, 'r');
plot(p_vec, power_vec_known_Q_Ridge, 'g');  
title(['Test Covariance, Power as function of p n=' num2str(n) ', p=' num2str(p) ', q=' num2str(q)]);
legend({'Unknown-Ridge', 'Known-Q',  'Known-Q-Ridge'}); legend('boxoff');
xlabel('p'); ylabel('Power');

