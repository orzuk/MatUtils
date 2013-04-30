% Try Ido's idea for transfer matrix with many HMM chains 
% p is the transition probability matrix S*S, where we have M different 
% markov chains and S = [M/2]+1
% Currently we use only odd chains for simplicity
%
% eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary
p = 0.2; min_eps = 0.0000004; max_eps = 0.000001; MAX_M = 2; 
eps_res=2;

phi_mat = zeros(MAX_M,eps_res);
for eps_iter=1:eps_res
    eps=min_eps + (eps_iter-1) * (max_eps-min_eps)/(eps_res-1);
    [eps_iter eps]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % First do the simplest 2*2 matrix, for M=2 :
    A = [   (p^2 + (1-p)^2)*(eps^2 + (1-eps)^2)    2*p*(1-p)*(eps^2 + (1-eps)^2); 2*p*(1-p)*2*eps*(1-eps)  (p^2 + (1-p)^2)*2*eps*(1-eps)]
    
    
    lambda = eig(A);
    log2(max(lambda));
    entropy([eps, 1-eps]);
    
    N = 1300; 
    phi = ( log2 ([1,1] * (A^N) * [1,1]') )/N;
    phiphi  = log2(max(lambda));
    
    Ido_Ent = entropy([eps, 1-eps]') - log2(max(lambda));
    lower_bound = H(p+eps-2*p*eps);
    upper_bound = H(p+ 2*eps*(1-eps) -1*p*2*eps*(1-eps));
    
    
    % Now compute the greatest eignvalue explicitly
    lambda = eig(A);
    (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)
    lambda_1 = 0.5* ( (p^2 + (1-p)^2) + sqrt ( (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)) )
    lambda_2 = 0.5* ( (p^2 + (1-p)^2) - sqrt ( (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)) )
    

    
    
    lamalm = 0.5 * ( A(1,1) + A(2,2) + sqrt( (A(1,1)-A(2,2))^2 + 4*A(1,2)*A(2,1) ) );
    A(1,1) + A(2,2) - (p^2 + (1-p)^2);
    
    log2(lambda_1);
    % End of simplest 2*2 matrix, for M=2 :
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    % Now do the rest
    
    % Prepare a binomial table 
    binom_tab = zeros(MAX_M+1, MAX_M+1);
    phi_vec = zeros(1,MAX_M);
    
    for i=1:MAX_M+1
        binom_tab(i,1:i) = binom(i-1, [0:i-1]);
    end
    binom_tab;
    eps_power_vec = eps .^ [0:MAX_M];
    one_minus_eps_power_vec = (1-eps) .^ [0:MAX_M];
    p_power_vec = p.^ [0:MAX_M];
    one_minus_p_power_vec = (1-p) .^ [0:MAX_M];
    phi_vec(2) = log2(lambda_1);
    
    
    % This is for independent
    ind_phi_vec = zeros(1,MAX_M);
    ind_phi_double_eps_vec = zeros(1,MAX_M);
    
    for M=2:MAX_M
        S = floor(M/2)+1;    % Set the matrix dimension
        
        A = zeros(S,S);
        
        % Compute the matrix elements one-by-one
        % Note : State i means that the minority has exactly i-1 elements !!
        for i=0:S-1
            for j=0:S-1
                Y_given_X = (eps*(1-eps))^j*(eps^(M-2*j)+(1-eps)^(M-2*j));
                eps^2 + (1-eps)^2;
                % First i->j
                k_vec = [max(0,i-j):min(i,M-j)];
                X_given_prev_X = sum(binom_tab(i+1,k_vec+1) .* binom_tab(M-i+1,j-i+k_vec+1) .* p_power_vec(j-i+2*k_vec+1) .* one_minus_p_power_vec(M+i-j-2*k_vec+1));
                (1-p)^2;
                % Now i->M-j
                k_vec = [max(0,i-M+j):min(i,j)];
                X_given_prev_X = X_given_prev_X + sum(binom_tab(i+1,k_vec+1) .* binom_tab(M-i+1,M-j-i+k_vec+1) .* p_power_vec(M-j-i+2*k_vec+1) .* one_minus_p_power_vec(i+j-2*k_vec+1));
                p^2;
                A(i+1,j+1) = Y_given_X * X_given_prev_X;
                (p^2 + (1-p)^2)*(eps^2 + (1-eps)^2);
                
                
            end
        end
        
        
        % Do boundary correction for S=M/2 if M is even
        if(mod(M,2)==0)
            A(:,S) = 0.5*A(:,S);        
        end
        
        A;
        lambda = eig(A);
        
        phi_vec(M) = log2(max(abs(lambda))); % Take the maximal eigenvalue
        
        phi_mat(:,eps_iter) = phi_vec';
        
        % Now compute the 'naive' phi, without hidden-markov : 
        ind_phi_vec(M) = log2( (p+eps-2*p*eps)^M + (1-(p+eps-2*p*eps))^M );
        ind_phi_double_eps_vec(M) = log2( (p+2*eps-2*p*2*eps)^M + (1-(p+2*eps-2*p*2*eps))^M );
    end
    
    
% % %     % Now compute the 'naive' phi, without hidden-markov : 
% % %     ind_phi = log2( (p+eps-2*p*eps)^2 + (1-(p+eps-2*p*eps))^2 );
% % %     
% % %     ind_phi_vec = repmat(ind_phi, 1,MAX_M);
% % %     
end


% % % OLD : 
% % % figure; hold on; plot(phi_vec, 'o');  plot(ind_phi_vec.* (max(1,[0:MAX_M-1])), 'r*'); plot(0.5*ind_phi_vec.*(max(1,[0:MAX_M-1])), 'g+');  title(['epsilon ' num2str(eps) ' p ' num2str(p) ]); xlabel('M : Chains Number'); ylabel('Phi : log prob. of equal Y'); 
% % % legend('HMM', 'Independent', 'Half Independent');
% % % figure; hold on; plot(phi_vec./  (max(1,[0:length(phi_vec)-1])), 'o');  plot(ind_phi_vec, 'r*'); plot(0.5*ind_phi_vec, 'g+'); title(['epsilon ' num2str(eps) ' p ' num2str(p) ]); xlabel('M : Chains Number'); ylabel('Phi : log prob. of equal Y NORMALIZED'); 
% % % legend('HMM', 'Independent', 'Half Independent');



figure; hold on; plot(phi_vec, 'o');  plot(ind_phi_vec, 'r*'); plot(ind_phi_double_eps_vec, 'g+');  
title(['epsilon ' num2str(eps) ' p ' num2str(p) ]); xlabel('M : Chains Number'); ylabel('Phi : log prob. of equal Y'); 
legend('HMM', 'Independent', 'Independent with double noise');
figure; hold on; plot(phi_vec./  (max(1,[0:length(phi_vec)-1])), 'o');  plot(ind_phi_vec ./ (max(1,[0:length(phi_vec)-1])), 'r*');  
plot(ind_phi_double_eps_vec ./ (max(1,[0:length(phi_vec)-1])), 'g+'); plot(repmat(log2(1-p-eps+2*p*eps), 1, MAX_M), 'm');
title(['epsilon ' num2str(eps) ' p ' num2str(p) ]); xlabel('M : Chains Number'); 
ylabel('Phi : log prob. of equal Y NORMALIZED'); 
legend('HMM', 'Independent', 'Independent with double noise', 'log-limit');


% % % figure; imagesc(phi_mat); colorbar;

