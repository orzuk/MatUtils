% Compute twp HMP's KL distance empirically
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary

ttt = cputime;

TOL = 0.0000000001;

% Discrete and Continunuos flags
DISCRETE = 0;
CONTINUOUS = 1;

BSC_flag = 0;

Entropy_flag = 0; % If this flag is on, than the other model is set to uniform !!! 

randomize_matrices = 0;

if(randomize_matrices)
    if(BSC_flag)
        S = 2;
        eps = 0.1;
        R = [eps, 1-eps; 1-eps, eps];
        T = [0.5, -0.5; -1, 1];
    else
        S = 3;
        % Create R, which should be ergodic
        R = abs(10.^rand(S));
        R_tag = abs(10.^rand(S));

        for i=1:S
            R(i,:) = R(i,:) ./ sum(R(i,:));
            R_tag(i,:) = R_tag(i,:) ./ sum(R_tag(i,:));
        end

        % Create T, with negatives on the diagonal and rows sum equal to zero
        T = rand(S); T_tag = rand(S);
        for i=1:S
            T(i,i)  = T(i,i)-sum(T(i,:));
            T_tag(i,i)  = T_tag(i,i)-sum(T_tag(i,:));
        end

        
        % Set T to zeros to see what happens
        %%%T = zeros(S);
        
    end
end



% crerate the vector of stationary distribution
pi = rand(1,S);
pi = pi ./ sum(pi);

pi = ones(1,S); pi = pi ./ S;

xi = ones(S,1);

U = ones(S); U = U ./ S;

U_minus_I = U - eye(S);
U_minus_I(:,end) = 1;
xi_T = xi' * T; xi_T(end) = 0;
phi = xi_T * inv(U_minus_I);
xi_T_tag = xi' * T; xi_T_tag(end) = 0;
phi_tag = xi_T_tag * inv(U_minus_I);


res = 0.00025;
delta_vec = [0.00:res:100*res]; %0.5;
KL_DIST = zeros(1, length(delta_vec));
C_2_vec = zeros(1, length(delta_vec));

HIGH_SNR=0; ALMOST_MEMORYLESS = 1; % regime flags.

save_R = R; save_R_tag = R_tag;


for regime_flag = 0:1
    if(regime_flag == HIGH_SNR)
            M = R; M_tag = R_tag;
    else
        delta_vec = [-50*res:res:50*res]; %0.5;
    end
    cur_delta_ind = 1;
    for delta = delta_vec

        if(regime_flag == HIGH_SNR)
            % Now do the same for the high-SNR regime :
        
            R = eye(S) + delta * T; R_tag = eye(S) + delta * T_tag;

        else % here almost memoryless
            R = save_R; R_tag = save_R_tag;
            M = U + delta * T; M_tag = U + delta * T_tag;
        end
        
        
        
        % We have to find pi :
        [vec lambda] = eig(M');
        ind = 1;
        for i=1:S
            if( abs(lambda(i,i) - 1) < TOL)
                ind = i;
            end
        end

        pi = vec(:,ind) ./ sum(vec(:,ind));

        [vec lambda] = eig(M_tag');
        ind = 1;
        for i=1:S
            if( abs(lambda(i,i) - 1) < TOL)
                ind = i;
            end
        end

        pi_tag = vec(:,ind) ./ sum(vec(:,ind));

        
        % Now find phi etc. 
        xi_T = (1/S) * xi' * T; xi_T(end) = 0;
        phi = -xi_T * inv(U_minus_I); % add a minus sign here !!!!
        xi_T_tag = (1/S) * xi' * T_tag; xi_T_tag(end) = 0;
        phi_tag = -xi_T_tag * inv(U_minus_I); % add a minus sign here !!!!
        
        % Now compute the exact entropy for C_2 (this is only for entropy now) : 
        F = R' * diag(pi) * M* R;         F_tag = R_tag' * diag(pi_tag) * M_tag* R_tag;
        if(Entropy_flag == 1)
            C_2_vec(cur_delta_ind) = (log(pi' * R) * F - xi' * (F .* log(F))) * xi;
        else
             C_2_vec(cur_delta_ind) = -[log(pi' * R) * F - xi' * (F .* log(F))] * xi + [log(pi_tag' * R_tag) * F - xi' * (F .* log(F_tag))] * xi;
        end
            
        if(regime_flag == HIGH_SNR)
            seq_len =1000; num_iters = 50;
        else
            seq_len = 1000; num_iters = 50;
        end

        
        
        
        % Now call the c function :
        if(Entropy_flag == 1)
            KL_DIST(cur_delta_ind) = ComputeHMMKLDistance(pi, M, R, 0, 1, ...
                (1/S) .* xi, U, U, 0, 1, seq_len, num_iters, DISCRETE);
            
        else
            KL_DIST(cur_delta_ind) = ComputeHMMKLDistance(pi, M, R, 0, 1, ...
                pi_tag, M_tag, R_tag, 0, 1, seq_len, num_iters, DISCRETE);
        end
        
        cur_delta_ind = cur_delta_ind+1;
    end % for loop on deltas


    % Comput the KL DIST when delta=0
    if(regime_flag == HIGH_SNR)
        
        if(Entropy_flag == 1)
            KL_zero_order =  sum(pi' * (M .* (log(M)-log(U))));
            KL_first_order = xi' * (    (diag(pi) * M * T + T' * diag(pi) * M) .* log(diag(pi) * M) - diag(log(pi)) * T' * diag(pi) * M     ) * xi;
            
            % Now compute the full order 
            Full_Zero(1) = xi'*diag(pi)*M*T*xi;
            Full_Zero(2) = xi'*[[diag(pi)*M*T].*[log(diag(pi)*M)]]*xi;
            Full_Zero(3) = xi'*T'*diag(pi)*M*xi;
            Full_Zero(4) = xi'* [[T' *diag(pi)*M].*[log(diag(pi)*M)]]*xi;
            Full_Zero(5) = -log(pi')*diag(pi)*M*T*xi; 
            Full_Zero(6) = -log(pi')*T'*diag(pi)*M*xi; 
            Full_Zero(7) = (pi'*T ./ pi')*diag(pi)*M*xi;
            
            Full_Zero
            sum(Full_Zero)
            KL_first_order
            
        else
            KL_zero_order =  sum(pi' * (M .* (log(M)-log(M_tag))));
            KL_first_order = xi' * (    (diag(pi) * M * T + T' * diag(pi) * M) .* log(diag(pi) * M) - diag(log(pi)) * T' * diag(pi) * M     ) * xi; % The contribution from 1st model
            %KL_first_order = KL_first_order - 1; %DUMMY !!! 
            
            
            
               % Now compute the full order 
            Full_Zero(1) = xi'*diag(pi)*M*T*xi;
            Full_Zero(2) = xi'*[[diag(pi)*M*T].*[log(diag(pi)*M)]]*xi;
            Full_Zero(3) = xi'*T'*diag(pi)*M*xi;
            Full_Zero(4) = xi'* [[T' *diag(pi)*M].*[log(diag(pi)*M)]]*xi;
            Full_Zero(5) = -log(pi')*diag(pi)*M*T*xi; 
            Full_Zero(6) = -log(pi')*T'*diag(pi)*M*xi; 
            Full_Zero(7) = -((pi'*T) ./ pi')*diag(pi)*M*xi;
            
            Full_Zero
            sum(Full_Zero)
            KL_first_order
        
            
               % Now compute the full order 
            Full_Zero_KL(1) = xi'*[[diag(pi_tag)*M_tag*T_tag]./[diag(pi_tag)*M_tag].*[diag(pi)*M]]*xi;
            Full_Zero_KL(2) = xi'*[[diag(pi)*M*T].*[log(diag(pi_tag)*M_tag)]]*xi;
            Full_Zero_KL(3) = xi'*[[T_tag'*diag(pi_tag)*M_tag]./[diag(pi_tag)*M_tag].*[diag(pi)*M]]*xi;
            Full_Zero_KL(4) = xi'* [[T' *diag(pi)*M].*[log(diag(pi_tag)*M_tag)]]*xi;
            Full_Zero_KL(5) = -log(pi_tag')*diag(pi)*M*T*xi;  % This is the only zero !!! 
            Full_Zero_KL(6) = -log(pi_tag')*T'*diag(pi)*M*xi; 
            Full_Zero_KL(7) = -((pi_tag'*T_tag) ./ pi_tag')*diag(pi)*M*xi;
            
            KL_first_order_Full = sum(Full_Zero)-sum( Full_Zero_KL)
            
            Full_Zero_KL
            sum(Full_Zero_KL)
            KL_first_order
            
        end
        
        
        
        % Now plot to see what we get :

        figure; hold on; plot(delta_vec, KL_DIST-KL_zero_order, '.');  plot(delta_vec, KL_first_order .* delta_vec, 'r');
        %        plot(delta_vec, log(S)-(C_2_vec+KL_zero_order), '--m'); 
        if(Entropy_flag == 1)
            plot(delta_vec, log(S)-(C_2_vec+KL_zero_order), '--m'); 
        else
            plot(delta_vec, C_2_vec-KL_zero_order, '--m'); 
        end
        plot(delta_vec, KL_first_order_Full .* delta_vec, 'c');
        
        
        legend('simulated', 'analytic first-order', 'C_2', 'first-order-full');
        xlabel('\delta'); ylabel('KL DIST'); title('KL DIST in the high-SNR regime');

    else  % Here regime is 'almost-memoryless'
        dumdum = 1;        
        if(Entropy_flag == 1)            
            KL_zero_order = sum( (1/S).*sum(R,1) .* (log((1/S).*sum(R,1)) - log((1/S).*sum(U,1))) )
            
            % Shortend - Is this correct ????????
            KL_first_order =  (   (1/S) *  xi' *(  (R' * T * R) * log((1/S) * R' * U * R)  )   * xi - log((1/S)*xi'*R) * R' * T' * R * xi )

            % Long - Maybe this is the correct one ????? 
            TTT = R' *[ (1/S).* diag(xi) * T + diag(phi)*U ] * R;
            KL_first_order_full = -( [-xi' * (TTT .* (S.* U + log((1/S)*R'*U*R) ))   + log((1/S) .* xi' * R) * TTT] * xi );
            KL_first_order_reduced = -(1/S) .* [xi' * ((R' * T * R) .* (log((1/S)*R'*U*R) ))] * xi;
        else
            % We need to change also the zero order :
            KL_zero_order = sum( (1/S).*sum(R,1) .* (log((1/S).*sum(R,1)) - log((1/S).*sum(R_tag,1))) )

            KL_first_order = (1/S) .* [xi' * ((R' * T * R) .* (log((1/S)*R'*U*R) ))] * xi;

            KL_first_order = KL_first_order + ...
                xi' * [  (1/S)* (R'*T*R).*(R'*U*R)./(R_tag'*U*R_tag)   + (R'*diag(phi)*U*R).*(R'*U*R)./(R_tag'*U*R_tag) + ...
                (R'*T*R).*log((1/S)*R_tag'*U*R_tag) ]*xi+ [((phi_tag*R_tag)./(xi'*R_tag))*(R'*U*R) ] * xi; % this is the p log q part

           

            % Check which of the six terms are zero !!!!
            Full_Zero  = zeros(1,6);
            TTT_1 = R' *[ (1/S) * T] * R;
            TTT_2 = R' *[diag(phi)*U ] * R;
            %Full_Zero(1) = [-xi' * (TTT_1 .* (S.* U))] * xi;
            %Full_Zero(2) = [-xi' * (TTT_2 .* (S.* U))] * xi;
            Full_Zero(1) = [-xi' * (TTT_1 )] * xi; % removed the multiplication by S*U
            Full_Zero(2) = [-xi' * (TTT_2 )] * xi; % removed the multiplication by S*U
            Full_Zero(3) = [-xi' * (TTT_1 .* (log((1/S)*R'*U*R) ))] * xi;
            Full_Zero(4) = [-xi' * (TTT_2 .* (log((1/S)*R'*U*R) ))] * xi;
            Full_Zero(5) = [log((1/S) .* xi' * R) * TTT_1] * xi;
            Full_Zero(6) = [log((1/S) .* xi' * R) * TTT_2] * xi;
            extra_term = [((phi*R)./(xi'*R)) * R'*U*R] * xi


            % Check the six terms for the KL :
            Full_Zero_KL  = zeros(1,6);
            TTT_tag_1 = R_tag' *[ (1/S) * T_tag] * R_tag;
            TTT_tag_2 = R_tag' *[diag(phi_tag)*U ] * R_tag;
            Full_Zero_KL(1) = [-xi' * (TTT_tag_1 .* (R'*U*R) ./ (R_tag'*U*R_tag))] * xi;
            Full_Zero_KL(2) = [-xi' * (TTT_tag_2 .* (R'*U*R) ./ (R_tag'*U*R_tag))] * xi;
            Full_Zero_KL(3) = [-xi' * (TTT_1 .* (log((1/S)*R_tag'*U*R_tag) ))] * xi;
            Full_Zero_KL(4) = [-xi' * (TTT_2 .* (log((1/S)*R_tag'*U*R_tag) ))] * xi;
            Full_Zero_KL(5) = [log((1/S) .* xi' * R_tag) * TTT_1] * xi; 
            Full_Zero_KL(6) = [log((1/S) .* xi' * R_tag) * TTT_2] * xi;
            extra_term_KL = [((phi_tag*R_tag)./(xi'*R_tag)) * R'*U*R] * xi;


            KL_first_order = [-xi' * (TTT_1 .* (log((1/S)*R'*U*R) ))] * xi + ...
                [ [xi' * (TTT_tag_1 .* (R'*U*R) ./ (R_tag'*U*R_tag))] * xi + [xi' * (TTT_tag_2 .* (R'*U*R) ./ (R_tag'*U*R_tag))] * xi + ...
                [xi' * (TTT_1 .* (log((1/S)*R_tag'*U*R_tag) ))] * xi - [((phi_tag*R_tag)./(xi'*R_tag)) * R'*U*R] * xi ];           
            KL_first_order_full = sum(Full_Zero)+extra_term - sum(Full_Zero_KL) - extra_term_KL; 
            
            sum(Full_Zero_KL) - extra_term_KL; 
            KL_first_order_reduced =  KL_first_order_full;
        end


        % Now plot to see what we get :
        figure; hold on; plot(delta_vec, KL_DIST-KL_zero_order, '.');  plot(delta_vec, -KL_first_order .* delta_vec, 'r');
        plot(delta_vec, KL_first_order_full .* delta_vec, 'c');  % try to catch the opposite
        
        if(Entropy_flag == 1)
            plot(delta_vec, log(S)-(C_2_vec+KL_zero_order), '--m'); 
        else
            plot(delta_vec, C_2_vec-KL_zero_order, '--m'); 
        end
        plot(delta_vec, -KL_first_order_reduced .* delta_vec, '+k');
        legend('simulated', 'analytic 1^{st}-order', 'Full analytic 1^{st}-order', 'C_2', 'Reduced Analytic 1^{st} order');
        xlabel('\delta'); ylabel('KL DIST'); title('KL DIST in the almost-memoryless regime');

        
    end
end % for loop on regimes



total_time = cputime - ttt








% Check which of the six terms are zero !!!! 
TTT_1 = R' *[ (1/S) * T] * R;
TTT_2 = R' *[diag(phi)*U ] * R;
%Full_Zero(1) = [-xi' * (TTT_1 .* (S.* U))] * xi;
%Full_Zero(2) = [-xi' * (TTT_2 .* (S.* U))] * xi;
Full_Zero(1) = [-xi' * (TTT_1 )] * xi; % removed the multiplication by S*U
Full_Zero(2) = [-xi' * (TTT_2 )] * xi; % removed the multiplication by S*U
Full_Zero(3) = [-xi' * (TTT_1 .* (log((1/S)*R'*U*R) ))] * xi;
Full_Zero(4) = [-xi' * (TTT_2 .* (log((1/S)*R'*U*R) ))] * xi;
Full_Zero(5) = [log((1/S) .* xi' * R) * TTT_1] * xi;
Full_Zero(6) = [log((1/S) .* xi' * R) * TTT_2] * xi;
extra_term = [((phi*R)./(xi'*R)) * R'*U*R] * xi

Full_Zero
sum(Full_Zero)
KL_first_order_full
[log((1/S) .* xi' * R) * R' *[diag(phi)*U ] * R] * xi
 [-xi' * ((R' *[ (1/S).* T] * R) .* (log((1/S)*R'*U*R) ))] * xi

% Check the six terms for the KL :  
TTT_tag_1 = R_tag' *[ (1/S) * T_tag] * R_tag;
TTT_tag_2 = R_tag' *[diag(phi_tag)*U ] * R_tag;
Full_Zero_KL(1) = [-xi' * (TTT_1 .* (R'*U*R) ./ (R_tag'*U*R_tag))] * xi;
Full_Zero_KL(2) = [-xi' * (TTT_2 .* (R'*U*R) ./ (R_tag'*U*R_tag))] * xi;
Full_Zero_KL(3) = [-xi' * (TTT_1 .* (log((1/S)*R_tag'*U*R_tag) ))] * xi;
Full_Zero_KL(4) = [-xi' * (TTT_2 .* (log((1/S)*R_tag'*U*R_tag) ))] * xi;
Full_Zero_KL(5) = [log((1/S) .* xi' * R_tag) * TTT_1] * xi;
Full_Zero_KL(6) = [log((1/S) .* xi' * R_tag) * TTT_2] * xi;
extra_term_KL = [((phi_tag*R_tag)./(xi'*R_tag)) * R'*U*R] * xi


Full_Zero_KL
sum(Full_Zero_KL)

 
 
 

  % Now plot to see what we get :
%         figure; hold on; plot(delta_vec, KL_DIST, '.');  plot(delta_vec,KL_zero_order + KL_first_order .* delta_vec, 'r'); 
%         plot(delta_vec, KL_zero_order-KL_first_order_full .* delta_vec, 'c');  plot(delta_vec, log(S)-(C_2_vec), '--m'); 
%         plot(delta_vec, KL_zero_order + KL_first_order_reduced .* delta_vec, '+k');
%         legend('simulated', 'analytic 1^{st}-order', 'Full analytic 1^{st}-order', 'C_2', 'Reduced Analytic 1_st order');
%         xlabel('\delta'); ylabel('KL DIST'); title('KL DIST in the almost-memoryless regime');




% B = R'* diag(pi) * M * R;
% 
% 
% 
% D  = R' * diag(pi) * U * R; % use D for simplicity
% 
% Exact_H = ( xi' * (B .* log(B)) - (log(pi*R)) * B) * xi
% 
% try_another_exact = sum(sum(B .* log(B))) - sum((log(pi*R)) * B)
% 
% 
% zero_order_should_be = xi' * ( D .* log(D))- log(pi * R) * D;
% zero_order_should_be = zero_order_should_be * xi
% log_S_is = log(S)
% 
% 
% try_try = (1/S) .* xi' * R * log(R' * xi) - log(S)
% 
% 
% 
% UIUI = eye(S)-U;
% UIUI(:,S) = 1;
% 
% piti = pi * T; piti(S)  = 0;
% 
% alpha = pi + delta * piti * inv(UIUI)
% 
% [vec lam] = eig(M');
% 
% ind = 1;
% for i=1:S
%     if( abs(lam(i,i) - 1) < TOL)
%         ind = i;
%     end
% end
% 
% alpha_correct = vec(:,ind) ./ sum(vec(:,ind))
% 
% alpha - alpha_correct'
% 
% % Now that we have alpha we can calculate the first order also :
% BIG_TEMP = R' * ( diag(pi) * T + diag (piti * inv(UIUI)) * U) * R;
% first_order_should_be = ( xi' * (BIG_TEMP .* (S .* U + log(D))) - log(pi * R) * BIG_TEMP ) * xi
% 
% 
% DDD = T' * diag(pi) * M + diag(pi) * M * T
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% p = 0.22; N = 80; iters = 4000;
% res = 100;
% 
% ttt = cputime;
% eps_vec = zeros(1,res);
% first_order = zeros(1,res);
% ent_diff = zeros(1,res);
% 
% for r=1:res
%     eps = 0.000000001 + (0.02*(r-1)/res);
%     eps_vec(r) = eps;
%     P = [1-p p; p 1-p];
%     EPS = [1-eps eps; eps 1-eps];
%     
%     sam_ent = 0;
%     
%     Y_Vec = sample_HMM_vec(P,EPS,N,iters);
%     
%     condprobs_vec = compute_HMM_condprob_vec(P, EPS, Y_Vec);
%     
%     EEEEEEE = entropy  ([condprobs_vec 1-condprobs_vec]');
%     
%     sam_ent = sum(  entropy  ([condprobs_vec 1-condprobs_vec]')) / iters;
%     
%     
% % % % % % % % %     for i=1:iters  % randomize and calculate entropy
% % % % % % % % %         
% % % % % % % % %    %%%     Y = sample_HMM(P, EPS, N); % sample the sequence
% % % % % % % % %         
% % % % % % % % %         Y = Y_Vec(i,:);
% % % % % % % % %         
% % % % % % % % %         condp = compute_HMM_condprob(P, EPS, Y);
% % % % % % % % %         
% % % % % % % % %         sam_ent = sam_ent + entropy([condp 1-condp]');
% % % % % % % % %         
% % % % % % % % %     end
% % % % % % % % %     
% % % % % % % % % % % %     r
% % % % % % % % % % % %     cputime-ttt
% % % % % % % % %     sam_ent = sam_ent / iters;
%     
%     
% %     (sam_ent - H(p)) / eps
% %     2*(2*p-1) * log2(p/(1-p))
%     
%     
%     first_order(r) = (sam_ent - entropy([p, 1-p]')) / eps;
%     ent_diff(r) = (sam_ent - entropy([p, 1-p]'));
%     
%     if(mod(r,10) == 0)
%         cur_iter = r
%     end
%    
%     did_res = r
% end
% 
% 
% cputime - ttt
% 
% % Comparison with what should be the correct answer
% eytan_val = repmat(2*(2*p-1) * log2(p/(1-p)), 1, res);
% 
% 
% RRR = entropy([ p+eps_vec - 2*p* eps_vec 1-(p+eps_vec - 2*p*eps_vec) ] )
% EEE = entropy([ p+eps_vec(1) - 2*p* eps_vec(1) 1-(p+eps_vec(1) - 2*p*eps_vec(1)) ]' );
% PPP = entropy([p 1-p]');
% DDD = (entropy([ p+eps_vec(1) - 2*p* eps_vec(1) 1-(p+eps_vec(1) - 2*p*eps_vec(1)) ]' ) - entropy([p 1-p]')) / eps_vec(1);
% 
% 
% new_eytan_val_half =  (entropy([ (p+eps_vec - 2*p* eps_vec)' (1-(p+eps_vec - 2*p*eps_vec))' ]') - entropy([p 1-p]')) ./ eps_vec
% 
% new_eps_vec = 2*eps_vec - 2*eps_vec .* eps_vec;
% 
% new_eytan_val =  (entropy([ (p+new_eps_vec - 2*p* new_eps_vec)' (1-(p+new_eps_vec - 2*p*new_eps_vec))' ]') - entropy([p 1-p]')) ./ eps_vec
% 
% 
% 
% 
% % Now insert also the 2nd term from taylor expansion
% eytan_val_2nd_term = eytan_val -   (-0.5   + 4*(1-2*p)*log((1-p)/p) + 0.5*(  (p*p+(1-p)*(1-p))/(p*(1-p))  )^2 )*eps_vec;
% 
% % Now insert also the 2nd term from taylor expansion
% eytan_val_2nd_term_corrected = eytan_val -   (-0.5   + 4*(1-2*p)*log((1-p)/p) + 0.5*(  1/(p*(1-p)) - 3  )^2 )*eps_vec;
% 
% Ido_2nd_term = eytan_val - (  2*(1-2*p)*log((1-p)/p)  + 0.5 * (1-2*p)*(1-2*p) /   (p*p*(1-p)*(1-p))  )*eps_vec;
% 
% figure; hold on; plot(eps_vec, first_order, '.');  plot(eps_vec, eytan_val, 'r');  %plot(eps_vec, eytan_val/2, 'g'); 
%  plot(eps_vec, eytan_val_2nd_term, 'm'); plot(eps_vec, eytan_val_2nd_term_corrected, 'k'); plot(eps_vec, Ido_2nd_term, 'g');
% 
%  %%%plot(eps_vec, new_eytan_val, 'm');  plot(eps_vec, new_eytan_val_half, 'k'); 
% % % % legend('points', '1st order', 'half 1st order', 'double_noise', 'one_noise');
% legend('points', '1st order',  '2nd order', '2nd order corrected', '2nd order Ido');
% title(['1st order coef. of epsilon for p='  num2str(p) ' N = ' num2str(N)]); xlabel('epsilon'); ylabel('1st order');
% 
% 
%     % Randomize Y according to the probabilities
%     % Now we have to calculate for a hidden markov process.  X_n --> Y_n
%     % This is done via viterby algorithm
%     % We need to calculate H(Y_n | Y_{n-1}, .., Y_1), or H(Y_n | Y_{n-1}, .., Y_1, X_1)
%     %
%     % For this, we will calculate the probability of Pr(Y_n = 1, Y_{n-1} = ?,
%     % .. Y_1 = ? ) and Pr(Y_n = 1, Y_{n-1} = ?, .. Y_1 = ?, X_1 = ?)
%     
% 
%     
%     
% % New plot :
% figure; hold on; plot(eps_vec, ent_diff, '.');  plot(eps_vec, 2 * (1-2*p) * log2((1-p)/p) .* eps_vec, 'r');
% legend('simulated', 'analytic first-order');
% 
% 
