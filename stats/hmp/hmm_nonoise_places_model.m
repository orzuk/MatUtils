% Plot the entropy of processes with different noises.
% The model is as follows : 
% p - the transition prob.
% eps - the total noise per-bit
% K - the ratio of not-noised places. K=2,3,...
% The model assumes : On places n = 0 mod K, we have no noise(Y_n = X_n)
% If n != 0 mod K, we have a noise of eps * (K-1) / K. So : 
%   Pr(Y_n = X_n) = 1 - eps*(k-1)/K

eps_tol = 0.0001;

eps_vec = [eps_tol:eps_tol:0.5];  % Vector of epsilons to consider
K_vec = [2:7];    % number of models. 

p = 0.1;   % Currently only one value of P !!!!

loc_ent_vec = zeros(K_vec(end),length(eps_vec));

% Now calculate the entropy (naively ....)
for K=K_vec
    
    NOWKKKKKKK = K 
    loc_eps_vec = eps_vec * (K/(K-1));
        
    % Initilize prob. (y_1|x_0=0) : 
    p_zero =  [1-loc_eps_vec ; loc_eps_vec] * (1-p);
    p_one = [loc_eps_vec; 1-loc_eps_vec] * p;
    
 %   begin_psize = size(p_zero)
 %   begin_psize = size(p_one)
    
    % add the contribution of all the K cases 
    for j=1:K-2  % Note : The last iteration is different !!!!!!
        % Compute the 'convulution' probability conv_prob, which is a
        % function of p, eps and the 'history length'.
        % p_zero + p_one = Pr(y_{1->k} | x_0=0) 
        new_p_zero = p_zero * (1-p) + p_one * p;
        new_p_zero = [new_p_zero .* repmat(1-loc_eps_vec, 2^j, 1); new_p_zero .* repmat(loc_eps_vec, 2^j, 1)];   
        
        new_p_one = p_zero * p + p_one * (1-p);
        new_p_one = [new_p_one .* repmat(loc_eps_vec, 2^j, 1); new_p_one .* repmat(1-loc_eps_vec, 2^j, 1)];   
        
        % update for next round
        p_zero = new_p_zero;
        p_one = new_p_one;
        
    end    
    
    % Now do the last iteration : Here X==Y
    %     # Note : The last iteration is different !!!  because X == Y !!!j := K-1;
    % > for ind from 0 to 2^(j)-1 do 
    % > new_p_zero[ind] := (p_zero[ind] * (1-p) + p_one[ind] * p) * 1;
    % > new_p_zero[2^j+ind] := 0;
    % > new_p_one[ind] := 0;
    % > new_p_one[2^j+ind] := (p_zero[ind] * p + p_one[ind] * (1-p)) * 1;
    % > end do: 
    % > # Now copy back to p_zero and p_one : 
    % > for ind from 0 to 2^(j+1)-1 do 
    % > p_zero[ind] := new_p_zero[ind];
    % > p_one[ind] := new_p_one[ind];
    % > end do:    
    % > # Ended last iteration
    
    j=K-1; 
    new_p_zero = p_zero * (1-p) + p_one * p;
    new_p_zero = [new_p_zero ; new_p_zero .* 0];   
    
    new_p_one = p_zero * p + p_one * (1-p);
    new_p_one = [new_p_one .* 0; new_p_one ];   
    
    % update for next round
    p_zero = new_p_zero;
    p_one = new_p_one;


    
    
            
    % Compute the entropy 
    loc_ent_size = size(loc_ent_vec(K,:))
    p_size = size(new_p_zero + new_p_one)
    entropy_size = size(entropy((new_p_zero + new_p_one)))
    loc_ent_vec(K,:) = (entropy((new_p_zero + new_p_one))) ./ K;    
    
end 


% Now plot the resulting entropies we've got : 
figure; imagesc(loc_ent_vec(2:end,:)); colorbar; xlabel('epsilon'); ylabel('K'); title(['Entropy rate for processes, p = ' num2str(p)]);

% Now plot for comparison the true entropy up to order 11 : First do orders 0-2
true_ent_first_orders = -(p*log(p) + (1-p)*log(1-p)) * (eps_vec .^ 0) + ...
                        2*(1-2*p) * log((1-p)/p) * eps_vec + ... 
                        (-2*(1-2*p) * log((1-p)/p) -  (1-2*p)^2/2*p^2*(1-p)^2) * (eps_vec .^ 2); 
true_ent_first_orders = true_ent_first_orders ./ log(2.0);  % transfer to bits !!!                     


% Now add the orders 3-11 : 
lambda = 1-2*p;

H_4 = 8*(109*lambda^8+20*lambda^6-114*lambda^4-140*lambda^2-3)*lambda^2 / (3*(1-lambda^2)^6);

H_5 = -128*(95*lambda^10+336*lambda^8+762*lambda^6-708*lambda^4-769*lambda^2-100)*lambda^4 / (15*(1-lambda^2)^8);

H_6 = 128*(125*lambda^14-321*lambda^12+9525*lambda^10+16511*lambda^8-7825*lambda^6- ...
17995*lambda^4-4001*lambda^2-115)*lambda^4/ (15*(1-lambda^2)^10);

H_7 = -256*(280*lambda^18-45941*lambda^16-110888*lambda^14+666580*lambda^12+1628568*lambda^10- ...
270014*lambda^8-1470296*lambda^6-524588*lambda^4-37296*lambda^2-245)*lambda^4/ (105*(1-lambda^2)^12);

H_8 = 64*(56*lambda^22-169169*lambda^20-2072958*lambda^18-5222301*lambda^16+12116328*lambda^14+ ...
35666574*lambda^12+3658284*lambda^10-29072946*lambda^8-14556080*lambda^6- ...
1872317*lambda^4-48286*lambda^2-49)*lambda^4/ (21*(1-lambda^2)^14);

H_9 = 2048*(37527*lambda^22+968829*lambda^20+8819501*lambda^18+20135431*lambda^16-23482698*lambda^14- ...
97554574*lambda^12-30319318*lambda^10+67137630*lambda^8+46641379*lambda^6+8950625*lambda^4+ ...
495993*lambda^2+4683)*lambda^6/ (63*(1-lambda^2)^16);

H_10 = -2048*(38757*lambda^26+1394199*lambda^24+31894966*lambda^22+243826482*lambda^20+ ...
571835031*lambda^18-326987427*lambda^16-2068579420*lambda^14-1054659252*lambda^12+ ...
1173787011*lambda^10+1120170657*lambda^8+296483526*lambda^6+26886370*lambda^4+ ...
684129*lambda^2+2187)*lambda^6/ (45*(1-lambda^2)^18);

H_11 = 8192*(98142*lambda^30-1899975*lambda^28+92425520*lambda^26+3095961215*lambda^24+ ...
25070557898*lambda^22+59810870313*lambda^20-11635283900*lambda^18-173686662185*lambda^16- ...
120533821070*lambda^14+74948247123*lambda^12+102982107048*lambda^10+35567469125*lambda^8+ ...
4673872550*lambda^6+217466315*lambda^4+2569380*lambda^2+2277)*lambda^6/ (495*(1-lambda^2)^20);


true_ent_first_orders = true_ent_first_orders + (-16*(5*lambda^4-10*lambda^2-3)*lambda^2 / (3*(1-lambda^2)^4)) * (eps_vec .^ 3) + ...
                  H_4 * (eps_vec .^ 4) + H_5 * (eps_vec .^ 5) + H_6 * (eps_vec .^ 6) + H_7 * (eps_vec .^ 7) + H_8 * (eps_vec .^ 8); % + ...
%                 H_9 * (eps_vec .^ 9) +  H_10 * (eps_vec .^ 10) +  H_11 * (eps_vec .^ 11);    

figure; hold on; plot(eps_vec, loc_ent_vec(2,:)); plot(eps_vec, loc_ent_vec(3,:), 'r'); 
plot(eps_vec, loc_ent_vec(4,:), 'g'); plot(eps_vec, loc_ent_vec(5,:), 'm'); plot(eps_vec, loc_ent_vec(6,:), 'k');  plot(eps_vec, loc_ent_vec(7,:), 'c'); 
% plot(eps_vec, true_ent_first_orders, '--');
legend('K=2', 'K=3', 'K=4', 'K=5', 'K=6','K=7', 2); % 'True Order 8', 2); 
title(['Entropy rate for p=' num2str(p)]); xlabel('eps'); ylabel('entropy rate');

% 
% figure; hold on; plot(eps_vec(1:500), loc_ent_vec(2,1:500)); plot(eps_vec(1:500), loc_ent_vec(3,1:500), 'r'); 
% plot(eps_vec(1:500), loc_ent_vec(4,1:500), 'g'); plot(eps_vec(1:500), loc_ent_vec(5,1:500), 'm'); 
% plot(eps_vec(1:500), loc_ent_vec(6,1:500), 'k');  plot(eps_vec(1:500), loc_ent_vec(7,1:500), 'c'); 
% plot(eps_vec(1:500), true_ent_first_orders(1:500), '--');
% legend('K=2', 'K=3', 'K=4', 'K=5', 'K=6','K=7', 'True Order 2', 2); title(['Entropy rate for p=' num2str(p)]); xlabel('eps'); ylabel('entropy rate');

                    