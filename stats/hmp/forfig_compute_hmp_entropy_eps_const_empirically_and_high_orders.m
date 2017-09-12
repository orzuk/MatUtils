% Generate HMP entropy plots for paper's figs (empirical and high-order) 
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary


res = 5; eps0 = 0.000000001; eps_jump = 0.01; 




eps_vec = [0.1]; % Set the epsilon we want to work with

figure; hold on; subplot(2,1,1);

i=1;

% Do to get two epsilons ...
for eps = eps_vec 
    
    tol=0.0025;
    p_vec = [0.34:tol:0.5-tol]; % Set the values of p we want

    if(eps == 0.01)
        p_vec = [0.08:tol:0.2-tol]; % Set the values of p we want
    end
    
    ttt = cputime;
    eps_vec = zeros(1,res);
    first_order = zeros(1,res);
    
    max_order_taken = 5;  % this says how many orders we use in the approximation
    
    
    % First plot the orders of the HMP in increasing value : 
    Hp = zeros(12, length(p_vec));
    
    
    
    Hp_zero = -(p_vec.*log(p_vec)+(1-p_vec).*log(1-p_vec));  % different since there is no zero in a matlab array
    Hp(1,:) = 2 *(1-2.*p_vec) .* log ( (1-p_vec)./p_vec );
    Hp(2,:) = -2 *(1-2.*p_vec) .* log ( (1-p_vec)./p_vec ) -(1-2.*p_vec).^2./(2.*p_vec.^2.*(1-p_vec).^2);
    
    % Now represent anything as a function of lambda = 1-2p
    lam_vec = 1-2.*p_vec; 
    
    Hp(3,:) = (-16.*(5.*lam_vec.^4-10.*lam_vec.^2-3) .* lam_vec.^2) ./ ( 3.*(1-lam_vec.^2).^4);
    
    Hp(4,:) = (8.*(109.*lam_vec.^8+20.*lam_vec.^6-114.*lam_vec.^4-140.*lam_vec.^2-3).*lam_vec.^2) ./ (3.*(1-lam_vec.^2).^6);
    
    Hp(5,:) = (-128.*(95.*lam_vec.^10+336.*lam_vec.^8+762.*lam_vec.^6-708.*lam_vec.^4-769.*lam_vec.^2-100).*lam_vec.^4) ./ (15.*(1-lam_vec.^2).^8);
    
    Hp(6,:) = 128.*(125.*lam_vec.^14-321.*lam_vec.^12+9525.*lam_vec.^10+16511.*lam_vec.^8-7825.*lam_vec.^6- ...
        17995.*lam_vec.^4-4001.*lam_vec.^2-115).*lam_vec.^4./   (15.*(1-lam_vec.^2).^10);
    
    Hp(7,:) = -256.*(280.*lam_vec.^18-45941.*lam_vec.^16-110888.*lam_vec.^14+666580.*lam_vec.^12+1628568.*lam_vec.^10- ...
        270014.*lam_vec.^8-1470296.*lam_vec.^6-524588.*lam_vec.^4-37296.*lam_vec.^2-245).*lam_vec.^4 ./ (105.*(1-lam_vec.^2).^12);
    
    Hp(8,:) = 64.*(56.*lam_vec.^22-169169.*lam_vec.^20-2072958.*lam_vec.^18-5222301.*lam_vec.^16+12116328.*lam_vec.^14+ ...
        35666574.*lam_vec.^12+3658284.*lam_vec.^10-29072946.*lam_vec.^8-14556080.*lam_vec.^6- ...
        1872317.*lam_vec.^4-48286.*lam_vec.^2-49).*lam_vec.^4 ./ (21.*(1-lam_vec.^2).^14);
    
    Hp(9,:) = 2048.*(37527.*lam_vec.^22+968829.*lam_vec.^20+8819501.*lam_vec.^18+20135431.*lam_vec.^16-23482698.*lam_vec.^14- ...
        97554574.*lam_vec.^12-30319318.*lam_vec.^10+67137630.*lam_vec.^8+46641379.*lam_vec.^6+8950625.*lam_vec.^4+ ...
        495993.*lam_vec.^2+4683).*lam_vec.^6 ./ (63.*(1-lam_vec.^2).^16);
    
    Hp(10,:) = -2048.*(38757.*lam_vec.^26+1394199.*lam_vec.^24+31894966.*lam_vec.^22+243826482.*lam_vec.^20+ ...
        571835031.*lam_vec.^18-326987427.*lam_vec.^16-2068579420.*lam_vec.^14-1054659252.*lam_vec.^12+1173787011.*lam_vec.^10+ ...
        1120170657.*lam_vec.^8+296483526.*lam_vec.^6+26886370.*lam_vec.^4+ 684129.*lam_vec.^2+2187).*lam_vec.^6./  (45.*(1-lam_vec.^2).^18);
    
    Hp(11,:) = 8192.*(98142.*lam_vec.^30-1899975.*lam_vec.^28+92425520.*lam_vec.^26+3095961215.*lam_vec.^24+ ...
        25070557898.*lam_vec.^22+59810870313.*lam_vec.^20-11635283900.*lam_vec.^18-173686662185.*lam_vec.^16- ...
        120533821070.*lam_vec.^14+74948247123.*lam_vec.^12+102982107048.*lam_vec.^10+35567469125.*lam_vec.^8+ ...
        4673872550.*lam_vec.^6+217466315.*lam_vec.^4+2569380.*lam_vec.^2+2277).*lam_vec.^6./ (495.*(1-lam_vec.^2).^20);
    
    
    
    % Compute Upper and Lower bounds 
    N = 2;  % How many bits to take in the upper/lower bound
    C_N = HMP_entropy_finite(eps, p_vec, N)-HMP_entropy_finite(eps, p_vec, N-1); % The upper bound
    cX_N = HMP_entropy_finite_X1(eps, p_vec, N)-HMP_entropy_finite_X1(eps, p_vec, N-1); % The lower bound
    
    
    % New plot : Constant epsilon. Verious values of p : 
    ent_orders_est_vec = Hp_zero + eps .* Hp(1,:) + (eps^2) .* Hp(2,:) + eps^3 .* Hp(3,:) + eps^4 .* Hp(4,:) + eps^5 .* Hp(5,:) + ...
        eps^6 .* Hp(6,:) + eps^7 .* Hp(7,:) + eps^8 .* Hp(8,:) + eps^9 .* Hp(9,:) + eps^10 .* Hp(10,:) + eps^11 .* Hp(11,:); 
    
    ent_orders_est_vec = ent_orders_est_vec ./ log(2.0); % Transfer from natural to binary logarithm
    
    % figure; subplot(1,2,1); hold on; 
    % subplot(1,2,1);  hold on; plot(p_vec, ent_orders_est_vec, 'r'); %plot([0:101], log(abs([Hp_zero Hp(1:101)'])), '+');  
    % plot(p_vec, C_N, '.b'); plot(p_vec, cX_N, '.m');
    % title(['HMP Order-Based Estimated Entropy for epsilon= ' num2str(eps)]); xlabel('p'); ylabel('Entropy Estimation'); %legend('I.I.D.', 'HMM');
    % legend('Expansion', 'Upper Bound', 'Lower Bound');                        
    % 
    % subplot(1,2,2); hold on; plot(p_vec, min((1-ent_orders_est_vec') ./ (1-0.5*(C_N+cX_N)), 3)   , '+r'); %plot([0:101], log(abs([Hp_zero Hp(1:101)'])), '+');  
    % plot(p_vec, min((1-C_N)./ (1-(C_N+cX_N)./2), 3), '.b'); plot(p_vec, min((1-cX_N) ./ (1-(C_N+cX_N)./2), 3), '.m');
    % title(['HMP Order-Based Estimated Entropy for epsilon= ' num2str(eps)]); xlabel('p'); ylabel('Entropy Estimation - (1-H)/(1-average(Upper,Lower))'); %legend('I.I.D.', 'HMM');
    % legend('1-Expansion', '1-UpperBound', '1-LowerBound');                        
    
    
    % Now the first plot in a seperate figure 
     hold on; 
    subplot(2,1,i); hold on; plot(p_vec, ent_orders_est_vec, 'r'); %plot([0:101], log(abs([Hp_zero Hp(1:101)'])), '+');  
    plot(p_vec, C_N, '.b'); plot(p_vec, cX_N, '-.m');
    title(['HMP Order-Based Estimated Entropy for epsilon= ' num2str(eps)]); xlabel('p'); ylabel('Entropy Estimation'); %legend('I.I.D.', 'HMM');
    legend('Expansion', 'Upper Bound', 'Lower Bound');                        
    
    i=i+1;
end % loop on epsilons 