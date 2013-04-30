% Generate HMP entropy plots for paper's fig for p constant (opposite)
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary


res = 5; eps0 = 0.000000001; eps_jump = 0.01; 




    tol=0.0035;

%  eps_vec = [tol:tol:0.5-tol]; % Set the values of epsilon we want

 p_vec = [tol:tol:0.5-tol];
 
 eps_vec = 0.2 * ones(1, length(p_vec)); % Set the epsilon we want to work with
 hold on; subplot(2,1,2);



i=1;

% Do to get two P's ...
%for p = p_vec 
    
   

%    if(p == 0.49)
%        eps_vec = [tol:tol:0.5-tol]; % Set the values of p we want
%    end
    
    ttt = cputime;
  
    
    max_order_taken = 5;  % this says how many orders we use in the approximation
    
    
    % First plot the orders of the HMP in increasing value : 
    Hp = zeros(12, length(eps_vec));
    
    
    mu = 1-2*eps_vec; 

    Hp_zero = log(2.0);  % different since there is no zero in a matlab array
    Hp(2,:) = 2;       
    Hp(4,:) = (4/3) .* (7.*mu.^4-12.*mu.^2+6);        
    Hp(6,:) = (32/15).*(46*mu.^8-120*mu.^6+120.*mu.^4-60.*mu.^2+15);        
    Hp(8,:) = (32/21).*(1137.*mu.^12-4088.*mu.^10+5964.*mu.^8-4536.*mu.^6+1946.*mu.^4-504.*mu.^2+84);        
    Hp(10,:) = (512/45).* (3346.*mu.^16-15120.*mu.^14+28800.*mu.^12-30120.*mu.^10+18990.*mu.^8-7560.*mu.^6+1980.*mu.^4-360.*mu.^2+45);
    Hp(12,:) = (1024/165).*(159230.*mu.^20-874632.*mu.^18+2091100.*mu.^16-2857360.*mu.^14+2465100.*mu.^12-1400960.*mu.^10+532312.*mu.^8-...
        135960.*mu.^6+24145.*mu.^4-3300.*mu.^2+330);
    
    Hp = repmat(-mu.^4,12,1) .* Hp;
    
    
    
    % Compute Upper and Lower bounds 
    N = 2;  % How many bits to take in the upper/lower bound
    
    ind=1;
    C_N = zeros(1,length(eps_vec));     cX_N = zeros(1,length(eps_vec));
    loc_eps = eps_vec(1); % Here all eps are the same
%         loc_eps
%         p
%         N
        C_N = HMP_entropy_finite(loc_eps, p_vec, N)-HMP_entropy_finite(loc_eps, p_vec, N-1); % The upper bound
        cX_N = HMP_entropy_finite_X1(loc_eps, p_vec, N)-HMP_entropy_finite_X1(loc_eps, p_vec, N-1); % The lower bound

    % New plot : Constant epsilon. Verious values of p : 
    ent_orders_est_vec = Hp_zero + (0.5-p_vec).^2 .* Hp(2,:) + (0.5-p_vec).^4 .* Hp(4,:) + ...
        (0.5-p_vec).^6 .* Hp(6,:) + (0.5-p_vec).^8 .* Hp(8,:) + (0.5-p_vec).^10 .* Hp(10,:) + (0.5-p_vec).^12 .* Hp(12,:);; 
    
    ent_orders_est_vec = ent_orders_est_vec ./ log(2.0); % Transfer from natural to binary logarithm
    
    ent_seven_orders_est_vec = Hp_zero + (0.5-p_vec).^2 .* Hp(2,:) + (0.5-p_vec).^4 .* Hp(4,:) + ...
        (0.5-p_vec).^6 .* Hp(6,:) + (0.5-p_vec).^8 .* Hp(8,:) + (0.5-p_vec).^10 .* Hp(10,:); %  + (0.5-p_vec).^12 .* Hp(12,:);; 
    
    ent_seven_orders_est_vec = ent_seven_orders_est_vec ./ log(2.0); % Transfer from natural to binary logarithm
    
    ent_five_orders_est_vec = Hp_zero + (0.5-p_vec).^2 .* Hp(2,:) + (0.5-p_vec).^4 .* Hp(4,:) + ...
        (0.5-p_vec).^6 .* Hp(6,:) + (0.5-p_vec).^8 .* Hp(8,:); %  + (0.5-p_vec).^10 .* Hp(10,:) + (0.5-p_vec).^12 .* Hp(12,:);; 
    
    ent_five_orders_est_vec = ent_five_orders_est_vec ./ log(2.0); % Transfer from natural to binary logarithm
    
    
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
  
     length(ent_orders_est_vec)
     length(C_N)
     length(cX_N)
     
    subplot(2,1,2); hold on; plot(p_vec, ent_orders_est_vec, 'g');  plot(p_vec, ent_seven_orders_est_vec, 'r');  plot(p_vec, ent_five_orders_est_vec, 'c');  %plot([0:101], log(abs([Hp_zero Hp(1:101)'])), '+');  
    plot(p_vec, C_N, 'b'); plot(p_vec, cX_N, '-m');
    title(['A-M Entropy-Rate Expansion for \epsilon= ' num2str(eps_vec(1))]); xlabel('p'); ylabel('Entropy-Rate'); %legend('I.I.D.', 'HMM');
    legend('12-term Expansion', '10-term Expansion', '8-term Expansion', 'Upper Bound', 'Lower Bound');        
%    legend('Expansion', 'Upper Bound', 'Lower Bound');                        
    
    i=i+1;
    %end % loop on epsilons 