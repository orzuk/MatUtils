% Compute empirical entropy and compare to high-order approximation 
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary

p = 0.3; N = 100; iters = 220;
res = 5; eps0 = 0.000000001; eps_jump = 0.02; 

ttt = cputime;
eps_vec = zeros(1,res);
first_order = zeros(1,res);

max_order_taken = 5;  % this says how many orders we use in the approximation

% First compute the orders for the i.i.d. case - Not markovian !
HpConv = zeros(110,1); 
HpConv_zero = -(p*log(p)+(1-p)*log(1-p));  % different since there is no zero in a matlab array
HpConv(1) = (1-2*p) * log ( (1-p)/p );  % Half from the markovian case

for n=2:110
    HpConv(n) = (-1) * ( (2*p-1)^n / p^(n-1) +  (1-2*p)^n / (1-p)^(n-1) ) / (n*(n-1)); 
end

% First plot the orders in increasing value : 
Hp = zeros(12, 1);

Hp_zero = -(p*log(p)+(1-p)*log(1-p));  % different since there is no zero in a matlab array
Hp(1) = 2 *(1-2*p) * log ( (1-p)/p );
Hp(2) = -2 *(1-2*p) * log ( (1-p)/p ) -(1-2*p)^2/(2*p^2*(1-p)^2);

% Now represent anything as a function of lambda = 1-2p
lam = 1-2*p; 

Hp(3) = (-16*(5*lam^4-10*lam^2-3) * lam^2) / ( 3*(1-lam^2)^4);

Hp(4) = (8*(109*lam^8+20*lam^6-114*lam^4-140*lam^2-3)*lam^2) / (3*(1-lam^2)^6)

Hp(5) = (-128*(95*lam^10+336*lam^8+762*lam^6-708*lam^4-769*lam^2-100)*lam^4) / (15*(1-lam^2)^8);

Hp(6) = 128*(125*lam^14-321*lam^12+9525*lam^10+16511*lam^8-7825*lam^6- ...
17995*lam^4-4001*lam^2-115)*lam^4/   (15*(1-lam^2)^10);

Hp(7) = -256*(280*lam^18-45941*lam^16-110888*lam^14+666580*lam^12+1628568*lam^10- ...
270014*lam^8-1470296*lam^6-524588*lam^4-37296*lam^2-245)*lam^4 / (105*(1-lam^2)^12);

Hp(8) = 64*(56*lam^22-169169*lam^20-2072958*lam^18-5222301*lam^16+12116328*lam^14+ ...
35666574*lam^12+3658284*lam^10-29072946*lam^8-14556080*lam^6- ...
1872317*lam^4-48286*lam^2-49)*lam^4 / (21*(1-lam^2)^14);

Hp(9) = 2048*(37527*lam^22+968829*lam^20+8819501*lam^18+20135431*lam^16-23482698*lam^14- ...
97554574*lam^12-30319318*lam^10+67137630*lam^8+46641379*lam^6+8950625*lam^4+ ...
495993*lam^2+4683)*lam^6 / (63*(1-lam^2)^16);

Hp(10) = -2048*(38757*lam^26+1394199*lam^24+31894966*lam^22+243826482*lam^20+ ...
571835031*lam^18-326987427*lam^16-2068579420*lam^14-1054659252*lam^12+1173787011*lam^10+ ...
1120170657*lam^8+296483526*lam^6+26886370*lam^4+ 684129*lam^2+2187)*lam^6/  (45*(1-lam^2)^18);

Hp(11) = 8192*(98142*lam^30-1899975*lam^28+92425520*lam^26+3095961215*lam^24+ ...
25070557898*lam^22+59810870313*lam^20-11635283900*lam^18-173686662185*lam^16- ...
120533821070*lam^14+74948247123*lam^12+102982107048*lam^10+35567469125*lam^8+ ...
4673872550*lam^6+217466315*lam^4+2569380*lam^2+2277)*lam^6/ (495*(1-lam^2)^20);



% First plot the orders for a given p
% figure; hold on; plot([0:130], log(abs([HpConv_zero HpConv(1:130)'])), '+'); plot([0:129], HpConv(1:130)' ./ [HpConv_zero HpConv(1:129)'], 'r*');
% title(['Ind. Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); legend('Log', 'Cosecutives Ratio');
% 
% figure; hold on; plot([0:11], log(abs([Hp_zero Hp(1:11)'])), '+'); plot([0:10], Hp(1:11)' ./ [Hp_zero Hp(1:10)'], 'r*');
% title(['HMM Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); legend('Log', 'Cosecutives Ratio');

figure; hold on; plot([0:101], log(abs([HpConv_zero HpConv(1:101)'])), '*r'); %plot([0:101], log(abs([Hp_zero Hp(1:101)'])), '+');  
title(['IID epsilon Log Abs Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); %legend('I.I.D.', 'HMM');

figure; hold on; plot([0:100], HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'], 'r*');  
plot([0:100], (-[3:103]./[1:101]) .* ((1-2*p)/p) .*    (( 1 + (p/(1-p)) .^[1:101]  )  ./ (( 1 + (p/(1-p)) .^[0:100]  )))  , 'm'); 
% plot([0:100], (-[3:103]./[1:101]) .* ((1-2*p)/p), 'm'); 
title(['epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');

figure; hold on; plot(log([1:101]), log(abs(HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'] - ((1-2*p)/p) )), 'r*'); 
title(['epsilon Log Ratio of Orders for p = ' num2str(p)]); xlabel('Log Order'); ylabel('Log Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');

figure; hold on; plot([0:100], log(abs(HpConv(1:101)' ./ [HpConv_zero HpConv(1:100)'] - ((1-2*p)/p) )), 'r*'); 
title(['epsilon Log Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Ratio of order to prev. order'); % legend('I.I.D.', 'HMM');

% figure; hold on; plot([0:11], log(abs([HpConv_zero HpConv(1:11)'])), '*r'); plot([0:11], log(abs([Hp_zero Hp(1:11)'])), '+');  
% title(['epsilon Log Abs Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Log Coefficient Value'); legend('I.I.D.', 'HMM');
% 
% figure; hold on; plot([0:10], HpConv(1:11)' ./ [HpConv_zero HpConv(1:10)'], 'r*'); plot([0:10], Hp(1:11)' ./ [Hp_zero Hp(1:10)'], '+'); 
% title(['epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); legend('I.I.D.', 'HMM');
% 
% figure; hold on; plot([0:10], log(HpConv(1:11)' ./ [HpConv_zero HpConv(1:10)']), 'r*'); plot([0:10], log(Hp(1:11)' ./ [Hp_zero Hp(1:10)']), '+'); 
% title(['Log epsilon Ratio of Orders for p = ' num2str(p)]); xlabel('Order'); ylabel('Ratio of order to prev. order'); legend('I.I.D.', 'HMM');

for r=1:res
    eps = eps0 + (eps_jump*(r-1)/res);
    eps_vec(r) = eps;
    P = [1-p p; p 1-p];
    EPS = [1-eps eps; eps 1-eps];
    
    sam_ent = 0;
    
    Y_Vec = sample_HMM_vec(P,EPS,N,iters);
    
    
    condprobs_vec = compute_HMM_condprob_vec(P, EPS, Y_Vec);
    
    EEEEEEE = entropy  ([condprobs_vec 1-condprobs_vec]');
    
    sam_ent = sum(  entropy  ([condprobs_vec 1-condprobs_vec]')) / iters;
    
    
% % % % % % % %     for i=1:iters  % randomize and calculate entropy
% % % % % % % %         
% % % % % % % %    %%%     Y = sample_HMM(P, EPS, N); % sample the sequence
% % % % % % % %         
% % % % % % % %         Y = Y_Vec(i,:);
% % % % % % % %         
% % % % % % % %         condp = compute_HMM_condprob(P, EPS, Y);
% % % % % % % %         
% % % % % % % %         sam_ent = sam_ent + entropy([condp 1-condp]');
% % % % % % % %         
% % % % % % % %     end
% % % % % % % %     
% % % % % % % % % % %     r
% % % % % % % % % % %     cputime-ttt
% % % % % % % %     sam_ent = sam_ent / iters;
    
    
%     (sam_ent - H(p)) / eps
%     2*(2*p-1) * log2(p/(1-p))
    
    
    first_order(r) = (sam_ent - H(p)) / eps;
    
    if(mod(r,10) == 0)
        cur_iter = r
    end
   
    did_res = r
end


cputime - ttt

% Comparison with what should be the correct answer
eytan_val = repmat(2*(2*p-1) * log2(p/(1-p)), 1, res);


RRR = entropy([ p+eps_vec - 2*p* eps_vec 1-(p+eps_vec - 2*p*eps_vec) ] )
EEE = entropy([ p+eps_vec(1) - 2*p* eps_vec(1) 1-(p+eps_vec(1) - 2*p*eps_vec(1)) ]' );
PPP = entropy([p 1-p]');
DDD = (entropy([ p+eps_vec(1) - 2*p* eps_vec(1) 1-(p+eps_vec(1) - 2*p*eps_vec(1)) ]' ) - entropy([p 1-p]')) / eps_vec(1);


new_eytan_val_half =  (entropy([ (p+eps_vec - 2*p* eps_vec)' (1-(p+eps_vec - 2*p*eps_vec))' ]') - entropy([p 1-p]')) ./ eps_vec

new_eps_vec = 2*eps_vec - 2*eps_vec .* eps_vec;

new_eytan_val =  (entropy([ (p+new_eps_vec - 2*p* new_eps_vec)' (1-(p+new_eps_vec - 2*p*new_eps_vec))' ]') - entropy([p 1-p]')) ./ eps_vec




% Now insert also the 2nd term from taylor expansion
eytan_val_2nd_term = eytan_val -   (-0.5   + 4*(1-2*p)*log((1-p)/p) + 0.5*(  (p*p+(1-p)*(1-p))/(p*(1-p))  )^2 )*eps_vec;

% Now insert also the 2nd term from taylor expansion
eytan_val_2nd_term_corrected = eytan_val -   (-0.5   + 4*(1-2*p)*log((1-p)/p) + 0.5*(  1/(p*(1-p)) - 3  )^2 )*eps_vec;

Ido_2nd_term = eytan_val - (  2*(1-2*p)*log((1-p)/p)  + 0.5 * (1-2*p)*(1-2*p) /   (p*p*(1-p)*(1-p))  )*eps_vec;

figure; hold on; plot(eps_vec, first_order, '.');  plot(eps_vec, eytan_val, 'r');  %plot(eps_vec, eytan_val/2, 'g'); 
 plot(eps_vec, eytan_val_2nd_term, 'm'); plot(eps_vec, eytan_val_2nd_term_corrected, 'k'); plot(eps_vec, Ido_2nd_term, 'g');

 %%%plot(eps_vec, new_eytan_val, 'm');  plot(eps_vec, new_eytan_val_half, 'k'); 
% % % legend('points', '1st order', 'half 1st order', 'double_noise', 'one_noise');
legend('points', '1st order',  '2nd order', '2nd order corrected', '2nd order Ido');
title(['1st order coef. of epsilon for p='  num2str(p) ' N = ' num2str(N)]); xlabel('epsilon'); ylabel('1st order');


    % Randomize Y according to the probabilities
    % Now we have to calculate for a hidden markov process.  X_n --> Y_n
    % This is done via viterby algorithm
    % We need to calculate H(Y_n | Y_{n-1}, .., Y_1), or H(Y_n | Y_{n-1}, .., Y_1, X_1)
    %
    % For this, we will calculate the probability of Pr(Y_n = 1, Y_{n-1} = ?,
    % .. Y_1 = ? ) and Pr(Y_n = 1, Y_{n-1} = ?, .. Y_1 = ?, X_1 = ?)
    
