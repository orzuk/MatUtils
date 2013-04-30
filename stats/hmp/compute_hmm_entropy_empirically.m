% Compute HMP entropy by sampling and computing the empirical entropy
% p is the transition probability matrix 2*2, eps is the probability of
% flipping the output, Y_vec is the vector of values Y_1 .. Y_N. We assume
% that the markov chain is stationary and binary

p = 0.22; N = 80; iters = 4000;
res = 100;

ttt = cputime;
eps_vec = zeros(1,res);
first_order = zeros(1,res);
ent_diff = zeros(1,res);

for r=1:res
    eps = 0.000000001 + (0.02*(r-1)/res);
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
    
    
    first_order(r) = (sam_ent - entropy([p, 1-p]')) / eps;
    ent_diff(r) = (sam_ent - entropy([p, 1-p]'));
    
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
    

    
    
% New plot :
figure; hold on; plot(eps_vec, ent_diff, '.');  plot(eps_vec, 2 * (1-2*p) * log2((1-p)/p) .* eps_vec, 'r');
legend('simulated', 'analytic first-order');


