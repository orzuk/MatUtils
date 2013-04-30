% Compute HMM entropy

p = 0.3;
eps = 0.1;

segs = 100;



p_vec = [0:0.5/(segs):0.5]';
eps_vec = [0:0.5/(segs):0.5]';
% p_vec = [0:0.005/(segs):0.005]';
% eps_vec = [0:0.005/(segs):0.005]';


tolerance = 0.000000000001;
mul_ent = zeros(segs+1,segs+1);
upper_ent = zeros(segs+1,segs+1);
sum_upper_ent = zeros(segs+1,segs+1);
ent_sum = zeros(segs+1,segs+1);
HP = zeros(segs+1,segs+1);
Heps = zeros(segs+1,segs+1);
Ido_ent = zeros(segs+1,segs+1);
Ido_ent_better_phi = zeros(segs+1,segs+1); % Here we estimate phi not based on intersection of two markov chains,
                                           % but of many.

% Try the new bounds and the 2nd order approximation
Ido_upper = zeros(segs+1,segs+1);
Ido_lower = zeros(segs+1,segs+1);
Eytan_2nd_order = zeros(segs+1,segs+1);
Eytan_2nd_order_corrected = zeros(segs+1,segs+1);

% Do small diriclet correction to avoid zero probabilities.
dirich = 0.000000000000001;
p_vec = (p_vec + dirich) / (1+2*dirich);

M=7;

for i=1:segs+1  % loop over p/eps
    i
    mul_ent(i,:) = entropy( [p_vec(i) + eps_vec - 2*p_vec(i)*eps_vec, 1 - (p_vec(i) + eps_vec - 2*p_vec(i)*eps_vec)]' );
    temp = p_vec(i) + 2*eps_vec - 4*p_vec(i)*eps_vec + 4*p_vec(i)*(eps_vec .^ 2) - 2* (eps_vec .^ 2);
    if(~isempty(find(temp > 1)))
        eps_vec
        ppp = p_vec(i)
        temp
    end
    
%     wrong = temp(find(temp > 1))
%     
%     
%     wrrrrong = find(temp < 0)
    upper_ent(i,:) = entropy( [temp, 1-temp]' );
    sum_upper_ent(i,:) = entropy([p_vec(i), 1-p_vec(i)]') + entropy([eps_vec, 1 - eps_vec]');  % A very naive sum upper-bound H(p) + H(eps)
    ent_sum(i,:) = entropy([ p_vec(i) + eps_vec, 1 - p_vec(i) - eps_vec]');    % Entropy,  H(p+eps)
    HP(:,i) = (entropy([ p_vec , 1 - p_vec ]'))';  % H(p)
    Heps(i,:) = entropy(  [ eps_vec , 1 - eps_vec ]' );  % H(eps)
    
    p = min(max(p_vec(i), tolerance), 1-tolerance);
    
    Ido_ent(i,:) = Heps(i,:) + 1 - ( log2( p^2 + (1-p)^2 + ...
        sqrt(       (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8.*eps_vec.*(1-eps_vec).*(eps_vec .^ 2 + (1-eps_vec).^2)  )   ) )';
    
%     (p^2 + (1-p)^2)^2 -8.*eps_vec.*(1-eps_vec).*(eps_vec .^ 2 + (1-eps_vec).^2) .* (p^2 - (1-p^2))^2   
% % %      (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8.*eps_vec.*(1-eps_vec).*(eps_vec .^ 2 + (1-eps_vec).^2)
% % %      for j=1:segs+1
% %          (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps_vec(j)*(1-eps_vec(j))*(eps_vec(j) ^ 2 + (1-eps_vec(j))^2)
% %      end
% % %     lambda_1 = 0.5* ( (p^2 + (1-p)^2) + sqrt ( (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)) );
% % %     lambda_2 = 0.5* ( (p^2 + (1-p)^2) - sqrt ( (p^2 + (1-p)^2)^2  - (p^2 - (1-p)^2)^2 * 8*eps*(1-eps)*(eps ^ 2 + (1-eps)^2)) );

    
    
%     badbad = (p^2 + (1-p)^2)^2 -8.*eps_vec.*(1-eps_vec).*(eps_vec .^ 2 + (1-eps_vec).^2) .* (p^2 - (1-p^2))^2 
    for j=1:segs+1
        Ido_ent_better_phi(i,j) = Heps(i,j) - Compute_Matrix_Many_Chains_Intersection(p, eps_vec(j), M)/(M-1);
    end
    % New bounds plus 2nd order 
% % %     Eytan_2nd_order(:,i) = HP(:,i) + eps_vec(i) * 2*(2*p_vec-1) .* log2( p_vec./(1-p_vec) ) + ... % 0 and 1st order 
% % %        eps_vec(i)^2 * ( 0.5 - 4*(2*p_vec-1).*log2(p_vec./(1-p_vec)) - 0.5 * ((p_vec.^2+(1-p_vec).^2)./(p_vec.*(1-p_vec))).^2 );  % 2nd order
        
   
    Eytan_2nd_order_corrected(:,i) = HP(:,i) + eps_vec(i) * 2*(2*p_vec-1) .* log2( p_vec./(1-p_vec) ) + ... % 0 and 1st order 
       eps_vec(i)^2 * ( 0.5/log(2) - 4*(2*p_vec-1).*log2(p_vec./(1-p_vec)) - (0.5/log(2)) * ( 1./(p_vec .* (1-p_vec)) - 3 ).^2 );  % 2nd order corrected

   % Do the upper bound
   IDO_PA = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.*(1-p_vec)*(eps_vec(i)^3 + (1-eps_vec(i))^3 + eps_vec(i)*(1-eps_vec(i))) + p_vec.^2*eps_vec(i)*(1-eps_vec(i));
   IDO_PB = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + 2*p_vec.*(1-p_vec)*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*(eps_vec(i)^3 + (1-eps_vec(i))^3);
   IDO_PC = (1-p_vec).^2*(eps_vec(i)^3 + (1-eps_vec(i))^3) + 2*p_vec.*(1-p_vec)*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*eps_vec(i)*(1-eps_vec(i));
   Ido_upper(:,i) = -( IDO_PA.*log2(IDO_PA) + 0.5*IDO_PB.*log2(IDO_PB) + 0.5*IDO_PC.*log2(IDO_PC));
   
   
   % Now do the lower bound
   IDO_PA = (1-p_vec).^2*(1-eps_vec(i))^2 + p_vec.^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.*(1-p_vec)*eps_vec(i);
   IDO_PB = (1-p_vec).^2*eps_vec(i)^2 + p_vec.^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.*(1-p_vec)*(1-eps_vec(i));
   IDO_PC = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*(1-eps_vec(i))^2 + p_vec.*(1-p_vec)*eps_vec(i);
   IDO_PD = (1-p_vec).^2*eps_vec(i)*(1-eps_vec(i)) + p_vec.^2*eps_vec(i)^2 + p_vec.*(1-p_vec)*(1-eps_vec(i));
   Ido_lower(:,i) = -0.5*( IDO_PA.*log2(IDO_PA) + IDO_PB.*log2(IDO_PB) + IDO_PC.*log2(IDO_PC) + IDO_PD.*log2(IDO_PD));
   
% % % % %     H(p,f) = H(p) + 2(2p-1)log(p/(1-p)) * f +
% % % % %  [ 1/2 - 4(2p-1)log(p/(1-p)) - 1/2 ((p^2+(1-p)^2)/p(1-p))^2 ] * f^2 + O(f^3)

    
end


figure; subplot(3,3,1); hold on;  imagesc(mul_ent); colorbar('vert'); title('HMM Lower Bound H(Y_n/X_{n-1})'); ylabel('p'); xlabel('eps');
subplot(3,3,2); hold on; imagesc(upper_ent); colorbar('vert'); title('HMM Upper Bound H(Y_n/Y_{n-1})'); ylabel('p'); xlabel('eps');
subplot(3,3,3); hold on; imagesc(sum_upper_ent); colorbar('vert'); title('HMM SUM Bound H(p) + H(eps)'); ylabel('p'); xlabel('eps');
subplot(3,3,4); hold on; imagesc(upper_ent - mul_ent); colorbar('vert'); title('H(Y_n/Y_{n-1}) - H(Y_n/X_{n-1})'); ylabel('p'); xlabel('eps');
subplot(3,3,5); hold on; imagesc(sum_upper_ent - upper_ent); colorbar('vert'); title('H(p) + H(eps) - H(Y_n/Y_{n-1})'); ylabel('p'); xlabel('eps');
subplot(3,3,6); hold on; imagesc(sum_upper_ent - ent_sum); colorbar('vert'); title('H(p) + H(eps) - H(p+eps)'); ylabel('p'); xlabel('eps');
subplot(3,3,7); hold on; imagesc(ent_sum); colorbar('vert'); title('H(p+eps)'); ylabel('p'); xlabel('eps');
subplot(3,3,8); hold on; imagesc(sign(sum_upper_ent - upper_ent)); colorbar('vert'); title('Sign of H(p) + H(eps) - H(Y_n/Y_{n-1})'); ylabel('p'); xlabel('eps');
%%subplot(3,3,9); hold on; imagesc(sign(sum_upper_ent - ent_sum)); colorbar('vert'); title('H(p) + H(eps) - H(p+eps) Sign '); ylabel('p (X flip probability)'); xlabel('eps (noise probability)');
subplot(3,3,9); hold on; imagesc(upper_ent-ent_sum); colorbar('vert'); title('H(Y_n/Y_{n-1}) - H(p+eps)'); ylabel('p'); xlabel('eps');



Eytan_2nd_order = min(Eytan_2nd_order,1); % avoid entropies above one !
Eytan_2nd_order = max(Eytan_2nd_order,0); % avoid entropies above one !
Eytan_2nd_order_corrected = min(Eytan_2nd_order_corrected,1); % avoid entropies above one !
Eytan_2nd_order_corrected = max(Eytan_2nd_order_corrected,0); % avoid entropies above one !


figure; subplot(3,3,1); hold on; 
% % % imagesc((Eytan_2nd_order)); colorbar('vert'); title('HMM 2nd order approximation'); ylabel('p'); xlabel('eps');
subplot(3,3,1); hold on; imagesc((Eytan_2nd_order_corrected)); colorbar('vert'); title('HMM 2nd order approximation Corrected'); ylabel('p'); xlabel('eps');
subplot(3,3,2); hold on; imagesc((Ido_upper)); colorbar('vert'); title('HMM Ido 2-step Upperbound'); ylabel('p'); xlabel('eps');
subplot(3,3,4); hold on; imagesc((Ido_upper-Eytan_2nd_order_corrected)); colorbar('vert'); title('Upperbound Minus 2nd order'); ylabel('p'); xlabel('eps');
subplot(3,3,3); hold on; imagesc((Ido_lower)); colorbar('vert'); title('HMM Ido 2-step Lowerbound'); ylabel('p'); xlabel('eps');
subplot(3,3,5); hold on; imagesc((Ido_upper-Ido_lower)); colorbar('vert'); title('Upperbound Minus Lowerbound'); ylabel('p'); xlabel('eps');
subplot(3,3,6); hold on; imagesc((Eytan_2nd_order_corrected-Ido_lower)); colorbar('vert'); title('2nd order Minus Lowerbound'); ylabel('p'); xlabel('eps');
% subplot(3,3,7); hold on; imagesc(sign(Ido_upper-Eytan_2nd_order_corrected)); colorbar('vert'); title('SIGN Upperbound Minus 2nd order'); ylabel('p'); xlabel('eps');
subplot(3,3,7); hold on; imagesc(Ido_ent); colorbar('vert'); title('Ido ent two-markovs phi'); ylabel('p'); xlabel('eps');
subplot(3,3,8); hold on; imagesc(2*(Ido_upper-Ido_lower)./(Ido_upper+Ido_lower)); colorbar('vert'); title('Upperbound Minus Lowerbound Relative Error'); ylabel('p'); xlabel('eps');
% subplot(3,3,9); hold on; imagesc(sign(Eytan_2nd_order_corrected-Ido_lower)); colorbar('vert'); title('SIGN 2nd order Minus Lowerbound'); ylabel('p'); xlabel('eps');
subplot(3,3,9); hold on; imagesc(Ido_ent_better_phi); colorbar('vert'); title('Ido ent with good phi'); ylabel('p'); xlabel('eps');




% % % % % % % % % figure; subplot(3, 2, 1); hold on; imagesc(sign(sum_upper_ent - upper_ent)); colorbar('vert'); title('HMM Upper Bounds Difference Sign'); ylabel('p'); xlabel('eps');
% % % % % % % % % subplot(3, 2, 2); hold on; imagesc( min(10, log(Heps) -log(max(ent_sum - HP, 0.0000000000000000001)))); colorbar('vert'); title('log H(eps) / (H(p+eps) - H(p))'); ylabel('p'); xlabel('eps');
% % % % % % % % % subplot(3, 2, 3); hold on; imagesc( min(10, log(Heps))); colorbar('vert'); title('log H(eps) '); ylabel('p'); xlabel('eps');
% % % % % % % % % subplot(3, 2, 4); hold on; imagesc( min(10, ent_sum - HP)); colorbar('vert'); title(' (H(p+eps) - H(p))'); ylabel('p'); xlabel('eps');
% % % % % % % % % subplot(3, 2, 5); hold on; imagesc( Ido_ent); colorbar('vert'); title(' Ido approx.'); ylabel('p'); xlabel('eps');
% % % % % % % % % %subplot(2, 2, 2); hold on; imagesc( HP); colorbar('vert'); title('H(p) direct'); ylabel('p'); xlabel('eps');
% % % % % % % % % %subplot(2, 2, 3); hold on; imagesc( Heps); colorbar('vert'); title('H(eps) direct'); ylabel('p'); xlabel('eps');
% % % % % % % % % %subplot(2, 2, 4); hold on; imagesc(sum_upper_ent - HP); colorbar('vert'); title('H(eps) diff'); ylabel('p'); xlabel('eps');



nsamples = 100; sample_len = 5000; 

P = [0.8, 0.2; 0.4, 0.6];
eps = [0.999 0.001; 0.001 0.999];
Y = (rand(1, 10) > 0.5) + 1;


% % % P = [1 0; 1 0];
% % % eps = [1 0; 0 1];


% Trying to guess the entropy
[V, D] = eig(P');

[val ind] = max(diag(D));

mew = V(:,ind) ./ sum(V(:,ind))  % The stationary distribution 

%%% guess_mul_ent = entropy([ mew(1) + eps(1,2) - 2 * mew(1) * eps(1,2), 1
%%% - mew(1) - eps(1,2) + 2 * mew(1) * eps(1,2) ]' )  %%% BADDDDD


% Give the entropy given X_{n-1}
Q = P * eps;

guess_mul_ent = mew(1) * entropy([Q(1,1), Q(1,2)]') + mew(2) * entropy([Q(2,1), Q(2,2)]')



Y = [2 2 2 1 1 1 2 1 2 2 2 1 1 1 2 2 2 2 2 1 1 2 1 2 2 1 2 1 2  1 1 2 2 2 1 1 2];

condp = compute_HMM_condprob(P, eps, Y);


sam_ent = 0;

for i=1:nsamples
   
    % Randomize Y according to the probabilities
    Y = sample_HMM(P, eps, sample_len);
    
    condp = compute_HMM_condprob(P, eps, Y);
    % Now we have to calculate for a hidden markov process.  X_n --> Y_n
    % This is done via viterby algorithm
    % We need to calculate H(Y_n | Y_{n-1}, .., Y_1), or H(Y_n | Y_{n-1}, .., Y_1, X_1)
    %
    % For this, we will calculate the probability of Pr(Y_n = 1, Y_{n-1} = ?,
    % .. Y_1 = ? ) and Pr(Y_n = 1, Y_{n-1} = ?, .. Y_1 = ?, X_1 = ?)
    
    
    sam_ent = sam_ent + entropy([condp 1-condp]');
    
end


sam_ent = sam_ent / nsamples