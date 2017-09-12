% Compute HMM entropy when x is fixed 

p = 0.3;
eps = 0.1;

segs = 500;



p_vec = [0:0.5/(segs):0.5]';
eps_vec = [0:0.5/(segs):0.5]';
% p_vec = [0:0.005/(segs):0.005]';
% eps_vec = [0:0.005/(segs):0.005]';



mul_ent = zeros(segs+1,segs+1);
upper_ent = zeros(segs+1,segs+1);
sum_upper_ent = zeros(segs+1,segs+1);
ent_sum = zeros(segs+1,segs+1);
HP = zeros(segs+1,segs+1);
Heps = zeros(segs+1,segs+1);
Ido_ent = zeros(segs+1,segs+1);

IID_X_eq_given_Y_eq = zeros(segs+1,segs+1);
HMM_X_eq_given_Y_eq = zeros(segs+1,segs+1);



for i=1:segs+1  % loop over p/eps
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
    
    
    p = p_vec(i);
    Ido_ent(i,:) = Heps(i,:) + 1 - ( log2( p^2 + (1-p)^2 + ...
        sqrt(    (p^2 + (1-p)^2)^2 -8.*eps_vec.*(1-eps_vec).*(eps_vec .^ 2 + (1-eps_vec).^2) .* (p^2 - (1-p^2))^2   )   ) )'; 
    
    
    IID_X_eq_given_Y_eq(i,:) = comp_prob_X_equal_given_Y_equal(p_vec(i), eps_vec, 0)';
    HMM_X_eq_given_Y_eq(i,:) = comp_prob_X_equal_given_Y_equal(p_vec(i), eps_vec, 1)';
end

figure; subplot(2,2,1); hold on; imagesc(IID_X_eq_given_Y_eq);colorbar('vert'); title('IID prob. X equal given Y equal'); ylabel('p'); xlabel('eps');
subplot(2,2,2); hold on; imagesc(HMM_X_eq_given_Y_eq);colorbar('vert'); title('HMM prob. X equal given Y equal'); ylabel('p'); xlabel('eps');
subplot(2,2,3); hold on; imagesc(IID_X_eq_given_Y_eq-HMM_X_eq_given_Y_eq);colorbar('vert'); title('Difference : IID prob. MINUS HMM prob. X equal given Y equal'); ylabel('p'); xlabel('eps');
subplot(2,2,4); hold on; imagesc(IID_X_eq_given_Y_eq./HMM_X_eq_given_Y_eq);colorbar('vert'); title('Fraction : IID prob. DIVIDED BY HMM prob. X equal given Y equal'); ylabel('p'); xlabel('eps');


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


figure; subplot(3, 2, 1); hold on; imagesc(sign(sum_upper_ent - upper_ent)); colorbar('vert'); title('HMM Upper Bounds Difference Sign'); ylabel('p'); xlabel('eps');
subplot(3, 2, 2); hold on; imagesc( min(10, log(Heps) -log(max(ent_sum - HP, 0.0000000000000000001)))); colorbar('vert'); title('log H(eps) / (H(p+eps) - H(p))'); ylabel('p'); xlabel('eps');
subplot(3, 2, 3); hold on; imagesc( min(10, log(Heps))); colorbar('vert'); title('log H(eps) '); ylabel('p'); xlabel('eps');
subplot(3, 2, 4); hold on; imagesc( min(10, ent_sum - HP)); colorbar('vert'); title(' (H(p+eps) - H(p))'); ylabel('p'); xlabel('eps');
subplot(3, 2, 5); hold on; imagesc( Ido_ent); colorbar('vert'); title(' Ido approx.'); ylabel('p'); xlabel('eps');
%subplot(2, 2, 2); hold on; imagesc( HP); colorbar('vert'); title('H(p) direct'); ylabel('p'); xlabel('eps');
%subplot(2, 2, 3); hold on; imagesc( Heps); colorbar('vert'); title('H(eps) direct'); ylabel('p'); xlabel('eps');
%subplot(2, 2, 4); hold on; imagesc(sum_upper_ent - HP); colorbar('vert'); title('H(eps) diff'); ylabel('p'); xlabel('eps');
nsamples = 100; sample_len = 5000; 

P = [0.8, 0.2; 0.4, 0.6];
eps = [0.999 0.001; 0.001 0.999];
Y = (rand(1, 10) > 0.5) + 1;


% % % P = [1 0; 1 0];
% % % eps = [1 0; 0 1];


% Trying to guess the entropy
[V, D] = eig(P');

[val ind] = max(diag(D));

mu = V(:,ind) ./ sum(V(:,ind))  % The stationary distribution 

%%% guess_mul_ent = entropy([ mu(1) + eps(1,2) - 2 * mu(1) * eps(1,2), 1
%%% - mu(1) - eps(1,2) + 2 * mu(1) * eps(1,2) ]' )  %%% BADDDDD


% Give the entropy given X_{n-1}
Q = P * eps;

guess_mul_ent = mu(1) * entropy([Q(1,1), Q(1,2)]') + mu(2) * entropy([Q(2,1), Q(2,2)]')



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