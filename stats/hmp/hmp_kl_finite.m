% Compute Exactly KL_N for a given N, epsilon and p (first model), and
% delta and q (second model). p,q are vectors.
% The method is like the maple one (which was symbollicaly). Enumerate all
% 2^N configurations and sum the P * log (P/Q)
% Expansion. 
function KL_N = HMP_KL_finite(eps, p_vec, delta, q_vec, N)

% We assume initial condition : X_1 = 1
p_len = length(p_vec); q_len = length(q_vec);
KL_N = zeros(p_len, 1);

P_Y = zeros(p_len, 2^N); % Make an array of Y probabilities given p.epsilon
Q_Y = zeros(q_len, 2^N); % Make an array of Y probabilities given q,delta

for X = 0:2^N-1  % Go over all the X's 
    neq_bonds = sum(bitget(bitxor(X, bitshift(X,-1)), 1:N-1));
    P_X = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)'; 
    Q_X = 0.5*((1-q_vec) .^ (N-1-neq_bonds) .* q_vec .^ neq_bonds)';    
    
    for Y = 0:2^(N-1)-1  % Go over all the Y's. We can save half here due to symmetry !!!!!
        neq_noise_bonds = sum(bitget(bitxor(X,Y), 1:N));
        P_Y(:,Y+1) = (1-eps) .^ (N-neq_noise_bonds) .* eps .^ neq_noise_bonds .* P_X + P_Y(:,Y+1);
        Q_Y(:,Y+1) = (1-delta) .^ (N-neq_noise_bonds) .* delta .^ neq_noise_bonds .* Q_X + Q_Y(:,Y+1);
    end
   
end

% Now compute the relative entropy 
for Y = 0:2^(N-1)-1
    KL_N = KL_N + P_Y(:,Y+1) .* ( log(P_Y(:,Y+1)) - log(Q_Y(:,Y+1)));
end

% Transfer to binary logarithm
KL_N = 2 .* KL_N ./ log(2.0);
    

    
    
    
    
% %  Maple needed operations    
% % > single_ent := proc(x) return -x*log(x); end proc;bin_ent := proc(x) return -x*log(x) - (1-x)*log(1-x); end proc;
% % > conv := proc(x,y); return x*(1-y) + y*(1-x); end proc;
% % > inv_conv := proc(x,y); return 1-x-y+2*x*y; end proc;
% % > Bit := proc(x, k); return floor(x/2^k) mod 2; end proc;
% % > Xor := proc(x,y); return (x+y) mod 2; end proc;




