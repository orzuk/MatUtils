% Compute KL(X1_N) for a given N, conditioned on X1.
% Everyhing is for given epsilon and p (first model), and
% delta and q (second model). p, q are vectors. Here everything is conditional 
% on X1 which is assumed to be known.
% The method is like the maple one (which was symbollicaly). Enumerate all
% 2^N configurations and sum the P * log (P/Q).
%
function KLX1_N = hmp_kl_finite_x1(eps, p_vec, delta, q_vec, N)


% We assume initial condition : X_1 = 1
p_len = length(p_vec);
q_len = length(q_vec);

KLX1_N = zeros(p_len, 1);

P_X1Y = zeros(p_len, 2^(N+1)); % Make an array of Y probabilities according to p,epsilon. Here we take also X1 !!!
Q_X1Y = zeros(p_len, 2^(N+1)); % Make an array of Y probabilities according to q,delta. Here we take also X1 !!!

for X1 = 0:0  % loop over the first X. From symmetry we can assume X=0 !!!
    for X = 0:2^(N-1)-1  % Go over all the other X's 
        XX = X + X1 * 2^(N-1);  % Make the whole X vector
        neq_bonds = sum(bitget(bitxor(XX, bitshift(XX,-1)), 1:N-1));
        P_XX = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)';         
        Q_XX = 0.5*((1-q_vec) .^ (N-1-neq_bonds) .* q_vec .^ neq_bonds)'; 

        
        for Y = 0:2^N-1  % Go over all the Y's
            neq_noise_bonds = sum(bitget(bitxor(XX,Y), 1:N));
            P_X1Y(:,X1*2^N+Y+1) = (1-eps) .^ (N-neq_noise_bonds) .* eps .^ neq_noise_bonds .* P_XX + P_X1Y(:,X1*2^N+Y+1);
            Q_X1Y(:,X1*2^N+Y+1) = (1-delta) .^ (N-neq_noise_bonds) .* delta .^ neq_noise_bonds .* Q_XX + Q_X1Y(:,X1*2^N+Y+1);
        end

    end
end

% Now compute the entropy 
for X1Y = 0:2^N-1
    KLX1_N = KLX1_N + P_X1Y(:,X1Y+1) .* ( log(P_X1Y(:,X1Y+1)) - log(Q_X1Y(:,X1Y+1)) );
end

% Transfer to binary logarithm
KLX1_N = 2 .* KLX1_N ./ log(2.0);
    

    
    
    
    
% %  Maple needed operations    
% % > single_ent := proc(x) return -x*log(x); end proc;bin_ent := proc(x) return -x*log(x) - (1-x)*log(1-x); end proc;
% % > conv := proc(x,y); return x*(1-y) + y*(1-x); end proc;
% % > inv_conv := proc(x,y); return 1-x-y+2*x*y; end proc;
% % > Bit := proc(x, k); return floor(x/2^k) mod 2; end proc;
% % > Xor := proc(x,y); return (x+y) mod 2; end proc;




