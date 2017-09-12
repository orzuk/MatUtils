% Compute finite HMP entropy H_N for a given N, epsilon and p and X1.
% The method is like the maple one (which was symbollicaly). Enumerate all
% 2^N configurations and sum the p * log p
%
function HX1_N = hmp_entropy_finite_X1(eps, p_vec, N)


% We assume initial condition : X_1 = 1
p_len = length(p_vec);
HX1_N = zeros(p_len, 1);

P_X1Y = zeros(p_len, 2^(N+1)); % Make an array of Y probabilities. Here we take also X1 !!!

for X1 = 0:0  % loop over the first X. From symmetry we can assume X=0 !!!
    for X = 0:2^(N-1)-1  % Go over all the other X's 
        XX = X + X1 * 2^(N-1);  % Make the whole X vector
        neq_bonds = sum(bitget(bitxor(XX, bitshift(XX,-1)), 1:N-1));
        P_XX = 0.5*((1-p_vec) .^ (N-1-neq_bonds) .* p_vec .^ neq_bonds)'; 
        
        for Y = 0:2^N-1  % Go over all the Y's
            neq_noise_bonds = sum(bitget(bitxor(XX,Y), 1:N));
            P_X1Y(:,X1*2^N+Y+1) = (1-eps) .^ (N-neq_noise_bonds) .* eps .^ neq_noise_bonds .* P_XX + P_X1Y(:,X1*2^N+Y+1);
        end
        
    end
end

% Now compute the entropy 
for X1Y = 0:2^N-1
    HX1_N = HX1_N - P_X1Y(:,X1Y+1) .* log(P_X1Y(:,X1Y+1));
end

% Transfer to binary logarithm
HX1_N = 2 .* HX1_N ./ log(2.0);
    

    
    
    
    
% %  Maple needed operations    
% % > single_ent := proc(x) return -x*log(x); end proc;bin_ent := proc(x) return -x*log(x) - (1-x)*log(1-x); end proc;
% % > conv := proc(x,y); return x*(1-y) + y*(1-x); end proc;
% % > inv_conv := proc(x,y); return 1-x-y+2*x*y; end proc;
% % > Bit := proc(x, k); return floor(x/2^k) mod 2; end proc;
% % > Xor := proc(x,y); return (x+y) mod 2; end proc;




