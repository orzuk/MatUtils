% Run flip-flop algorithm with affine constraints
% 
% Input: 
% Z - data matrix 
% P0 - initial guess for P 
% Q0 - initial guess for Q
% epsilon - tolerance in stopping criteria 
% options - add restrictions on Q or P 
% 
% Output: 
% P - final value for P part
% Q - final value for Q part
% 
function [P, Q] = KroneckerFlipFlop(Z, P0, Q0, epsilon, options)

[q, p] = size(Z);
 
prev_P = P0; prev_Q = Q0; % Initialize

diff_mat = 10 * epsilon;
while(diff_mat > epsilon)
    P = Z*inv(prev_Q)*Z ./ q;
    
    if(isfeild('regularized', options))
       diff_mat_P = 10 * epsilon;
       while (diff_mat_P > epsilon)
          P =  Z*inv(prev_Q)*Z ./ q - ((4*options.lambda)/q)*inv(P); 
          diff_mat_P = norm(P - P2, 2) / norm(P2, 2); P2=P; 
          
       end
    end
    
    Q = Z'*inv(prev_P)*Z ./ p;

    
    
    
    diff_mat = norm(P - prev_P, 2) / norm(prev_P, 2) + ...
        norm(Q - prev_Q, 2) / norm(prev_Q, 2);
    prev_P = P; prev_Q = Q; 
end
    
