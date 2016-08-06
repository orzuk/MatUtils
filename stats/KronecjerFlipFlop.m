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
    P = Z*inv(prev_Q)*Z ./ p;
    
    
    
    Q = Z'*inv(prev_P)*Z ./ q;

    
    
    
    diff_mat = norm(Q - prev_P, 2) / norm(prev_P, 2) + ...
        norm(Q - prev_Q, 2) / norm(prev_Q, 2);
    prev_P = P; prev_Q = Q; 
end
    
