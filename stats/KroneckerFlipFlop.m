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
%    P = Z'*inv(prev_Q)*Z ./ q;
%    Q = Z*inv(prev_P)*Z' ./ p;

    Q = KroneckerFlip_internal(Z, prev_P, options);  % Q = Q .* (norm(prev_Q, 'fro') / norm(Q, 'fro')); % Normalize. Should accelerate convergence 
    P = KroneckerFlop_internal(Z, prev_Q, options); 
% %     if(isfield(options, 'regularized'))
% %        diff_mat_P = 10 * epsilon;
% %        P2 = Z'*inv(prev_Q + options.lambda.*eye(q))*Z ./ q;
% %        while (diff_mat_P > epsilon)
% %           P =  Z'*inv(prev_Q)*Z ./ q - ((4*options.lambda)/q)*inv(P2); 
% %           P = fsolve(@(x) x^2 - (4*options.lambda/q)*eye(p) - x*(Z'*inv(prev_Q)*Z ./ q), P2)
% %           P = fminsearch(@(x) norm(x^2 - (4*options.lambda/q)*eye(p) - x*(Z'*inv(prev_Q)*Z ./ q), 'fro'), P2)
% %           diff_mat_P = norm(P - P2, 2) / norm(P2, 2)
% %           P2=P; 
% %           
% %        end
% %     end
% %     

    diff_mat = norm(P - prev_P, 2) / norm(prev_P, 2) + ...
        norm(Q - prev_Q, 2) / norm(prev_Q, 2);
    prev_P = P; prev_Q = Q; 
end
    
% Internal function for performing one iteration of flip-flop 
function new_Q = KroneckerFlip_internal(Z, P, options)

[q, p] = size(Z); r = rank(P); 

if(isfield(options, 'regularized')) % solve regularized version
   ZZ = Z*inv(P)*Z'; ZZ = 0.5*(ZZ+ZZ'); 
   [W, Gamma] = eig(ZZ); % perform eigendecomposition       
   gamma = abs(diag(Gamma));  gamma(1:(q-r) ) = 0; % take eigenvalues
   delta = (gamma + sqrt( gamma.^2 +4*p*options.lambda)) ./ p;    
   new_Q = W*diag(delta)*W'; % get new Q 
end


% Internal function for performing one iteration of flip-flop 
function new_P = KroneckerFlop_internal(Z, Q, options)

[q, p] = size(Z); r = rank(Q); 
p1 = options.p1; % p2 = options.p2; 

% compute under H0
if(ismember(options.h , {'0', 'h0'})) % under H0 we fit both parts seperately
    options.h = 'h1';
    new_P = zeros(p); 
    new_P(1:p1,1:p1) = ...
        KroneckerFlop_internal(Z(:,1:p1), Q, options);
    new_P((p1+1):end, (p1+1):end) = ...
        KroneckerFlop_internal(Z(:,(p1+1):end), Q, options);
    options.h = 'h0';
else    % use standard fitting under H1
    if(isfield(options, 'regularized')) % solve regularized version
        ZZ = Z'*inv(Q)*Z; ZZ = 0.5*(ZZ+ZZ'); 
        [W, Gamma] = eig(ZZ); % perform eigendecomposition
        gamma = abs(diag(Gamma));  gamma(1:(q-r) ) = 0; % take eigenvalues
        delta = (gamma + sqrt( gamma.^2 +4*p*options.lambda)) ./ q;    
%        delta = (diag(Gamma) + sqrt( diag(Gamma).^2 +4*q*options.lambda)) ./ q;
        new_P = W*diag(delta)*W'; % get new Q
    end
    
end
