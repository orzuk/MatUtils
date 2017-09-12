% Internal function for performing one iteration of flip-flop 
function new_P = KroneckerFlop(Z, Q, options)

[q, p] = size(Z); r = rank(Q); 
p1 = options.p1; % p2 = options.p2; 

% compute under H0. Here P has two blocks 
if(ismember(options.h , {'0', 'h0'})) % under H0 we fit both parts seperately
    options.h = 'h1';
    new_P = zeros(p); 
    new_P(1:p1,1:p1) = ...
        KroneckerFlop(Z(:,1:p1), Q, options);
    new_P((p1+1):end, (p1+1):end) = ...
        KroneckerFlop(Z(:,(p1+1):end), Q, options);
    options.h = 'h0';
else    % use standard fitting under H1
    if(isfield(options, 'regularized')) % solve regularized version
        ZZ = Z'*inv(Q)*Z; ZZ = 0.5*(ZZ+ZZ'); 
        [W, Gamma] = eig(ZZ); % perform eigendecomposition
        gamma = abs(diag(Gamma));  gamma(1:(q-r) ) = 0; % take eigenvalues
        delta = (gamma + sqrt( gamma.^2 +16*q*options.lambda)) ./ (2*q);    
%        delta = (diag(Gamma) + sqrt( diag(Gamma).^2 +4*q*options.lambda)) ./ q;
        new_P = W*diag(delta)*W'; % get new Q
    end
    
end
