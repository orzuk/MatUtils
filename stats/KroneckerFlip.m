% Internal function for performing one iteration of flip-flop 
function new_Q = KroneckerFlip(Z, P, options)

[q, p] = size(Z); r = rank(P); 

if(isfield(options, 'regularized')) % solve regularized version
   ZZ = Z*inv(P)*Z'; ZZ = 0.5*(ZZ+ZZ'); 
   [W, Gamma] = eig(ZZ); % perform eigendecomposition       
   gamma = abs(diag(Gamma));  gamma(1:(q-r) ) = 0; % take eigenvalues
   delta = (gamma + sqrt( gamma.^2 +4*p*options.lambda)) ./ p;    
   new_Q = W*diag(delta)*W'; % get new Q 
end
