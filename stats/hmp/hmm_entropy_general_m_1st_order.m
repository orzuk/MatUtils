% Compute HMM entropy for general alphabet

p = 0.3;
eps = 0.1;

segs = 100;


M = [1-p p; p 1-p];

Noise = [1-eps, eps; eps, 1-eps];
delta = [-1 1; 1 -1];  % The noise matrix difference

PI = [0.5, 0.5];   % The stationary vector

Ones = [1, 1]';   % Ones vector to sum results

pi_mat = [0.5 0; 0 0.5];
log_pi_mat = [log(0.5) 0; 0 log(0.5)];


D = delta'*pi_mat*M + pi_mat*M*delta

D + D*log(pi_mat*M) - log_pi_mat*D - pi_mat*delta*M
D + D*log(pi_mat*M) - log_pi_mat*D - pi_mat*delta*M
new_order=-Ones'*(D + D.*log(pi_mat*M))*Ones  - (log(Ones'*pi_mat)*D +  ((Ones'*pi_mat*delta)./(Ones'*pi_mat))*pi_mat*M)*Ones
new_order_2=-Ones'*(D + D.*log(pi_mat*M))*Ones  - (log(Ones'*pi_mat)*D +  Ones'*pi_mat*delta*M)*Ones

%first_order = PI * (delta * M); %% + M * delta'); %+ (2*delta*M+M*delta') .* log(M)) * Ones
%first_order = -PI * ((delta*M + M*delta') + (2*delta*M+M*delta') .* log(M)) * Ones
first_order_old_luck = -PI * ((delta*M + M*delta') + (delta*M+M*delta') .* log(M) + delta * (M .* log(M))) * Ones


% Now try to compute the total entropy
B = Noise'*pi_mat*M*Noise

exact_ent = Ones'*(B .* log(B))*Ones  - log(Ones'*pi_mat*Noise)*B*Ones 
exact_ent = sum(sum(B .* log(B)))  - sum(log(PI*Noise)*B)

%exact_ent = - sum(sum(log(pi_mat*Noise)*B))

temp = p + 2*eps - 4*p*eps + 4*p*eps^2 - 2*eps^2;
exact_ent_with_p = entropy( [temp, 1-temp]' )




correct_first = 2*(1-2*p) * log((1-p)/p)