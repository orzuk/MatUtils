% Test first order of entropy for HMM for general alphabet
TOL = 0.0000000001;

BSC_flag = 0;

if(BSC_flag)
    S = 2;
    eps = 0.1;
    R = [eps, 1-eps; 1-eps, eps];
    T = [0.5, -0.5; -1, 1];
else
    S = 5;
    % Create R, which should be ergodic
    R = abs(rand(S));
    R_tag = abs(rand(S));
    
    for i=1:S
        R(i,:) = R(i,:) ./ sum(R(i,:));
        R_tag(i,:) = R_tag(i,:) ./ sum(R_tag(i,:));
    end

    % Create T, with negatives on the diagonal and rows sum equal to zero
    T = rand(S); T_tag = rand(S);
    for i=1:S
        T(i,i)  = T(i,i)-sum(T(i,:));
        T_tag(i,i)  = T_tag(i,i)-sum(T_tag(i,:));
    end

end


% crerate the vector of stationary distribution
pi = rand(1,S);
pi = pi ./ sum(pi);

pi = ones(1,S); pi = pi ./ S;

xi = ones(S,1);

U = ones(S); U = U ./ S;

delta = 0.00; %0.5;
M = U + delta * T;

B = R'* diag(pi) * M * R;



D  = R' * diag(pi) * U * R; % use D for simplicity

Exact_H = ( xi' * (B .* log(B)) - (log(pi*R)) * B) * xi

try_another_exact = sum(sum(B .* log(B))) - sum((log(pi*R)) * B)


zero_order_should_be = xi' * ( D .* log(D))- log(pi * R) * D;
zero_order_should_be = zero_order_should_be * xi
log_S_is = log(S)


try_try = (1/S) .* xi' * R * log(R' * xi) - log(S)



UIUI = eye(S)-U;
UIUI(:,S) = 1;

piti = pi * T; piti(S)  = 0;

alpha = pi + delta * piti * inv(UIUI)

[vec lam] = eig(M');

ind = 1;
for i=1:S
    if( abs(lam(i,i) - 1) < TOL)
        ind = i;
    end
end

alpha_correct = vec(:,ind) ./ sum(vec(:,ind))

alpha - alpha_correct'

% Now that we have alpha we can calculate the first order also :
BIG_TEMP = R' * ( diag(pi) * T + diag (piti * inv(UIUI)) * U) * R;
first_order_should_be = ( xi' * (BIG_TEMP .* (S .* U + log(D))) - log(pi * R) * BIG_TEMP ) * xi


DDD = T' * diag(pi) * M + diag(pi) * M * T






% Try to compute the KL orders

% Now we want to check what is zero for zero-order. Nothing is zero here
% ...
is_zero = xi' * ( (R'*U*R) .* log((1/S)*R'*U*R) ) * xi;
is_zero = xi' * ( (R'*U*R) .* log((1/S)*R_tag'*U*R_tag) ) * xi;
is_zero = log((1/S) * xi' * R) * (R'*U*R) * xi;



KL_zero_order = sum( (1/S).*sum(R,1) .* (log((1/S).*sum(R,1)) - log((1/S).*sum(R_tag,1))) )


KL_zero_order_matrix_way = -(1/S) .* ( log((1/S) * xi' * R) - log((1/S) * xi' * R_tag) ) * (R'*U*R) * xi;
KL_zero_order_matrix_way = 1*KL_zero_order_matrix_way + 0 .* xi' * ((1/S) .* (R'*U*R) .* (log((1/S)*R'*U*R)-log((1/S)*R_tag'*U*R_tag)) ) * xi

F = (1/S) .* R' * U * R; F_tag = (1/S) .* R_tag' * U * R_tag;
KL_Full = (xi' * (F .* log(F)) - log((1/S).* xi' * R) * F) * xi - (xi' * (F .* log(F_tag)) - log((1/S).* xi' * R_tag) * F) * xi 




% Now go to the first order

% First, calculate phi : 
U_minus_I = U - eye(S);
U_minus_I(:,end) = 1;
xi_T = xi' * T; xi_T(end) = 0;
phi = xi_T * inv(U_minus_I);
xi_T_tag = xi' * T; xi_T_tag(end) = 0;
phi_tag = xi_T_tag * inv(U_minus_I);


% Now we want to check what is zero for first-order
is_zero = xi' * [R' * ((1/S)*T + diag(phi) * U) * R] * xi;
is_zero = xi' * [R' * ((1/S)*T + diag(phi) * U) * R] * log((1/S)*R'*U*R) * xi;
is_zero = xi' * [R' * ((1/S)*T + diag(phi) * U) * R] * log((1/S)*R_tag'*U*R_tag) * xi;
is_zero = xi' * [R' * ((1/S)*T) * R] * log((1/S)*R'*U*R) * xi;
is_zero = xi' * [R' * ( diag(phi) * U) * R] * log((1/S)*R'*U*R) * xi;
is_zero = xi' * [R' * ((1/S)*T) * R] * log((1/S)*R_tag'*U*R_tag) * xi;
is_zero = xi' * [R' * ( diag(phi) * U) * R] * log((1/S)*R_tag'*U*R_tag) * xi;



is_zero = xi' * [(R_tag' * ((1/S)*T_tag + diag(phi_tag) * U) * R_tag) .* ((R'*U*R) ./ (R_tag'*U*R))] * xi;
is_zero = xi' * [(R_tag' * ((1/S)*T_tag ) * R_tag) .* ((R'*U*R) ./ (R_tag'*U*R))] * xi;
is_zero = xi' * [(R_tag' * ( diag(phi_tag) * U) * R_tag) .* ((R'*U*R) ./ (R_tag'*U*R))] * xi;

is_zero = (log((1/S)*xi'*R) - log((1/S)*xi'*R_tag)) * [R' * ((1/S)*T + diag(phi) * U) * R] * xi
is_zero = (log((1/S)*xi'*R) ) * [R' * ((1/S)*T + diag(phi) * U) * R] * xi
is_zero = ( - log((1/S)*xi'*R_tag)) * [R' * ((1/S)*T + diag(phi) * U) * R] * xi
is_zero = (log((1/S)*xi'*R) - log((1/S)*xi'*R_tag)) * [R' * ((1/S)*T ) * R] * xi
is_zero = (log((1/S)*xi'*R) - log((1/S)*xi'*R_tag)) * [R' * (diag(phi) * U) * R] * xi
is_zero = (log((1/S)*xi'*R) ) * [R' * ((1/S)*T ) * R] * xi
is_zero = (log((1/S)*xi'*R) ) * [R' * (diag(phi) * U) * R] * xi

is_zero = ( (phi * R) ./ ((1/S)*xi'*R) - (phi_tag * R_tag) ./ ((1/S)*xi'*R_tag) ) * ((1/S)*R'*U*R) * xi
is_zero = ( (phi * R) ./ ((1/S)*xi'*R)  ) * ((1/S)*R'*U*R) * xi


KL_first_order_matrix_way(1) = xi' * [(1/S) * R' * T * R] * (log((1/S)*R'*U*R)-log((1/S)*R_tag'*U*R_tag)) * xi 
KL_first_order_matrix_way(2) = xi' * [(R_tag' * diag(phi_tag) * U * R_tag) .* ((R'*U*R) ./ (R_tag'*U*R))] * xi 
KL_first_order_matrix_way(3) = (log((1/S)*xi'*R) - log((1/S)*xi'*R_tag)) * [R' * diag(phi) * U * R] * xi
KL_first_order_matrix_way(4) = ( (phi_tag * R_tag) ./ ((1/S)*xi'*R_tag)  ) * ((1/S)*R'*U*R) * xi


KL_first_order_matrix_way(1) = xi' * [(1/S) * R' * T * R] * (log((1/S)*R'*U*R)) * xi 
KL_first_order_matrix_way(2) = xi' * [(R_tag' * diag(phi_tag) * U * R_tag) .* ((R'*U*R) ./ (R_tag'*U*R))] * xi 
KL_first_order_matrix_way(3) = (log((1/S)*xi'*R) ) * [R' * diag(phi) * U * R] * xi
KL_first_order_matrix_way(4) = ( (phi_tag * R_tag) ./ ((1/S)*xi'*R_tag)  ) * ((1/S)*R'*U*R) * xi

