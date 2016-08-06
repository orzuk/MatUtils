% Compute density of Normalized Gaussian vector 
Sigma= [1 8; 8 100];
n = 100000;

[Z, T] = mvnrnd([0, 0]' , Sigma, n);

figure; plot(Z(:,1), Z(:,2), '.');

Z_norm = sqrt(sum(Z.^2, 2));

Z_u = Z ./ repmat(Z_norm, 1, 2);
figure; plot(Z_u(:,1), Z_u(:,2), '.');

Sigma_empir = cov(Z)
Sigma_u_empir = cov(Z_u)

[U, Lambda] = eig(Sigma)
[U_empir, Lambda_empir] = eig(Sigma_empir)

[U_u_empir, Lambda_u_empir] = eig(Sigma_u_empir)


W = ( inv(sqrt(Lambda)) * U' * Z')'

W_u = ( inv(sqrt(Lambda_u_empir)) * U_u_empir' * Z_u')'

cov(W_u)

Theta = atan(Z(:,2) ./  Z(:,1)); figure; hist(Theta, 500); xlabel('Theta'); 
Theta_W = atan(W(:,2) ./  W(:,1)); figure; hist(Theta_W, 500); xlabel('Theta_W'); 
Theta_u = atan(Z_u(:,2) ./  Z_u(:,1)); figure; hist(Theta_u, 500); xlabel('Theta_u'); 
Theta_W_u = atan(W_u(:,2) ./  W_u(:,1)); figure; hist(Theta_W_u, 500); xlabel('Theta_W_u'); 

[Z_s] = mvnrnd([0, 0]' , eye(2), n);
Z_s_norm = sqrt(sum(Z_s.^2, 2));

W_s = Z_s ./ repmat(Z_s_norm, 1, 2); 

cov(W_s)
Theta_W_s = atan(W_s(:,2) ./  W_s(:,1)); figure; hist(Theta_W_s, 500); xlabel('Theta_W'); 

