% Improving Zigi's code - better readability and performance.
% The function returns the the estimated rotation matrix A_bar and the
% 'original' (estimated) indepndent data Y_bar
function [A_bar Y_bar] = spectral_ICA_Zuk(N, epsilon);

% Y1 ~ U[-sqrt(3),sqrt(3)], Y2 ~ N(0,1]. Y1, Y2 independent random variables
Y = randn(N,2); Y(:,1) = 2*sqrt(3)*(rand(N,1)-1/2);

% Orthogonal transformation: rotation
A = [1 1; 1 -1]./sqrt(2); X = Y * A;

figure; scatter(X(:,1),X(:,2),5); title('X (mixed) Data');

% Construct Laplacian matrix. W == weight matrix
W = exp(-((repmat(X(:,1), 1, N)-repmat(X(:,1)', N, 1)).^2 + (repmat(X(:,2), 1, N)-repmat(X(:,2)', N, 1)).^2) ./ (2.*epsilon)); 
s = sum(W,2); % s == proportional to kernel approximation of the density 

% find how many elements of s are small and remove them
small_inds = find(s <= 20); big_inds = setdiff([1:N], small_inds);
if(isempty(big_inds)) % all points are far away from each other
    sprintf('Data is too sparse. Try a larger epsilon ..')
    A_bar = []; Y_bar = []; return;
end
N_new = length(big_inds); Y_new = Y(big_inds,:); X_new = X(big_inds,:);  
figure; scatter(X_new(:,1),X_new(:,2),5); title('Data after removing isolated points');

% W_new == new weight matrix
W_new = W(big_inds, big_inds);
W_new = W_new ./ repmat(sum(W_new,2), 1, N_new);

% Laplacian: L = D^{-1}*W - I
% largest eigenvalues of D^{-1}*W <=> lowest eigenvalues of -L - WHY ? 
[eigenvectors, lambda] = eigs(W_new,5);
eigenvectors = sqrt(N)*eigenvectors ./ repmat(sum(eigenvectors.^2), N_new, 1); % normalize eigenvectors ||v||^2=N
(1-diag(lambda)) ./ (epsilon/2) % Print eigenvalues (?)

% Get the 2nd and 3rd eigenvectors (First one is trivial - not important)
vec2 = eigenvectors(:,2); vec3 = eigenvectors(:,3);

% The Z-statistics
Z = X_new' * vec2; Z = Z ./ norm(Z);  A_bar = zeros(2); A_bar(:,1) = Z; A_bar(1,2) = -Z(2); A_bar(2,2) = Z(1);
Z3 = X_new' * vec3; Z3 = Z3 ./ norm(Z3);
A_bar = X_new' * eigenvectors(:,2:3); A_bar(:,1) = A_bar(:,1)./norm(A_bar(:,1)); A_bar(:,2) = A_bar(:,2)./norm(A_bar(:,2));

Y_bar = X_new * A_bar'; % The recovered data

figure; hold on; scatter(X_new(:,1),X_new(:,2),5,vec2); title('Data colored according to first component');
plot(Z(1)*X_new(:,1), Z(2)*X_new(:,1), 'r'); plot(Z3(1)*X_new(:,1), Z3(2)*X_new(:,1), 'r'); 
figure; hold on; scatter(X_new(:,1),X_new(:,2),5,vec3); title('Data colored according to second component');
plot(Z(1)*X_new(:,1), Z(2)*X_new(:,1), 'r'); plot(Z3(1)*X_new(:,1), Z3(2)*X_new(:,1), 'r'); 
figure; hold on; scatter(Y_bar(:,1),Y_bar(:,2),5,vec2); title('Recovered (Ind. Comp.) Data ');
