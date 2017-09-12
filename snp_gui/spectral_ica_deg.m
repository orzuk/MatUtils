function result = spectral_ICA_deg(N, epsilon);

% The degenarated spectral ICA algorithm for two uniform distributions

% Y1,Y2 ~ U[-sqrt(3),sqrt(3)]
% Y1, Y2 independent random variables
Y1 = 2*sqrt(3)*(rand(N,1)-1/2);
Y2 = 2*sqrt(3)*(rand(N,1)-1/2);

% Orthogonal transformation: rotation
X1 = (Y1 + Y2)/sqrt(2);
X2 = (Y1 - Y2)/sqrt(2);

figure;
scatter(X1,X2,5);

% Construct Laplacian matrix

% W == weight matrix
W = zeros(N);

for i=1:N;
    for j=1:N;
        W(i,j) = exp(-((X1(i)-X1(j))^2+(X2(i)-X2(j))^2)/(2*epsilon));   
    end;
    s = sum(W(i,:));
    % Row normalization of W: D^{-1}*W
    W(i,:) = W(i,:) / s;
end;

% Laplacian: L = D^{-1}*W - I
% largest eigenvalues of D^{-1}*W <=> lowest eigenvalues of -L
[eigenvectors, lambda] = eigs(W,4);

% normalize eigenvectors ||v||^2=1
for j=1:4;
    eigenvectors(:,j) = eigenvectors(:,j)/norm(eigenvectors(:,j));
end;

for i=1:4;
    (1-lambda(i,i))/(epsilon/2)
end;

% MATLAB eigs function has a little bug -- the eigenvalues are not necessarily ordered
eigenvalues = zeros(4,1);
for i=1:4;
    eigenvalues(i) = (1-lambda(i,i))/(epsilon/2);
end;

[sorted_eigenvals, idx_eigenvals] = sort(eigenvalues);

% The 2nd, 3rd and 4th eigenvectors
vec1 = eigenvectors(:,idx_eigenvals(2));
vec2 = eigenvectors(:,idx_eigenvals(3));
vec3 = eigenvectors(:,idx_eigenvals(4));

% Gram-Schmidt
vec1 = vec1 / norm(vec1);
vec2 = vec2 - dot(vec1,vec2)*vec1;
vec2 = vec2 / norm(vec2);

tmp = 1/2 * dot(vec3, (vec1 .* vec1 - vec2 .* vec2))/dot(vec3, vec1 .* vec2);
theta = atan(tmp)/2;

v1 = cos(theta) * vec1 - sin(theta) * vec2;
v2 = sin(theta) * vec1 + cos(theta) * vec2;

% The Z-statistics
x = [0;0];
for i=1:N;
    x(1) = x(1) + X1(i)*v1(i);
    x(2) = x(2) + X2(i)*v1(i);
end;


figure;
scatter(X1,X2,5,v1);

figure;
scatter(X1,X2,5,v2);

result = x/norm(x);
