function result = spectral_ICA(N, epsilon);

% Y1 ~ U[-sqrt(3),sqrt(3)]
% Y2 ~ N(0,1]
% Y1, Y2 independent random variables
Y1 = 2*sqrt(3)*(rand(N,1)-1/2);
Y2 = randn(N,1);

% Orthogonal transformation: rotation
X1 = (Y1 + Y2)/sqrt(2);
X2 = (Y1 - Y2)/sqrt(2);

figure;
scatter(X1,X2,5);

% Construct Laplacian matrix

% W == weight matrix
W = zeros(N);

% s == proportional to kernel approximation of the density 
s = zeros(N,1);
for i=1:N;
    for j=1:N;
        W(i,j) = exp(-((X1(i)-X1(j))^2+(X2(i)-X2(j))^2)/(2*epsilon));   
    end;
    s(i) = sum(W(i,:));
end;

% Remove isolated points
[s_sort, s_idx] = sort(s);

% find how many elements of s are small
small = 0;
i = 1;
while (s_sort(i) <= 20)
    small  = small + 1;
    i = i+1;
end;

N_new = N - small

Y1_new = zeros(N_new,1);
Y2_new = zeros(N_new,1);
X1_new = zeros(N_new,1);
X2_new = zeros(N_new,1);

i=1;
k=0;
for i=1:N; 
    if (s(i) > 20)
        k = k + 1;
        X1_new(k) = X1(i);
        X2_new(k) = X2(i);
        Y1_new(k) = Y1(i);
        Y2_new(k) = Y2(i);
    end;
end;

figure;
scatter(X1_new,X2_new,5);

% W_new == weight matrix
W_new = zeros(N_new);
 
s_new = zeros(N_new,1);
for i=1:N_new;
    for j=1:N_new;
        W_new(i,j) = exp(-((X1_new(i)-X1_new(j))^2+(X2_new(i)-X2_new(j))^2)/(2*epsilon));   
    end;
    % Row normalization of W: D^{-1}*W
    s_new(i) = sum(W_new(i,:));
    W_new(i,:) = W_new(i,:) / s_new(i);
end;


% Laplacian: L = D^{-1}*W - I
% largest eigenvalues of D^{-1}*W <=> lowest eigenvalues of -L
[eigenvectors, lambda] = eigs(W_new,5);


% normalize eigenvectors ||v||^2=N
for j=1:5;
    eigenvectors(:,j) = sqrt(N)*eigenvectors(:,j)/norm(eigenvectors(:,j));
end;

for i=1:5;
    (1-lambda(i,i))/(epsilon/2)
end;

% The 2nd eigenvector
vec2 = eigenvectors(:,2);

% The 3rd eigenvector
vec3 = eigenvectors(:,3);

% The Z-statistics
x = [0;0];
for i=1:N_new;
    x(1) = x(1) + X1_new(i)*vec2(i);
    x(2) = x(2) + X2_new(i)*vec2(i);
end;


figure;
scatter(X1_new,X2_new,5,vec2);

figure;
scatter(X1_new,X2_new,5,vec3);

result = x/norm(x);
