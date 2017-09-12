function result = spectral_ICA3(N, epsilon);

rand('state',sum(100*clock))

% Y1 ~ U[-sqrt(3),sqrt(3)]
% Y2 ~ (U + U) / sqrt(2)
% Y3 ~ U^3
% Y1, Y2, Y3 independent random variables with E(Y)=0, E(Y^2)=1.
Y1 = 2*sqrt(3)*(rand(N,1)-1/2);

Z1 = 2*sqrt(3)*(rand(N,1)-1/2);
Z2 = 2*sqrt(3)*(rand(N,1)-1/2);
Y2 = (Z1 + Z2)/sqrt(2);

Y3 = (2*(7^(1/6))*(rand(N,1)-1/2)).^3;

std(Y1)
std(Y2)
std(Y3)

% Orthogonal transformation: rotation
% The output matrix B (see below) 
% should be B = [1/sqrt(3), 1/sqrt(3), 1/sqrt(3) ; 1/sqrt(2), -1/sqrt(2), 0 ; 1/sqrt(6), 1/sqrt(6), -2/sqrt(6)] 
X1 = (Y1 + Y2 + Y3) / sqrt(3);
X2 = (Y1 - Y2) / sqrt(2);
X3 = (Y1 + Y2 - 2*Y3) / sqrt(6);

% Construct Laplacian matrix

% W == weight matrix
W = zeros(N);

% s == proportional to kernel approximation of the density 
s = zeros(N,1);
for i=1:N;
    for j=1:N;
        W(i,j) = exp(-((X1(i)-X1(j))^2+(X2(i)-X2(j))^2 + (X3(i)-X3(j))^2)/(2*epsilon));   
    end;
    s(i) = sum(W(i,:));
end;

% Laplacian: L = D^{-1}*W - I
% largest eigenvalues of D^{-1}*W <=> lowest eigenvalues of -L
[eigenvectors, lambda] = eigs(W,5);

% normalize eigenvectors ||v||^2=N
for j=1:5;
    eigenvectors(:,j) = sqrt(N)*eigenvectors(:,j)/norm(eigenvectors(:,j));
end;

for i=1:5;
    (lambda(i,i)-1)/(epsilon/2)
end;

% The 2nd eigenvector
phi2 = eigenvectors(:,2);
phi3 = eigenvectors(:,3);
phi4 = eigenvectors(:,4);

% The Z-statistics
vec2 = [0;0;0];
vec3 = [0;0;0];
vec4 = [0;0;0];

for i=1:N;
    vec2(1) = vec2(1) + X1(i)*phi2(i);
    vec2(2) = vec2(2) + X2(i)*phi2(i);
    vec2(3) = vec2(3) + X3(i)*phi2(i);
    
    vec3(1) = vec3(1) + X1(i)*phi3(i);
    vec3(2) = vec3(2) + X2(i)*phi3(i);
    vec3(3) = vec3(3) + X3(i)*phi3(i);
    
    vec4(1) = vec4(1) + X1(i)*phi4(i);
    vec4(2) = vec4(2) + X2(i)*phi4(i);
    vec4(3) = vec4(3) + X3(i)*phi4(i);
end;

% Normalize vec2, vec3, vec4
vec2 = vec2 / norm(vec2);
vec3 = vec3 / norm(vec3);
vec4 = vec4 / norm(vec4);

% Find best orthogonal matrix approximation
A = [vec2 , vec3 , vec4];
[U,D,V] = svd(A);
Q = U * V'
%Q = [vec2, vec3, vec4]

b1 = [1,1,1]/sqrt(3);
b2 = [1,-1,0]/sqrt(2);
b3 = [1,1,-2]/sqrt(6);

B = [b1;b2;b3]
