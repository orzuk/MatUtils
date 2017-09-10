% Test the PCA function

nsamples = 1000;
nsnps = 120; 

x = 0.11*randn(nsamples, nsnps); % x  = floor(x*4);x(x == 2) = 3; x=x'; % snps example
x(:,101) = 120*rand(nsamples,1)-1; 
x(:,102) = 120*rand(nsamples,1)-1;
%x(:,102) = 2*sqrt(1-x(:,101).^2) .* sign(randn(nsamples,1)); % 10*rand(nsamples,1);;

 
figure; 
plot(x(:,101), x(:,102), '.'); title('Original Data');

[V Y] = PachterPolytopeAnalysis(x', 1);

L = 2; % How many top components to take
figure; plot(Y(1,:), Y(2,:), '.'); title('First Components Projection');

