function [BETA, DEV, STAT, probit_data, p_x_x_z] =  test_probit(N)

% Run a test of matlab probit regression. Determine coefficients, simulate
% data, fit parameters and compare
beta = sqrt(0.05/0.09); %    0.75;
f = 0.1; 
mu = 0.01; 
%N = 2; % number of predictors;

h = beta^2*f*(1-f);

x_mu = norminv(1-mu); 

iters = 10000000;
x = (rand(iters,N) < f);  
y = randn(iters,1) ./ sqrt(1-N*h); 
z = (sum(beta*(x-f),2)+y) > x_mu; % make x unbiased (zero mean)

switch N
    case 1
        x_probit_matrix = [0; 1]
        probit_data = [sum((1-x).*z) sum(x.*z); iters-sum(x) sum(x)]'
    case 2
        x_probit_matrix = [0 0; 0 1; 1 0; 1 1];
        probit_data = ...
            [sum(z.*(1-x(:,1)).*(1-x(:,2))) sum(z.*(1-x(:,1)).*x(:,2)) sum(z.*x(:,1).*(1-x(:,2))) sum(z.*x(:,1).*x(:,2)); ...
            sum((1-x(:,1)).*(1-x(:,2))) sum((1-x(:,1)).*x(:,2)) sum(x(:,1).*(1-x(:,2))) sum(x(:,1).*x(:,2))]'
end

[BETA, DEV, STAT] = glmfit( x_probit_matrix, probit_data, 'binomial', 'link', 'probit')
p_x_x_z = mat2vec([probit_data(:,2) - probit_data(:,1) probit_data(:,1)]')' ./ iters; % output 2x2x2 empirical table 


% function LL = probit_loglike(x_probit_matrix, probit_data, beta) 
% LL=0;
% for i=1:size(probit_data,1) % loop on all cells

    
