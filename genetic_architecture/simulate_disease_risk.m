% Simulate disease risk according to additive model
% This part is heavy as it needs many simulations
% 
% Input:
% f - genotypes vector
% d - disease contribution vector (additive)
% a - parameter of the disease prob. logistic function 1 / (1 +a*exp(-x))
% iters - number of simulations
% generations - maximal number of generations
%
% Output:
% S - distribution of disease risk score
%
function [S SW] = simulate_disease_risk(f, d, iters, generations)

N = length(d); % number of loci

ff = repmat(single(f), 1, iters); % watch that N*iters isn't too big (no memory)
dd = repmat(single(d), 1, iters); % watch that N*iters isn't too big (no memory)
x = single(rand(N,iters, 'single') < ff); y = single(rand(N,iters, 'single') < ff); % simulate genotypes
z = rand(N,iters,generations, 'single'); %  < ff); % Mother's chosen genotypes. Needed for relatives

w = single(rand(N,iters,generations, 'single') < 0.5); % Choice of each genotype. Needed for relatives
S = disease_risk(x,y,ff,dd);  % risk function
SW = cell(generations,1); 
for i=1:generations % loop over generations
    z(:,:,i) = single(z(:,:,i) < ff); % Mother's chosen genotypes. Needed for relatives
    if(i == 1)
        w(:,:,i) = x.*w(:,:,i) + y.*(1-w(:,:,i));
    else
        w(:,:,i) = w(:,:,i-1).*w(:,:,i) + z(:,:,i-1).*(1-w(:,:,i));
    end
    SW{i} = disease_risk(w(:,:,i),z(:,:,i),ff,dd); % risk function
end


