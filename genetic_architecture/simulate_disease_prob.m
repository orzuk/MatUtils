% Simulate disease probability according to additive model
%
% Input:
% f - genotypes vector
% d - disease contribution vector (additive)
% a - parameter of the disease prob. logistic function 1 / (1 +a*exp(-x))
% iters - number of simulations
% generations - maximal number of generations
% S - ??? 
% SW - ??? 
%
% Output:
% P - distribution of disease prob.
% S - distribution of disease risk score
% rho - overall disease prob.
% h - heretability content. Disease prob. given a twin has disease
% r - relative heretability content. Disease prob. given a father has disease
%
function [P rho h r] = simulate_disease_prob(f, d, a, b, c, generations, S, SW)

N = length(d); % number of loci

% %S = zeros(N,1, 'single'); SW = zeros(N,1, 'single');
% ff = repmat(single(f), 1, iters); % watch that N*iters isn't too big (no memory)
% dd = repmat(single(d), 1, iters); % watch that N*iters isn't too big (no memory)
% x = single(rand(N,iters, 'single') < ff); y = single(rand(N,iters, 'single') < ff); % simulate genotypes
% z = rand(N,iters,generations, 'single'); %  < ff); % Mother's chosen genotypes. Needed for relatives
% w = single(rand(N,iters,generations, 'single') < 0.5); % Choice of each genotype. Needed for relatives
% S = disease_risk(x,y,ff,dd);  % risk function

P = 1 ./ (a+b.*exp(-c.*S)); % risk probability (sigmoid function)

rho = mean(P); % overall disease probability
h = mean(P.^2) / rho; % heretability content (prob. disease | twin has disease)
r = zeros(generations,1);
for i=1:generations % loop over generations
    PW = 1 ./ (a+b.*exp(-c.*SW{i}));
    r(i) = mean(P.*PW) / rho; % heretability content for relative (prob. disease | fater has disease)
end

% % % Old: loop (good for large values which don't fit into memory)
% % % for i=1:iters
% % %     x = double(rand(N,1) < f); y = double(rand(N,1) < f); % simulate genotypes
% % %     z = double(rand(N,1) < f); w = double(rand(N,1) < 0.5); % needed for relatives
% % %     S(i) = sum( ((-f).^(1-x) .* (1-f).^x + (-f).^(1-y) .* (1-f).^y) .* d ); % risk function
% % %     w = x.*w + y.*(1-w);
% % %     SW(i) = sum( ((-f).^(1-w) .* (1-f).^w + (-f).^(1-z) .* (1-f).^z) .* d ); % risk function
% % % end


