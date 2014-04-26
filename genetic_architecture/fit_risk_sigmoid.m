% Fit the sigmoid risk function 1/(a+b*exp(-c*x)) to 
% match a certain disease prob. rho and twin risk h
%
% Input:
% S - distribution of disease risk score
% rho - overall disease prob.
% h - heretability content. Disease prob. given a twin has disease
% iters - number of simulations
%
% Output:
% P - distribution of disease prob.
% a - parameter of the disease prob. logistic function 1 / (1 +a*exp(-x))
% b - parameter of the disease prob. logistic function 1 / (1 +a*exp(-x))
% c - parameter of the disease prob. logistic function 1 / (a +b*exp(-c*x))
%
% Where is the middle? 
% We need to solve: 1 / (a + b*exp(-c*x)) = 1/2a
% So this gives: exp(-c*x) = a/b
% -c*x = log(a/b)
% x = -log(a/b) / c
%
function [a b c] = fit_risk_sigmoid(S, rho, h)

%num_bins = 100; 
% [h_freq x_vec] = hist_density(S, num_bins); % get histogram of risk function 

% a = 1; % fzero(@(x) hist_sigmoid(x, S, rho), 1); % fit rho to get a. The maximum disease prob is 1/a, therefore we must have a>1
% b = fzero(@(x) hist_sigmoid(1, x, 1, S, rho), 1); % fit rho to get b. Assume a=c=1
sigmoid_mu_0 = quantile(S, 1-rho); % an initial guess for mean of the sigmoid

% Try a grid from one to 10
min_err = 1000000; 
for i=-5:10 % loop on b
    try_opt_i = i
    for j=-5:10 % loop on c
        abc_vec = fminsearch(@(x) hist_sigmoid(x, S, rho, h), [exp(i) exp(j)]); %  log(2) / double(sigmoid_mu_0)]);
        cur_err = hist_sigmoid(abc_vec, S, rho, h);
        if(cur_err < min_err)
            min_err = cur_err; 
            best_abc_vec = abc_vec;
        end
    end
end
abc_vec = best_abc_vec
a = 1; % abc_vec(1); % Temp ( don't have data to fit it ..) 
b = abc_vec(1); c = abc_vec(2); 

% Compute the rho and h for disease
% Input: 
% abc_vec - a,b and c parameters
% S - ?? 
% rho_0 - ?? 
% h_0 - ?? 
% 
% Output: 
% rmsf - L2 error
% rho - ?? 
% h - ?? 
% 
function [rmsd rho h] = hist_sigmoid(abc_vec, S, rho_0, h_0) 

b = abc_vec(1); c = abc_vec(2); 
a=1; % a = abc_vec(1); % Temp. (don't have data to fit it ..) 
P = 1 ./ (a+b.*exp(-c.*S)); % risk probability (sigmoid function)
rho = mean(P); 
h = (mean(P.^2)) / mean(P);

rmsd = (rho - rho_0)^2/rho_0^2 + (h - h_0)^2/h_0^2 + ...
    exp(-100*b) + exp(-100*c); % heretability content (prob. disease | twin has disease); % overall disease probability

%a = binary_search(min_a, max_a, density(

