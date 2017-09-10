% Compute: E [ X_i / (\sum_{j=1}^n X_j) ] where X_j ~ exp(lambda_j) independent
% 
% Input: 
% lambda_vec - vector of rates 
% 
% Output: 
% r_vec - vector of expectations
function r_vec = mean_exp_ratio_sum(lambda_vec)

n = length(lambda_vec); r_vec = zeros(n,1); 

for i=1:n 
    cur_lambda = vec2column(lambda_vec([1:(i-1) (i+1):n])); 
    
    r_vec(i) = quad2d( @(x,t) (t ./ (t+x)) .* lambda_vec(i).*exp(-lambda_vec(i)*t) .* expsumpdf(x, cur_lambda), ...
        0, 20 .* sum(1./lambda_vec) , 0, 20 ./ lambda_vec(i) ); % set limits 
end

