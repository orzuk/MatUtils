% Probability density function of sum of independent exponential random variables
% 
% Input: 
% x - where to evaluate density
% lambda_vec - vector of rates (assume column vector!!!)
% 
% Output
% f - density 
% 
function f = expsumpdf(x, lambda_vec)

if(size(x,1)>1) % here assume multiple rows 
    f = zeros(size(x)); 
    for i=1:size(x,1)
        f(i,:) = expsumpdf(x(i,:), lambda_vec); 
    end
    return; 
end

n=length(lambda_vec); m = length(x); % get # exponents and # points 


f = repmat(lambda_vec, 1, m) .* exp(-x .* lambda_vec); 

for i=1:n
    f(i,:) = f(i,:) .* prod( lambda_vec([1:i-1 i+1:n]) ./ (lambda_vec([1:i-1 i+1:n]) - lambda_vec(i)) ); 
end
if(n>1)
    f = sum(f); 
end

