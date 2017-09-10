% Gerenate an instance of chinese restaurant process
% 
% Input: 
% n - number of people 
% alpha - process param (Poisson rate - default is one)
% 
% Output: 
% s - matrix of generated partitions of n 
% 
function s = indian_buffet_process(n, alpha)

if(~exist('alpha', 'var') || isempty(alpha)) % set default parameters 
    alpha = 0; 
end

s = zeros(n, n*alpha); 

new_dishes = poissrnd(alpha); s(1,1:new_dishes) = 1; num_dishes = new_dishes; 

for i=2:n % loop and add person
    dish_weights = sum(s(:,1:num_dishes)) ./ i; 
    s(i,1:num_dishes) = rand(1,num_dishes) < dish_weights; % sample each existing dish
    
    new_dishes = poissrnd(alpha / i); s(i,num_dishes+1:num_dishes+new_dishes) = 1; 
    num_dishes = num_dishes + new_dishes;
end
s = s(:,1:num_dishes); 

