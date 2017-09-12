% Compute properties of coalescent with analytic formulas
% 
% Input: 
% n - number of individuals in sample
% theta - population mutation rate (currently not used) 
% Output: 
% Mut - structure with mutations information 
% Trees - structure with trees information
%
function [Mut, Trees] = AnalyticCoalescent(n, theta)

har_vec = 1:(n-1); har_sum = sum(1./har_vec); 

Trees = []; 

Trees.T_Total_mean = 2*har_sum; 
Trees.T_MRCA_mean = 2*(1-1/n); 

Mut = []; 
Mut.age_mean = 0;  Mut.age_var = 0; 
for i=1:(n-1)
    Mut.age_mean = Mut.age_mean + sum( 1 ./ ( i .* (i:(n-1)) .* ((i+1):n) ));
    Mut.age_var = Mut.age_var + (1/i) * ...
        ( sum( 1 ./ ((i:(n-1)).^2 .* ((i+1):n).^2 )) + sum( 1 ./ ((i:(n-1)) .* ((i+1):n) )).^2 );
end
Mut.age_mean = Mut.age_mean * 2 / har_sum; % Normalize by H_{n-1}
Mut.age_var = Mut.age_var * 4 / har_sum - Mut.age_mean^2; % subtract squared mean 
