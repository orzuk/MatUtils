% Compute a simple estimator for selection coefficients from Ka and Ks values
% Use formula from [Nielsen and Yang, MBE 2003] for Ka/Ks as function of S 
% 
% Input: 
% Ka - count of non-synonymous variants 
% Ks - count of synonymous variants 
% 
% Output: 
% w - ratiof counts 
% s - corresponding selection coefficient 
% 
function [s, w] = KaKs_to_selection_coefficient(Ka, Ks)

AssignRVASConstants; % Set N_eff 

n = length(Ka); s = zeros(n, 1); 
w = Ka/Ks; % compute ratio 
for i=1:n
    if(mod(i, 100) == 0)
        compute_s = i
    end
    s(i) = fzero(@(S) S - (1-exp(-S))*max(10^(-6),w(i)), 2*log(max(10^(-6),w(i)))); % solve with close-by initial condition
end
s = s ./ (4*N_eff); % go from Big S to little s
