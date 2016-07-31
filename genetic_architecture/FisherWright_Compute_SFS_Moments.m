% Compute moments of the site-frequency spectrum for variable population
% size using the method in Evans, Shvets&Slatkin [2007]
% 
% Input: 
% N_vec - vector of population sizes 
% s - selection coefficients (currently works only for s=0). s shoudl be NEGATIVE for deletirious alleles
% max_k - maximum moment to compute 
% 
% Output: 
% mu_vec - vector of moments values 
% mu_vec_equilibrium - vector of moments values at equilibrium 
% 
function [mu_vec, mu_vec_equilibrium] = FisherWright_Compute_SFS_Moments(N_vec, s, max_k)

mu_vec = zeros(max_k, 1); 

N = N_vec(1); 
rho_vec = N_vec ./ N; 

if(~exist('s', 'var') || isempty(s)) % Assume s=0
    s = 0;
end

% Set t = ??? (scaling of time!) 
% First compute moments at equilirbium 
mu_vec_equilibrium =  zeros(max_k, 1); % need  a recurrence formula here !! 
    
int_one_over_rho = sum( 1./rho_vec) / (2*N); % take integral of 1/rho. Time unit is 1/2N with N=N_0 initial time
int_one_over_rho_cum_vec = cumsum( 1./rho_vec) ./ (2*N); % take integral of 1/rho

theta = 1; 
for k=1:max_k
    mu_vec_equilibrium(k) = absorption_time_by_selection(-abs(s), theta, N, 0, 1, -k-1);
end

% Compute non-equilibrium vec (for s=0)
mu_vec(1) = exp(-int_one_over_rho) * (mu_vec_equilibrium(1) + 0.5*theta * sum ( exp (int_one_over_rho_cum_vec) ) / (2*N));
for k=2:max_k
   mu_vec(k) = exp(-nchoosek(k+1,2)*int_one_over_rho) * ...
       ( mu_vec_equilibrium(k) + nchoosek(k,2) * sum( ( 1./rho_vec) .* mu_vec(k-1) .* exp (nchoosek(k+1,2) .* int_one_over_rho_cum_vec) ) /(2*N) ); 
end







