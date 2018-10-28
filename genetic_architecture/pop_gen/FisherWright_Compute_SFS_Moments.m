% Compute moments of the site-frequency spectrum for variable population size using the method in Evans, Shvets&Slatkin [2007]
% Here moment k is: theta * \int_{0}^1 \psi_0(f) f^k (1-f) df, where k starts from 1    
%
% Input: 
% N_vec - vector of population sizes 
% s - selection coefficients (currently works only for s=0). s should be NEGATIVE for deleterious alleles
% max_k - maximum moment to compute 
% init_str - initial conditions 
% theta - effective mutation rate (2*N*mu for initial N) 
% 
% Output: 
% mu_vec - matrix of moments values. mu_vec(k, j) is k-th moment at generation j.  
% mu_vec_equilibrium - vector of moments values at equilibrium 
% 
function [mu_vec, mu_vec_equilibrium] = FisherWright_Compute_SFS_Moments(N_vec, s, max_k, init_str, theta)

if(~exist('init_str', 'var') || isempty(init_str))
    init_str = 'equilibrium'; % 'newly_born'; % default: start at equilibrium
end
if(~exist('theta', 'var') || isempty(theta))
    theta = 1; 
end
if(~exist('s', 'var') || isempty(s)) % Assume s=0
    s = 0;
end
num_generations = length(N_vec); 
mu_vec = zeros(max_k, num_generations); 
N = N_vec(1); 
rho_vec = N_vec ./ N; % relative population size 

    % int_one_over_rho = sum( 1./rho_vec) / (2*N); % take integral of 1/rho. Time unit is 1/2N with N=N_0 initial time
int_one_over_rho_cum_vec = cumsum( 1./rho_vec) ./ (2*N); % take integral of 1/rho
exp_minus_int_one_over_rho_cum_vec = exp(-int_one_over_rho_cum_vec); 
exp_int_one_over_rho_cum_vec = exp(int_one_over_rho_cum_vec); 
exp_int_one_over_rho_cum_cum_vec = cumsum(exp_int_one_over_rho_cum_vec); 

% First compute moments at equilirbium 
mu_vec_equilibrium =  zeros(max_k, 1); % need a recurrence formula here !! 
for k=1:max_k
    mu_vec_equilibrium(k) = absorption_time_by_selection(s, theta, N, 0, 1, -k-1);  % at equilibrium compute: int_{f=0}^1   f(1-f) \psi_s(f) df 
end
switch init_str
    case 'equilibrium'
        mu_vec_init = mu_vec_equilibrium;
    case {'newly_born', 'empty'}
        mu_vec_init = zeros(size(mu_vec_equilibrium));         
end

% Compute non-equilibrium vec (only for neutral! s=0) - use eqs. (36-37) in in Evans et al.:  
mu_vec(1,:) = exp_minus_int_one_over_rho_cum_vec .* (mu_vec_init(1) + 0.5*theta * exp_int_one_over_rho_cum_cum_vec ./ (2*N)); % first moment 
for j=1:num_generations
%     if(j == 1) % alternative: exact, discrete calculation  
%        mu_vec2(j) = theta; % 4*N*mu; 
%     else
%        mu_vec2(j) = (mu_vec2(j-1) + theta / (2*N) ) .* (1- 1/(2*N_vec(j))); 
%     end        
    for k=2:max_k
        mu_vec(k,j) = exp(-0.5*k*(k+1)*int_one_over_rho_cum_vec(j)) * ...
            ( mu_vec_init(k) + 0.5*k*(k-1) * sum( ( 1./rho_vec(1:j)) .* mu_vec(k-1,j) .* ...
            exp ( 0.5*k*(k+1) .* int_one_over_rho_cum_vec(1:j)) ) /(2*N) );
    end
end
% figure; plot(mu_vec(1,:), mu_vec2, '*'); xlabel('Ewans Het.'); ylabel('Exact Het.'); 
% figure; plot(2*mu_vec(1,:) - mu_vec2, '*'); xlabel('Ewans Het.'); ylabel('Exact Het.'); 

% compute special case: expansion solution (eq. (39) in Evans et al.) . Here initial conditions are zero !! 
%R = 2*N_vec(1)*log(N_vec(2)/N_vec(1)); % get expansion factor 
%t_vec = (1:num_generations) ./ (2*N); 
%mu_vec_expansion = (theta/(2*R)) .* exp ( exp (-R.*t_vec) ./ R) .* ( Ei(-1/R) - Ei(-exp(-R.*t_vec) ./ R) );
%mu_vec_expansion_old = real( (theta/(2*R)) .* exp ( exp (-R.*t_vec) ./ R) .* Ei(-1/R) ); % seperate old and new alleles
%mu_vec_expansion_new = real( -(theta/(2*R)) .* exp ( exp (-R.*t_vec) ./ R) .* Ei(-exp(-R.*t_vec) ./ R) );

