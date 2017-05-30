% Apply Dalal procedure
%
% Input:
% N0 - # samples to draw at first stage
% PCS - desired Prob. of correct selection
% Delta - guarantee on differences between true means
% mu_vec - vector of true means
% sigma_vec - vector of true st.d.s
%
% Output:
% CS - yes or now

function [max_I, PCS_res, h] = dalal_procedure(N0, PCS, Delta, mu_vec, sigma_vec, iters, method_str)

n_pop = length(mu_vec);

q_p = gevinv(PCS, N0-1, 1,0); % what should we take here? how many deg. freedom? what is v-Frechet? 

% Set init condition
gamma_v = (gamma(N0/2) / ((N0-1)*gamma((N0-1)/2)*sqrt(pi))) ^ (1/(N0-1)) * (n_pop-1)^(1/(N0-1)) * q_p; 

switch lower(method_str)
    case 'rinott'
        h = fsolve(@(x) two_stage_integral_rinott(x, n_pop-1, N0-1)-PCS, gamma_v);
    case 'dalal'
        h = fsolve(@(x) two_stage_integral_dalal(x, n_pop-1, N0-1)-PCS,  gamma_v);
end

for t=1:iters
    run_iter = t
    X = zeros(N0, n_pop);
    for i=1:n_pop % new sampling
        X(:,i) = normrnd(mu_vec(i), sigma_vec(i), 1, N0);
    end
    S_vec = var(X);
    N_vec = ceil(max(N0+1, (h/Delta)^2 .* S_vec)); % Choose new sample sizes
    
    for i=1:n_pop % new sampling
        new_X{i} = normrnd(mu_vec(i), sigma_vec(i), N_vec(i), 1);
        
        switch lower(method_str)
            case 'rinott'
                a_vec{i} = repmat(1/N_vec(i), N_vec(i), 1); % simple mean
            case 'dalal'
                C2 = (Delta/h)^2 / S_vec(i);
                % set a vec
                a = ( 2*N0 + sqrt( 4*N0^2 - 4*N0*N_vec(i)*( 1 - C2*(N_vec(i)-N0) ) ) ) / (2*N0*N_vec(i));
                b = (1-N0*a)^2 / (N_vec(i)-N0);
                a_vec{i} = repmat(a, N_vec(i), 1); a_vec{i}(N0+1:end) = b;
        end % switch
        new_mu(i) = sum(a_vec{i}.*new_X{i});
    end
    [~, max_I(t)] = max(new_mu);
end % loop on iters
[~, true_max_I] = max(mu_vec);
PCS_res = mean(max_I == true_max_I);
