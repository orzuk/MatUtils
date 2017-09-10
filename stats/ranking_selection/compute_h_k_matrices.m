% Compute constants for two-stage procedures 
function [h_k_rinott, h_k_dalal, h_k_rinott_approx, h_k_dalal_approx] = compute_h_k_matrices(pcs_vec, nu_vec, k_vec)

num_p = length(pcs_vec); num_nu = length(nu_vec); num_k = length(k_vec); 

[h_k_rinott, h_k_dalal, h_k_rinott_approx, h_k_dalal_approx] = deal(zeros(num_p, num_nu, num_k+1));

%h_k_rinott = []; h_k_dalal = []; h_k_rinott_approx = []; h_k_dalal_approx = [];
for i_pcs = 1:length(pcs_vec)
    for i_nu = 1:length(nu_vec)
        p = pcs_vec(i_pcs), nu = nu_vec(i_nu)
        q_p = (-1/log(p))^(1/nu); % gevinv(p, nu, 1, 0); % the p-th quantile of nu-Frechet distribution with c.d.f.: e^{-x^{-nu}}
%        [h_k_rinott{i_pcs, i_nu}, h_k_rinott_approx{i_pcs, i_nu}, h_k_dalal{i_pcs, i_nu}, h_k_dalal_approx{i_pcs, i_nu}] = deal(zeros(1, num_k+1));
        
        for i_k=1:length(k_vec)
            k=k_vec(i_k);
            if(nu == Inf)
                h_k_dalal_approx(i_pcs, i_nu, i_k+1) = sqrt(2*log(k));
                h_k_rinott_approx(i_pcs, i_nu, i_k+1) = 2*sqrt(log(k));
            else
                h_k_dalal_approx(i_pcs, i_nu, i_k+1) = dalal_h_k_1(nu, k, p, 0);
                h_k_rinott_approx(i_pcs, i_nu, i_k+1) = dalal_h_k_1(nu, k, p, 0) * 2^(1/nu);
%                 
%                 
%                 (gamma((nu+1)/2) / (nu^(1-nu/2)*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p; % Compute approximations
%                 h_k_rinott_approx(i_pcs, i_nu, i_k+1) = (gamma((nu+1)/2) / (nu^(1-nu/2)*gamma(nu/2)*sqrt(pi))) ^ (1/nu) * k^(1/nu) * q_p * 2^(1/nu); % sqrt(2); % Compute approximations
            end
            h_k_dalal(i_pcs, i_nu, i_k+1) = fzero(@(x) two_stage_integral_dalal(x, k, nu)-p, ...
                [-0.01 3] .* h_k_rinott_approx(i_pcs, i_nu, i_k+1));
            h_k_rinott(i_pcs, i_nu, i_k+1) = fzero(@(x) two_stage_integral_rinott(x, nu)-(1-p^(1/k)), ... %   p, ... % , [tinv(epsilon, nu)-x -tinv(epsilon, nu)])-p, ...
                [-0.01 3] .* h_k_dalal_approx(i_pcs, i_nu, i_k+1));
        end
    end % loop on nu
end % loop on pcs

h_k_dalal = h_k_dalal(:, :, 2:end);
h_k_rinott = h_k_rinott(:, :, 2:end);
h_k_dalal_approx = h_k_dalal_approx(:, :, 2:end);
h_k_rinott_approx = h_k_rinott_approx(:, :, 2:end);



