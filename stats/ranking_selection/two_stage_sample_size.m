% Compute sample size 
% Input:
% nu - degrees of freedom (N0-1 where N0 is initial sample size)
% k - # of populations
% p - desired pcs
% proc_str - which procedure
% mode_str - how to compute: asymptotic or exact (numeric)
% 
function N = two_stage_sample_size(pcs, nu, k, proc_str, mode_str, Delta, sigma2_vec)

% Set default
if(~exist('proc_str', 'var') || isempty(proc_str))
    proc_str = 'dalal';
end
if(~exist('mode_str', 'var') || isempty(mode_str))
    mode_str = 'asymptotic';
end
if(isscalar(sigma2_vec))
    sigma2_vec = repmat(sigma2_vec, k, 1); 
end

N = max(nu+2, two_stage_compute_h_k(pcs, nu, k, proc_str, mode_str).^2 .* sum(sigma2_vec) ./ Delta^2); 


