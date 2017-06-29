% Compute h
% Input:
% nu - degrees of freedom (N0-1 where N0 is initial sample size)
% k - # of populations
% p - desired pcs
% proc_str - which procedure
% mode_str - how to compute: asymptotic or exact (numeric)
% 
% Output: 
% h - optimal constant h 
%
function h = two_stage_compute_h_k(pcs, nu, k, proc_str, mode_str)

% Set default
if(~exist('proc_str', 'var') || isempty(proc_str))
    proc_str = 'dalal';
end
if(~exist('mode_str', 'var') || isempty(mode_str))
    mode_str = 'asymptotic';
end

% Start with approx to get initial condition 
h = dalal_h_k_1(nu, k, pcs, 0);
if(strcmp(proc_str, 'rinott'))
    h = 2.^(1./nu) .* h;
end
switch mode_str
    case 'asymptotic'
      
    case 'numeric' % here solve integral exactly 
        if(strcmp(proc_str, 'rinott'))
            h = fzero(@(x) two_stage_integral_rinott(x, nu)-(1-pcs.^(1./k)), ...
            1 .* h ./ 2.^(1./nu)); % [-0.01 3] .* h ./ 2.^(1./nu));
        else
            h = fzero(@(x) two_stage_integral_dalal(x, k, nu)-pcs, ...
                1 .*  h .* 2.^(1./nu)); %[-0.01 3] .*  h .* 2.^(1./nu));
        end
end



