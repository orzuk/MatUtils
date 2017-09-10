% Compute enrichment p-value for pwms on regions
% CURRENTLY THIS FUNCTION IS CRAP AND NOT USED !!! 
function regions_pval_mat = regions_pwm_pval_from_file(regions_file, pwm_ind, do_log, strand, beta, ...
    in_matlab_flag, shuffle_pwm_flag, iters, varargin)

Assign24MammalsGlobalConstants;
if(machine == PC)
    in_matlab_flag = 1;
end

load(regions_file); num_pwms = size(pwms,1); num_regions = length(regions.packed_seqs);
if(in_matlab_flag) % here perform everything in matlab
    regions_affinity_mat = regions_pwm_affinity_from_file(regions_file, pwm_ind, do_log, strand, beta, ...
        in_matlab_flag, shuffle_pwm_flag, iters, varargin)
else
    in_matlab_flag = 0;
    for p=1:num_pwms % loop over pwms
        regions_pval_mat = regions_pwm_pval_from_file(regions_file, pwm_ind, do_log, strand, beta, ...
            in_matlab_flag, shuffle_pwm_flag, iters, varargin)
    end
end

