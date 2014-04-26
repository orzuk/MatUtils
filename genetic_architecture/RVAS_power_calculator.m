% New: stand-alone power calculator for rare-variants association studies 
function power_output  = RVAS_power_calculator(power_parameters_file, varargin)


% Read parameters from file:
params_struct = ReadParametersFile(power_parameters_file); 

args = varargin;

for i=1:nargin % read all input variables 
    
end


%function [power_mat p_vals_vec stat_vec non_centrality_parameter] = ...
%    compute_association_power(p_mat, n_cases_vec, n_controls_vec, alpha_vec, ...
%    iters, test_type, test_stat, sampling_type, const_effect_flag, model_params, varargin)

compute_association_power(); %this actually computers power
