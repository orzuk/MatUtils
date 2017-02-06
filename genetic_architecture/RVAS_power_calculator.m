% New: stand-alone power calculator for rare-variants association studies 
% 
% Input: 
% power_parameters_file - file with all power parameters
% Output: 
% power_output - power 
function power_output  = RVAS_power_calculator(power_parameters_file, varargin)


% Read parameters from file:
params_struct = ReadParametersFile(power_parameters_file); 

args = varargin;

for i=1:nargin % read all input variables 
    
end


%function [power_mat p_vals_vec stat_vec non_centrality_parameter] = ...
%    compute_association_power(p_mat, n_cases_vec, n_controls_vec, alpha_vec, ...
%    iters, test_type, test_stat, sampling_type, const_effect_flag, model_params, varargin)

%compute_association_power(); %this actually computers power


[power_mat, p_vals_vec, stat_vec, non_centrality_parameter] = ...
    compute_association_power(params_struct.p_mat, params_struct.n_cases_vec, params_struct.n_controls_vec, params_struct.alpha_vec, ...
    params_struct.iters, params_struct.test_type, params_struct.test_stat, params_struct.sampling_type, params_struct.const_effect_flag, params_struct.model_params);


% Organize output: 
power_output = power_mat; 

