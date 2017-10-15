% New: stand-alone power calculator for rare-variants association studies
%
% Input:
% power_parameters_file - file with all power parameters
% Output:
% power_output - power
function power_output  = RVAS_power_calculator(power_parameters_file, varargin)


% Read parameters from file:
params_struct = ReadParametersFile(power_parameters_file);

if(~isfield(params_struct, 'output_file'))
    params_struct.output_file = [power_parameters_file '.out'];
end



% Preprocessing
power_struct = [];
num_OR = length(params_struct.OR_vec);
num_prev = length(params_struct.prevalence);
num_freq = length(params_struct.f_vec);
num_alpha = length(params_struct.alpha_vec);

if(~isfield(params_struct, 'iters'))
    params_struct.iters = 1000; % set default iterations for sampling
end
if(~isfield(params_struct, 'test_stat')) % set default test statistic
    params_struct.test_stat = 'chi_square_analytic';
end
params_struct.const_effect_flag = 0;
params_struct.model_parms = [];
params_struct.sampling_type = 'case_control';
ctr=1;
for i=1:num_OR
    for j=1:num_prev
        for k=1:num_freq
            params_struct.GRR_vec = odds_ratio_to_genetic_relative_risk(params_struct.OR_vec(i), params_struct.f_vec(k), params_struct.prevalence(j));
            p_z_x_marginal = genetic_relative_risk_to_p_z_x_marginal(params_struct.f_vec(k), params_struct.GRR_vec, params_struct.prevalence(j));
            for l=1:num_alpha
                power_mat = ...
                    compute_association_power(p_z_x_marginal, params_struct.n_cases, params_struct.n_controls, params_struct.alpha_vec(l), ...
                    params_struct.iters, params_struct.test_type, params_struct.test_stat); %
                %, params_struct.sampling_type,    params_struct.const_effect_flag, params_struct.model_params);
                
                power_struct.OR(ctr) = params_struct.OR_vec(i);
                power_struct.prev(ctr) = params_struct.prevalence(j);
                power_struct.freq(ctr) = params_struct.f_vec(k);
                power_struct.n_cases(ctr) = params_struct.n_cases;
                power_struct.n_controls(ctr) = params_struct.n_controls;
                power_struct.alpha(ctr) = params_struct.alpha_vec(l);
                power_struct.power(ctr) = round(power_mat, 3);
                ctr=ctr+1;
            end
        end
    end
end
% Organize output and write to file
cell_to_mat = 0;
WriteDataFile(power_struct, params_struct.output_file, cell_to_mat); %  skip_lines)

power_output = power_mat;

