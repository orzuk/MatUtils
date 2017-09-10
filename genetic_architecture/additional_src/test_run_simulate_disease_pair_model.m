% Script for running simulations to find good architectures
queue_str = 'broad'; % 'priority';
architecture_str_vec = {'and-of-sigmoids'}; %    'DNF', 'CNF', 'sigmoid' ... % }; % 'and-of-sigmoids',
%    'or-of-sigmoids', 'sum-of-ors', 'sum-of-ands'}; % , 'DNF',  }; % , ...
%       'sigmoid'  ,'random', 'monotone', 'and', 'xor', 'or', 'additive'}; % , 'monotone'};

% action_mode = 1; in_matlab_flag = 0;
if(~exist('in_matlab_flag', 'var'))
    in_matlab_flag = 1;
end
if(~exist('action_mode', 'var'))
    action_mode = 0;
end
N_vec =   [2:2:100]; % a really long run .. (sent jobs)
try_inds = 25; % set one example 1:length(N_vec);
N_vec = N_vec(try_inds); % look at the N=28 example .. 1:length(N_vec2); % 1:5; % 1:length(N_vec2); N_vec2 = N_vec2(try_inds);
f_vec = 0.1; % [0.01 0.05 0.1 0.2 0.5]; % try different f_vecs

% intervals_struct.ratio_interval = [0 0.20]; % highlight any ratios below this threshold
% intervals_struct.h_interval = [0.25 0.8]; % allowable genetic heretability content
% intervals_struct.freq_interval = [0.01 0.10]; % allowable frequency of disease in the population
% intervals_struct.penetrance_interval = [0.2 0.8]; % allowable penetrance (risk for a twin)
% intervals_struct.h_add_interval = [0.00002 0.92]; % allowable fraction of variance explained by additive effects
% intervals_struct.lods_interval = [1 1.5]; % allowable lods-ratio for disease between 0 and 1 markers

disease_type_str = {'rare', 'common', 'high-heret', 'high-ratio'};
for i=1:1%4 % 2:length(disease_type_str)
    run_simulate_disease_pair_model(architecture_str_vec, ...
        action_mode, in_matlab_flag, N_vec, f_vec, disease_type_str{i}, queue_str);
end


