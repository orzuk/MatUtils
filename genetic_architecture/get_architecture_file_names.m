% Get standard names to files for saving architectures
%
% Input:
% good_architectures_dir
% architecture_str_vec
% N
% 
% Output:
% good_architectures_file 
% good_architectures_plot_file
% good_architectures_latex_file
%
function [good_architectures_file good_architectures_plot_file good_architectures_latex_file] = ...
    get_architecture_file_names(good_architectures_dir, architecture_str_vec, N)

good_architectures_file = fullfile(good_architectures_dir, ...
    'data', ['good_architectures_' cell2vec(architecture_str_vec, '-') ...
    '_N_' num2str(N) '.mat']);
good_architectures_plot_file = fullfile(good_architectures_dir, ...
    'figures', ['architecture_figures_' cell2vec(architecture_str_vec, '-')  ...
    '_N_' num2str(N)]);
good_architectures_latex_file = fullfile(good_architectures_dir, ...
    'docs', 'architectures', ['architecture_statistics_', ...
    cell2vec(architecture_str_vec, '_') '_N_' num2str(N) '.tex']);
